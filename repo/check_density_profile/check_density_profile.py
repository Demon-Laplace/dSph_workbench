import argparse
import contextlib
import gc
import io
import os
import re
import shutil
import subprocess
import sys
import time
from functools import lru_cache
from multiprocessing import get_all_start_methods, get_context
from pathlib import Path
from typing import Optional

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig-dsph")

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None

REPO_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = REPO_ROOT.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from basefunc import Analysis, DataProcessor, GalaxySimulation  # noqa: E402

CODE_VERSION = "v0.3.10"
SHAPE_SAMPLE_SIZE = 50_000
BOUNDARY_SAMPLE_SIZE = 50_000
MAX_PANEL_IMAGE_CELLS = 1_000_000


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Density-profile analysis for isolated dwarf snapshots. "
            "Supports both legacy batch mode and explicit single-snapshot mode."
        )
    )
    parser.add_argument(
        "--model-dir",
        help=(
            "Directory containing snapshot_XXX.hdf5 and dwarf ini file. "
            "If omitted, legacy batch mode uses the current working directory."
        ),
    )
    parser.add_argument("--snapshot", type=int, help="Snapshot number to analyze in explicit single-snapshot mode.")
    parser.add_argument(
        "--range",
        nargs="+",
        help=(
            'Snapshot selection for legacy batch mode. Use "0-300" for an inclusive range, '
            'or "0,300" to process only the listed snapshots. Defaults to all detected snapshots.'
        ),
    )
    parser.add_argument("--processes", type=int, default=1, help="Number of worker processes for legacy batch mode.")
    parser.add_argument("--radius-max", type=float, default=10.0, help="3D radial cut in kpc for analysis samples.")
    parser.add_argument("--nbin", type=int, default=100, help="Number of radial bins.")
    parser.add_argument("--sigma-radius", type=float, default=1.0, help="Projected radius in kpc for velocity-dispersion annotations.")
    parser.add_argument(
        "--output-dir",
        help=(
            "Directory for generated figure and summary files. "
            "Defaults to ./density_output in legacy batch mode."
        ),
    )
    parser.add_argument("--replace", action="store_true", help="Replace existing snapshot outputs instead of skipping them.")
    parser.add_argument(
        "--no-video",
        action="store_true",
        help="Do not assemble generated density-profile PNGs into an MP4 in legacy batch mode.",
    )
    parser.add_argument(
        "--video",
        action="store_true",
        help="Force MP4 assembly even for comma-separated snapshot lists.",
    )
    parser.add_argument(
        "--video-only",
        action="store_true",
        help="Assemble existing density-profile PNGs into an MP4 and exit without analyzing snapshots.",
    )
    parser.add_argument(
        "--video-scale",
        default="3416:1886",
        help='ffmpeg scale filter for the legacy batch MP4. Defaults to "3416:1886".',
    )
    parser.add_argument("--video-crf", type=int, default=18, help="ffmpeg CRF value for the legacy batch MP4.")
    return parser.parse_args()


def parse_snapshot_selection(range_text: str):
    if isinstance(range_text, (list, tuple)):
        range_text = " ".join(str(value) for value in range_text)
    range_text = range_text.strip()
    if not range_text:
        raise ValueError("Invalid --range value: empty selection.")

    if "," in range_text:
        if "-" in range_text:
            raise ValueError(f'Invalid --range value "{range_text}". Do not mix "," and "-".')
        try:
            selected = [int(value.strip()) for value in range_text.split(",") if value.strip()]
        except ValueError as exc:
            raise ValueError(f'Invalid --range value "{range_text}". Expected comma-separated integers.') from exc
        if not selected:
            raise ValueError(f'Invalid --range value "{range_text}". Expected at least one snapshot number.')
        return "list", selected

    if "-" in range_text and "," not in range_text:
        try:
            start, end = (int(value.strip()) for value in range_text.split("-", maxsplit=1))
        except ValueError as exc:
            raise ValueError(f'Invalid --range value "{range_text}". Expected format "start-end".') from exc
        if start > end:
            raise ValueError(f'Invalid --range value "{range_text}": start must be <= end.')
        return "interval", list(range(start, end + 1))

    try:
        return "list", [int(range_text)]
    except ValueError as exc:
        raise ValueError(
            f'Invalid --range value "{range_text}". Expected "start-end" or comma-separated integers.'
        ) from exc


def detect_batch_snapshots(model_dir: Path, range_text: Optional[str]):
    snapshot_numbers = DataProcessor.list_snapshot_numbers(model_dir)
    if not snapshot_numbers:
        raise FileNotFoundError(f"No snapshots found under {model_dir}")

    if range_text is None:
        return snapshot_numbers, "all"

    selection_mode, selected_numbers = parse_snapshot_selection(range_text)
    if selection_mode == "interval":
        start = selected_numbers[0]
        end = selected_numbers[-1]
        requested = [num for num in snapshot_numbers if start <= num <= end]
    else:
        available = set(snapshot_numbers)
        missing = [num for num in selected_numbers if num not in available]
        if missing:
            missing_text = ", ".join(f"{num:03d}" for num in missing)
            raise FileNotFoundError(f"Requested snapshot(s) not found under {model_dir}: {missing_text}")
        seen = set()
        requested = []
        for num in selected_numbers:
            if num not in seen:
                requested.append(num)
                seen.add(num)

    if not requested:
        raise FileNotFoundError(f"No snapshots in requested selection under {model_dir}: {range_text}")
    return requested, selection_mode


def resolve_legacy_model_dir(model_dir_arg):
    if model_dir_arg:
        return Path(model_dir_arg).resolve()

    legacy_output_dir = (Path.cwd() / "output").resolve()
    if legacy_output_dir.exists():
        return legacy_output_dir
    return Path.cwd().resolve()


def resolve_output_dir(output_dir_arg, explicit_mode: bool):
    if output_dir_arg:
        return Path(output_dir_arg).resolve()
    if explicit_mode:
        return (WORKSPACE_ROOT / "sandbox_runs" / "check_density_profile").resolve()
    return (Path.cwd() / "density_output").resolve()


class ProgressReporter:
    def __init__(self, total, desc="CheckDensity", unit="snap", mininterval=0.5):
        self.total = total
        self.desc = desc
        self.unit = unit
        self.mininterval = mininterval
        self.count = 0
        self._last_print = 0.0
        self._use_tqdm = total > 0 and tqdm is not None and sys.stdout.isatty()
        self._progress = None
        if self._use_tqdm:
            self._progress = tqdm(
                total=total,
                desc=desc,
                unit=unit,
                dynamic_ncols=True,
                mininterval=mininterval,
                file=sys.stdout,
            )
        elif total > 0:
            print(f"[{self.desc}Progress] 0/{self.total} {self.unit}", flush=True)

    def update(self, n=1):
        if self.total <= 0:
            return
        self.count += n
        if self._progress is not None:
            self._progress.update(n)
            return

        now = time.perf_counter()
        if self.count >= self.total or now - self._last_print >= self.mininterval:
            print(f"[{self.desc}Progress] {self.count}/{self.total} {self.unit}", flush=True)
            self._last_print = now

    def close(self):
        if self._progress is not None:
            self._progress.close()


def get_output_paths(model_dir: Path, snapshot: int, output_dir: Path):
    stem = f"snapshot_{snapshot:03d}_density_profile"
    return output_dir / f"{stem}.png"


def detect_frame_snapshots(output_dir: Path):
    snapshots = []
    for frame_path in output_dir.glob("snapshot_*_density_profile.png"):
        match = re.fullmatch(r"snapshot_(\d+)_density_profile\.png", frame_path.name)
        if match:
            snapshots.append(int(match.group(1)))
    return sorted(set(snapshots))


def snapshot_output_exists(model_dir: Path, snapshot: int, output_dir: Path):
    fig_path = get_output_paths(model_dir, snapshot, output_dir)
    return fig_path.exists()


def sanitize_filename_part(text: str):
    clean = re.sub(r"[^A-Za-z0-9_.-]+", "", str(text))
    return clean or "unknown"


def detect_model_number(model_dir: Path):
    candidates = [model_dir.name, model_dir.parent.name]
    if model_dir.name == "output":
        candidates.insert(0, model_dir.parent.name)
    for candidate in candidates:
        match = re.search(r"(\d+)$", candidate)
        if match:
            return match.group(1)
    return ""


def build_video_output_path(output_dir: Path, dwarf_name: str, model_dir: Path):
    dwarf_part = sanitize_filename_part(dwarf_name)
    number_part = detect_model_number(model_dir)
    return output_dir / f"IC_density_{dwarf_part}{number_part}_{CODE_VERSION}.mp4"


def create_density_video(
    output_dir: Path,
    model_dir: Path,
    dwarf_name: str,
    snapshots,
    scale: str,
    crf: int,
):
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg is None:
        print("Skipped MP4 assembly: ffmpeg was not found on PATH.")
        return None

    frame_paths = [get_output_paths(model_dir, snapshot, output_dir) for snapshot in sorted(snapshots)]
    frame_paths = [path for path in frame_paths if path.exists()]
    if not frame_paths:
        print("Skipped MP4 assembly: no density-profile PNG frames were found.")
        return None

    output_path = build_video_output_path(output_dir, dwarf_name=dwarf_name, model_dir=model_dir)
    concat_path = output_dir / f".{output_path.stem}_frames.txt"
    with concat_path.open("w") as f:
        for frame_path in frame_paths:
            escaped = str(frame_path.resolve()).replace("'", r"'\''")
            f.write(f"file '{escaped}'\n")

    cmd = [
        ffmpeg,
        "-loglevel",
        "error",
        "-f",
        "concat",
        "-safe",
        "0",
        "-i",
        str(concat_path),
        "-vf",
        f"scale={scale}",
        "-crf",
        str(crf),
        "-pix_fmt",
        "yuv420p",
        "-y",
        str(output_path),
    ]
    try:
        subprocess.run(cmd, cwd=output_dir, check=True)
    except subprocess.CalledProcessError as exc:
        print(f"Skipped MP4 assembly: ffmpeg failed with exit code {exc.returncode}.")
        return None
    finally:
        try:
            concat_path.unlink()
        except OSError:
            pass

    return output_path


def summarize_snapshot_numbers(numbers, limit=12):
    values = sorted(numbers)
    if not values:
        return "none"
    shown = ", ".join(f"{value:03d}" for value in values[:limit])
    if len(values) > limit:
        shown = f"{shown}, ... (+{len(values) - limit} more)"
    return shown


def detect_dwarf_name_for_video(model_dir: Path, snapshots, results=None):
    if results:
        for result in sorted(results, key=lambda item: item["snapshot"]):
            dwarf_name = result.get("dwarf_name")
            if dwarf_name:
                return dwarf_name

    if snapshots:
        _, dwarf_name, _ = resolve_snapshot_context(model_dir, snapshots[0])
        return dwarf_name

    return detect_model_name(model_dir)


def detect_dwarf_name_for_video_safe(model_dir: Path, snapshots, results=None):
    try:
        return detect_dwarf_name_for_video(model_dir, snapshots, results=results)
    except Exception as exc:
        print(f"Video assembly: using model directory name for MP4 name after dwarf-name detection failed: {exc}")
        return detect_model_name(model_dir)


def assemble_density_video_with_log(
    output_dir: Path,
    model_dir: Path,
    snapshots,
    scale: str,
    crf: int,
    results=None,
):
    frame_count = sum(1 for snapshot in snapshots if get_output_paths(model_dir, snapshot, output_dir).exists())
    print(
        f"Video assembly: found {frame_count}/{len(snapshots)} density-profile PNG frame(s) "
        f"in {output_dir}"
    )
    dwarf_name = detect_dwarf_name_for_video_safe(model_dir, snapshots, results=results)
    output_path = build_video_output_path(output_dir, dwarf_name=dwarf_name, model_dir=model_dir)
    print(f"Video assembly: target MP4: {output_path}")
    video_path = create_density_video(
        output_dir=output_dir,
        model_dir=model_dir,
        dwarf_name=dwarf_name,
        snapshots=snapshots,
        scale=scale,
        crf=crf,
    )
    if video_path is not None:
        print(f"Saved video: {video_path}")
    return video_path


def check_video_batch_result(
    output_dir: Path,
    model_dir: Path,
    snapshots,
    make_video: bool,
    selection_mode: str,
    no_video: bool,
    forced_video: bool,
    video_path: Optional[Path] = None,
    results=None,
):
    snapshots = sorted(snapshots)
    frame_paths = [get_output_paths(model_dir, snapshot, output_dir) for snapshot in snapshots]
    existing_frames = [path for path in frame_paths if path.exists()]
    missing_snapshots = [
        snapshot for snapshot, frame_path in zip(snapshots, frame_paths) if not frame_path.exists()
    ]

    if no_video:
        reason = "--no-video"
    elif selection_mode == "list" and not forced_video:
        reason = 'comma-separated --range without --video'
    elif make_video:
        reason = "enabled"
    else:
        reason = "disabled"

    dwarf_name = detect_dwarf_name_for_video_safe(model_dir, snapshots, results=results)
    expected_video = build_video_output_path(output_dir, dwarf_name=dwarf_name, model_dir=model_dir)
    actual_video = Path(video_path) if video_path is not None else expected_video
    video_exists = actual_video.exists()
    video_size = actual_video.stat().st_size if video_exists else 0
    ffmpeg = shutil.which("ffmpeg")
    mp4_candidates = sorted(output_dir.glob("*.mp4")) if output_dir.exists() else []

    print(
        f"[VideoCheck] enabled={make_video}, reason={reason}, ffmpeg={ffmpeg or 'not found'}, "
        f"output_dir={output_dir}"
    )
    print(
        f"[VideoCheck] frames={len(existing_frames)}/{len(snapshots)}, "
        f"missing={summarize_snapshot_numbers(missing_snapshots)}"
    )
    print(f"[VideoCheck] expected_mp4={expected_video}, exists={video_exists}, size_bytes={video_size}")
    if not video_exists and mp4_candidates:
        candidates = ", ".join(str(path) for path in mp4_candidates[-5:])
        print(f"[VideoCheck] existing_mp4_candidates={candidates}")
    if make_video and not video_exists:
        print("[VideoCheck] problem: video assembly was enabled but the expected MP4 was not created.")


def mass_weighted_center(x, y, z, m):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    m = np.asarray(m, dtype=float)
    if x.size == 0 or m.size == 0 or np.sum(m) <= 0:
        return 0.0, 0.0, 0.0
    return np.average(x, weights=m), np.average(y, weights=m), np.average(z, weights=m)


def safe_std(values):
    values = np.asarray(values, dtype=float)
    return float(np.std(values)) if values.size > 0 else float("nan")


def projected_half_mass_radius(x, z, m):
    x = np.asarray(x, dtype=float)
    z = np.asarray(z, dtype=float)
    m = np.asarray(m, dtype=float)
    if x.size == 0 or z.size == 0 or m.size == 0:
        return float("nan")
    r_proj = np.sqrt(x**2 + z**2)
    sort_idx = np.argsort(r_proj)
    r_sorted = r_proj[sort_idx]
    m_sorted = m[sort_idx]
    m_cum = np.cumsum(m_sorted)
    if m_cum[-1] <= 0:
        return float("nan")
    m_half = m_cum[-1] / 2.0
    idx_half = np.searchsorted(m_cum, m_half)
    return float(r_sorted[idx_half] if idx_half < len(r_sorted) else np.nan)


def find_ini_file(search_dirs):
    for search_dir in search_dirs:
        if search_dir is None or not search_dir.exists() or not search_dir.is_dir():
            continue
        ini_files = sorted(p for p in search_dir.glob("*.ini") if not p.name.startswith("IC_"))
        if not ini_files:
            ini_files = sorted(search_dir.glob("*.ini"))
        if ini_files:
            ini_file = ini_files[0]
            return ini_file.stem, ini_file
    searched = ", ".join(str(path) for path in search_dirs if path is not None)
    raise FileNotFoundError(f"No .ini file found in any of: {searched}")


def infer_dwarf_name_without_ini(model_dir: Path, snapshot_dir: Path):
    for candidate in (snapshot_dir, model_dir, snapshot_dir.parent, model_dir.parent):
        if candidate.name and candidate.name != "output":
            return candidate.name
    return "dwarf"


def resolve_snapshot_dir(model_dir: Path, snapshot_num: int):
    snapshot_name = f"snapshot_{snapshot_num:03d}.hdf5"
    if (model_dir / snapshot_name).exists():
        return model_dir

    candidates = []
    for subdir in sorted(p for p in model_dir.iterdir() if p.is_dir()):
        if (subdir / snapshot_name).exists():
            candidates.append(subdir)

    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        raise ValueError(f"Multiple model subdirectories under {model_dir} contain {snapshot_name}: {candidates}")
    raise FileNotFoundError(f"Could not find a snapshot directory under {model_dir} containing {snapshot_name}")


def resolve_snapshot_context(model_dir: Path, snapshot_num: int):
    snapshot_dir = resolve_snapshot_dir(model_dir, snapshot_num)
    try:
        dwarf_name, ini_file = find_ini_file(
            [
                snapshot_dir,
                snapshot_dir.parent,
                model_dir,
                model_dir.parent,
            ]
        )
    except FileNotFoundError:
        dwarf_name = infer_dwarf_name_without_ini(model_dir, snapshot_dir)
        ini_file = None
    return snapshot_dir, dwarf_name, ini_file


def select_dwarf_mask(df: pd.DataFrame, simulation: GalaxySimulation, dw_particles: int):
    if dw_particles <= 0 or dw_particles >= len(df):
        return np.ones(len(df), dtype=bool)

    dwarf_mask = simulation.find_dwarf_particles(dw_particles_num=dw_particles).to_numpy()
    if dwarf_mask.sum() == 0:
        return np.ones(len(df), dtype=bool)

    tp = df["tp"].to_numpy()
    dwarf_star_count = int(np.count_nonzero(dwarf_mask & ((tp == 2) | (tp == 4))))
    if dwarf_star_count < 128:
        return np.ones(len(df), dtype=bool)

    return dwarf_mask


def load_centered_dwarf_dataframe(model_dir: Path, snapshot_num: int):
    snapshot_dir, dwarf_name, ini_file = resolve_snapshot_context(model_dir, snapshot_num)
    simulation = GalaxySimulation(folder_path=str(snapshot_dir) + "/", snapshot_num=snapshot_num)
    # Keep legacy CheckDensityProfile behavior: isolated-model snapshots here should not apply box centering.
    with contextlib.redirect_stdout(io.StringIO()):
        df, tsnap = simulation.load_snapshot(
            dftype=2,
            box=False,
            gas_fields=("temp",),
            star_fields=(),
            particle_types=(0, 2, 3, 4),
        )

    dwarf_df = df
    tp = dwarf_df["tp"].to_numpy(copy=False)
    x = dwarf_df["x"].to_numpy(copy=False)
    y = dwarf_df["y"].to_numpy(copy=False)
    z = dwarf_df["z"].to_numpy(copy=False)
    m = dwarf_df["m"].to_numpy(copy=False)
    dwarf_star_mask = (tp == 2) | (tp == 4)
    dwarf_gas_mask = tp == 0
    dwarf_total_mask = dwarf_star_mask | dwarf_gas_mask

    xc, yc, zc = mass_weighted_center(
        x[dwarf_star_mask],
        y[dwarf_star_mask],
        z[dwarf_star_mask],
        m[dwarf_star_mask],
    )

    x -= np.float32(xc)
    y -= np.float32(yc)
    z -= np.float32(zc)

    for velocity_axis in ("vx", "vy", "vz"):
        velocity = dwarf_df[velocity_axis].to_numpy(copy=False)
        if np.any(dwarf_total_mask):
            velocity -= np.float32(velocity[dwarf_total_mask].mean(dtype=np.float64))

    r3d = np.empty_like(x)
    scratch = np.empty_like(x)
    np.multiply(x, x, out=r3d)
    np.multiply(y, y, out=scratch)
    r3d += scratch
    np.multiply(z, z, out=scratch)
    r3d += scratch
    np.sqrt(r3d, out=r3d)
    dwarf_df.loc[:, "r3d"] = r3d

    return dwarf_name, snapshot_dir, tsnap, dwarf_df


def radial_surface_density(radius, mass, bins):
    radius = np.asarray(radius, dtype=float)
    mass = np.asarray(mass, dtype=float)
    mass_in_bin, _ = np.histogram(radius, bins=bins, weights=mass)
    area = np.pi * (bins[1:] ** 2 - bins[:-1] ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        sigma = np.where(area > 0, mass_in_bin / area, np.nan)
    sigma_h = sigma * 1.25e-6
    return sigma, sigma_h


def safe_boundary_path(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 3 or y.size < 3:
        return None
    if np.unique(x).size < 2 or np.unique(y).size < 2:
        return None
    try:
        return Analysis.estimate_galaxy_boundary(x, y)
    except Exception:
        return None


def sample_arrays(arrays, max_points, *, seed=0, weight_array=None):
    if not arrays:
        return tuple()

    size = len(arrays[0])
    if size <= max_points:
        return tuple(np.asarray(arr) for arr in arrays)

    rng = np.random.default_rng(seed)
    if weight_array is not None:
        weights = np.asarray(weight_array, dtype=float)
        valid = np.isfinite(weights) & (weights > 0)
        if np.any(valid):
            prob = np.zeros(size, dtype=float)
            prob[valid] = weights[valid]
            prob /= prob.sum()
            indices = rng.choice(size, size=max_points, replace=False, p=prob)
        else:
            indices = rng.choice(size, size=max_points, replace=False)
    else:
        indices = rng.choice(size, size=max_points, replace=False)

    indices.sort()
    return tuple(np.asarray(arr)[indices] for arr in arrays)


def estimate_stellar_shape(x, y, mass, n_neighbors=30):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mass = np.asarray(mass, dtype=float)
    if x.size == 0 or y.size == 0 or mass.size == 0:
        return np.nan, np.nan

    x_eval, y_eval, mass_eval = sample_arrays(
        (x, y, mass),
        SHAPE_SAMPLE_SIZE,
        seed=0,
        weight_array=mass,
    )
    return Analysis.calculate_ellipticity(x_eval, y_eval, mass=mass_eval, n_neighbors=n_neighbors)


def safe_boundary_path_sampled(x, y, mass=None):
    arrays = (x, y) if mass is None else (x, y, mass)
    sampled = sample_arrays(
        arrays,
        BOUNDARY_SAMPLE_SIZE,
        seed=1,
        weight_array=mass,
    )
    return safe_boundary_path(sampled[0], sampled[1])


def clamp_panel_bins(bins):
    bins = int(bins)
    if bins <= 1:
        return 1
    if bins * bins <= MAX_PANEL_IMAGE_CELLS:
        return bins
    return max(32, int(np.sqrt(MAX_PANEL_IMAGE_CELLS)))


def smooth_histogram2d(hist, passes=2):
    kernel = np.array([1.0, 4.0, 6.0, 4.0, 1.0], dtype=float)
    kernel /= kernel.sum()
    smoothed = np.asarray(hist, dtype=float)

    for _ in range(passes):
        padded_x = np.pad(smoothed, ((0, 0), (2, 2)), mode="edge")
        smoothed = (
            kernel[0] * padded_x[:, :-4]
            + kernel[1] * padded_x[:, 1:-3]
            + kernel[2] * padded_x[:, 2:-2]
            + kernel[3] * padded_x[:, 3:-1]
            + kernel[4] * padded_x[:, 4:]
        )
        padded_y = np.pad(smoothed, ((2, 2), (0, 0)), mode="edge")
        smoothed = (
            kernel[0] * padded_y[:-4, :]
            + kernel[1] * padded_y[1:-3, :]
            + kernel[2] * padded_y[2:-2, :]
            + kernel[3] * padded_y[3:-1, :]
            + kernel[4] * padded_y[4:, :]
        )

    return smoothed


PANEL_EXTENT_KPC = 9.4
GAS_PANEL_BINS = 600
STELLAR_PANEL_BINS = 1500
GAS_PANEL_VMIN = 3.6
GAS_PANEL_VMAX = 5.0
STELLAR_PANEL_VMIN = -2.5
STELLAR_PANEL_VMAX = 5.0

ASTRO_GAS_CMAP = LinearSegmentedColormap.from_list(
    "astro_gas",
    [
        (0.00, "#1d234f"),
        (0.10, "#383071"),
        (0.22, "#6b3c86"),
        (0.36, "#a24a78"),
        (0.50, "#d06055"),
        (0.64, "#ee9148"),
        (0.78, "#f6c85f"),
        (0.90, "#fff0a8"),
        (0.97, "#fffbe5"),
        (1.00, "#ffffff"),
    ],
)

ASTRO_STELLAR_CMAP = ListedColormap(
    [
        "#070707",
        "#0d1018",
        "#111a2d",
        "#18345c",
        "#2b5f92",
        "#5d8fb6",
        "#99bfd3",
        "#c5d8dd",
        "#dbe2df",
        "#ede8d8",
        "#faf4e8",
        "#ffffff",
    ],
    "astro_stellar",
)

ASTRO_GAS_CMAP.set_bad((0.0, 0.0, 0.0, 0.0))
ASTRO_STELLAR_CMAP.set_bad((0.0, 0.0, 0.0, 0.0))


def plot_mass_density_panel(
    ax,
    x,
    y,
    mass,
    title,
    xlabel,
    ylabel,
    cmap,
    colorbar_label,
    no_mass_text,
    no_view_text,
    extent=PANEL_EXTENT_KPC,
    bins=GAS_PANEL_BINS,
    vmin=None,
    vmax=None,
    smooth_passes=2,
    show_colorbar=False,
):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mass = np.asarray(mass, dtype=float)
    bins = clamp_panel_bins(bins)

    ax.set_facecolor("black")
    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="box")
    ax.tick_params(axis="both", which="both", direction="in", color="white", labelcolor="black")
    ax.xaxis.label.set_color("black")
    ax.yaxis.label.set_color("black")
    ax.title.set_color("black")
    for spine in ax.spines.values():
        spine.set_color("white")

    if x.size == 0 or y.size == 0 or mass.size == 0 or np.sum(mass) <= 0:
        ax.text(
            0.5,
            0.5,
            no_mass_text,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=9,
            color="white",
            bbox=dict(facecolor="black", alpha=0.85, edgecolor="white"),
        )
        return

    view_mask = (np.abs(x) <= extent) & (np.abs(y) <= extent)
    x_view = x[view_mask]
    y_view = y[view_mask]
    mass_view = mass[view_mask]
    if x_view.size == 0 or np.sum(mass_view) <= 0:
        ax.text(
            0.5,
            0.5,
            no_view_text,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=9,
            color="white",
            bbox=dict(facecolor="black", alpha=0.85, edgecolor="white"),
        )
        return

    x_edges = np.linspace(-extent, extent, bins + 1)
    y_edges = np.linspace(-extent, extent, bins + 1)
    hist, _, _ = np.histogram2d(x_view, y_view, bins=[x_edges, y_edges], weights=mass_view)
    hist = smooth_histogram2d(hist.T, passes=smooth_passes)
    positive = hist > 0

    if not np.any(positive):
        ax.text(
            0.5,
            0.5,
            no_view_text,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=9,
            color="white",
            bbox=dict(facecolor="black", alpha=0.85, edgecolor="white"),
        )
        return

    log_hist = np.full_like(hist, np.nan, dtype=float)
    log_hist[positive] = np.log10(hist[positive])
    values = log_hist[positive]
    plot_vmin = vmin
    plot_vmax = vmax
    if plot_vmin is None or plot_vmax is None:
        auto_vmin = float(np.percentile(values, 8.0))
        auto_vmax = float(np.percentile(values, 99.5))
        if not np.isfinite(auto_vmin) or not np.isfinite(auto_vmax) or auto_vmin >= auto_vmax:
            auto_vmin = float(np.nanmin(values))
            auto_vmax = float(np.nanmax(values))
        if not np.isfinite(auto_vmin) or not np.isfinite(auto_vmax) or auto_vmin >= auto_vmax:
            auto_vmax = float(values[0])
            auto_vmin = auto_vmax - 1.0
        plot_vmin = auto_vmin if plot_vmin is None else plot_vmin
        plot_vmax = auto_vmax if plot_vmax is None else plot_vmax

    image = ax.imshow(
        log_hist,
        origin="lower",
        extent=[-extent, extent, -extent, extent],
        cmap=cmap,
        vmin=plot_vmin,
        vmax=plot_vmax,
        interpolation="bicubic",
    )

    if show_colorbar:
        colorbar = plt.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
        colorbar.outline.set_edgecolor("black")
        colorbar.ax.yaxis.set_tick_params(color="black")
        plt.setp(colorbar.ax.get_yticklabels(), color="black")
        colorbar.set_label(colorbar_label)
        colorbar.ax.yaxis.label.set_color("black")


def plot_gas_density_panel(ax, x, y, mass, title, xlabel, ylabel):
    plot_mass_density_panel(
        ax,
        x,
        y,
        mass,
        title,
        xlabel,
        ylabel,
        cmap=ASTRO_GAS_CMAP,
        colorbar_label=r"$\log_{10}(\mathrm{smoothed\ gas\ mass\ per\ pixel})$",
        no_mass_text="No cold gas",
        no_view_text="No gas mass in view",
        vmin=GAS_PANEL_VMIN,
        vmax=GAS_PANEL_VMAX,
    )


def plot_stellar_density_panel(ax, x, y, mass, title, xlabel, ylabel):
    plot_mass_density_panel(
        ax,
        x,
        y,
        mass,
        title,
        xlabel,
        ylabel,
        cmap=ASTRO_STELLAR_CMAP,
        colorbar_label=r"$\log_{10}(\mathrm{smoothed\ stellar\ mass\ per\ pixel})$",
        no_mass_text="No stars",
        no_view_text="No stellar mass in view",
        bins=STELLAR_PANEL_BINS,
        vmin=STELLAR_PANEL_VMIN,
        vmax=STELLAR_PANEL_VMAX,
        smooth_passes=3,
    )


def compute_summary(dwarf_df: pd.DataFrame, radius_max: float, nbin: int, sigma_radius: float):
    x = dwarf_df["x"].to_numpy()
    y = dwarf_df["y"].to_numpy()
    z = dwarf_df["z"].to_numpy()
    vx = dwarf_df["vx"].to_numpy()
    vy = dwarf_df["vy"].to_numpy()
    vz = dwarf_df["vz"].to_numpy()
    m = dwarf_df["m"].to_numpy()
    tp = dwarf_df["tp"].to_numpy()
    temp = dwarf_df["temp"].to_numpy()
    r3d = dwarf_df["r3d"].to_numpy()

    star_mask = (tp == 2) | (tp == 4)
    gas_mask = tp == 0
    cold_gas_mask = gas_mask & (temp < 20000.0)
    hot_gas_mask = gas_mask & ~cold_gas_mask

    star_within = star_mask & (r3d <= radius_max)
    cold_within = cold_gas_mask & (r3d <= radius_max)

    xs = x[star_within]
    ys = y[star_within]
    zs = z[star_within]
    ms = m[star_within]

    xcg = x[cold_within]
    ycg = y[cold_within]
    mcg = m[cold_within]

    r_star = np.hypot(xs, ys)
    r_cold = np.hypot(xcg, ycg)

    bins = np.linspace(0.0, radius_max, nbin + 1)
    r_cent = 0.5 * (bins[:-1] + bins[1:])
    sigma_star, sigma_h_star = radial_surface_density(r_star, ms, bins)
    sigma_cold, sigma_h_cold = radial_surface_density(r_cold, mcg, bins)

    ellipticity, pa = estimate_stellar_shape(xs, zs, ms, n_neighbors=30)
    rotate = np.array([[np.cos(-pa), -np.sin(-pa)], [np.sin(-pa), np.cos(-pa)]])
    coords = np.vstack((xs, zs))
    x_rot, z_rot = (rotate @ coords)
    r_half = Analysis.half_light_radius(x_rot, z_rot, ms, ep=ellipticity)

    inner_mask = r3d <= 5.0
    subset_star_mask = star_mask & inner_mask
    subset_cold_mask = cold_gas_mask & inner_mask
    subset_hot_mask = hot_gas_mask & inner_mask

    subset_star_x = x[subset_star_mask]
    subset_star_y = y[subset_star_mask]
    subset_star_z = z[subset_star_mask]
    subset_star_m = m[subset_star_mask]
    subset_cold_x = x[subset_cold_mask]
    subset_cold_z = z[subset_cold_mask]
    subset_cold_m = m[subset_cold_mask]
    subset_hot_x = x[subset_hot_mask]
    subset_hot_z = z[subset_hot_mask]
    subset_hot_m = m[subset_hot_mask]

    stellar_region_star_x = xs
    stellar_region_star_z = zs
    stellar_region_star_m = ms
    stellar_region_gas_mask = gas_mask & (r3d <= radius_max)
    stellar_region_gas_x = x[stellar_region_gas_mask]
    stellar_region_gas_z = z[stellar_region_gas_mask]
    stellar_region_gas_m = m[stellar_region_gas_mask]

    boundary_path = safe_boundary_path_sampled(xs, zs, mass=ms)
    if boundary_path is not None and stellar_region_star_x.size > 0:
        star_in = boundary_path.contains_points(np.column_stack([stellar_region_star_x, stellar_region_star_z]))
    else:
        star_in = np.ones(stellar_region_star_x.size, dtype=bool)

    m_star_boundary = float(stellar_region_star_m[star_in].sum())
    m_cold_gas = float(subset_cold_m.sum())
    m_hot_gas = float(subset_hot_m.sum())

    q = 1.0 - ellipticity
    if not np.isfinite(q) or q <= 0:
        q = np.finfo(float).eps

    def stellar_ellipse_mask(x_values, z_values):
        if x_values.size == 0 or not np.isfinite(r_half):
            return np.zeros(x_values.size, dtype=bool)
        x_ell, z_ell = rotate @ np.vstack((x_values, z_values))
        return np.sqrt(x_ell**2 + (z_ell / q) ** 2) <= r_half

    star_rhalf_mask = stellar_ellipse_mask(subset_star_x, subset_star_z)
    cold_rhalf_mask = stellar_ellipse_mask(subset_cold_x, subset_cold_z)
    hot_rhalf_mask = stellar_ellipse_mask(subset_hot_x, subset_hot_z)

    m_star_rhalf = float(subset_star_m[star_rhalf_mask].sum())
    m_cold_rhalf = float(subset_cold_m[cold_rhalf_mask].sum())
    m_hot_rhalf = float(subset_hot_m[hot_rhalf_mask].sum())

    if boundary_path is not None:
        if stellar_region_gas_x.size > 0:
            gas_stellar_region_mask = boundary_path.contains_points(
                np.column_stack([stellar_region_gas_x, stellar_region_gas_z])
            )
            m_gas_stellar_region = float(stellar_region_gas_m[gas_stellar_region_mask].sum())
        else:
            m_gas_stellar_region = 0.0
        m_star_stellar_region = m_star_boundary
    else:
        gas_stellar_region_mask = stellar_ellipse_mask(stellar_region_gas_x, stellar_region_gas_z)
        star_stellar_region_mask = stellar_ellipse_mask(stellar_region_star_x, stellar_region_star_z)
        m_gas_stellar_region = float(stellar_region_gas_m[gas_stellar_region_mask].sum())
        m_star_stellar_region = float(stellar_region_star_m[star_stellar_region_mask].sum())

    stellar_region_denom = m_star_stellar_region + m_gas_stellar_region
    gas_fraction_stellar_region = (
        float(m_gas_stellar_region / stellar_region_denom) if stellar_region_denom > 0 else 0.0
    )

    all_star_x = x[star_mask]
    all_star_y = y[star_mask]
    all_star_z = z[star_mask]
    all_star_vx = vx[star_mask]
    all_star_vy = vy[star_mask]
    all_star_vz = vz[star_mask]

    sigma_z = safe_std(all_star_vz[np.hypot(all_star_x, all_star_y) < sigma_radius])
    sigma_y = safe_std(all_star_vy[np.hypot(all_star_x, all_star_z) < sigma_radius])
    sigma_x = safe_std(all_star_vx[np.hypot(all_star_y, all_star_z) < sigma_radius])

    stellar_mass_5kpc = float(m[star_mask & (r3d < 5.0)].sum())
    total_gas = float(m[gas_mask].sum())
    n_gas = int(np.count_nonzero(gas_mask))
    n_star = int(np.count_nonzero(tp >= 2))

    gas_fraction_denom = m_cold_gas + m_hot_gas + stellar_mass_5kpc
    gas_fraction = float(m_cold_gas / gas_fraction_denom) if gas_fraction_denom > 0 else 0.0
    gas_fraction_rhalf_denom = m_cold_rhalf + m_hot_rhalf + m_star_rhalf
    gas_fraction_rhalf = (
        float((m_cold_rhalf + m_hot_rhalf) / gas_fraction_rhalf_denom) if gas_fraction_rhalf_denom > 0 else 0.0
    )

    area_inner = np.pi * (1.0**2)
    sigma_h_inner_gas = float(mcg[r_cold <= 1.0].sum() / area_inner * 1.25e-6) if mcg.size else 0.0
    sigma_h_inner_star = float(ms[r_star <= 1.0].sum() / area_inner * 1.25e-6) if ms.size else 0.0

    return {
        "r_cent": r_cent,
        "sigma_h_cold": sigma_h_cold,
        "sigma_h_star": sigma_h_star,
        "ellipticity": ellipticity,
        "pa": pa,
        "r_half": r_half,
        "m_star_boundary": m_star_boundary,
        "m_cold_gas": m_cold_gas,
        "m_hot_gas": m_hot_gas,
        "m_star_rhalf": m_star_rhalf,
        "m_cold_rhalf": m_cold_rhalf,
        "m_hot_rhalf": m_hot_rhalf,
        "m_star_stellar_region": m_star_stellar_region,
        "m_gas_stellar_region": m_gas_stellar_region,
        "gas_fraction_stellar_region": gas_fraction_stellar_region,
        "sigma_x": sigma_x,
        "sigma_y": sigma_y,
        "sigma_z": sigma_z,
        "stellar_mass_5kpc": stellar_mass_5kpc,
        "total_gas": total_gas,
        "n_gas": n_gas,
        "n_star": n_star,
        "gas_fraction": gas_fraction,
        "gas_fraction_rhalf": gas_fraction_rhalf,
        "sigma_h_inner_gas": sigma_h_inner_gas,
        "sigma_h_inner_star": sigma_h_inner_star,
    }


def detect_model_name(model_dir: Path):
    if model_dir.name == "output":
        return model_dir.parent.name
    return model_dir.name


@lru_cache(maxsize=1)
def load_run_metadata():
    fd_value = None
    sfr_eff = None
    try:
        fd_value = DataProcessor.read_feedback_value()
    except Exception:
        pass
    try:
        sfr_eff = DataProcessor.read_sfreff_value()
    except Exception:
        pass
    return fd_value, sfr_eff


@lru_cache(maxsize=32)
def load_wlm_reference(model_dir: Path):
    search_dirs = [
        Path.cwd().parent,
        model_dir.parent,
        model_dir,
    ]
    for search_dir in search_dirs:
        radiusbin_path = search_dir / "radiusbin.npy"
        flux_path = search_dir / "flux.npy"
        if radiusbin_path.exists() and flux_path.exists():
            return np.load(radiusbin_path), np.load(flux_path)
    return None, None


def write_density_snapshot(output_dir: Path, snapshot: int, tsnap: float, summary: dict):
    density_file = output_dir / f"density_snapshot_{snapshot:03d}.txt"
    np.savetxt(
        density_file,
        [[snapshot, tsnap, summary["sigma_h_inner_gas"], summary["sigma_h_inner_star"]]],
        header="snapshot time gas_density stellar_density",
        fmt="%d %.6f %.6e %.6e",
    )


def load_density_history(output_dir: Path, max_snapshot: Optional[int] = None, max_time: Optional[float] = None):
    time_tol = 1e-9

    def record_within_limits(snapshot_value, time_value):
        if max_snapshot is not None and snapshot_value is not None and snapshot_value > max_snapshot:
            return False
        if max_time is not None and time_value is not None and time_value > max_time + time_tol:
            return False
        return True

    records = {}
    evo_path = output_dir / "density_evolution.csv"
    if evo_path.exists():
        try:
            df_evo = pd.read_csv(evo_path, comment="#")
            for _, row in df_evo.iterrows():
                snapshot_value = None
                if "snapshot" in df_evo.columns and not pd.isna(row["snapshot"]):
                    snapshot_value = int(row["snapshot"])
                time_value = float(row["time"])
                if not record_within_limits(snapshot_value, time_value):
                    continue
                record_key = (
                    ("snapshot", snapshot_value)
                    if snapshot_value is not None
                    else ("time", round(time_value, 10))
                )
                records[record_key] = (
                    snapshot_value,
                    time_value,
                    float(row["gas_density"]),
                    float(row["stellar_density"]),
                )
        except Exception:
            pass

    for density_file in sorted(output_dir.glob("density_snapshot_*.txt")):
        match = re.fullmatch(r"density_snapshot_(\d+)\.txt", density_file.name)
        file_snapshot = int(match.group(1)) if match else None
        if max_snapshot is not None and file_snapshot is not None and file_snapshot > max_snapshot:
            continue
        try:
            with density_file.open() as f:
                non_empty_lines = 0
                for line in f:
                    if line.strip():
                        non_empty_lines += 1
                    if non_empty_lines > 1:
                        break
            if non_empty_lines <= 1:
                continue
        except OSError:
            continue

        try:
            data = np.loadtxt(density_file, skiprows=1)
        except Exception:
            continue
        if np.size(data) == 0:
            continue
        data = np.atleast_2d(data)
        for row in data:
            if row.size < 3:
                continue
            try:
                if row.size >= 4:
                    snapshot_value = int(row[0])
                    time_value = float(row[1])
                    gas_density = float(row[2])
                    stellar_density = float(row[3])
                else:
                    snapshot_value = file_snapshot
                    time_value = float(row[0])
                    gas_density = float(row[1])
                    stellar_density = float(row[2])
                if not record_within_limits(snapshot_value, time_value):
                    continue
                record_key = (
                    ("snapshot", snapshot_value)
                    if snapshot_value is not None
                    else ("time", round(time_value, 10))
                )
                records[record_key] = (snapshot_value, time_value, gas_density, stellar_density)
            except (TypeError, ValueError):
                continue

    if not records:
        return None

    all_data = np.array(sorted(records.values(), key=lambda item: item[1]), dtype=float)
    return {
        "snapshot": all_data[:, 0],
        "time": all_data[:, 1],
        "gas_density": all_data[:, 2],
        "stellar_density": all_data[:, 3],
    }


def create_density_evolution_csv(output_dir: Path, fd_value=None, sfr_eff=None):
    history = load_density_history(output_dir)
    if history is None:
        return

    evo_path = output_dir / "density_evolution.csv"
    with evo_path.open("w") as f:
        if fd_value is not None:
            f.write(f"#Timescale_fd: {fd_value}\n")
        if sfr_eff is not None:
            f.write(f"#sfr_effiency: {sfr_eff}\n")
        f.write("snapshot,time,gas_density,stellar_density\n")
        data_to_save = np.column_stack(
            [history["snapshot"], history["time"], history["gas_density"], history["stellar_density"]]
        )
        np.savetxt(f, data_to_save, delimiter=",", fmt="%.0f,%.6f,%.6e,%.6e")


def cleanup_density_snapshots(output_dir: Path):
    for density_file in output_dir.glob("density_snapshot_*.txt"):
        try:
            density_file.unlink()
        except Exception:
            pass


def make_figure(dwarf_df: pd.DataFrame, tsnap: float, summary: dict, density_history=None, wlm_r=None, wlm_f=None, title_text=None):
    x = dwarf_df["x"].to_numpy()
    y = dwarf_df["y"].to_numpy()
    z = dwarf_df["z"].to_numpy()
    tp = dwarf_df["tp"].to_numpy()
    temp = dwarf_df["temp"].to_numpy()
    star_mask = (tp == 2) | (tp == 4)
    cold_gas_mask = (tp == 0) & (temp < 20000.0)
    mass = dwarf_df["m"].to_numpy()

    star_x = x[star_mask]
    star_y = y[star_mask]
    star_z = z[star_mask]
    star_mass = mass[star_mask]
    cold_x = x[cold_gas_mask]
    cold_y = y[cold_gas_mask]
    cold_z = z[cold_gas_mask]
    cold_mass = mass[cold_gas_mask]

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(16.0, 9.2),
        dpi=180,
        gridspec_kw={"width_ratios": [1.0, 1.0, 1.0], "wspace": 0.14, "hspace": 0.26},
    )

    ax = axes[0, 0]
    plot_gas_density_panel(
        ax,
        cold_x,
        cold_y,
        cold_mass,
        title="Cold Gas Distribution (x-y)",
        xlabel="x [kpc]",
        ylabel="y [kpc]",
    )

    ax = axes[1, 0]
    plot_gas_density_panel(
        ax,
        cold_x,
        cold_z,
        cold_mass,
        title="Cold Gas Distribution (x-z)",
        xlabel="x [kpc]",
        ylabel="z [kpc]",
    )

    ax = axes[0, 1]
    plot_stellar_density_panel(
        ax,
        star_x,
        star_y,
        star_mass,
        title="Stellar Distribution (x-y)",
        xlabel="x [kpc]",
        ylabel="y [kpc]",
    )

    ax = axes[0, 2]
    ax.plot(summary["r_cent"], summary["sigma_h_cold"], lw=2.0, c="cyan", label="Gas")
    ax.plot(summary["r_cent"], summary["sigma_h_star"], lw=2.0, c="black", label="Stars")
    if wlm_r is not None and wlm_f is not None:
        ax.plot(wlm_r, wlm_f, label="WLM (Neel et al.)", color="orange", lw=2, ls="--")
    ax.set_yscale("log")
    ax.set_xlabel(r"$R\ [{\rm kpc}]$")
    ax.set_ylabel(r"$\Sigma(R)\ [10^{20}\ {\rm H\ atoms\ cm^{-2}}]$")
    ax.set_xlim(0, 10)
    ax.set_ylim(0.01, 55)
    ax.set_title("Radial Surface Density Profile")
    ax.grid(True, ls="--", alpha=0.4)
    ax.legend()
    ax.tick_params(axis="both", which="both", direction="in")

    ax = axes[1, 1]
    plot_stellar_density_panel(
        ax,
        star_x,
        star_z,
        star_mass,
        title="Stellar Distribution (x-z)",
        xlabel="x [kpc]",
        ylabel="z [kpc]",
    )
    ax.text(
        0.05, 0.95,
        (
            f"Stellar Mass (<5 kpc): {summary['stellar_mass_5kpc']:.2e} M⊙\n"
            f"Stellar Half-light Radius: {summary['r_half']:.2f} kpc\n"
            f"Gas Mass: {summary['m_gas_stellar_region']:.2e} M⊙\n"
            f"Gas fraction: {summary['gas_fraction_stellar_region']:.2f}"
        ),
        transform=ax.transAxes,
        verticalalignment="top",
        color="white",
        bbox=dict(facecolor="black", alpha=0.65, edgecolor="white", linewidth=0.4),
        fontsize=8,
    )
    ax.text(
        0.05,
        0.05,
        (
            f"σx (yz-proj) = {summary['sigma_x']:.1f} km/s\n"
            f"σy (xz-proj) = {summary['sigma_y']:.1f} km/s\n"
            f"σz (xy-proj) = {summary['sigma_z']:.1f} km/s"
        ),
        transform=ax.transAxes,
        verticalalignment="bottom",
        color="white",
        bbox=dict(facecolor="black", alpha=0.65, edgecolor="white", linewidth=0.4),
        fontsize=8,
    )

    ax = axes[1, 2]
    if density_history is not None:
        gas_dens_plot = density_history["gas_density"]
        star_dens_plot = density_history["stellar_density"]
        time_plot = density_history["time"]
        pos_vals = np.concatenate([gas_dens_plot[gas_dens_plot > 0], star_dens_plot[star_dens_plot > 0]])
        if pos_vals.size > 0:
            eps = pos_vals.min() * 1e-3 if pos_vals.min() > 0 else 1e-12
            ax.plot(time_plot, np.where(gas_dens_plot <= 0, eps, gas_dens_plot), lw=2.0, color="cyan", label="Gas")
            ax.plot(time_plot, np.where(star_dens_plot <= 0, eps, star_dens_plot), lw=2.0, color="black", label="Stars")
            ax.set_yscale("log")
        else:
            ax.plot(time_plot, gas_dens_plot, lw=2.0, color="cyan", label="Gas")
            ax.plot(time_plot, star_dens_plot, lw=2.0, color="black", label="Stars")
        ax.legend(fontsize=8, loc="lower right")

    ax.set_xlabel("Time [Gyr]")
    ax.set_ylim(0.01, 55)
    ax.set_xlim(0, 8)
    ax.set_ylabel(r"$\Sigma(R)\ [10^{20}\ {\rm H\ atoms\ cm^{-2}}]$")
    ax.grid(True, ls="--", alpha=0.4)
    ax.tick_params(axis="both", which="both", direction="in")

    time_title = f"T = {tsnap:.2f} Gyr"
    suptitle = time_title
    if title_text is not None:
        title_lines = str(title_text).splitlines() or [str(title_text)]
        title_lines[0] = f"{title_lines[0]}    {time_title}"
        suptitle = "\n".join(title_lines)
    fig.suptitle(suptitle, y=0.99)
    fig.subplots_adjust(left=0.055, right=0.985, bottom=0.065, top=0.92)
    return fig


def run_single_analysis(
    model_dir: Path,
    snapshot: int,
    output_dir: Path,
    radius_max: float,
    nbin: int,
    sigma_radius: float,
    run_metadata=None,
):
    output_dir.mkdir(parents=True, exist_ok=True)

    dwarf_name, snapshot_dir, tsnap, dwarf_df = load_centered_dwarf_dataframe(model_dir, snapshot)
    fig = None
    try:
        summary = compute_summary(dwarf_df, radius_max=radius_max, nbin=nbin, sigma_radius=sigma_radius)
        write_density_snapshot(output_dir, snapshot, tsnap, summary)

        wlm_r, wlm_f = load_wlm_reference(snapshot_dir)
        model_name = detect_model_name(model_dir)
        fd_value, sfr_eff = run_metadata if run_metadata is not None else load_run_metadata()
        title_text = None
        if fd_value is not None and sfr_eff is not None:
            title_text = f"IC_{model_name}\nfeedback = {fd_value / 8.3:.2f} median, sfr_eff = {sfr_eff:.3f}"
        elif model_name:
            title_text = f"IC_{model_name}"

        fig = make_figure(
            dwarf_df,
            tsnap,
            summary,
            density_history=load_density_history(output_dir, max_snapshot=snapshot, max_time=tsnap),
            wlm_r=wlm_r,
            wlm_f=wlm_f,
            title_text=title_text,
        )

        stem = f"snapshot_{snapshot:03d}_density_profile"
        fig_path = output_dir / f"{stem}.png"
        fig.savefig(fig_path, dpi=240)
    finally:
        if fig is not None:
            plt.close(fig)
        del dwarf_df
        gc.collect()

    return {
        "snapshot": snapshot,
        "figure_path": fig_path,
        "dwarf_name": dwarf_name,
    }


def run_single_analysis_task(task):
    model_dir, snapshot, output_dir, radius_max, nbin, sigma_radius, run_metadata = task
    result = run_single_analysis(
        model_dir=Path(model_dir),
        snapshot=snapshot,
        output_dir=Path(output_dir),
        radius_max=radius_max,
        nbin=nbin,
        sigma_radius=sigma_radius,
        run_metadata=run_metadata,
    )
    return {
        "snapshot": result["snapshot"],
        "figure_path": str(result["figure_path"]),
        "dwarf_name": result["dwarf_name"],
    }


def select_pool_start_method():
    available = set(get_all_start_methods())
    if "fork" in available:
        return "fork"
    if "spawn" in available:
        return "spawn"
    return next(iter(available))


def main():
    args = parse_args()
    if args.no_video and args.video:
        raise SystemExit("--video and --no-video cannot be used together.")
    if args.no_video and args.video_only:
        raise SystemExit("--video-only and --no-video cannot be used together.")
    if args.snapshot is not None and args.video_only:
        raise SystemExit("--video-only cannot be used with --snapshot.")

    explicit_mode = args.snapshot is not None
    output_dir = resolve_output_dir(args.output_dir, explicit_mode=explicit_mode)

    if args.video_only:
        model_dir = resolve_legacy_model_dir(args.model_dir)
        try:
            snapshots, selection_mode = detect_batch_snapshots(model_dir, args.range)
        except FileNotFoundError:
            if args.range is not None:
                raise
            snapshots = detect_frame_snapshots(output_dir)
            selection_mode = "frames"

        if not snapshots:
            raise SystemExit(f"--video-only found no snapshots or density-profile PNG frames under {output_dir}")

        selection_text = (
            f"{snapshots[0]}-{snapshots[-1]}"
            if selection_mode != "list"
            else ",".join(str(num) for num in snapshots)
        )
        print(f"Video-only mode: model_dir={model_dir}, snapshots={selection_text}, count={len(snapshots)}")
        video_path = assemble_density_video_with_log(
            output_dir=output_dir,
            model_dir=model_dir,
            snapshots=snapshots,
            scale=args.video_scale,
            crf=args.video_crf,
        )
        check_video_batch_result(
            output_dir=output_dir,
            model_dir=model_dir,
            snapshots=snapshots,
            make_video=True,
            selection_mode=selection_mode,
            no_video=False,
            forced_video=True,
            video_path=video_path,
        )
        if video_path is None:
            raise SystemExit("Video assembly did not produce an MP4.")
        return

    if explicit_mode:
        if args.model_dir is None:
            raise SystemExit("--snapshot requires --model-dir in explicit single-snapshot mode.")

        model_dir = Path(args.model_dir).resolve()
        run_metadata = load_run_metadata()
        if not args.replace and snapshot_output_exists(model_dir, args.snapshot, output_dir):
            fig_path = get_output_paths(model_dir, args.snapshot, output_dir)
            print(f"Skipped snapshot {args.snapshot:03d}: existing outputs found")
            print(f"Existing figure: {fig_path}")
            return

        result = run_single_analysis(
            model_dir=model_dir,
            snapshot=args.snapshot,
            output_dir=output_dir,
            radius_max=args.radius_max,
            nbin=args.nbin,
            sigma_radius=args.sigma_radius,
            run_metadata=run_metadata,
        )
        fd_value, sfr_eff = run_metadata
        create_density_evolution_csv(output_dir, fd_value=fd_value, sfr_eff=sfr_eff)
        cleanup_density_snapshots(output_dir)
        print(f"Saved figure: {result['figure_path']}")
        return

    model_dir = resolve_legacy_model_dir(args.model_dir)
    run_metadata = load_run_metadata()
    snapshots, selection_mode = detect_batch_snapshots(model_dir, args.range)
    make_video = (not args.no_video) and (selection_mode != "list" or args.video)
    selection_text = f"{snapshots[0]}-{snapshots[-1]}" if selection_mode != "list" else ",".join(str(num) for num in snapshots)
    print(f"Legacy batch mode: model_dir={model_dir}, snapshots={selection_text}, count={len(snapshots)}")
    if selection_mode == "list" and not make_video:
        print("Video assembly disabled by default for comma-separated --range selections.")

    if args.replace:
        pending_snapshots = snapshots
        skipped_snapshots = []
    else:
        pending_snapshots = []
        skipped_snapshots = []
        for snapshot in snapshots:
            if snapshot_output_exists(model_dir, snapshot, output_dir):
                skipped_snapshots.append(snapshot)
            else:
                pending_snapshots.append(snapshot)

    print(
        f"Resume status: total={len(snapshots)}, pending={len(pending_snapshots)}, "
        f"skipped_existing={len(skipped_snapshots)}, replace={args.replace}"
    )

    if not pending_snapshots:
        print("No pending snapshots to process.")
        video_path = None
        if make_video:
            video_path = assemble_density_video_with_log(
                output_dir=output_dir,
                model_dir=model_dir,
                snapshots=snapshots,
                scale=args.video_scale,
                crf=args.video_crf,
            )
        check_video_batch_result(
            output_dir=output_dir,
            model_dir=model_dir,
            snapshots=snapshots,
            make_video=make_video,
            selection_mode=selection_mode,
            no_video=args.no_video,
            forced_video=args.video,
            video_path=video_path,
        )
        return

    tasks = [
        (str(model_dir), snapshot, str(output_dir), args.radius_max, args.nbin, args.sigma_radius, run_metadata)
        for snapshot in pending_snapshots
    ]

    results = []
    progress = ProgressReporter(len(tasks), desc="CheckDensity", unit="snap")
    try:
        if args.processes <= 1:
            for task in tasks:
                results.append(run_single_analysis_task(task))
                progress.update(1)
        else:
            start_method = select_pool_start_method()
            print(
                f"Parallel batch mode requested with {args.processes} workers; "
                f"memory use scales roughly with worker count. start_method={start_method}"
            )
            with get_context(start_method).Pool(processes=args.processes, maxtasksperchild=1) as pool:
                for result in pool.imap_unordered(run_single_analysis_task, tasks):
                    results.append(result)
                    progress.update(1)
    finally:
        progress.close()

    fd_value, sfr_eff = run_metadata
    create_density_evolution_csv(output_dir, fd_value=fd_value, sfr_eff=sfr_eff)
    cleanup_density_snapshots(output_dir)

    video_path = None
    if make_video:
        video_path = assemble_density_video_with_log(
            output_dir=output_dir,
            model_dir=model_dir,
            snapshots=snapshots,
            scale=args.video_scale,
            crf=args.video_crf,
            results=results,
        )
    check_video_batch_result(
        output_dir=output_dir,
        model_dir=model_dir,
        snapshots=snapshots,
        make_video=make_video,
        selection_mode=selection_mode,
        no_video=args.no_video,
        forced_video=args.video,
        video_path=video_path,
        results=results,
    )

    print(f"Completed {len(results)} snapshot(s). Output directory: {output_dir}")
    for result in sorted(results, key=lambda item: item["snapshot"]):
        print(f"[{result['snapshot']:03d}] {result['figure_path']}")


if __name__ == "__main__":
    main()
