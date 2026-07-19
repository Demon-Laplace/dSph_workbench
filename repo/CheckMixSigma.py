import argparse
import importlib
import importlib.util
import os
import re
import sys
import warnings
from pathlib import Path
from typing import Optional

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig-dsph")
warnings.filterwarnings(
    "ignore",
    message=r"mergecube\.par not found under .*; falling back to legacy cube_path=.*",
    category=RuntimeWarning,
)

import matplotlib.pyplot as plt
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
WORKSPACE_ROOT = SCRIPT_DIR.parent
MODEL_LOCAL_MODULES = [
    "PlotFig",
    "snapshot_context",
    "variable",
    "basefunc",
    "analysis_core",
    "data_processing",
    "mysubs",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compare projected r-sigma from PlotFig panel (1,0) with a true 3D r-sigma "
            "profile for one model snapshot."
        )
    )
    parser.add_argument(
        "model_snapshot",
        help="Model and snapshot spec, e.g. 1035-277 or /path/to/Fornax1035-277.",
    )
    parser.add_argument(
        "--dwarf-name",
        default="Fornax",
        help="Dwarf prefix used when the model token is numeric. Default: Fornax",
    )
    parser.add_argument(
        "--model-root",
        help=(
            "Optional base directory for model lookup. "
            "Useful when models are not under the default fixture/current-workdir search roots."
        ),
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=50,
        help="Number of radial-bin edges for the 3D sigma profile. Default: 50",
    )
    parser.add_argument(
        "--output-dir",
        help="Directory for the output figure. Defaults to ../sandbox_runs/check_mix_sigma",
    )
    parser.add_argument(
        "--output-name",
        help="Optional custom filename for the generated figure.",
    )
    return parser.parse_args()


def parse_model_snapshot_spec(raw_spec: str):
    path_spec = Path(raw_spec)
    if path_spec.exists():
        model_dir = path_spec if path_spec.is_dir() else path_spec.parent
        match = re.match(r"^(?P<model>.+)-(?P<snapshot>\d+)$", model_dir.name)
        if match:
            return {
                "model_spec": match.group("model"),
                "snapshot": int(match.group("snapshot")),
                "raw_spec": raw_spec,
                "path_override": model_dir.resolve(),
            }
        raise ValueError(
            f"Could not infer snapshot from path '{raw_spec}'. Expected directory name like Fornax1035-277."
        )

    match = re.match(r"^(?P<model>.+)-(?P<snapshot>\d+)$", raw_spec)
    if not match:
        raise ValueError(
            f"Invalid model-snapshot spec '{raw_spec}'. Expected format like 1035-277."
        )
    return {
        "model_spec": match.group("model"),
        "snapshot": int(match.group("snapshot")),
        "raw_spec": raw_spec,
        "path_override": None,
    }


def build_model_label(model_spec: str, dwarf_name: str):
    if Path(model_spec).exists():
        return Path(model_spec).resolve().name
    if model_spec.isdigit():
        return f"{dwarf_name}{int(model_spec)}"
    return model_spec


def candidate_search_roots(dwarf_name: str, model_root: Optional[str]):
    roots = []
    if model_root:
        root = Path(model_root).resolve()
        roots.extend([root, root / dwarf_name])
    else:
        cwd = Path.cwd().resolve()
        roots.extend(
            [
                cwd,
                cwd / dwarf_name,
                cwd / "IsolationModel" / dwarf_name,
                WORKSPACE_ROOT / "fixtures_realistic" / "local_pc" / "IsolationModel" / dwarf_name,
                WORKSPACE_ROOT / "fixtures_realistic" / "local_pc" / dwarf_name,
                WORKSPACE_ROOT / "fixtures_realistic" / "mock_server",
            ]
        )
    return [root for root in roots if root.exists()]


def resolve_model_dir(
    model_spec: str,
    dwarf_name: str,
    model_root: Optional[str],
    path_override: Optional[Path],
):
    if path_override is not None:
        if not path_override.is_dir():
            raise NotADirectoryError(f"Model path is not a directory: {path_override}")
        return path_override.resolve()

    model_path = Path(model_spec)
    if model_path.exists():
        if not model_path.is_dir():
            raise NotADirectoryError(f"Model path is not a directory: {model_path}")
        return model_path.resolve()

    model_label = build_model_label(model_spec, dwarf_name)
    candidates = []
    for root in candidate_search_roots(dwarf_name, model_root):
        candidate = root / model_label
        if candidate.exists():
            if not candidate.is_dir():
                raise NotADirectoryError(f"Resolved model path is not a directory: {candidate}")
            candidates.append(candidate.resolve())

    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        return candidates

    tried_text = "\n".join(f"  - {root / model_label}" for root in candidate_search_roots(dwarf_name, model_root))
    if not tried_text:
        tried_text = "  - <no search roots>"
    raise FileNotFoundError(
        f'Could not resolve model "{model_spec}" (label={model_label}). Tried:\n{tried_text}'
    )

def select_model_dir_for_snapshot(model_dir_or_dirs, snapshot: int):
    if isinstance(model_dir_or_dirs, Path):
        ensure_snapshot_exists(model_dir_or_dirs, snapshot)
        return model_dir_or_dirs

    matched = []
    for candidate in model_dir_or_dirs:
        if snapshot in detect_available_snapshots(candidate):
            matched.append(candidate)

    if len(matched) == 1:
        return matched[0]
    if len(matched) > 1:
        return matched[0]

    available_text = "\n".join(
        f"  - {candidate}: {detect_available_snapshots(candidate)}" for candidate in model_dir_or_dirs
    )
    raise FileNotFoundError(
        f"Snapshot {snapshot:03d} not found in any resolved model directory. Candidates:\n{available_text}"
    )


def list_snapshot_numbers(folder: Path):
    numbers = []
    for file_path in folder.glob("snapshot_*.hdf5"):
        match = re.match(r"snapshot_(\d+)\.hdf5", file_path.name)
        if match:
            numbers.append(int(match.group(1)))
    return sorted(numbers)


def detect_available_snapshots(model_dir: Path):
    snapshot_numbers = set(list_snapshot_numbers(model_dir))
    for subdir in sorted(path for path in model_dir.iterdir() if path.is_dir()):
        snapshot_numbers.update(list_snapshot_numbers(subdir))
    return sorted(snapshot_numbers)


def ensure_snapshot_exists(model_dir: Path, snapshot: int):
    snapshots = detect_available_snapshots(model_dir)
    if snapshot not in snapshots:
        raise FileNotFoundError(
            f"Snapshot {snapshot:03d} not found under {model_dir}. "
            f"Available snapshots: {snapshots}"
        )


def make_output_dir(output_dir_arg: Optional[str]):
    if output_dir_arg:
        return Path(output_dir_arg).resolve()
    sandbox_root = WORKSPACE_ROOT / "sandbox_runs"
    if sandbox_root.exists():
        return (sandbox_root / "check_mix_sigma").resolve()
    return (Path.cwd().resolve() / "check_mix_sigma_output").resolve()


def make_output_name(model_label: str, snapshot: int, custom_name: Optional[str]):
    if custom_name:
        return custom_name
    return f"CheckMixSigma_{model_label}_{snapshot:03d}.png"


def compute_3d_sigma_profile(snapshot, bins: int):
    df = snapshot["df"]
    star_mask = np.asarray(snapshot["total_dw_star_mask"], dtype=bool)
    if not np.any(star_mask):
        return np.array([], dtype=float), np.array([], dtype=float)

    star_df = df.loc[star_mask]
    return compute_3d_sigma_profile_from_arrays(
        star_df["x"].to_numpy(dtype=float) - float(snapshot["dw_xc"]),
        star_df["y"].to_numpy(dtype=float) - float(snapshot["dw_yc"]),
        star_df["z"].to_numpy(dtype=float) - float(snapshot["dw_zc"]),
        star_df["vx"].to_numpy(dtype=float),
        star_df["vy"].to_numpy(dtype=float),
        star_df["vz"].to_numpy(dtype=float),
        bins=bins,
    )


def rotate_xy(x, y, pa):
    cos_pa = np.cos(pa)
    sin_pa = np.sin(pa)
    x_rot = cos_pa * x + sin_pa * y
    y_rot = -sin_pa * x + cos_pa * y
    return x_rot, y_rot


def compute_3d_sigma_profile_from_arrays(x, y, z, vx, vy, vz, bins: int):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    vx = np.asarray(vx, dtype=float)
    vy = np.asarray(vy, dtype=float)
    vz = np.asarray(vz, dtype=float)

    r_3d = np.sqrt(x**2 + y**2 + z**2)
    finite = np.isfinite(r_3d) & np.isfinite(vx) & np.isfinite(vy) & np.isfinite(vz)
    if not np.any(finite):
        return np.array([], dtype=float), np.array([], dtype=float)

    r_3d = r_3d[finite]
    vx = vx[finite]
    vy = vy[finite]
    vz = vz[finite]

    positive = r_3d > 0
    if not np.any(positive):
        return np.array([], dtype=float), np.array([], dtype=float)

    r_max = r_3d[positive].max()
    if not np.isfinite(r_max) or r_max <= 0:
        return np.array([], dtype=float), np.array([], dtype=float)

    bin_count = max(int(bins), 2)
    r_bins = np.linspace(0.0, r_max, bin_count)
    bin_centers = 0.5 * (r_bins[:-1] + r_bins[1:])

    sigma_xyz = np.full(bin_centers.shape, np.nan, dtype=float)
    for idx in range(len(bin_centers)):
        if idx == len(bin_centers) - 1:
            in_bin = (r_3d >= r_bins[idx]) & (r_3d <= r_bins[idx + 1])
        else:
            in_bin = (r_3d >= r_bins[idx]) & (r_3d < r_bins[idx + 1])
        if np.count_nonzero(in_bin) < 2:
            continue
        sigma_x = np.std(vx[in_bin])
        sigma_y = np.std(vy[in_bin])
        sigma_z = np.std(vz[in_bin])
        sigma_xyz[idx] = np.sqrt((sigma_x**2 + sigma_y**2 + sigma_z**2) / 3.0)

    valid = np.isfinite(sigma_xyz)
    return bin_centers[valid], sigma_xyz[valid]


def compute_3d_ellipsoidal_sigma_profile(snapshot, eps: float, pa: float, bins: int):
    df = snapshot["df"]
    star_mask = np.asarray(snapshot["total_dw_star_mask"], dtype=bool)
    if not np.any(star_mask):
        return np.array([], dtype=float), np.array([], dtype=float)

    star_df = df.loc[star_mask]
    x = star_df["x"].to_numpy(dtype=float) - float(snapshot["dw_xc"])
    y = star_df["y"].to_numpy(dtype=float) - float(snapshot["dw_yc"])
    z = star_df["z"].to_numpy(dtype=float) - float(snapshot["dw_zc"])
    vx = star_df["vx"].to_numpy(dtype=float)
    vy = star_df["vy"].to_numpy(dtype=float)
    vz = star_df["vz"].to_numpy(dtype=float)

    x_rot, y_rot = rotate_xy(x, y, pa)
    q = 1.0 - float(eps)
    if not np.isfinite(q) or q <= 0:
        q = np.finfo(float).eps

    r_ell = np.sqrt(x_rot**2 + (y_rot / q) ** 2 + (z / q) ** 2)
    return compute_sigma_profile_from_radius(r_ell, vx, vy, vz, bins=bins)


def compute_sigma_profile_from_radius(radius, vx, vy, vz, bins: int):
    radius = np.asarray(radius, dtype=float)
    vx = np.asarray(vx, dtype=float)
    vy = np.asarray(vy, dtype=float)
    vz = np.asarray(vz, dtype=float)

    finite = np.isfinite(radius) & np.isfinite(vx) & np.isfinite(vy) & np.isfinite(vz)
    if not np.any(finite):
        return np.array([], dtype=float), np.array([], dtype=float)

    radius = radius[finite]
    vx = vx[finite]
    vy = vy[finite]
    vz = vz[finite]

    positive = radius > 0
    if not np.any(positive):
        return np.array([], dtype=float), np.array([], dtype=float)

    r_max = radius[positive].max()
    if not np.isfinite(r_max) or r_max <= 0:
        return np.array([], dtype=float), np.array([], dtype=float)

    bin_count = max(int(bins), 2)
    r_bins = np.linspace(0.0, r_max, bin_count)
    bin_centers = 0.5 * (r_bins[:-1] + r_bins[1:])

    sigma_xyz = np.full(bin_centers.shape, np.nan, dtype=float)
    for idx in range(len(bin_centers)):
        if idx == len(bin_centers) - 1:
            in_bin = (radius >= r_bins[idx]) & (radius <= r_bins[idx + 1])
        else:
            in_bin = (radius >= r_bins[idx]) & (radius < r_bins[idx + 1])
        if np.count_nonzero(in_bin) < 2:
            continue
        sigma_x = np.std(vx[in_bin])
        sigma_y = np.std(vy[in_bin])
        sigma_z = np.std(vz[in_bin])
        sigma_xyz[idx] = np.sqrt((sigma_x**2 + sigma_y**2 + sigma_z**2) / 3.0)

    valid = np.isfinite(sigma_xyz)
    return bin_centers[valid], sigma_xyz[valid]


def remove_velocity_dispersion_markers(ax):
    red_rgba = np.array([1.0, 0.0, 0.0, 1.0])
    gray_rgba = np.array([0.50196078, 0.50196078, 0.50196078, 1.0])

    for collection in list(ax.collections):
        if hasattr(collection, "get_offsets"):
            offsets = collection.get_offsets()
            facecolors = collection.get_facecolors()
            if len(offsets) == 1 and len(facecolors) > 0 and np.allclose(facecolors[0], red_rgba):
                collection.remove()
                continue

        if hasattr(collection, "get_segments"):
            segments = collection.get_segments()
            colors = collection.get_colors()
            if len(segments) == 1 and len(colors) > 0 and np.allclose(colors[0], gray_rgba):
                collection.remove()


def unload_model_local_modules():
    for module_name in MODEL_LOCAL_MODULES:
        sys.modules.pop(module_name, None)


def load_model_local_plot_modules(model_dir: Path):
    plotfig_path = model_dir / "PlotFig.py"
    if not plotfig_path.exists():
        raise FileNotFoundError(f"PlotFig.py not found under model directory: {plotfig_path}")

    unload_model_local_modules()
    sys.path.insert(0, str(model_dir))
    try:
        spec = importlib.util.spec_from_file_location("PlotFig", plotfig_path)
        if spec is None or spec.loader is None:
            raise ImportError(f"Could not create import spec for {plotfig_path}")

        plotfig_module = importlib.util.module_from_spec(spec)
        sys.modules["PlotFig"] = plotfig_module
        spec.loader.exec_module(plotfig_module)
        snapshot_context_module = importlib.import_module("snapshot_context")
        return plotfig_module, snapshot_context_module
    finally:
        try:
            sys.path.remove(str(model_dir))
        except ValueError:
            pass


def build_sigma_payload(PlotFig, snapshot_context, snapshot_num: int, bins: int):
    context = PlotFig.get_plot_context()
    snapshot = snapshot_context.prepare_snapshot_context(
        folder_path=PlotFig.folder_path,
        snapshot_num=snapshot_num,
        core_radius=context["core_radius"],
    )
    metrics = PlotFig._compute_plot_snapshot_metrics(snapshot, context, snapshot_num)
    row = context["elinfo_by_numsp"].loc[snapshot_num]
    sigma_3d_r, sigma_3d = compute_3d_sigma_profile(snapshot, bins=bins)
    sigma_3d_ell_r, sigma_3d_ell = compute_3d_ellipsoidal_sigma_profile(
        snapshot,
        eps=metrics["eps"],
        pa=metrics["pa"],
        bins=bins,
    )
    return {
        "mode": "modern",
        "modelname": context["modelname"],
        "tsnap": snapshot["tsnap"],
        "bin_centers": metrics["bin_centers"],
        "vlos_dispersion": metrics["vlos_dispersion"],
        "tsigma": row["tsigma"],
        "r_half": metrics["r_half_circularized"],
        "sigma_3d_r": sigma_3d_r,
        "sigma_3d": sigma_3d,
        "sigma_3d_ell_r": sigma_3d_ell_r,
        "sigma_3d_ell": sigma_3d_ell,
    }

def draw_projected_sigma_panel(ax, PlotFig, payload):
    PlotFig._plot_velocity_dispersion_panel(
        ax,
        payload["bin_centers"],
        payload["vlos_dispersion"],
        payload["tsigma"],
        payload["r_half"],
    )


def render_mix_sigma_figure(model_dir: Path, snapshot_num: int, bins: int, output_path: Path):
    original_cwd = Path.cwd()
    try:
        os.chdir(model_dir)

        PlotFig, snapshot_context = load_model_local_plot_modules(model_dir)
        payload = build_sigma_payload(PlotFig, snapshot_context, snapshot_num, bins)

        fig, ax = plt.subplots(1, 1, figsize=(7.0, 5.4), dpi=150)
        draw_projected_sigma_panel(ax, PlotFig, payload)
        remove_velocity_dispersion_markers(ax)
        if payload["sigma_3d_r"].size > 0:
            ax.plot(
                payload["sigma_3d_r"],
                payload["sigma_3d"],
                color="tab:orange",
                lw=2.0,
                label=r"3D spherical $\sigma_{xyz}(r)$",
                zorder=4,
            )
        if payload["sigma_3d_ell_r"].size > 0:
            ax.plot(
                payload["sigma_3d_ell_r"],
                payload["sigma_3d_ell"],
                color="tab:green",
                lw=2.0,
                ls="--",
                label=r"3D ellipsoidal $\sigma_{xyz}(r_{\rm ell})$",
                zorder=4,
            )
        ax.set_title("Projected vs 3D r-sigma", fontsize=13)
        ax.legend(fontsize=11, loc="upper right")
        fig.suptitle(
            f"{payload['modelname']} snapshot={snapshot_num:03d}  T={payload['tsnap']:.2f} Gyr",
            fontsize=13,
        )
        fig.subplots_adjust(left=0.11, right=0.98, bottom=0.12, top=0.88)
        fig.savefig(output_path, bbox_inches="tight", pad_inches=0.12)
        plt.close(fig)

        return payload["modelname"], payload["tsnap"]
    finally:
        os.chdir(original_cwd)


def main():
    args = parse_args()
    request = parse_model_snapshot_spec(args.model_snapshot)
    model_dir_candidates = resolve_model_dir(
        request["model_spec"],
        args.dwarf_name,
        args.model_root,
        request["path_override"],
    )
    model_dir = select_model_dir_for_snapshot(model_dir_candidates, request["snapshot"])

    output_dir = make_output_dir(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    model_label = build_model_label(request["model_spec"], args.dwarf_name)
    output_path = output_dir / make_output_name(model_label, request["snapshot"], args.output_name)

    modelname, tsnap = render_mix_sigma_figure(
        model_dir=model_dir,
        snapshot_num=request["snapshot"],
        bins=args.bins,
        output_path=output_path,
    )
    print(f"Saved figure: {output_path}")
    print(f"Model: {modelname}")
    print(f"Snapshot: {request['snapshot']:03d}")
    print(f"Time: {tsnap:.3f} Gyr")


if __name__ == "__main__":
    main()
