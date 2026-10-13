import argparse
import contextlib
import io
import os
import sys
import time
from multiprocessing import get_context
from pathlib import Path
from model_paths import model_id
from typing import Optional

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig-dsph")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm = None

REPO_ROOT = Path(__file__).resolve().parent
WORKSPACE_ROOT = Path.cwd()
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from basefunc import Analysis, DataProcessor, GalaxySimulation  # noqa: E402


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
        type=str,
        help='Snapshot range for legacy batch mode, e.g. "0,400". Defaults to all detected snapshots.',
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
    return parser.parse_args()


def parse_snapshot_range(range_text: str):
    try:
        start, end = map(int, range_text.split(","))
    except ValueError as exc:
        raise ValueError(f'Invalid --range value "{range_text}". Expected format "start,end".') from exc
    if start > end:
        raise ValueError(f"Invalid --range value \"{range_text}\": start must be <= end.")
    return start, end


def detect_batch_snapshots(model_dir: Path, range_text: Optional[str]):
    snapshot_numbers = DataProcessor.list_snapshot_numbers(model_dir)
    if not snapshot_numbers:
        raise FileNotFoundError(f"No snapshots found under {model_dir}")

    if range_text is None:
        return snapshot_numbers

    start, end = parse_snapshot_range(range_text)
    requested = [num for num in snapshot_numbers if start <= num <= end]
    if not requested:
        raise FileNotFoundError(f"No snapshots in requested range {start},{end} under {model_dir}")
    return requested


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


def snapshot_output_exists(model_dir: Path, snapshot: int, output_dir: Path):
    fig_path = get_output_paths(model_dir, snapshot, output_dir)
    return fig_path.exists()


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
    dwarf_name, ini_file = find_ini_file(
        [
            snapshot_dir,
            snapshot_dir.parent,
            model_dir,
            model_dir.parent,
        ]
    )
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
            gas_fields=("temp", "nh"),
            star_fields=(),
            particle_types=(0, 2, 3, 4),
        )

    dwarf_df = df.copy()
    tp = dwarf_df["tp"].to_numpy()
    dwarf_star_mask = (tp == 2) | (tp == 4)
    dwarf_gas_mask = tp == 0
    dwarf_total_mask = dwarf_star_mask | dwarf_gas_mask

    xc, yc, zc = mass_weighted_center(
        dwarf_df.loc[dwarf_star_mask, "x"].to_numpy(),
        dwarf_df.loc[dwarf_star_mask, "y"].to_numpy(),
        dwarf_df.loc[dwarf_star_mask, "z"].to_numpy(),
        dwarf_df.loc[dwarf_star_mask, "m"].to_numpy(),
    )

    x_centered = dwarf_df["x"].to_numpy(dtype=np.float64, copy=True) - xc
    y_centered = dwarf_df["y"].to_numpy(dtype=np.float64, copy=True) - yc
    z_centered = dwarf_df["z"].to_numpy(dtype=np.float64, copy=True) - zc
    dwarf_df.loc[:, "x"] = x_centered.astype(dwarf_df["x"].dtype, copy=False)
    dwarf_df.loc[:, "y"] = y_centered.astype(dwarf_df["y"].dtype, copy=False)
    dwarf_df.loc[:, "z"] = z_centered.astype(dwarf_df["z"].dtype, copy=False)

    for velocity_axis in ("vx", "vy", "vz"):
        v_centered = dwarf_df[velocity_axis].to_numpy(dtype=np.float64, copy=True)
        v_centered -= dwarf_df.loc[dwarf_total_mask, velocity_axis].to_numpy(dtype=np.float64).mean()
        dwarf_df.loc[:, velocity_axis] = v_centered.astype(dwarf_df[velocity_axis].dtype, copy=False)

    dwarf_df.loc[:, "r3d"] = np.sqrt(
        dwarf_df["x"].to_numpy() ** 2 + dwarf_df["y"].to_numpy() ** 2 + dwarf_df["z"].to_numpy() ** 2
    )

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
    cold_gas_mask = (tp == 0) & (temp < 20000.0)
    hot_gas_mask = (tp == 0) & ~cold_gas_mask

    star_within = star_mask & (r3d <= radius_max)
    cold_within = cold_gas_mask & (r3d <= radius_max)
    hot_within = hot_gas_mask & (r3d <= radius_max)

    xs = x[star_within]
    ys = y[star_within]
    zs = z[star_within]
    ms = m[star_within]

    xcg = x[cold_within]
    ycg = y[cold_within]
    zcg = z[cold_within]
    mcg = m[cold_within]

    xhg = x[hot_within]
    yhg = y[hot_within]
    zhg = z[hot_within]
    mhg = m[hot_within]

    r_star = np.sqrt(xs**2 + ys**2)
    r_cold = np.sqrt(xcg**2 + ycg**2)

    bins = np.linspace(0.0, radius_max, nbin + 1)
    r_cent = 0.5 * (bins[:-1] + bins[1:])
    sigma_star, sigma_h_star = radial_surface_density(r_star, ms, bins)
    sigma_cold, sigma_h_cold = radial_surface_density(r_cold, mcg, bins)

    ellipticity, pa = Analysis.calculate_ellipticity(xs, zs, mass=ms, n_neighbors=30)
    rotate = np.array([[np.cos(-pa), -np.sin(-pa)], [np.sin(-pa), np.cos(-pa)]])
    coords = np.vstack((xs, zs))
    x_rot, z_rot = (rotate @ coords)
    r_half = Analysis.half_light_radius(x_rot, z_rot, ms, ep=ellipticity)

    subset_star = dwarf_df.loc[star_mask & (r3d <= 5.0)].copy()
    subset_cold = dwarf_df.loc[cold_gas_mask & (r3d <= 5.0)].copy()
    subset_hot = dwarf_df.loc[hot_gas_mask & (r3d <= 5.0)].copy()

    boundary_path = safe_boundary_path(xs, ys)
    if boundary_path is not None and not subset_star.empty:
        star_in = boundary_path.contains_points(np.column_stack([subset_star["x"], subset_star["y"]]))
    else:
        star_in = np.ones(len(subset_star), dtype=bool)

    m_star_boundary = float(subset_star.loc[star_in, "m"].sum())
    m_cold_gas = float(subset_cold["m"].sum())
    m_hot_gas = float(subset_hot["m"].sum())

    r_proj_star = np.sqrt(subset_star["x"].to_numpy() ** 2 + subset_star["z"].to_numpy() ** 2)
    r_proj_cold = np.sqrt(subset_cold["x"].to_numpy() ** 2 + subset_cold["z"].to_numpy() ** 2)
    r_proj_hot = np.sqrt(subset_hot["x"].to_numpy() ** 2 + subset_hot["z"].to_numpy() ** 2)
    m_star_rhalf = float(subset_star.loc[r_proj_star <= r_half, "m"].sum()) if len(subset_star) else 0.0
    m_cold_rhalf = float(subset_cold.loc[r_proj_cold <= r_half, "m"].sum()) if len(subset_cold) else 0.0
    m_hot_rhalf = float(subset_hot.loc[r_proj_hot <= r_half, "m"].sum()) if len(subset_hot) else 0.0

    sigma_z = safe_std(vz[star_mask & (np.sqrt(x**2 + y**2) < sigma_radius)])
    sigma_y = safe_std(vy[star_mask & (np.sqrt(x**2 + z**2) < sigma_radius)])
    sigma_x = safe_std(vx[star_mask & (np.sqrt(y**2 + z**2) < sigma_radius)])

    stellar_mass_5kpc = float(m[star_mask & (r3d < 5.0)].sum())
    total_gas = float(m[tp == 0].sum())
    n_gas = int(np.count_nonzero(tp == 0))
    n_star = int(np.count_nonzero(tp >= 2))

    gas_fraction_denom = m_cold_gas + m_hot_gas + stellar_mass_5kpc
    gas_fraction = float(m_cold_gas / gas_fraction_denom) if gas_fraction_denom > 0 else 0.0
    gas_fraction_rhalf_denom = m_cold_rhalf + m_hot_rhalf + m_star_rhalf
    gas_fraction_rhalf = float(m_cold_rhalf / gas_fraction_rhalf_denom) if gas_fraction_rhalf_denom > 0 else 0.0

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
    return model_id(model_dir.parent if model_dir.name == "output" else model_dir)


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
        [[tsnap, summary["sigma_h_inner_gas"], summary["sigma_h_inner_star"]]],
        header="time gas_density stellar_density",
        fmt="%.6f %.6e %.6e",
    )


def load_density_history(output_dir: Path):
    records = {}
    evo_path = output_dir / "density_evolution.csv"
    if evo_path.exists():
        try:
            df_evo = pd.read_csv(evo_path, comment="#")
            for _, row in df_evo.iterrows():
                time_key = round(float(row["time"]), 10)
                records[time_key] = (
                    float(row["time"]),
                    float(row["gas_density"]),
                    float(row["stellar_density"]),
                )
        except Exception:
            pass

    for density_file in sorted(output_dir.glob("density_snapshot_*.txt")):
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
                time_key = round(float(row[0]), 10)
                records[time_key] = (float(row[0]), float(row[1]), float(row[2]))
            except (TypeError, ValueError):
                continue

    if not records:
        return None

    all_data = np.array(sorted(records.values(), key=lambda item: item[0]), dtype=float)
    return {
        "time": all_data[:, 0],
        "gas_density": all_data[:, 1],
        "stellar_density": all_data[:, 2],
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
        f.write("time,gas_density,stellar_density\n")
        data_to_save = np.column_stack([history["time"], history["gas_density"], history["stellar_density"]])
        np.savetxt(f, data_to_save, delimiter=",", fmt="%.6f,%.6e,%.6e")


def cleanup_density_snapshots(output_dir: Path):
    for density_file in output_dir.glob("density_snapshot_*.txt"):
        try:
            density_file.unlink()
        except Exception:
            pass


def make_figure(dwarf_df: pd.DataFrame, tsnap: float, summary: dict, density_history=None, wlm_r=None, wlm_f=None, title_text=None):
    tp = dwarf_df["tp"].to_numpy()
    star_mask = (tp == 2) | (tp == 4)
    cold_gas_mask = (tp == 0) & (dwarf_df["temp"].to_numpy() < 20000.0)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=150)

    ax = axes[0, 0]
    ax.scatter(dwarf_df.loc[star_mask, "x"], dwarf_df.loc[star_mask, "y"], s=0.1, c="black", label="Stars")
    ax.scatter(dwarf_df.loc[cold_gas_mask, "x"], dwarf_df.loc[cold_gas_mask, "y"], s=0.2, c="cyan", alpha=0.6, zorder=0, label="Cold gas")
    disp_text = (
        f"r < 1 kpc:\n"
        f"σx (yz-proj) = {summary['sigma_x']:.1f} km/s\n"
        f"σy (xz-proj) = {summary['sigma_y']:.1f} km/s\n"
        f"σz (xy-proj) = {summary['sigma_z']:.1f} km/s"
    )
    ax.text(
        0.05, 0.95, disp_text,
        transform=ax.transAxes,
        verticalalignment="top",
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
        fontsize=8,
    )
    ax.set_xlim(-20, 20)
    ax.set_ylim(-20, 20)
    ax.set_xlabel("x [kpc]")
    ax.set_ylabel("y [kpc]")
    ax.set_title(f"T = {tsnap:.2f} Gyr")
    ax.grid(True, ls="--", alpha=0.4)
    ax.legend(markerscale=5, fontsize=7, loc="lower right")
    ax.tick_params(axis="both", which="both", direction="in")

    ax = axes[0, 1]
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

    ax = axes[1, 0]
    ax.scatter(dwarf_df.loc[star_mask, "x"], dwarf_df.loc[star_mask, "z"], s=0.1, c="black", label="Stars")
    ax.scatter(dwarf_df.loc[cold_gas_mask, "x"], dwarf_df.loc[cold_gas_mask, "z"], s=0.2, c="cyan", alpha=0.6, zorder=0, label="Cold gas")
    ax.text(
        0.05, 0.95,
        (
            f"Stellar Mass (<5kpc):{summary['stellar_mass_5kpc']:.2e} M⊙\n"
            f"Stellar Half-light Radius: {summary['r_half']:.2f} kpc\n"
            f"Stellar Mass within R_half: {summary['m_star_rhalf']:.2e} M⊙\n"
            f"Cold Gas Mass within R_half: {summary['m_cold_rhalf']:.2e} M⊙\n"
            f"Total Gas Mass (tp==0): {summary['total_gas']:.2e} M⊙\n"
            f"Gas Particles (tp==0): {summary['n_gas']}\n"
            f"Star Particles (tp>=2): {summary['n_star']}\n"
            f"Gas fraction: {summary['gas_fraction']:.2f}\n"
            f"Gas fraction within R_half: {summary['gas_fraction_rhalf']:.2f}"
        ),
        transform=ax.transAxes,
        verticalalignment="top",
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
        fontsize=8,
    )
    ax.set_xlim(-20, 20)
    ax.set_ylim(-20, 20)
    ax.set_xlabel("x [kpc]")
    ax.set_ylabel("z [kpc]")
    ax.grid(True, ls="--", alpha=0.4)
    ax.legend(markerscale=5, fontsize=7, loc="lower right")
    ax.tick_params(axis="both", which="both", direction="in")

    ax = axes[1, 1]
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

    if title_text is not None:
        fig.suptitle(title_text)
    fig.tight_layout()
    return fig


def run_single_analysis(model_dir: Path, snapshot: int, output_dir: Path, radius_max: float, nbin: int, sigma_radius: float):
    output_dir.mkdir(parents=True, exist_ok=True)

    _, snapshot_dir, tsnap, dwarf_df = load_centered_dwarf_dataframe(model_dir, snapshot)
    summary = compute_summary(dwarf_df, radius_max=radius_max, nbin=nbin, sigma_radius=sigma_radius)
    write_density_snapshot(output_dir, snapshot, tsnap, summary)

    wlm_r, wlm_f = load_wlm_reference(snapshot_dir)
    model_name = detect_model_name(model_dir)
    fd_value, sfr_eff = load_run_metadata()
    title_text = None
    if fd_value is not None and sfr_eff is not None:
        title_text = f"IC_{model_name}\nfeedback = {fd_value / 8.3:.2f} median, sfr_eff = {sfr_eff:.3f}"
    elif model_name:
        title_text = f"IC_{model_name}"

    fig = make_figure(
        dwarf_df,
        tsnap,
        summary,
        density_history=load_density_history(output_dir),
        wlm_r=wlm_r,
        wlm_f=wlm_f,
        title_text=title_text,
    )

    stem = f"snapshot_{snapshot:03d}_density_profile"
    fig_path = output_dir / f"{stem}.png"
    fig.savefig(fig_path, dpi=200)
    plt.close(fig)

    return {
        "snapshot": snapshot,
        "figure_path": fig_path,
    }


def run_single_analysis_task(task):
    model_dir, snapshot, output_dir, radius_max, nbin, sigma_radius = task
    result = run_single_analysis(
        model_dir=Path(model_dir),
        snapshot=snapshot,
        output_dir=Path(output_dir),
        radius_max=radius_max,
        nbin=nbin,
        sigma_radius=sigma_radius,
    )
    return {
        "snapshot": result["snapshot"],
        "figure_path": str(result["figure_path"]),
    }


def main():
    args = parse_args()
    explicit_mode = args.snapshot is not None
    output_dir = resolve_output_dir(args.output_dir, explicit_mode=explicit_mode)

    if explicit_mode:
        if args.model_dir is None:
            raise SystemExit("--snapshot requires --model-dir in explicit single-snapshot mode.")

        model_dir = Path(args.model_dir).resolve()
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
        )
        fd_value, sfr_eff = load_run_metadata()
        create_density_evolution_csv(output_dir, fd_value=fd_value, sfr_eff=sfr_eff)
        cleanup_density_snapshots(output_dir)
        print(f"Saved figure: {result['figure_path']}")
        return

    model_dir = resolve_legacy_model_dir(args.model_dir)
    snapshots = detect_batch_snapshots(model_dir, args.range)
    print(f"Legacy batch mode: model_dir={model_dir}, snapshots={snapshots[0]}-{snapshots[-1]}, count={len(snapshots)}")

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
        return

    tasks = [
        (str(model_dir), snapshot, str(output_dir), args.radius_max, args.nbin, args.sigma_radius)
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
            with get_context("fork").Pool(processes=args.processes) as pool:
                for result in pool.imap_unordered(run_single_analysis_task, tasks):
                    results.append(result)
                    progress.update(1)
    finally:
        progress.close()

    fd_value, sfr_eff = load_run_metadata()
    create_density_evolution_csv(output_dir, fd_value=fd_value, sfr_eff=sfr_eff)
    cleanup_density_snapshots(output_dir)

    print(f"Completed {len(results)} snapshot(s). Output directory: {output_dir}")
    for result in sorted(results, key=lambda item: item["snapshot"]):
        print(f"[{result['snapshot']:03d}] {result['figure_path']}")


if __name__ == "__main__":
    main()
