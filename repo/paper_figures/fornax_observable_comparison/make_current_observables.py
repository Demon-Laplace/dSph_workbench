"""Build the current-epoch three-panel observable comparison for the paper.

The script must be run with a Fornax model directory as its working directory.
By default it selects the elinfo row whose heliocentric distance is closest to
``variable.d_today`` and processes only that snapshot.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Optional


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "fornax_figure2_mplconfig")
)

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from basefunc import Analysis
from data_processing import DataProcessor
from snapshot_context import prepare_snapshot_context
from snapshot_metrics import old_star_projected_kinematics
from variable import d_today, fornax_core_radius, r_pc, sigma, err_sigma


SIMULATION_COLOR = "#171717"
OBSERVATION_COLOR = "#7b3294"
SFH_OBSERVATION_COLOR = "#d95f02"
MASS_TO_LIGHT_V = 2.6
SFH_SFR_UNIT_MSUN_PER_YEAR = 1.0e-4
SFH_APERTURE_DEG = 0.8
SFH_BIN_WIDTH_GYR = 0.5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "model_dir",
        nargs="?",
        type=Path,
        default=Path.cwd(),
        help="Model directory containing output/ and elinfo_<model>.csv.",
    )
    parser.add_argument("--snapshot", type=int, default=None)
    parser.add_argument("--distance", type=float, default=d_today)
    parser.add_argument("--mass-to-light", type=float, default=MASS_TO_LIGHT_V)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def load_elinfo(model_dir: Path) -> tuple[pd.DataFrame, Path]:
    model_name = model_dir.name
    elinfo_path = model_dir / f"elinfo_{model_name}.csv"
    if not elinfo_path.exists():
        matches = sorted(model_dir.glob("elinfo_*.csv"))
        if len(matches) != 1:
            raise FileNotFoundError(
                f"Expected one elinfo file under {model_dir}; found {len(matches)}."
            )
        elinfo_path = matches[0]
    return DataProcessor.read_csv_with_comments(elinfo_path), elinfo_path


def select_current_row(
    elinfo: pd.DataFrame,
    target_distance_kpc: float,
    snapshot_override: Optional[int],
) -> pd.Series:
    if snapshot_override is not None:
        selected = elinfo.loc[elinfo["numsp"].astype(int) == int(snapshot_override)]
        if len(selected) != 1:
            raise ValueError(f"Snapshot {snapshot_override} is not unique in elinfo.")
        return selected.iloc[0]
    finite = np.isfinite(elinfo["distance"].to_numpy(dtype=float))
    if not np.any(finite):
        raise ValueError("elinfo contains no finite heliocentric distances.")
    candidates = elinfo.loc[finite]
    index = (candidates["distance"] - target_distance_kpc).abs().idxmin()
    return elinfo.loc[index]


def observation_paths(model_dir: Path) -> tuple[Path, Path]:
    data_root = model_dir.parent
    profile_path = data_root / "Fornax_surface_brightness_profile.csv"
    sfh_path = data_root / "SFH_Fornax.csv"
    for path in (profile_path, sfh_path):
        if not path.exists():
            raise FileNotFoundError(path)
    return profile_path, sfh_path


def corrected_surface_brightness_intervals(
    observation: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    required = {
        "r_outer_kpc",
        "r_mid_kpc",
        "mu_v_mag_arcsec2",
        "mu_v_err_total",
    }
    missing = required - set(observation.columns)
    if missing:
        raise ValueError(f"Surface-brightness table is missing {sorted(missing)}")

    outer = observation["r_outer_kpc"].to_numpy(dtype=float)
    inner = np.r_[0.0, outer[:-1]]
    mid = observation["r_mid_kpc"].to_numpy(dtype=float)
    valid = (
        np.isfinite(inner)
        & np.isfinite(outer)
        & np.isfinite(mid)
        & (outer > inner)
        & (mid >= inner)
        & (mid <= outer)
    )
    intervals = np.column_stack([inner[valid], outer[valid]])
    return intervals, mid[valid], valid


def surface_brightness_profile(
    snapshot: dict,
    row: pd.Series,
    observation: pd.DataFrame,
    mass_to_light: float,
) -> dict[str, np.ndarray]:
    kinematics = old_star_projected_kinematics(snapshot)
    old_mask = np.asarray(kinematics["old_local_mask"], dtype=bool)
    star_mass = snapshot["df"].loc[
        snapshot["total_dw_star_mask"], "m"
    ].to_numpy(dtype=float)[old_mask]
    intervals, observed_mid, interval_mask = corrected_surface_brightness_intervals(
        observation
    )
    _, simulated_mu = Analysis.radial_magnitude_profile(
        kinematics["x_kpc"],
        kinematics["y_kpc"],
        mass=star_mass,
        d_kpc=float(snapshot["d_mean"]),
        ep=float(row["eps"]),
        pa=float(row["pa"]),
        center_x=float(row["shape_center_x_kpc"]),
        center_y=float(row["shape_center_y_kpc"]),
        r_intervals=intervals,
        m_l_relation=mass_to_light,
        mag_sys="vega",
    )
    return {
        "radius_kpc": observed_mid,
        "simulated_mu": simulated_mu,
        "observed_mu": observation.loc[
            interval_mask, "mu_v_mag_arcsec2"
        ].to_numpy(dtype=float),
        "observed_mu_error": observation.loc[
            interval_mask, "mu_v_err_total"
        ].to_numpy(dtype=float),
    }


def radial_edges_from_centres(centres: np.ndarray) -> np.ndarray:
    centres = np.asarray(centres, dtype=float)
    if centres.ndim != 1 or centres.size < 2 or np.any(np.diff(centres) <= 0):
        raise ValueError("Radial centres must be a strictly increasing 1D array.")
    midpoints = 0.5 * (centres[:-1] + centres[1:])
    last = centres[-1] + 0.5 * (centres[-1] - centres[-2])
    return np.r_[0.0, midpoints, last]


def velocity_dispersion_profile(snapshot: dict) -> dict[str, np.ndarray]:
    kinematics = old_star_projected_kinematics(snapshot)
    radius = np.hypot(kinematics["x_kpc"], kinematics["y_kpc"])
    velocity = np.asarray(kinematics["vlos_detrended"], dtype=float)
    centres = np.asarray(r_pc, dtype=float) / 1000.0
    edges = radial_edges_from_centres(centres)
    simulated = np.full(centres.shape, np.nan, dtype=float)
    simulated_error = np.full(centres.shape, np.nan, dtype=float)
    counts = np.zeros(centres.shape, dtype=int)
    for index, (left, right) in enumerate(zip(edges[:-1], edges[1:])):
        selected = (
            np.isfinite(radius)
            & np.isfinite(velocity)
            & (radius >= left)
            & (radius < right)
        )
        counts[index] = int(np.sum(selected))
        if counts[index] < 3:
            continue
        value = float(np.std(velocity[selected], ddof=1))
        simulated[index] = value
        simulated_error[index] = value / np.sqrt(2.0 * (counts[index] - 1.0))
    return {
        "radius_kpc": centres,
        "simulated_sigma": simulated,
        "simulated_sigma_error": simulated_error,
        "particle_count": counts,
        "observed_sigma": np.asarray(sigma, dtype=float),
        "observed_sigma_error": np.asarray(err_sigma, dtype=float),
        "velocity_gradient_amplitude": np.array(
            [kinematics["velocity_gradient"]["grad_amp"]], dtype=float
        ),
    }


def fixed_lookback_edges(current_time_gyr: float) -> np.ndarray:
    upper_edge = np.ceil(current_time_gyr / SFH_BIN_WIDTH_GYR) * SFH_BIN_WIDTH_GYR
    return np.arange(0.0, upper_edge + SFH_BIN_WIDTH_GYR, SFH_BIN_WIDTH_GYR)


def rebin_observed_sfh(
    observation: pd.DataFrame,
    edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    age = observation["age"].to_numpy(dtype=float)
    sfr = observation["sfr"].to_numpy(dtype=float)
    error = observation["esfr"].to_numpy(dtype=float)
    valid = np.isfinite(age) & np.isfinite(sfr)
    sfr_sum, _ = np.histogram(age[valid], bins=edges, weights=sfr[valid])
    count, _ = np.histogram(age[valid], bins=edges)
    mean_sfr = np.divide(
        sfr_sum,
        count,
        out=np.zeros_like(sfr_sum, dtype=float),
        where=count > 0,
    )
    valid_error = valid & np.isfinite(error)
    error_sum, _ = np.histogram(
        age[valid_error], bins=edges, weights=error[valid_error]
    )
    error_count, _ = np.histogram(age[valid_error], bins=edges)
    mean_error = np.divide(
        error_sum,
        error_count,
        out=np.zeros_like(error_sum, dtype=float),
        where=error_count > 0,
    )
    return mean_sfr, mean_error


def recent_star_formation_history(
    snapshot: dict,
    current_time_gyr: float,
    observation: pd.DataFrame,
    row: pd.Series,
) -> dict[str, np.ndarray]:
    required = {"age", "sfr", "esfr"}
    missing = required - set(observation.columns)
    if missing:
        raise ValueError(f"Observed SFH table is missing {sorted(missing)}")
    edges = fixed_lookback_edges(current_time_gyr)
    dwarf_stars = snapshot["df"].loc[
        snapshot["total_dw_star_mask"], ["birth", "m"]
    ]
    birth = dwarf_stars["birth"].to_numpy(dtype=float)
    mass = dwarf_stars["m"].to_numpy(dtype=float)
    elliptical_radius = Analysis.projected_elliptical_radius(
        np.asarray(snapshot["x_kpc"], dtype=float),
        np.asarray(snapshot["y_kpc"], dtype=float),
        ep=float(row["eps"]),
        pa=float(row["pa"]),
        center_x=float(row["shape_center_x_kpc"]),
        center_y=float(row["shape_center_y_kpc"]),
    )
    aperture_radius_kpc = np.radians(SFH_APERTURE_DEG) * float(snapshot["d_mean"])
    inside_observed_aperture = elliptical_radius <= aperture_radius_kpc
    new_stars = (
        np.isfinite(birth)
        & np.isfinite(mass)
        & (birth > 0.0)
        & (birth <= current_time_gyr)
        & (mass > 0.0)
        & (mass < 1.0e4)
        & inside_observed_aperture
    )
    lookback = current_time_gyr - birth[new_stars]
    formed_mass, _ = np.histogram(lookback, bins=edges, weights=mass[new_stars])
    bin_width_gyr = np.diff(edges)
    simulated_sfr = (
        formed_mass
        / (bin_width_gyr * 1.0e9)
        / SFH_SFR_UNIT_MSUN_PER_YEAR
    )
    observed_sfr, observed_sfr_error = rebin_observed_sfh(observation, edges)
    return {
        "lookback_gyr": edges[:-1],
        "bin_edges_gyr": edges,
        "simulated_sfr": simulated_sfr,
        "observed_sfr": observed_sfr,
        "observed_sfr_error": observed_sfr_error,
        "formed_mass_msun": formed_mass,
        "new_star_particle_count": np.array([np.sum(new_stars)], dtype=int),
        "aperture_deg": np.array([SFH_APERTURE_DEG], dtype=float),
        "aperture_radius_kpc": np.array([aperture_radius_kpc], dtype=float),
    }


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 7.2,
            "axes.labelsize": 7.5,
            "axes.titlesize": 8.0,
            "legend.fontsize": 6.2,
            "xtick.labelsize": 6.6,
            "ytick.labelsize": 6.6,
            "axes.linewidth": 0.7,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )


def style_axis(
    axis: plt.Axes,
    panel_label: str,
    title: str,
    panel_label_x: float = 0.975,
    panel_label_ha: str = "right",
) -> None:
    axis.minorticks_on()
    axis.tick_params(which="major", length=3.2, width=0.65)
    axis.tick_params(which="minor", length=1.8, width=0.5)
    axis.text(
        panel_label_x,
        0.97,
        panel_label,
        transform=axis.transAxes,
        ha=panel_label_ha,
        va="top",
        fontweight="bold",
        fontsize=8.0,
        zorder=20,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 0.6},
    )
    axis.set_title(title, pad=3.0)


def draw_figure(
    surface: dict[str, np.ndarray],
    velocity: dict[str, np.ndarray],
    sfh: dict[str, np.ndarray],
    row: pd.Series,
    mass_to_light: float,
    output_path: Path,
    dpi: int,
) -> None:
    configure_matplotlib()
    fig, axes = plt.subplots(1, 3, figsize=(7.05, 2.35))

    axis = axes[0]
    obs_good = (
        np.isfinite(surface["radius_kpc"])
        & np.isfinite(surface["observed_mu"])
        & np.isfinite(surface["observed_mu_error"])
        & (surface["radius_kpc"] <= 5.1)
    )
    axis.errorbar(
        surface["radius_kpc"][obs_good],
        surface["observed_mu"][obs_good],
        yerr=surface["observed_mu_error"][obs_good],
        fmt="o",
        ms=2.2,
        mew=0,
        color=OBSERVATION_COLOR,
        ecolor=OBSERVATION_COLOR,
        elinewidth=0.55,
        alpha=0.55,
        label="Yang et al. (2022)",
        zorder=1,
    )
    sim_radius = surface["radius_kpc"].copy()
    sim_mu = surface["simulated_mu"].copy()
    sim_mu[sim_radius > 5.1] = np.nan
    axis.plot(
        sim_radius,
        sim_mu,
        color=SIMULATION_COLOR,
        lw=1.35,
        label=rf"Simulation ($M_\star/L_V={mass_to_light:g}$)",
        zorder=3,
    )
    finite_sim = np.flatnonzero(np.isfinite(sim_mu))
    marker_indices = finite_sim[:: max(1, finite_sim.size // 18)]
    axis.plot(
        sim_radius[marker_indices],
        sim_mu[marker_indices],
        linestyle="none",
        marker="o",
        ms=2.0,
        color=SIMULATION_COLOR,
        zorder=3.1,
    )
    axis.set(xlabel=r"$R_{\rm ell}$ (kpc)", ylabel=r"$\mu_V$ (mag arcsec$^{-2}$)")
    axis.set_xlim(0.0, 5.1)
    axis.set_ylim(36.3, 22.0)
    axis.legend(loc="lower left", frameon=False, handlelength=1.5)
    style_axis(axis, "(a)", "Surface brightness")

    axis = axes[1]
    axis.errorbar(
        velocity["radius_kpc"],
        velocity["observed_sigma"],
        yerr=velocity["observed_sigma_error"],
        fmt="o",
        ms=2.5,
        color=OBSERVATION_COLOR,
        ecolor=OBSERVATION_COLOR,
        elinewidth=0.65,
        capsize=1.2,
        alpha=0.55,
        label="Walker et al. (2009)",
        zorder=1,
    )
    sigma_good = (
        np.isfinite(velocity["simulated_sigma"])
        & (velocity["particle_count"] >= 3)
    )
    sim_sigma = velocity["simulated_sigma"].copy()
    sim_sigma[~sigma_good] = np.nan
    sim_error = velocity["simulated_sigma_error"].copy()
    lower = sim_sigma - sim_error
    upper = sim_sigma + sim_error
    axis.fill_between(
        velocity["radius_kpc"],
        lower,
        upper,
        where=np.isfinite(lower) & np.isfinite(upper),
        color="#bdbdbd",
        alpha=0.55,
        linewidth=0,
        zorder=2,
    )
    axis.plot(
        velocity["radius_kpc"],
        sim_sigma,
        color=SIMULATION_COLOR,
        marker="o",
        ms=2.1,
        lw=1.25,
        label="Simulation",
        zorder=3,
    )
    axis.set(xlabel=r"$R$ (kpc)", ylabel=r"$\sigma_{\rm los}$ (km s$^{-1}$)")
    axis.set_xlim(0.0, 1.8)
    axis.set_ylim(5.0, 13.2)
    axis.legend(loc="lower left", frameon=False, handlelength=1.5)
    style_axis(axis, "(b)", "Velocity dispersion")

    axis = axes[2]
    axis.step(
        sfh["bin_edges_gyr"],
        np.r_[sfh["observed_sfr"], sfh["observed_sfr"][-1]],
        where="post",
        color=SFH_OBSERVATION_COLOR,
        lw=1.35,
        label="de Boer et al. (2012)",
        zorder=2,
    )
    axis.step(
        sfh["bin_edges_gyr"],
        np.r_[sfh["simulated_sfr"], sfh["simulated_sfr"][-1]],
        where="post",
        color=SIMULATION_COLOR,
        lw=1.35,
        label="Simulation",
        zorder=3,
    )
    axis.set(
        xlabel="Lookback time (Gyr)",
        ylabel=r"SFR ($10^{-4}\,M_\odot\,{\rm yr}^{-1}$)",
    )
    axis.set_xlim(0.0, float(row["age"]))
    ymax = max(
        1.0,
        float(np.nanmax(sfh["observed_sfr"])),
        float(np.nanmax(sfh["simulated_sfr"])),
    )
    axis.set_ylim(0.0, 1.18 * ymax)
    axis.legend(loc="upper center", frameon=False, handlelength=1.5)
    style_axis(axis, "(c)", "Recent star formation")

    fig.subplots_adjust(left=0.075, right=0.995, bottom=0.22, top=0.89, wspace=0.34)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight", pad_inches=0.015)
    plt.close(fig)


def serializable_array(values: np.ndarray) -> list:
    output = []
    for value in np.asarray(values).ravel():
        if isinstance(value, (np.integer, int)):
            output.append(int(value))
        elif np.isfinite(value):
            output.append(float(value))
        else:
            output.append(None)
    return output


def main() -> None:
    args = parse_args()
    model_dir = args.model_dir.resolve()
    elinfo, elinfo_path = load_elinfo(model_dir)
    row = select_current_row(elinfo, args.distance, args.snapshot)
    snapshot_number = int(row["numsp"])
    profile_path, sfh_path = observation_paths(model_dir)

    snapshot = prepare_snapshot_context(
        folder_path=str(model_dir / "output"),
        snapshot_num=snapshot_number,
        core_radius=fornax_core_radius,
        include_star_birth=True,
    )
    surface_observation = pd.read_csv(profile_path)
    sfh_observation = pd.read_csv(sfh_path)
    surface = surface_brightness_profile(
        snapshot, row, surface_observation, args.mass_to_light
    )
    velocity = velocity_dispersion_profile(snapshot)
    sfh = recent_star_formation_history(
        snapshot, float(row["age"]), sfh_observation, row
    )

    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else model_dir / "paper_figure2"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"fornax_current_observables_snapshot_{snapshot_number:03d}"
    output_png = output_dir / f"{stem}.png"
    metadata_path = output_dir / f"{stem}_metadata.json"
    draw_figure(
        surface,
        velocity,
        sfh,
        row,
        args.mass_to_light,
        output_png,
        args.dpi,
    )

    metadata = {
        "model": model_dir.name,
        "snapshot": snapshot_number,
        "target_heliocentric_distance_kpc": float(args.distance),
        "elinfo_heliocentric_distance_kpc": float(row["distance"]),
        "transformed_snapshot_distance_kpc": float(snapshot["d_mean"]),
        "simulation_time_gyr": float(row["age"]),
        "mass_to_light_v": float(args.mass_to_light),
        "elinfo_path": str(elinfo_path),
        "surface_brightness_observation": str(profile_path),
        "sfh_observation": str(sfh_path),
        "surface_brightness": {
            key: serializable_array(value) for key, value in surface.items()
        },
        "velocity_dispersion": {
            key: serializable_array(value) for key, value in velocity.items()
        },
        "recent_sfh": {key: serializable_array(value) for key, value in sfh.items()},
        "sfh_sfr_unit_msun_per_year": SFH_SFR_UNIT_MSUN_PER_YEAR,
    }
    metadata_path.write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(output_png)
    print(metadata_path)


if __name__ == "__main__":
    main()
