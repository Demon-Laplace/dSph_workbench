"""Measure the radial CGM structure of the MW52 isolated MW+CGM run.

The script reads GIZMO HDF5 snapshots and writes non-destructive profile
tables.  It can process one snapshot directly or several snapshots with an
explicit multiprocessing worker count.
"""

from __future__ import annotations

import argparse
import json
from multiprocessing import get_context
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


KPC_CM = 3.085678e21
MSUN_G = 1.989e33
PROTON_G = 1.67262192369e-24
BOLTZMANN = 1.380649e-16
HYDROGEN_FRACTION = 0.76
MEAN_MOLECULAR_WEIGHT_HOT = 0.59
GAMMA_MINUS_ONE = 2.0 / 3.0

UNIT_LENGTH_CM = KPC_CM
UNIT_MASS_G = 1.0e10 * MSUN_G
UNIT_VELOCITY_CM_S = 1.0e5
UNIT_DENSITY_CGS = UNIT_MASS_G / UNIT_LENGTH_CM**3

RADIAL_EDGES_KPC = np.array(
    [0, 10, 20, 30, 50, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500],
    dtype=float,
)


def temperature_from_internal_energy(internal_energy: np.ndarray) -> np.ndarray:
    return (
        GAMMA_MINUS_ONE
        * np.asarray(internal_energy, dtype=np.float64)
        * UNIT_VELOCITY_CM_S**2
        * MEAN_MOLECULAR_WEIGHT_HOT
        * PROTON_G
        / BOLTZMANN
    )


def weighted_center(group: h5py.Group) -> np.ndarray:
    coordinates = np.asarray(group["Coordinates"], dtype=np.float64)
    if "Masses" in group:
        masses = np.asarray(group["Masses"], dtype=np.float64)
        return np.average(coordinates, axis=0, weights=masses)
    return np.mean(coordinates, axis=0)


def shell_profile(
    radius: np.ndarray,
    density_code: np.ndarray,
    temperature: np.ndarray,
    masses_code: np.ndarray,
    radial_velocity: np.ndarray,
) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    hot = np.isfinite(temperature) & (temperature > 2.0e4)
    density_nh = (
        np.asarray(density_code, dtype=np.float64)
        * UNIT_DENSITY_CGS
        * HYDROGEN_FRACTION
        / PROTON_G
    )
    masses_msun = np.asarray(masses_code, dtype=np.float64) * 1.0e10

    for inner, outer in zip(RADIAL_EDGES_KPC[:-1], RADIAL_EDGES_KPC[1:]):
        selection = hot & (radius >= inner) & (radius < outer)
        count = int(np.count_nonzero(selection))
        shell_volume_kpc3 = 4.0 * np.pi / 3.0 * (outer**3 - inner**3)
        shell_mass_msun = float(np.sum(masses_msun[selection]))
        shell_density_nh = (
            shell_mass_msun
            * MSUN_G
            * HYDROGEN_FRACTION
            / (shell_volume_kpc3 * KPC_CM**3 * PROTON_G)
        )
        if count:
            density_values = density_nh[selection]
            temperature_values = temperature[selection]
            radial_values = radial_velocity[selection]
            p16, median, p84 = np.percentile(density_values, [16, 50, 84])
            t16, tmedian, t84 = np.percentile(temperature_values, [16, 50, 84])
            vr16, vrmedian, vr84 = np.percentile(radial_values, [16, 50, 84])
        else:
            p16 = median = p84 = np.nan
            t16 = tmedian = t84 = np.nan
            vr16 = vrmedian = vr84 = np.nan
        rows.append(
            {
                "r_inner_kpc": inner,
                "r_outer_kpc": outer,
                "r_mid_kpc": 0.5 * (inner + outer),
                "particle_count": count,
                "median_nH_cm3": median,
                "p16_nH_cm3": p16,
                "p84_nH_cm3": p84,
                "shell_mass_density_nH_cm3": shell_density_nh,
                "hot_mass_shell_msun": shell_mass_msun,
                "median_temperature_K": tmedian,
                "p16_temperature_K": t16,
                "p84_temperature_K": t84,
                "median_radial_velocity_km_s": vrmedian,
                "p16_radial_velocity_km_s": vr16,
                "p84_radial_velocity_km_s": vr84,
            }
        )
    return rows


def analyze_snapshot(arguments: tuple[str, int, str]) -> tuple[int, str, dict]:
    model_dir_text, snapshot_number, output_dir_text = arguments
    model_dir = Path(model_dir_text)
    output_dir = Path(output_dir_text)
    snapshot_path = model_dir / "output" / f"snapshot_{snapshot_number:03d}.hdf5"
    if not snapshot_path.exists():
        raise FileNotFoundError(snapshot_path)

    with h5py.File(snapshot_path, "r") as handle:
        header = handle["Header"].attrs
        time_gyr = float(header["Time"])
        if "PartType2" in handle:
            center = weighted_center(handle["PartType2"])
        elif "PartType3" in handle:
            center = weighted_center(handle["PartType3"])
        elif "PartType1" in handle:
            center = weighted_center(handle["PartType1"])
        else:
            center = np.array([500.0, 500.0, 500.0])

        gas = handle["PartType0"]
        coordinates = np.asarray(gas["Coordinates"], dtype=np.float64)
        velocities = np.asarray(gas["Velocities"], dtype=np.float64)
        density = np.asarray(gas["Density"], dtype=np.float64)
        internal_energy = np.asarray(gas["InternalEnergy"], dtype=np.float64)
        masses = np.asarray(gas["Masses"], dtype=np.float64)

        displacement = coordinates - center
        radius = np.linalg.norm(displacement, axis=1)
        safe_radius = np.where(radius > 0.0, radius, 1.0)
        radial_velocity = np.sum(velocities * displacement, axis=1) / safe_radius
        temperature = temperature_from_internal_energy(internal_energy)
        rows = shell_profile(
            radius, density, temperature, masses, radial_velocity
        )

        hot = np.isfinite(temperature) & (temperature > 2.0e4)
        mass_msun = masses * 1.0e10
        enclosed = {}
        for radial_limit in [50, 100, 150, 200, 250, 300, 400, 500]:
            enclosed[str(radial_limit)] = float(
                np.sum(mass_msun[hot & (radius < radial_limit)])
            )

        unique_masses, unique_counts = np.unique(masses, return_counts=True)
        mass_groups = sorted(
            zip(unique_masses.tolist(), unique_counts.tolist()),
            key=lambda item: item[1],
            reverse=True,
        )[:10]
        metadata = {
            "snapshot": snapshot_number,
            "time_gyr": time_gyr,
            "center_kpc": center.tolist(),
            "number_of_gas_particles": int(len(radius)),
            "number_of_hot_particles": int(np.count_nonzero(hot)),
            "hot_mass_within_500kpc_msun": enclosed["500"],
            "hot_mass_enclosed_msun": enclosed,
            "most_common_particle_masses_code": [
                {"mass_code": float(mass), "count": int(count)}
                for mass, count in mass_groups
            ],
            "header_num_part_total": np.asarray(
                header["NumPart_Total"], dtype=np.int64
            ).tolist(),
        }

    profile = pd.DataFrame(rows)
    profile.insert(0, "time_gyr", time_gyr)
    profile.insert(0, "snapshot", snapshot_number)
    csv_path = output_dir / f"profile_snapshot_{snapshot_number:03d}.csv"
    json_path = output_dir / f"metadata_snapshot_{snapshot_number:03d}.json"
    profile.to_csv(csv_path, index=False)
    json_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return snapshot_number, str(csv_path), metadata


def gamma0_nh(radius_kpc: np.ndarray, mass_1e10_msun: float, scale_kpc: float) -> np.ndarray:
    mass_msun = mass_1e10_msun * 1.0e10
    density_msun_kpc3 = (
        3.0
        * mass_msun
        * scale_kpc
        / (4.0 * np.pi * (radius_kpc + scale_kpc) ** 4)
    )
    return (
        density_msun_kpc3
        * MSUN_G
        / KPC_CM**3
        * HYDROGEN_FRACTION
        / PROTON_G
    )


def make_combined_figure(combined: pd.DataFrame, output_dir: Path) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(13.2, 4.1), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, combined["snapshot"].nunique()))
    for color, (snapshot, profile) in zip(colors, combined.groupby("snapshot")):
        label = f"{profile['time_gyr'].iloc[0]:.2f} Gyr"
        radius = profile["r_mid_kpc"].to_numpy()
        axes[0].plot(
            radius,
            profile["shell_mass_density_nH_cm3"].to_numpy(),
            color=color,
            lw=2,
            label=label,
        )
        axes[0].scatter(
            radius,
            profile["median_nH_cm3"].to_numpy(),
            color=color,
            s=14,
            alpha=0.7,
        )
        axes[1].plot(
            radius,
            profile["median_temperature_K"].to_numpy(),
            color=color,
            lw=2,
            label=label,
        )
        axes[2].plot(
            radius,
            profile["median_radial_velocity_km_s"].to_numpy(),
            color=color,
            lw=2,
            label=label,
        )

    analytic_radius = np.geomspace(10, 500, 300)
    axes[0].plot(
        analytic_radius,
        gamma0_nh(analytic_radius, 68.0, 130.0),
        color="black",
        ls="--",
        lw=1.5,
        label=r"input $\gamma=0$: $M=68$, $a=130$",
    )
    axes[0].set_yscale("log")
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Galactocentric radius [kpc]")
    axes[0].set_ylabel(r"$n_{\rm H}$ [cm$^{-3}$]")
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Galactocentric radius [kpc]")
    axes[1].set_ylabel("Median gas temperature [K]")

    axes[2].set_xscale("log")
    axes[2].axhline(0.0, color="0.5", lw=1)
    axes[2].set_xlabel("Galactocentric radius [kpc]")
    axes[2].set_ylabel("Median radial gas velocity [km/s]")

    for axis in axes:
        axis.set_xlim(10, 500)
        axis.tick_params(direction="in", which="both", top=True, right=True)
    figure.savefig(output_dir / "mw52_cgm_evolution_profiles.png", dpi=220)
    plt.close(figure)


def parse_snapshot_list(value: str) -> list[int]:
    return [int(item.strip()) for item in value.split(",") if item.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model_dir", type=Path)
    parser.add_argument("--snapshots", default="0")
    parser.add_argument("--processes", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    snapshots = parse_snapshot_list(args.snapshots)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tasks = [
        (str(args.model_dir.resolve()), snapshot, str(args.output_dir.resolve()))
        for snapshot in snapshots
    ]
    if args.processes == 1:
        results = [analyze_snapshot(task) for task in tasks]
    else:
        with get_context("spawn").Pool(processes=args.processes) as pool:
            results = pool.map(analyze_snapshot, tasks)

    profiles = [pd.read_csv(csv_path) for _, csv_path, _ in sorted(results)]
    combined = pd.concat(profiles, ignore_index=True)
    combined.to_csv(args.output_dir / "mw52_cgm_profiles.csv", index=False)
    make_combined_figure(combined, args.output_dir)
    print(json.dumps([metadata for _, _, metadata in sorted(results)], indent=2))


if __name__ == "__main__":
    main()
