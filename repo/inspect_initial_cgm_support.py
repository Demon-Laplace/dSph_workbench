#!/usr/bin/env python3
"""Inspect thermal and kinetic support of the MW52 CGM in one snapshot."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np


PROTON_G = 1.67262192369e-24
BOLTZMANN = 1.380649e-16
MEAN_MOLECULAR_WEIGHT_HOT = 0.59
UNIT_VELOCITY_CM_S = 1.0e5
GAMMA_MINUS_ONE = 2.0 / 3.0
MW_CGM_PARTICLE_MASS_CODE = 3.666669954327517e-6
RADIAL_EDGES_KPC = np.array(
    [0, 25, 50, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500],
    dtype=float,
)


def temperature_from_internal_energy(internal_energy: np.ndarray) -> np.ndarray:
    return (
        GAMMA_MINUS_ONE
        * internal_energy
        * UNIT_VELOCITY_CM_S**2
        * MEAN_MOLECULAR_WEIGHT_HOT
        * PROTON_G
        / BOLTZMANN
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--processes", type=int, default=1)
    args = parser.parse_args()
    if args.processes != 1:
        raise ValueError("This single-snapshot diagnostic must use --processes 1")

    with h5py.File(args.snapshot, "r") as handle:
        star = handle["PartType2"]
        star_mass = np.asarray(star["Masses"], dtype=np.float64)
        center = np.average(
            np.asarray(star["Coordinates"], dtype=np.float64),
            axis=0,
            weights=star_mass,
        )

        gas = handle["PartType0"]
        coordinates = np.asarray(gas["Coordinates"], dtype=np.float64)
        velocities = np.asarray(gas["Velocities"], dtype=np.float64)
        masses = np.asarray(gas["Masses"], dtype=np.float64)
        internal_energy = np.asarray(gas["InternalEnergy"], dtype=np.float64)
        time_gyr = float(handle["Header"].attrs["Time"])

    displacement = coordinates - center
    radius = np.linalg.norm(displacement, axis=1)
    safe_radius = np.where(radius > 0, radius, 1.0)
    radial_velocity = np.sum(velocities * displacement, axis=1) / safe_radius
    mw_cgm = np.isclose(masses, MW_CGM_PARTICLE_MASS_CODE, rtol=2e-5, atol=0.0)
    temperature = temperature_from_internal_energy(internal_energy)

    rows: list[dict[str, float | int]] = []
    for inner, outer in zip(RADIAL_EDGES_KPC[:-1], RADIAL_EDGES_KPC[1:]):
        selected = mw_cgm & (radius >= inner) & (radius < outer)
        count = int(np.count_nonzero(selected))
        selected_velocity = velocities[selected]
        selected_vr = radial_velocity[selected]
        selected_u = internal_energy[selected]
        selected_temperature = temperature[selected]
        mean_velocity = np.mean(selected_velocity, axis=0)
        variance_components = np.var(selected_velocity, axis=0)
        sigma_1d = float(np.sqrt(np.sum(variance_components) / 3.0))
        equivalent_temperature = (
            MEAN_MOLECULAR_WEIGHT_HOT
            * PROTON_G
            * (sigma_1d * UNIT_VELOCITY_CM_S) ** 2
            / BOLTZMANN
        )
        vr16, vr50, vr84 = np.percentile(selected_vr, [16, 50, 84])
        rows.append(
            {
                "r_inner_kpc": inner,
                "r_outer_kpc": outer,
                "r_mid_kpc": 0.5 * (inner + outer),
                "particle_count": count,
                "mean_vx_km_s": mean_velocity[0],
                "mean_vy_km_s": mean_velocity[1],
                "mean_vz_km_s": mean_velocity[2],
                "sigma_1d_km_s": sigma_1d,
                "kinetic_equivalent_temperature_K": equivalent_temperature,
                "median_internal_energy_code": float(np.median(selected_u)),
                "median_thermal_temperature_K": float(np.median(selected_temperature)),
                "p16_radial_velocity_km_s": vr16,
                "median_radial_velocity_km_s": vr50,
                "p84_radial_velocity_km_s": vr84,
            }
        )

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_prefix.with_suffix(".csv")
    json_path = args.output_prefix.with_suffix(".json")
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "snapshot": str(args.snapshot),
        "time_gyr": time_gyr,
        "center_kpc": center.tolist(),
        "mw_cgm_particle_count": int(np.count_nonzero(mw_cgm)),
        "mw_cgm_internal_energy_min_code": float(np.min(internal_energy[mw_cgm])),
        "mw_cgm_internal_energy_median_code": float(np.median(internal_energy[mw_cgm])),
        "mw_cgm_internal_energy_max_code": float(np.max(internal_energy[mw_cgm])),
    }
    json_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
