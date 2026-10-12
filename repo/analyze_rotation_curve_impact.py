#!/usr/bin/env python3
"""Estimate how alternative MW52 CGM profiles change the circular-speed curve."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


G_KPC_KMS2_MSUN = 4.30091e-6
MASS_UNIT_MSUN = 1.0e10
CURRENT_CGM_MASS_PARAMETER = 68.0
CURRENT_CGM_SCALE_KPC = 130.0
MW_CGM_PARTICLE_MASS_CODE = 3.666669954327517e-6
RADII_KPC = np.array(
    [5, 8, 10, 15, 20, 25, 30, 40, 50, 75, 100, 125, 140, 150, 175,
     200, 250, 300, 338, 400, 450, 500],
    dtype=float,
)
CANDIDATES = [
    ("current_68_130", 68.0, 130.0),
    ("strict_rotation_65_120", 65.0, 120.0),
    ("rotation_preserving_62_110", 62.0, 110.0),
    ("intermediate_61_107", 61.0, 107.0),
    ("formal_nearby_65_100", 65.0, 100.0),
    ("balanced_61_100", 61.0, 100.0),
    ("mass500_fixed_59_100", 59.0, 100.0),
    ("formal_fixed_68_100", 68.0, 100.0),
]


def enclosed_dehnen_mass_msun(radius_kpc: np.ndarray, mass_parameter: float, scale_kpc: float) -> np.ndarray:
    return mass_parameter * MASS_UNIT_MSUN * (radius_kpc / (radius_kpc + scale_kpc)) ** 3


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--processes", type=int, default=1)
    args = parser.parse_args()
    if args.processes != 1:
        raise ValueError("This single-snapshot analysis must use --processes 1")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with h5py.File(args.snapshot, "r") as handle:
        star = handle["PartType2"]
        star_mass = np.asarray(star["Masses"], dtype=np.float64)
        center = np.average(
            np.asarray(star["Coordinates"], dtype=np.float64),
            axis=0,
            weights=star_mass,
        )
        components: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for particle_type in range(6):
            key = f"PartType{particle_type}"
            if key not in handle:
                continue
            group = handle[key]
            coordinates = np.asarray(group["Coordinates"], dtype=np.float64)
            masses = np.asarray(group["Masses"], dtype=np.float64) * MASS_UNIT_MSUN
            radius = np.linalg.norm(coordinates - center, axis=1)
            components[key] = (radius, masses)

        gas = handle["PartType0"]
        gas_coordinates = np.asarray(gas["Coordinates"], dtype=np.float64)
        gas_masses_code = np.asarray(gas["Masses"], dtype=np.float64)
        gas_radius = np.linalg.norm(gas_coordinates - center, axis=1)
        mw_cgm_mask = np.isclose(
            gas_masses_code,
            MW_CGM_PARTICLE_MASS_CODE,
            rtol=2e-5,
            atol=0.0,
        )
        mw_cgm_radius = gas_radius[mw_cgm_mask]
        mw_cgm_mass = gas_masses_code[mw_cgm_mask] * MASS_UNIT_MSUN

    component_enclosed: dict[str, np.ndarray] = {}
    for name, (radius, masses) in components.items():
        order = np.argsort(radius)
        sorted_radius = radius[order]
        cumulative_mass = np.cumsum(masses[order])
        indices = np.searchsorted(sorted_radius, RADII_KPC, side="right") - 1
        enclosed = np.where(indices >= 0, cumulative_mass[np.maximum(indices, 0)], 0.0)
        component_enclosed[name] = enclosed

    cgm_order = np.argsort(mw_cgm_radius)
    sorted_cgm_radius = mw_cgm_radius[cgm_order]
    cumulative_cgm_mass = np.cumsum(mw_cgm_mass[cgm_order])
    cgm_indices = np.searchsorted(sorted_cgm_radius, RADII_KPC, side="right") - 1
    current_cgm_sampled = np.where(
        cgm_indices >= 0,
        cumulative_cgm_mass[np.maximum(cgm_indices, 0)],
        0.0,
    )
    current_total_mass = np.sum(np.vstack(list(component_enclosed.values())), axis=0)
    current_vc = np.sqrt(G_KPC_KMS2_MSUN * current_total_mass / RADII_KPC)

    rows: list[dict[str, float | str]] = []
    curves: dict[str, np.ndarray] = {}
    summaries: list[dict[str, float | str]] = []
    for name, mass_parameter, scale_kpc in CANDIDATES:
        candidate_cgm = enclosed_dehnen_mass_msun(RADII_KPC, mass_parameter, scale_kpc)
        candidate_total_mass = current_total_mass - current_cgm_sampled + candidate_cgm
        candidate_vc = np.sqrt(G_KPC_KMS2_MSUN * candidate_total_mass / RADII_KPC)
        curves[name] = candidate_vc
        fractional_change = candidate_vc / current_vc - 1.0

        for radius, total_mass, cgm_mass, vc, delta in zip(
            RADII_KPC,
            candidate_total_mass,
            candidate_cgm,
            candidate_vc,
            fractional_change,
        ):
            rows.append(
                {
                    "option": name,
                    "mass_parameter": mass_parameter,
                    "scale_kpc": scale_kpc,
                    "radius_kpc": radius,
                    "candidate_cgm_enclosed_msun": cgm_mass,
                    "candidate_total_enclosed_msun": total_mass,
                    "candidate_spherical_vc_km_s": vc,
                    "fractional_vc_change": delta,
                    "percent_vc_change": 100.0 * delta,
                }
            )

        def maximum_change(inner: float, outer: float) -> float:
            selected = (RADII_KPC >= inner) & (RADII_KPC <= outer)
            return 100.0 * float(np.max(np.abs(fractional_change[selected])))

        summaries.append(
            {
                "option": name,
                "mass_parameter": mass_parameter,
                "scale_kpc": scale_kpc,
                "max_abs_percent_vc_change_5_25kpc": maximum_change(5, 25),
                "max_abs_percent_vc_change_5_30kpc": maximum_change(5, 30),
                "max_abs_percent_vc_change_50_150kpc": maximum_change(50, 150),
                "max_abs_percent_vc_change_100_400kpc": maximum_change(100, 400),
                "percent_vc_change_140kpc": float(
                    100.0 * fractional_change[np.where(RADII_KPC == 140)[0][0]]
                ),
                "percent_vc_change_338kpc": float(
                    100.0 * fractional_change[np.where(RADII_KPC == 338)[0][0]]
                ),
            }
        )

    with (args.output_dir / "rotation_curve_candidates.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (args.output_dir / "rotation_curve_candidate_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)

    metadata = {
        "snapshot": str(args.snapshot),
        "center_kpc": center.tolist(),
        "method": (
            "Spherical enclosed-mass estimate. The CGM change is exact for a spherical "
            "Dehnen component; the absolute disk contribution is approximated spherically."
        ),
        "summaries": summaries,
    }
    (args.output_dir / "rotation_curve_candidate_summary.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )

    figure, (axis_curve, axis_delta) = plt.subplots(1, 2, figsize=(11.0, 4.2))
    axis_curve.plot(RADII_KPC, current_vc, color="black", lw=2.4, label="current sampled MW52")
    colors = plt.cm.tab10(np.linspace(0, 0.8, len(CANDIDATES) - 1))
    for color, (name, _, _) in zip(colors, CANDIDATES[1:]):
        axis_curve.plot(RADII_KPC, curves[name], color=color, lw=1.7, label=name)
        axis_delta.plot(
            RADII_KPC,
            100.0 * (curves[name] / current_vc - 1.0),
            color=color,
            lw=1.9,
            label=name,
        )
    axis_curve.set(
        xscale="log",
        xlabel=r"$r$ [kpc]",
        ylabel=r"spherical-equivalent $V_c$ [km s$^{-1}$]",
        xlim=(5, 500),
    )
    axis_delta.axhline(0.0, color="0.4", lw=1.0)
    axis_delta.axhspan(-1, 1, color="0.5", alpha=0.10)
    axis_delta.set(
        xscale="log",
        xlabel=r"$r$ [kpc]",
        ylabel=r"$\Delta V_c/V_c$ [per cent]",
        xlim=(5, 500),
    )
    axis_curve.legend(frameon=False, fontsize=7.5)
    axis_delta.legend(frameon=False, fontsize=7.5)
    for axis in (axis_curve, axis_delta):
        axis.tick_params(direction="in", which="both", top=True, right=True)
        axis.grid(alpha=0.16)
    figure.tight_layout()
    figure.savefig(args.output_dir / "rotation_curve_candidate_comparison.png", dpi=220)
    plt.close(figure)
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
