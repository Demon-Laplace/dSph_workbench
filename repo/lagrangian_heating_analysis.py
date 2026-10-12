#!/usr/bin/env python3
"""Track the same old-star populations in a formal and an isolation run.

The analysis is deliberately based on intrinsic, one-dimensional velocity
dispersion.  This removes viewing-angle effects while comparing external
heating with the matched isolation baseline.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


G = 4.30091e-6  # kpc (km/s)^2 Msun^-1
OLD_STAR_COUNT = 600_000
RADIAL_EDGES = np.array([0.0, 0.10, 0.20, 0.40, 0.70, 1.00, 1.50, 3.00, np.inf])
RADIAL_LABELS = [
    "0.00-0.10", "0.10-0.20", "0.20-0.40", "0.40-0.70",
    "0.70-1.00", "1.00-1.50", "1.50-3.00", ">3.00",
]
DENSITY_EDGES = np.array([0.0, 0.005, 0.010, 0.020, 0.040, 0.080,
                          0.120, 0.170, 0.250, 0.350, 0.500, 1.000, np.inf])
DENSITY_LABELS = [
    "<0.005", "0.005-0.010", "0.010-0.020", "0.020-0.040",
    "0.040-0.080", "0.080-0.120", "0.120-0.170", "0.170-0.250",
    "0.250-0.350", "0.350-0.500", "0.500-1.000", ">1.000",
]
CURRENT_RADIAL_EDGES = np.array([0.0, 0.10, 0.20, 0.40, 0.70, 1.00, 1.50, 3.00, np.inf])
DENSITY_PROFILE_EDGES = np.concatenate([
    np.arange(0.0, 2.0001, 0.025),
    np.arange(2.10, 5.0001, 0.10),
    np.arange(5.50, 20.0001, 0.50),
    np.array([30.0, 50.0, 100.0, np.inf]),
])


def weighted_mean(values, weights):
    return np.average(values, axis=0, weights=weights)


def weighted_std(values, weights):
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if values.size < 2 or weights.sum() <= 0:
        return float("nan")
    mean = np.average(values, weights=weights)
    return float(np.sqrt(np.average((values - mean) ** 2, weights=weights)))


def weighted_quantile(values, weights, quantile):
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(good):
        return float("nan")
    order = np.argsort(values[good])
    sorted_values = values[good][order]
    cumulative = np.cumsum(weights[good][order])
    return float(sorted_values[np.searchsorted(cumulative, quantile * cumulative[-1])])


def shrinking_center(pos, mass, radius_start=6.0, radius_stop=0.10, shrink=0.80):
    center = weighted_mean(pos, mass)
    radius = radius_start
    while radius > radius_stop:
        distance = np.linalg.norm(pos - center, axis=1)
        keep = distance < radius
        if keep.sum() < 500:
            break
        center = weighted_mean(pos[keep], mass[keep])
        radius *= shrink
    return center


def load_old_stars(handle):
    group = handle["PartType2"]
    ids = group["ParticleIDs"][:]
    if len(ids) < OLD_STAR_COUNT:
        raise RuntimeError(f"PartType2 contains only {len(ids)} particles")
    selected = np.sort(np.argpartition(ids, -OLD_STAR_COUNT)[-OLD_STAR_COUNT:])
    box_center = float(handle["Header"].attrs["BoxSize"]) / 2.0
    ids = ids[selected].astype(np.int64)
    order = np.argsort(ids)
    pos = group["Coordinates"][selected].astype(np.float64)[order] - box_center
    vel = group["Velocities"][selected].astype(np.float64)[order]
    mass = group["Masses"][selected].astype(np.float64)[order] * 1.0e10
    return ids[order], pos, vel, mass


def select_initial_gas_ids(handle, center, radius=10.0, chunk=500_000):
    group = handle["PartType0"]
    box_center = float(handle["Header"].attrs["BoxSize"]) / 2.0
    selected = []
    for start in range(0, len(group["Coordinates"]), chunk):
        stop = min(start + chunk, len(group["Coordinates"]))
        pos = group["Coordinates"][start:stop].astype(np.float64) - box_center
        keep = np.linalg.norm(pos - center, axis=1) < radius
        if np.any(keep):
            selected.append(group["ParticleIDs"][start:stop][keep].astype(np.int64))
    if not selected:
        return np.empty(0, dtype=np.int64)
    return np.unique(np.concatenate(selected))


def load_tracked_gas(handle, tracked_ids, chunk=500_000):
    if len(tracked_ids) == 0 or "PartType0" not in handle:
        return np.empty((0, 3)), np.empty(0)
    group = handle["PartType0"]
    box_center = float(handle["Header"].attrs["BoxSize"]) / 2.0
    pos_parts, mass_parts = [], []
    for start in range(0, len(group["ParticleIDs"]), chunk):
        stop = min(start + chunk, len(group["ParticleIDs"]))
        ids = group["ParticleIDs"][start:stop].astype(np.int64)
        keep = np.isin(ids, tracked_ids, assume_unique=False)
        if not np.any(keep):
            continue
        pos_parts.append(group["Coordinates"][start:stop][keep].astype(np.float64) - box_center)
        mass_parts.append(group["Masses"][start:stop][keep].astype(np.float64) * 1.0e10)
    return (
        np.concatenate(pos_parts) if pos_parts else np.empty((0, 3)),
        np.concatenate(mass_parts) if mass_parts else np.empty(0),
    )


def load_new_stars_near_center(handle, center, radius=20.0, chunk=500_000):
    if "PartType4" not in handle:
        return np.empty((0, 3)), np.empty(0)
    group = handle["PartType4"]
    box_center = float(handle["Header"].attrs["BoxSize"]) / 2.0
    pos_parts, mass_parts = [], []
    for start in range(0, len(group["Coordinates"]), chunk):
        stop = min(start + chunk, len(group["Coordinates"]))
        pos = group["Coordinates"][start:stop].astype(np.float64) - box_center
        keep = np.linalg.norm(pos - center, axis=1) < radius
        if not np.any(keep):
            continue
        pos_parts.append(pos[keep])
        mass_parts.append(group["Masses"][start:stop][keep].astype(np.float64) * 1.0e10)
    return (
        np.concatenate(pos_parts) if pos_parts else np.empty((0, 3)),
        np.concatenate(mass_parts) if mass_parts else np.empty(0),
    )


def density_profile(source_radius, source_mass, edges=DENSITY_PROFILE_EDGES):
    finite_edges = edges.copy()
    finite_edges[-1] = max(1000.0, float(np.nanmax(source_radius)) + 1.0)
    mass, _ = np.histogram(source_radius, bins=finite_edges, weights=source_mass)
    volume = 4.0 * math.pi / 3.0 * (finite_edges[1:] ** 3 - finite_edges[:-1] ** 3)
    return mass / volume, finite_edges


def sample_shell_density(radius, density, edges):
    index = np.searchsorted(edges, radius, side="right") - 1
    index = np.clip(index, 0, len(density) - 1)
    return density[index]


def spherical_potential(target_radius, source_radius, source_mass):
    order = np.argsort(source_radius)
    radius = np.maximum(source_radius[order], 1.0e-4)
    mass = source_mass[order]
    cumulative_mass = np.cumsum(mass)
    outer_m_over_r = np.cumsum((mass / radius)[::-1])[::-1]
    index = np.searchsorted(radius, target_radius, side="right")
    inner = np.zeros_like(target_radius, dtype=float)
    has_inner = index > 0
    inner[has_inner] = cumulative_mass[index[has_inner] - 1] / np.maximum(target_radius[has_inner], 1.0e-4)
    outer = np.zeros_like(target_radius, dtype=float)
    has_outer = index < len(radius)
    outer[has_outer] = outer_m_over_r[index[has_outer]]
    return -G * (inner + outer)


def spherical_velocity_components(relative_pos, relative_vel):
    radius = np.linalg.norm(relative_pos, axis=1)
    rhat = relative_pos / np.maximum(radius[:, None], 1.0e-12)
    cylindrical = np.hypot(relative_pos[:, 0], relative_pos[:, 1])
    phihat = np.column_stack([
        -relative_pos[:, 1] / np.maximum(cylindrical, 1.0e-12),
        relative_pos[:, 0] / np.maximum(cylindrical, 1.0e-12),
        np.zeros(len(radius)),
    ])
    polar = cylindrical < 1.0e-10
    phihat[polar] = np.array([0.0, 1.0, 0.0])
    thetahat = np.cross(phihat, rhat)
    return relative_vel @ rhat.T if False else (
        np.einsum("ij,ij->i", relative_vel, rhat),
        np.einsum("ij,ij->i", relative_vel, thetahat),
        np.einsum("ij,ij->i", relative_vel, phihat),
    )


def group_metrics(mask, mass, radius, velocity, vr, vtheta, vphi, energy,
                  delta_energy, angular_momentum, delta_angular_momentum,
                  initial_radius, initial_density, current_density):
    count = int(np.count_nonzero(mask))
    if count < 10:
        return {"count": count}
    weights = mass[mask]
    values_r = radius[mask]
    sigma_xyz = [weighted_std(velocity[mask, axis], weights) for axis in range(3)]
    sigma_1d = math.sqrt(sum(value * value for value in sigma_xyz) / 3.0)
    sigma_r = weighted_std(vr[mask], weights)
    sigma_theta = weighted_std(vtheta[mask], weights)
    sigma_phi = weighted_std(vphi[mask], weights)
    beta = 1.0 - (sigma_theta**2 + sigma_phi**2) / (2.0 * sigma_r**2) if sigma_r > 0 else float("nan")
    return {
        "count": count,
        "mass_msun": float(weights.sum()),
        "initial_radius_median_kpc": weighted_quantile(initial_radius[mask], weights, 0.5),
        "initial_density_median_msun_pc3": weighted_quantile(initial_density[mask], weights, 0.5),
        "current_density_median_msun_pc3": weighted_quantile(current_density[mask], weights, 0.5),
        "current_radius_p16_kpc": weighted_quantile(values_r, weights, 0.16),
        "current_radius_median_kpc": weighted_quantile(values_r, weights, 0.50),
        "current_radius_p84_kpc": weighted_quantile(values_r, weights, 0.84),
        "fraction_within_0p2": float(weights[values_r < 0.2].sum() / weights.sum()),
        "fraction_within_0p4": float(weights[values_r < 0.4].sum() / weights.sum()),
        "fraction_within_0p7": float(weights[values_r < 0.7].sum() / weights.sum()),
        "fraction_beyond_3": float(weights[values_r >= 3.0].sum() / weights.sum()),
        "sigma_1d_km_s": sigma_1d,
        "sigma_x_km_s": sigma_xyz[0],
        "sigma_y_km_s": sigma_xyz[1],
        "sigma_z_km_s": sigma_xyz[2],
        "sigma_radial_km_s": sigma_r,
        "sigma_theta_km_s": sigma_theta,
        "sigma_phi_km_s": sigma_phi,
        "beta": beta,
        "mean_radial_velocity_km_s": float(np.average(vr[mask], weights=weights)),
        "energy_median_km2_s2": weighted_quantile(energy[mask], weights, 0.5),
        "delta_energy_median_km2_s2": weighted_quantile(delta_energy[mask], weights, 0.5),
        "unbound_proxy_fraction": float(weights[energy[mask] > 0].sum() / weights.sum()),
        "angular_momentum_median_kpc_km_s": weighted_quantile(angular_momentum[mask], weights, 0.5),
        "delta_angular_momentum_median_kpc_km_s": weighted_quantile(delta_angular_momentum[mask], weights, 0.5),
    }


def snapshot_path(model_dir, number):
    return Path(model_dir) / "output" / f"snapshot_{number:03d}.hdf5"


def prepare_reference(model_dir, output_path, sequence, reference_number):
    path = snapshot_path(model_dir, reference_number)
    with h5py.File(path, "r") as handle:
        time_gyr = float(handle["Header"].attrs["Time"])
        ids, pos, vel, mass = load_old_stars(handle)
        center = shrinking_center(pos, mass)
        r = np.linalg.norm(pos - center, axis=1)
        core = r < 0.5
        center_velocity = weighted_mean(vel[core], mass[core])
        tracked_gas_ids = select_initial_gas_ids(handle, center)
        gas_pos, gas_mass = load_tracked_gas(handle, tracked_gas_ids)
        new_pos, new_mass = load_new_stars_near_center(handle, center)

    relative_pos = pos - center
    relative_vel = vel - center_velocity
    gas_r = np.linalg.norm(gas_pos - center, axis=1) if len(gas_mass) else np.empty(0)
    new_r = np.linalg.norm(new_pos - center, axis=1) if len(new_mass) else np.empty(0)
    source_r = np.concatenate([r, gas_r, new_r])
    source_m = np.concatenate([mass, gas_mass, new_mass])
    baryon_density, profile_edges = density_profile(source_r, source_m)
    old_density, _ = density_profile(r, mass)
    initial_density = sample_shell_density(r, baryon_density, profile_edges) / 1.0e9
    initial_old_density = sample_shell_density(r, old_density, profile_edges) / 1.0e9
    potential = spherical_potential(r, source_r, source_m)
    energy = 0.5 * np.sum(relative_vel**2, axis=1) + potential
    angular_momentum = np.linalg.norm(np.cross(relative_pos, relative_vel), axis=1)
    radial_group = np.searchsorted(RADIAL_EDGES, r, side="right") - 1
    radial_group = np.clip(radial_group, 0, len(RADIAL_LABELS) - 1)
    density_group = np.searchsorted(DENSITY_EDGES, initial_density, side="right") - 1
    density_group = np.clip(density_group, 0, len(DENSITY_LABELS) - 1)
    np.savez_compressed(
        output_path,
        ids=ids,
        mass=mass,
        initial_radius=r,
        initial_density=initial_density,
        initial_old_density=initial_old_density,
        initial_energy=energy,
        initial_angular_momentum=angular_momentum,
        radial_group=radial_group,
        density_group=density_group,
        tracked_gas_ids=tracked_gas_ids,
        center=center,
        center_velocity=center_velocity,
        time_gyr=np.array([time_gyr]),
    )
    return {
        "sequence": sequence,
        "snapshot": 0,
        "source_snapshot": reference_number,
        "time_gyr": time_gyr,
        "center_kpc": center.tolist(),
        "old_star_count": int(len(ids)),
        "old_star_mass_msun": float(mass.sum()),
        "tracked_initial_gas_count": int(len(tracked_gas_ids)),
        "tracked_initial_gas_mass_msun": float(gas_mass.sum()),
    }


def process_snapshot(task):
    sequence, model_dir, number, source_number, reference_path = task
    reference = np.load(reference_path)
    ref_ids = reference["ids"]
    tracked_gas_ids = reference["tracked_gas_ids"]
    reference_time = float(reference["time_gyr"][0])
    path = snapshot_path(model_dir, source_number)
    with h5py.File(path, "r") as handle:
        time_gyr = float(handle["Header"].attrs["Time"])
        ids, pos, vel, mass = load_old_stars(handle)
        center = shrinking_center(pos, mass)
        radius = np.linalg.norm(pos - center, axis=1)
        core = radius < 0.5
        center_velocity = weighted_mean(vel[core], mass[core])
        gas_pos, gas_mass = load_tracked_gas(handle, tracked_gas_ids)
        new_pos, new_mass = load_new_stars_near_center(handle, center)

    ref_index = np.searchsorted(ref_ids, ids)
    valid = (ref_index >= 0) & (ref_index < len(ref_ids)) & (ref_ids[np.clip(ref_index, 0, len(ref_ids)-1)] == ids)
    if not np.all(valid):
        ids, pos, vel, mass, radius, ref_index = (
            ids[valid], pos[valid], vel[valid], mass[valid], radius[valid], ref_index[valid]
        )
    initial_radius = reference["initial_radius"][ref_index]
    initial_density = reference["initial_density"][ref_index]
    radial_group = reference["radial_group"][ref_index]
    density_group = reference["density_group"][ref_index]
    initial_energy = reference["initial_energy"][ref_index]
    initial_angular_momentum = reference["initial_angular_momentum"][ref_index]

    relative_pos = pos - center
    relative_vel = vel - center_velocity
    gas_r = np.linalg.norm(gas_pos - center, axis=1) if len(gas_mass) else np.empty(0)
    new_r = np.linalg.norm(new_pos - center, axis=1) if len(new_mass) else np.empty(0)
    source_r = np.concatenate([radius, gas_r, new_r])
    source_m = np.concatenate([mass, gas_mass, new_mass])
    baryon_density, profile_edges = density_profile(source_r, source_m)
    current_density = sample_shell_density(radius, baryon_density, profile_edges) / 1.0e9
    potential = spherical_potential(radius, source_r, source_m)
    energy = 0.5 * np.sum(relative_vel**2, axis=1) + potential
    delta_energy = energy - initial_energy
    angular_momentum = np.linalg.norm(np.cross(relative_pos, relative_vel), axis=1)
    delta_angular_momentum = angular_momentum - initial_angular_momentum
    vr, vtheta, vphi = spherical_velocity_components(relative_pos, relative_vel)

    radial_rows = []
    for group_index, label in enumerate(RADIAL_LABELS):
        row = group_metrics(
            radial_group == group_index, mass, radius, relative_vel, vr, vtheta, vphi,
            energy, delta_energy, angular_momentum, delta_angular_momentum,
            initial_radius, initial_density, current_density,
        )
        row.update({"sequence": sequence, "snapshot": number, "source_snapshot": source_number,
                    "time_gyr": time_gyr - reference_time, "source_time_gyr": time_gyr,
                    "group_kind": "initial_radius", "group": label, "group_index": group_index})
        radial_rows.append(row)

    density_rows = []
    for group_index, label in enumerate(DENSITY_LABELS):
        row = group_metrics(
            density_group == group_index, mass, radius, relative_vel, vr, vtheta, vphi,
            energy, delta_energy, angular_momentum, delta_angular_momentum,
            initial_radius, initial_density, current_density,
        )
        row.update({"sequence": sequence, "snapshot": number, "source_snapshot": source_number,
                    "time_gyr": time_gyr - reference_time, "source_time_gyr": time_gyr,
                    "group_kind": "initial_density", "group": label, "group_index": group_index})
        density_rows.append(row)

    eulerian_group = np.searchsorted(CURRENT_RADIAL_EDGES, radius, side="right") - 1
    eulerian_group = np.clip(eulerian_group, 0, len(RADIAL_LABELS) - 1)
    eulerian_rows = []
    transition = np.zeros((len(RADIAL_LABELS), len(RADIAL_LABELS)), dtype=float)
    for group_index, label in enumerate(RADIAL_LABELS):
        mask = eulerian_group == group_index
        row = group_metrics(
            mask, mass, radius, relative_vel, vr, vtheta, vphi,
            energy, delta_energy, angular_momentum, delta_angular_momentum,
            initial_radius, initial_density, current_density,
        )
        row.update({"sequence": sequence, "snapshot": number, "source_snapshot": source_number,
                    "time_gyr": time_gyr - reference_time, "source_time_gyr": time_gyr,
                    "group_kind": "current_radius", "group": label, "group_index": group_index})
        if np.any(mask):
            weights = mass[mask]
            row["origin_fraction_initial_lt_0p2"] = float(weights[initial_radius[mask] < 0.2].sum() / weights.sum())
            row["origin_fraction_initial_0p2_0p7"] = float(weights[(initial_radius[mask] >= 0.2) & (initial_radius[mask] < 0.7)].sum() / weights.sum())
        eulerian_rows.append(row)
    for initial_index in range(len(RADIAL_LABELS)):
        mask_initial = radial_group == initial_index
        denominator = mass[mask_initial].sum()
        if denominator <= 0:
            continue
        for current_index in range(len(RADIAL_LABELS)):
            transition[initial_index, current_index] = mass[mask_initial & (eulerian_group == current_index)].sum() / denominator

    summary = {
        "sequence": sequence,
        "snapshot": number,
        "source_snapshot": source_number,
        "time_gyr": time_gyr - reference_time,
        "source_time_gyr": time_gyr,
        "path": str(path),
        "center_kpc": center.tolist(),
        "center_velocity_km_s": center_velocity.tolist(),
        "matched_old_star_count": int(len(ids)),
        "tracked_gas_count": int(len(gas_mass)),
        "tracked_gas_mass_msun": float(gas_mass.sum()),
        "new_star_mass_within_20kpc_msun": float(new_mass.sum()),
        "old_star_mass_msun": float(mass.sum()),
    }
    return {
        "summary": summary,
        "radial_rows": radial_rows,
        "density_rows": density_rows,
        "eulerian_rows": eulerian_rows,
        "transition": transition,
    }


def write_csv(path, rows):
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def comparison_rows(formal_rows, isolation_rows):
    isolation_lookup = {(row["snapshot"], row["group_index"]): row for row in isolation_rows}
    rows = []
    for formal in formal_rows:
        isolation = isolation_lookup.get((formal["snapshot"], formal["group_index"]))
        if isolation is None or formal.get("count", 0) < 10 or isolation.get("count", 0) < 10:
            continue
        sf = float(formal["sigma_1d_km_s"])
        si = float(isolation["sigma_1d_km_s"])
        delta2 = sf**2 - si**2
        rows.append({
            "snapshot": formal["snapshot"],
            "time_gyr": formal["time_gyr"],
            "group_kind": formal["group_kind"],
            "group": formal["group"],
            "group_index": formal["group_index"],
            "count_formal": formal["count"],
            "count_isolation": isolation["count"],
            "initial_density_formal_median_msun_pc3": formal.get("initial_density_median_msun_pc3"),
            "initial_density_isolation_median_msun_pc3": isolation.get("initial_density_median_msun_pc3"),
            "formal_sigma_1d_km_s": sf,
            "isolation_sigma_1d_km_s": si,
            "external_delta_sigma2_km2_s2": delta2,
            "external_equivalent_sigma_km_s": math.sqrt(max(delta2, 0.0)),
            "external_fraction_of_isolation_sigma2": delta2 / si**2 if si > 0 else float("nan"),
            "formal_delta_energy_median_km2_s2": formal.get("delta_energy_median_km2_s2"),
            "isolation_delta_energy_median_km2_s2": isolation.get("delta_energy_median_km2_s2"),
            "external_delta_energy_median_km2_s2": (
                float(formal.get("delta_energy_median_km2_s2", np.nan))
                - float(isolation.get("delta_energy_median_km2_s2", np.nan))
            ),
            "formal_current_radius_median_kpc": formal.get("current_radius_median_kpc"),
            "isolation_current_radius_median_kpc": isolation.get("current_radius_median_kpc"),
            "formal_beta": formal.get("beta"),
            "isolation_beta": isolation.get("beta"),
        })
    return rows


def add_initial_offset_correction(rows):
    baseline = {
        (row["group_kind"], row["group_index"]): row["external_delta_sigma2_km2_s2"]
        for row in rows if row["snapshot"] == 0
    }
    for row in rows:
        offset = baseline.get((row["group_kind"], row["group_index"]), float("nan"))
        row["initial_variance_offset_km2_s2"] = offset
        row["change_corrected_external_delta_sigma2_km2_s2"] = (
            row["external_delta_sigma2_km2_s2"] - offset
        )
        row["change_corrected_external_equivalent_sigma_km_s"] = math.copysign(
            math.sqrt(abs(row["change_corrected_external_delta_sigma2_km2_s2"])),
            row["change_corrected_external_delta_sigma2_km2_s2"],
        )


def read_orbit_csv(path):
    if not path or not Path(path).exists():
        return []
    with Path(path).open(encoding="utf-8") as handle:
        lines = [line for line in handle if not line.startswith("#")]
    return list(csv.DictReader(lines))


def make_evolution_figure(radial_rows, radial_comparison, orbit_rows, output):
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.3), constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, 6))
    chosen = list(range(6))
    for color, index in zip(colors, chosen):
        for sequence, style in (("formal", "-"), ("isolation", "--")):
            rows = sorted(
                [row for row in radial_rows if row["sequence"] == sequence and row["group_index"] == index and row.get("count", 0) >= 10],
                key=lambda row: row["time_gyr"],
            )
            axes[0, 0].plot([row["time_gyr"] for row in rows], [row["sigma_1d_km_s"] for row in rows],
                            color=color, ls=style, marker="o", ms=3, lw=1.4,
                            label=f"{RADIAL_LABELS[index]} formal" if sequence == "formal" else None)
            axes[1, 0].plot([row["time_gyr"] for row in rows], [row["current_radius_median_kpc"] for row in rows],
                            color=color, ls=style, marker="o", ms=3, lw=1.3)
            axes[1, 1].plot([row["time_gyr"] for row in rows], [row["delta_energy_median_km2_s2"] for row in rows],
                            color=color, ls=style, marker="o", ms=3, lw=1.3)
        compare = sorted([row for row in radial_comparison if row["group_index"] == index], key=lambda row: row["time_gyr"])
        axes[0, 1].plot([row["time_gyr"] for row in compare], [row["change_corrected_external_delta_sigma2_km2_s2"] for row in compare],
                        color=color, marker="o", ms=3, lw=1.4, label=RADIAL_LABELS[index])

    if orbit_rows:
        peri = min(orbit_rows, key=lambda row: float(row["distance_gal"]))
        peri_time = float(peri["age"])
        for axis in axes.flat:
            axis.axvline(peri_time, color="0.25", ls=":", lw=1.0)
    axes[0, 0].set(ylabel=r"Intrinsic $\sigma_{1D}$ [km s$^{-1}$]", title="Same initial-radius populations")
    axes[0, 1].axhline(0, color="0.4", lw=0.8)
    axes[0, 1].set(ylabel=r"$\sigma_{formal}^2-\sigma_{iso}^2$ [km$^2$ s$^{-2}$]", title="Net external response")
    axes[1, 0].set(xlabel="Time [Gyr]", ylabel="Median current radius [kpc]", title="Lagrangian migration")
    axes[1, 1].axhline(0, color="0.4", lw=0.8)
    axes[1, 1].set(xlabel="Time [Gyr]", ylabel=r"Median $\Delta E$ [km$^2$ s$^{-2}$]", title="Internal-energy proxy")
    axes[0, 0].text(0.02, 0.03, "solid: formal; dashed: isolation", transform=axes[0, 0].transAxes, fontsize=9)
    for axis in axes.flat:
        axis.grid(alpha=0.2)
    axes[0, 0].legend(frameon=False, fontsize=7, ncol=2)
    axes[0, 1].legend(frameon=False, fontsize=8, ncol=2)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def make_density_figure(density_rows, density_comparison, target_snapshot, output):
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.2), constrained_layout=True)
    target_compare = sorted(
        [row for row in density_comparison if row["snapshot"] == target_snapshot and row["count_formal"] >= 100],
        key=lambda row: row["initial_density_formal_median_msun_pc3"],
    )
    density = np.array([row["initial_density_formal_median_msun_pc3"] for row in target_compare], dtype=float)
    formal_sigma = np.array([row["formal_sigma_1d_km_s"] for row in target_compare], dtype=float)
    iso_sigma = np.array([row["isolation_sigma_1d_km_s"] for row in target_compare], dtype=float)
    delta2 = np.array([row["change_corrected_external_delta_sigma2_km2_s2"] for row in target_compare], dtype=float)
    axes[0, 0].plot(density, iso_sigma, "o-", label="Isolation baseline")
    axes[0, 0].plot(density, formal_sigma, "o-", label="Formal total")
    axes[0, 0].set(ylabel=r"$\sigma_{1D}$ [km s$^{-1}$]", title=f"Snapshot {target_snapshot:03d}: baseline and total")
    axes[0, 1].axhline(0, color="0.4", lw=0.8)
    axes[0, 1].plot(density, delta2, "o-", color="#d62728")
    axes[0, 1].set(ylabel=r"External $\Delta\sigma^2$ [km$^2$ s$^{-2}$]", title="Net heating or cooling")
    axes[1, 0].plot(density, np.maximum(delta2, 0), "o-", color="#2ca02c", label="External contribution")
    axes[1, 0].plot(density, iso_sigma**2, "o-", color="#1f77b4", label="Isolation contribution")
    axes[1, 0].plot(density, formal_sigma**2, "o-", color="k", label="Total")
    axes[1, 0].set(xlabel=r"Initial baryonic density [$M_\odot$ pc$^{-3}$]", ylabel=r"Variance [km$^2$ s$^{-2}$]", title="Quadrature decomposition")

    snapshots = sorted({row["snapshot"] for row in density_comparison})
    group_indices = sorted({row["group_index"] for row in density_comparison})
    matrix = np.full((len(group_indices), len(snapshots)), np.nan)
    density_axis = np.full(len(group_indices), np.nan)
    lookup_group = {value: index for index, value in enumerate(group_indices)}
    lookup_snapshot = {value: index for index, value in enumerate(snapshots)}
    for row in density_comparison:
        i = lookup_group[row["group_index"]]
        j = lookup_snapshot[row["snapshot"]]
        matrix[i, j] = row["change_corrected_external_delta_sigma2_km2_s2"]
        density_axis[i] = row["initial_density_formal_median_msun_pc3"]
    image = axes[1, 1].imshow(matrix, origin="lower", aspect="auto", cmap="coolwarm",
                              extent=[min(snapshots) / 100, max(snapshots) / 100, -0.5, len(group_indices)-0.5])
    axes[1, 1].set(xlabel="Time [Gyr]", ylabel="Initial density bin", title=r"External $\Delta\sigma^2$ evolution")
    axes[1, 1].set_yticks(range(len(group_indices)))
    axes[1, 1].set_yticklabels([f"{density_axis[i]:.3g}" for i in range(len(group_indices))], fontsize=8)
    fig.colorbar(image, ax=axes[1, 1], label=r"km$^2$ s$^{-2}$")
    for axis in (axes[0, 0], axes[0, 1], axes[1, 0]):
        axis.set_xscale("log")
        axis.set_xlabel(r"Initial baryonic density [$M_\odot$ pc$^{-3}$]")
        axis.grid(alpha=0.2)
    axes[0, 0].legend(frameon=False)
    axes[1, 0].legend(frameon=False)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def make_transition_figure(transitions, eulerian_rows, target_snapshot, output):
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.2), constrained_layout=True)
    for axis, sequence in zip(axes[:2], ("formal", "isolation")):
        matrix = transitions[(sequence, target_snapshot)]
        image = axis.imshow(matrix, origin="lower", vmin=0, vmax=max(0.5, np.nanmax(matrix)), cmap="magma")
        axis.set(title=sequence.capitalize(), xlabel="Current radial bin", ylabel="Initial radial bin")
        axis.set_xticks(range(len(RADIAL_LABELS)))
        axis.set_yticks(range(len(RADIAL_LABELS)))
        axis.set_xticklabels(RADIAL_LABELS, rotation=60, ha="right", fontsize=7)
        axis.set_yticklabels(RADIAL_LABELS, fontsize=7)
        fig.colorbar(image, ax=axis, label="Mass fraction")
    for sequence, style, color in (("formal", "-", "#1f77b4"), ("isolation", "--", "#d62728")):
        rows = sorted([row for row in eulerian_rows if row["sequence"] == sequence and row["snapshot"] == target_snapshot and row.get("count", 0) >= 10], key=lambda row: row["group_index"])
        axes[2].plot([row["group_index"] for row in rows], [row["sigma_1d_km_s"] for row in rows], style, marker="o", color=color, label=sequence)
    axes[2].set(title="Current intrinsic profile", xlabel="Current radial bin", ylabel=r"$\sigma_{1D}$ [km s$^{-1}$]")
    axes[2].set_xticks(range(len(RADIAL_LABELS)))
    axes[2].set_xticklabels(RADIAL_LABELS, rotation=60, ha="right", fontsize=7)
    axes[2].grid(alpha=0.2)
    axes[2].legend(frameon=False)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--formal-dir", required=True)
    parser.add_argument("--isolation-dir", required=True)
    parser.add_argument("--snapshots", default="0,1,7,25,50,100,150,200,225,235,246,250,270,300")
    parser.add_argument("--workers", type=int, default=max(1, min(12, int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))))
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--orbit-csv")
    parser.add_argument("--target-snapshot", type=int, default=235)
    parser.add_argument("--isolation-offset", type=int, default=0,
                        help="Isolation source snapshot corresponding to formal snapshot 000")
    args = parser.parse_args()

    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    snapshots = sorted({int(value) for value in args.snapshots.split(",")})
    for number in snapshots:
        if not snapshot_path(args.formal_dir, number).exists():
            raise FileNotFoundError(snapshot_path(args.formal_dir, number))
        isolation_number = number + args.isolation_offset
        if not snapshot_path(args.isolation_dir, isolation_number).exists():
            raise FileNotFoundError(snapshot_path(args.isolation_dir, isolation_number))

    reference_paths = {
        "formal": output / "reference_formal.npz",
        "isolation": output / "reference_isolation.npz",
    }
    reference_summary = [
        prepare_reference(args.formal_dir, reference_paths["formal"], "formal", 0),
        prepare_reference(args.isolation_dir, reference_paths["isolation"], "isolation", args.isolation_offset),
    ]
    tasks = []
    for number in snapshots:
        tasks.append(("formal", args.formal_dir, number, number, str(reference_paths["formal"])))
        tasks.append(("isolation", args.isolation_dir, number, number + args.isolation_offset,
                      str(reference_paths["isolation"])))

    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        future_map = {pool.submit(process_snapshot, task): task for task in tasks}
        for future in as_completed(future_map):
            task = future_map[future]
            result = future.result()
            results.append(result)
            print(f"completed {task[0]} analysis {task[2]:03d} source {task[3]:03d}", flush=True)

    results.sort(key=lambda result: (result["summary"]["sequence"], result["summary"]["snapshot"]))
    summaries = [result["summary"] for result in results]
    radial_rows = [row for result in results for row in result["radial_rows"]]
    density_rows = [row for result in results for row in result["density_rows"]]
    eulerian_rows = [row for result in results for row in result["eulerian_rows"]]
    transitions = {(result["summary"]["sequence"], result["summary"]["snapshot"]): result["transition"] for result in results}

    radial_comparison = comparison_rows(
        [row for row in radial_rows if row["sequence"] == "formal"],
        [row for row in radial_rows if row["sequence"] == "isolation"],
    )
    density_comparison = comparison_rows(
        [row for row in density_rows if row["sequence"] == "formal"],
        [row for row in density_rows if row["sequence"] == "isolation"],
    )
    add_initial_offset_correction(radial_comparison)
    add_initial_offset_correction(density_comparison)
    write_csv(output / "lagrangian_radial_groups.csv", radial_rows)
    write_csv(output / "initial_density_groups.csv", density_rows)
    write_csv(output / "eulerian_radial_profiles.csv", eulerian_rows)
    write_csv(output / "formal_isolation_radial_comparison.csv", radial_comparison)
    write_csv(output / "formal_isolation_density_comparison.csv", density_comparison)
    np.savez_compressed(output / "radial_transition_matrices.npz", **{
        f"{sequence}_{number:03d}": matrix for (sequence, number), matrix in transitions.items()
    })
    metadata = {
        "formal_dir": args.formal_dir,
        "isolation_dir": args.isolation_dir,
        "snapshots": snapshots,
        "workers": args.workers,
        "target_snapshot": args.target_snapshot,
        "isolation_offset": args.isolation_offset,
        "radial_edges_kpc": RADIAL_EDGES.tolist(),
        "density_edges_msun_pc3": DENSITY_EDGES.tolist(),
        "reference_summary": reference_summary,
        "snapshot_summary": summaries,
    }
    (output / "lagrangian_heating_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    orbit_rows = read_orbit_csv(args.orbit_csv)
    make_evolution_figure(radial_rows, radial_comparison, orbit_rows, output / "lagrangian_heating_evolution.png")
    make_density_figure(density_rows, density_comparison, args.target_snapshot, output / "density_heating_response.png")
    make_transition_figure(transitions, eulerian_rows, args.target_snapshot, output / "radial_origin_mixing.png")
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == "__main__":
    main()
