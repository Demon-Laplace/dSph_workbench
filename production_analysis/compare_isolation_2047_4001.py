#!/usr/bin/env python3
"""Compare intrinsic structure and kinematics of two isolated Fornax ICs."""

import argparse
import csv
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


OLD_STAR_COUNT = 600_000
RADIAL_EDGES = np.array([0.0, 0.05, 0.10, 0.20, 0.40, 0.70, 1.00, 1.50, 3.00, 5.00])
RADIAL_LABELS = [f"{left:.2f}-{right:.2f}" for left, right in zip(RADIAL_EDGES[:-1], RADIAL_EDGES[1:])]


def weighted_mean(values, weights):
    return np.average(values, axis=0, weights=weights)


def weighted_std(values, weights):
    mean = np.average(values, weights=weights)
    return float(np.sqrt(np.average((values - mean) ** 2, weights=weights)))


def weighted_quantile(values, weights, quantile):
    order = np.argsort(values)
    sorted_values = values[order]
    cumulative = np.cumsum(weights[order])
    if cumulative[-1] <= 0:
        return float("nan")
    return float(np.interp(quantile * cumulative[-1], cumulative, sorted_values))


def shrinking_center(pos, mass, radius_start=6.0, radius_stop=0.10, shrink=0.8):
    center = weighted_mean(pos, mass)
    radius = radius_start
    while radius > radius_stop:
        keep = np.linalg.norm(pos - center, axis=1) < radius
        if np.count_nonzero(keep) < 500:
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
    pos = group["Coordinates"][selected].astype(np.float64) - box_center
    vel = group["Velocities"][selected].astype(np.float64)
    mass = group["Masses"][selected].astype(np.float64) * 1.0e10
    return pos, vel, mass


def load_nearby(handle, particle_type, center, radius=10.0, chunk=500_000):
    if particle_type not in handle:
        return {"pos": np.empty((0, 3)), "vel": np.empty((0, 3)), "mass": np.empty(0), "sfr": np.empty(0), "density": np.empty(0), "internal_energy": np.empty(0)}
    group = handle[particle_type]
    box_center = float(handle["Header"].attrs["BoxSize"]) / 2.0
    pos_parts, vel_parts, mass_parts = [], [], []
    extras = {"sfr": [], "density": [], "internal_energy": []}
    datasets = {"sfr": "StarFormationRate", "density": "Density", "internal_energy": "InternalEnergy"}
    for start in range(0, len(group["Coordinates"]), chunk):
        stop = min(start + chunk, len(group["Coordinates"]))
        pos = group["Coordinates"][start:stop].astype(np.float64) - box_center
        keep = np.linalg.norm(pos - center, axis=1) < radius
        if np.any(keep):
            pos_parts.append(pos[keep])
            vel_parts.append(group["Velocities"][start:stop][keep].astype(np.float64))
            mass_parts.append(group["Masses"][start:stop][keep].astype(np.float64) * 1.0e10)
            for key, dataset in datasets.items():
                if dataset in group:
                    extras[key].append(group[dataset][start:stop][keep].astype(np.float64))
                else:
                    extras[key].append(np.full(np.count_nonzero(keep), np.nan))
    return {
        "pos": np.concatenate(pos_parts) if pos_parts else np.empty((0, 3)),
        "vel": np.concatenate(vel_parts) if vel_parts else np.empty((0, 3)),
        "mass": np.concatenate(mass_parts) if mass_parts else np.empty(0),
        **{key: np.concatenate(parts) if parts else np.empty(0) for key, parts in extras.items()},
    }


def intrinsic_axis(relative_pos, relative_vel, mass, radius):
    keep = radius < 2.0
    covariance = np.einsum("ni,nj,n->ij", relative_pos[keep], relative_pos[keep], mass[keep]) / mass[keep].sum()
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    z_axis = eigenvectors[:, np.argmin(eigenvalues)]
    angular_momentum = np.sum(np.cross(relative_pos[keep], relative_vel[keep]) * mass[keep, None], axis=0)
    if np.dot(z_axis, angular_momentum) < 0:
        z_axis *= -1.0
    return z_axis


def cylindrical_components(relative_pos, relative_vel, z_axis):
    z = relative_pos @ z_axis
    plane = relative_pos - z[:, None] * z_axis
    cylindrical_radius = np.linalg.norm(plane, axis=1)
    e_radius = plane / np.maximum(cylindrical_radius[:, None], 1.0e-12)
    e_phi = np.cross(z_axis[None, :], e_radius)
    v_radius = np.einsum("ij,ij->i", relative_vel, e_radius)
    v_phi = np.einsum("ij,ij->i", relative_vel, e_phi)
    v_z = relative_vel @ z_axis
    return cylindrical_radius, z, v_radius, v_phi, v_z


def shape_ratios(relative_pos, mass, radius, aperture):
    keep = radius < aperture
    covariance = np.einsum("ni,nj,n->ij", relative_pos[keep], relative_pos[keep], mass[keep]) / mass[keep].sum()
    eigenvalues = np.sort(np.linalg.eigvalsh(covariance))[::-1]
    return float(np.sqrt(eigenvalues[1] / eigenvalues[0])), float(np.sqrt(eigenvalues[2] / eigenvalues[0]))


def group_metrics(mask, mass, radius, cylindrical_radius, z, v_radius, v_phi, v_z, gas_radius, gas_mass, new_radius, new_mass):
    count = int(np.count_nonzero(mask))
    if count < 10:
        return {"count": count}
    weights = mass[mask]
    sigma_r = weighted_std(v_radius[mask], weights)
    sigma_phi = weighted_std(v_phi[mask], weights)
    sigma_z = weighted_std(v_z[mask], weights)
    sigma_random = math.sqrt((sigma_r**2 + sigma_phi**2 + sigma_z**2) / 3.0)
    mean_vphi = float(np.average(v_phi[mask], weights=weights))
    mean_vr = float(np.average(v_radius[mask], weights=weights))
    mean_vz = float(np.average(v_z[mask], weights=weights))
    return {
        "count": count,
        "old_star_mass_msun": float(weights.sum()),
        "old_star_density_msun_pc3": float(weights.sum() / (4.0 * math.pi / 3.0 * (RADIAL_EDGES[np.searchsorted(RADIAL_EDGES, np.max(radius[mask]), side='left')] ** 3 if False else 1.0))),
        "radius_median_kpc": weighted_quantile(radius[mask], weights, 0.5),
        "cylindrical_radius_median_kpc": weighted_quantile(cylindrical_radius[mask], weights, 0.5),
        "abs_z_median_kpc": weighted_quantile(np.abs(z[mask]), weights, 0.5),
        "sigma_random_1d_km_s": sigma_random,
        "sigma_R_km_s": sigma_r,
        "sigma_phi_km_s": sigma_phi,
        "sigma_z_km_s": sigma_z,
        "mean_vphi_km_s": mean_vphi,
        "mean_vR_km_s": mean_vr,
        "mean_vz_km_s": mean_vz,
        "vrot_to_sigma": mean_vphi / sigma_random if sigma_random > 0 else float("nan"),
        "ordered_phi_fraction": mean_vphi**2 / (mean_vphi**2 + sigma_phi**2) if mean_vphi**2 + sigma_phi**2 > 0 else float("nan"),
        "sigma_z_to_sigma_R": sigma_z / sigma_r if sigma_r > 0 else float("nan"),
        "sigma_phi_to_sigma_R": sigma_phi / sigma_r if sigma_r > 0 else float("nan"),
    }


def analyze_snapshot(task):
    model, model_dir, snapshot = task
    path = Path(model_dir) / "output" / f"snapshot_{snapshot:03d}.hdf5"
    with h5py.File(path, "r") as handle:
        time_gyr = float(handle["Header"].attrs["Time"])
        pos, vel, mass = load_old_stars(handle)
        center = shrinking_center(pos, mass)
        relative_pos = pos - center
        radius = np.linalg.norm(relative_pos, axis=1)
        center_velocity = weighted_mean(vel[radius < 0.5], mass[radius < 0.5])
        relative_vel = vel - center_velocity
        gas = load_nearby(handle, "PartType0", center)
        new = load_nearby(handle, "PartType4", center)

    gas_pos, gas_vel, gas_mass = gas["pos"], gas["vel"], gas["mass"]
    new_pos, new_mass = new["pos"], new["mass"]
    gas_radius = np.linalg.norm(gas_pos - center, axis=1) if len(gas_mass) else np.empty(0)
    new_radius = np.linalg.norm(new_pos - center, axis=1) if len(new_mass) else np.empty(0)
    z_axis = intrinsic_axis(relative_pos, relative_vel, mass, radius)
    cylindrical_radius, z, v_radius, v_phi, v_z = cylindrical_components(relative_pos, relative_vel, z_axis)
    if len(gas_mass):
        gas_relative_pos = gas_pos - center
        gas_relative_vel = gas_vel - center_velocity
        gas_cylindrical_radius, gas_z, gas_vR, gas_vphi, gas_vz = cylindrical_components(gas_relative_pos, gas_relative_vel, z_axis)
    else:
        gas_cylindrical_radius = gas_z = gas_vR = gas_vphi = gas_vz = np.empty(0)
    half_mass_radius = weighted_quantile(radius, mass, 0.5)

    rows = []
    for index, label in enumerate(RADIAL_LABELS):
        left, right = RADIAL_EDGES[index], RADIAL_EDGES[index + 1]
        mask = (radius >= left) & (radius < right)
        row = group_metrics(mask, mass, radius, cylindrical_radius, z, v_radius, v_phi, v_z, gas_radius, gas_mass, new_radius, new_mass)
        volume = 4.0 * math.pi / 3.0 * (right**3 - left**3)
        gas_shell = float(gas_mass[(gas_radius >= left) & (gas_radius < right)].sum()) if len(gas_mass) else 0.0
        new_shell = float(new_mass[(new_radius >= left) & (new_radius < right)].sum()) if len(new_mass) else 0.0
        row["old_star_density_msun_pc3"] = row.get("old_star_mass_msun", 0.0) / volume / 1.0e9
        row["baryon_density_msun_pc3"] = (row.get("old_star_mass_msun", 0.0) + gas_shell + new_shell) / volume / 1.0e9
        row["gas_density_msun_pc3"] = gas_shell / volume / 1.0e9
        gas_mask = (gas_radius >= left) & (gas_radius < right) if len(gas_mass) else np.zeros(0, dtype=bool)
        if np.count_nonzero(gas_mask) >= 10:
            gas_weights = gas_mass[gas_mask]
            row["gas_sigma_R_km_s"] = weighted_std(gas_vR[gas_mask], gas_weights)
            row["gas_sigma_phi_km_s"] = weighted_std(gas_vphi[gas_mask], gas_weights)
            row["gas_sigma_z_km_s"] = weighted_std(gas_vz[gas_mask], gas_weights)
            row["gas_sigma_random_1d_km_s"] = math.sqrt((row["gas_sigma_R_km_s"]**2 + row["gas_sigma_phi_km_s"]**2 + row["gas_sigma_z_km_s"]**2) / 3.0)
            row["gas_mean_vphi_km_s"] = float(np.average(gas_vphi[gas_mask], weights=gas_weights))
        else:
            row.update({key: float("nan") for key in ("gas_sigma_R_km_s", "gas_sigma_phi_km_s", "gas_sigma_z_km_s", "gas_sigma_random_1d_km_s", "gas_mean_vphi_km_s")})
        row.update({
            "model": model,
            "snapshot": snapshot,
            "time_gyr": time_gyr,
            "radial_bin": label,
            "radial_bin_index": index,
            "r_left_kpc": left,
            "r_right_kpc": right,
            "gas_mass_shell_msun": gas_shell,
            "new_star_mass_shell_msun": new_shell,
        })
        rows.append(row)

    central = (radius < 0.2)
    plateau = (radius >= 0.4) & (radius < 0.7)
    central_metrics = group_metrics(central, mass, radius, cylindrical_radius, z, v_radius, v_phi, v_z, gas_radius, gas_mass, new_radius, new_mass)
    plateau_metrics = group_metrics(plateau, mass, radius, cylindrical_radius, z, v_radius, v_phi, v_z, gas_radius, gas_mass, new_radius, new_mass)
    p1, q1 = shape_ratios(relative_pos, mass, radius, 1.0)
    p2, q2 = shape_ratios(relative_pos, mass, radius, 2.0)
    central_cylinder = cylindrical_radius < 0.2
    central_sphere = radius < 0.2
    star_mass_r02 = float(mass[central_sphere].sum())
    gas_mass_r02 = float(gas_mass[gas_radius < 0.2].sum()) if len(gas_mass) else 0.0
    new_mass_r02 = float(new_mass[new_radius < 0.2].sum()) if len(new_mass) else 0.0
    gas_inner = gas_radius < 0.7 if len(gas_mass) else np.zeros(0, dtype=bool)
    if np.count_nonzero(gas_inner) >= 10:
        gas_inner_sigma_R = weighted_std(gas_vR[gas_inner], gas_mass[gas_inner])
        gas_inner_sigma_phi = weighted_std(gas_vphi[gas_inner], gas_mass[gas_inner])
        gas_inner_sigma_z = weighted_std(gas_vz[gas_inner], gas_mass[gas_inner])
        gas_inner_sigma_random = math.sqrt((gas_inner_sigma_R**2 + gas_inner_sigma_phi**2 + gas_inner_sigma_z**2) / 3.0)
        gas_inner_mean_vphi = float(np.average(gas_vphi[gas_inner], weights=gas_mass[gas_inner]))
    else:
        gas_inner_sigma_R = gas_inner_sigma_phi = gas_inner_sigma_z = gas_inner_sigma_random = gas_inner_mean_vphi = float("nan")
    gas_half_mass_radius = weighted_quantile(gas_radius, gas_mass, 0.5) if len(gas_mass) else float("nan")
    gas_half_mass_cylindrical_radius = weighted_quantile(gas_cylindrical_radius, gas_mass, 0.5) if len(gas_mass) else float("nan")
    gas_abs_z_median = weighted_quantile(np.abs(gas_z), gas_mass, 0.5) if len(gas_mass) else float("nan")
    summary = {
        "model": model,
        "snapshot": snapshot,
        "time_gyr": time_gyr,
        "path": str(path),
        "old_star_mass_msun": float(mass.sum()),
        "new_star_mass_within_10kpc_msun": float(new_mass.sum()),
        "gas_mass_within_10kpc_msun": float(gas_mass.sum()),
        "gas_mass_within_0p2kpc_msun": gas_mass_r02,
        "gas_mass_within_0p7kpc_msun": float(gas_mass[gas_radius < 0.7].sum()) if len(gas_mass) else 0.0,
        "gas_mass_within_1kpc_msun": float(gas_mass[gas_radius < 1.0].sum()) if len(gas_mass) else 0.0,
        "gas_mass_within_2kpc_msun": float(gas_mass[gas_radius < 2.0].sum()) if len(gas_mass) else 0.0,
        "gas_mass_within_5kpc_msun": float(gas_mass[gas_radius < 5.0].sum()) if len(gas_mass) else 0.0,
        "gas_half_mass_radius_3d_kpc": gas_half_mass_radius,
        "gas_half_mass_radius_cylindrical_kpc": gas_half_mass_cylindrical_radius,
        "gas_abs_z_median_kpc": gas_abs_z_median,
        "gas_inner_sigma_random_km_s": gas_inner_sigma_random,
        "gas_inner_sigma_R_km_s": gas_inner_sigma_R,
        "gas_inner_sigma_phi_km_s": gas_inner_sigma_phi,
        "gas_inner_sigma_z_km_s": gas_inner_sigma_z,
        "gas_inner_mean_vphi_km_s": gas_inner_mean_vphi,
        "gas_sfr_within_10kpc_msun_yr": float(np.nansum(gas["sfr"])) if len(gas_mass) else 0.0,
        "gas_density_field_mass_weighted": float(np.average(gas["density"], weights=gas_mass)) if len(gas_mass) and np.any(np.isfinite(gas["density"])) else float("nan"),
        "gas_internal_energy_mass_weighted": float(np.average(gas["internal_energy"], weights=gas_mass)) if len(gas_mass) and np.any(np.isfinite(gas["internal_energy"])) else float("nan"),
        "half_mass_radius_3d_kpc": half_mass_radius,
        "shape_p_within_1kpc": p1,
        "shape_q_within_1kpc": q1,
        "shape_p_within_2kpc": p2,
        "shape_q_within_2kpc": q2,
        "central_sigma_random_km_s": central_metrics["sigma_random_1d_km_s"],
        "plateau_sigma_random_km_s": plateau_metrics["sigma_random_1d_km_s"],
        "central_to_plateau_sigma_ratio": central_metrics["sigma_random_1d_km_s"] / plateau_metrics["sigma_random_1d_km_s"],
        "central_sigma_R_km_s": central_metrics["sigma_R_km_s"],
        "central_sigma_phi_km_s": central_metrics["sigma_phi_km_s"],
        "central_sigma_z_km_s": central_metrics["sigma_z_km_s"],
        "central_mean_vphi_km_s": central_metrics["mean_vphi_km_s"],
        "central_vrot_to_sigma": central_metrics["vrot_to_sigma"],
        "central_ordered_phi_fraction": central_metrics["ordered_phi_fraction"],
        "plateau_sigma_R_km_s": plateau_metrics["sigma_R_km_s"],
        "plateau_sigma_phi_km_s": plateau_metrics["sigma_phi_km_s"],
        "plateau_sigma_z_km_s": plateau_metrics["sigma_z_km_s"],
        "plateau_mean_vphi_km_s": plateau_metrics["mean_vphi_km_s"],
        "plateau_vrot_to_sigma": plateau_metrics["vrot_to_sigma"],
        "central_abs_z_median_kpc": weighted_quantile(np.abs(z[central_cylinder]), mass[central_cylinder], 0.5),
        "central_surface_density_msun_pc2": float(mass[central_cylinder].sum() / (math.pi * 0.2**2) / 1.0e6),
        "central_old_density_msun_pc3": star_mass_r02 / (4.0 * math.pi / 3.0 * 0.2**3) / 1.0e9,
        "central_baryon_density_msun_pc3": (star_mass_r02 + gas_mass_r02 + new_mass_r02) / (4.0 * math.pi / 3.0 * 0.2**3) / 1.0e9,
        "central_old_mass_msun": star_mass_r02,
        "central_gas_mass_msun": gas_mass_r02,
        "central_new_star_mass_msun": new_mass_r02,
        "intrinsic_axis_x": float(z_axis[0]),
        "intrinsic_axis_y": float(z_axis[1]),
        "intrinsic_axis_z": float(z_axis[2]),
    }
    return {"summary": summary, "profiles": rows}


def write_csv(path, rows):
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path):
    rows = []
    with Path(path).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            converted = {}
            for key, value in row.items():
                if value == "":
                    converted[key] = float("nan")
                    continue
                if key in {"model", "path", "radial_bin"}:
                    converted[key] = value
                    continue
                try:
                    number = float(value)
                    converted[key] = int(number) if key in {"snapshot", "radial_bin_index", "count"} else number
                except ValueError:
                    converted[key] = value
            rows.append(converted)
    return rows


def select_rows(rows, snapshot):
    return sorted([row for row in rows if row["snapshot"] == snapshot], key=lambda row: row["radial_bin_index"])


def make_figure(summaries, profiles, output):
    model_names = sorted({row["model"] for row in summaries})
    palette = ("#1f77b4", "#d62728", "#2ca02c", "#9467bd")
    colors = {model: palette[index % len(palette)] for index, model in enumerate(model_names)}
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 9.5), constrained_layout=True)
    comparison_snapshots = (15, 386)
    for model in colors:
        model_profiles = [row for row in profiles if row["model"] == model]
        for snapshot, linestyle in zip(comparison_snapshots, ("-", "--")):
            rows = select_rows(model_profiles, snapshot)
            radius = [(row["r_left_kpc"] + row["r_right_kpc"]) / 2.0 for row in rows]
            label = f"{model} {snapshot / 100:.2f} Gyr"
            axes[0, 0].plot(radius, [row.get("sigma_random_1d_km_s", np.nan) for row in rows], color=colors[model], ls=linestyle, marker="o", ms=3, label=label)
            axes[0, 1].plot(radius, [row["baryon_density_msun_pc3"] for row in rows], color=colors[model], ls=linestyle, marker="o", ms=3, label=label)
            axes[0, 2].plot(radius, [row.get("mean_vphi_km_s", np.nan) for row in rows], color=colors[model], ls=linestyle, marker="o", ms=3, label=label)

        rows = select_rows(model_profiles, 15)
        radius = [(row["r_left_kpc"] + row["r_right_kpc"]) / 2.0 for row in rows]
        for key, ls, component in (("sigma_R_km_s", "-", "R"), ("sigma_phi_km_s", "--", "phi"), ("sigma_z_km_s", ":", "z")):
            axes[1, 0].plot(radius, [row.get(key, np.nan) for row in rows], color=colors[model], ls=ls, label=f"{model} sigma_{component}")

        model_summary = sorted([row for row in summaries if row["model"] == model], key=lambda row: row["time_gyr"])
        axes[1, 1].plot([row["time_gyr"] for row in model_summary], [row["central_to_plateau_sigma_ratio"] for row in model_summary], color=colors[model], marker="o", label=model)
        axes[1, 2].plot([row["time_gyr"] for row in model_summary], [row["shape_q_within_1kpc"] for row in model_summary], color=colors[model], marker="o", ls="-", label=f"{model} q")
        axes[1, 2].plot([row["time_gyr"] for row in model_summary], [row["half_mass_radius_3d_kpc"] for row in model_summary], color=colors[model], marker="s", ls="--", label=f"{model} r_half")

    axes[0, 0].set(xlim=(0, 1.6), xlabel="3-D radius [kpc]", ylabel="Random 1-D dispersion [km/s]", title="Intrinsic random-dispersion profiles")
    axes[0, 1].set(xlim=(0.02, 1.6), yscale="log", xlabel="3-D radius [kpc]", ylabel="Baryonic density [Msun/pc3]", title="Spherical-shell density")
    axes[0, 2].set(xlim=(0, 1.6), xlabel="3-D radius [kpc]", ylabel="Mean azimuthal velocity [km/s]", title="Ordered rotation")
    axes[1, 0].set(xlim=(0, 1.6), xlabel="3-D radius [kpc]", ylabel="Component dispersion [km/s]", title="Velocity ellipsoid at 0.15 Gyr")
    axes[1, 1].axhline(1.0, color="0.4", lw=1.0, ls=":")
    axes[1, 1].set(xlabel="Isolation time [Gyr]", ylabel="sigma(<0.2) / sigma(0.4-0.7)", title="Persistence of the central dip")
    axes[1, 2].set(xlabel="Isolation time [Gyr]", ylabel="Axis ratio q or r_half [kpc]", title="Shape and size evolution")
    for axis in axes.flat:
        axis.grid(alpha=0.2)
        axis.legend(frameon=False, fontsize=7)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def make_gas_figure(summaries, profiles, output):
    model_names = sorted({row["model"] for row in summaries})
    palette = ("#1f77b4", "#d62728", "#2ca02c", "#9467bd")
    colors = {model: palette[index % len(palette)] for index, model in enumerate(model_names)}
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 9.5), constrained_layout=True)
    for model in colors:
        model_profiles = [row for row in profiles if row["model"] == model]
        for snapshot, linestyle in ((15, "-"), (386, "--")):
            rows = select_rows(model_profiles, snapshot)
            radius = np.array([(row["r_left_kpc"] + row["r_right_kpc"]) / 2.0 for row in rows])
            gas_shell = np.array([row["gas_mass_shell_msun"] for row in rows])
            label = f"{model} {snapshot / 100:.2f} Gyr"
            axes[0, 0].plot(radius, np.cumsum(gas_shell), color=colors[model], ls=linestyle, marker="o", ms=3, label=label)
            axes[0, 1].plot(radius, [row["gas_density_msun_pc3"] for row in rows], color=colors[model], ls=linestyle, marker="o", ms=3, label=label)
            axes[0, 2].plot(radius, [row.get("gas_mean_vphi_km_s", np.nan) for row in rows], color=colors[model], ls=linestyle, marker="o", ms=3, label=label)
        model_summary = sorted([row for row in summaries if row["model"] == model], key=lambda row: row["time_gyr"])
        time = [row["time_gyr"] for row in model_summary]
        axes[1, 0].plot(time, [row["gas_mass_within_0p7kpc_msun"] for row in model_summary], color=colors[model], marker="o", ls="-", label=f"{model} Mgas(<0.7)")
        axes[1, 0].plot(time, [row["gas_mass_within_2kpc_msun"] for row in model_summary], color=colors[model], marker="s", ls="--", label=f"{model} Mgas(<2)")
        axes[1, 1].plot(time, [row["gas_half_mass_radius_cylindrical_kpc"] for row in model_summary], color=colors[model], marker="o", ls="-", label=f"{model} Rhalf")
        axes[1, 1].plot(time, [row["gas_abs_z_median_kpc"] for row in model_summary], color=colors[model], marker="s", ls="--", label=f"{model} |z|med")
        axes[1, 2].plot(time, [row["gas_sfr_within_10kpc_msun_yr"] for row in model_summary], color=colors[model], marker="o", ls="-", label=f"{model} SFR")
    axes[0, 0].set(xlim=(0, 3.1), xlabel="3-D radius [kpc]", ylabel="Cumulative gas mass [Msun]", title="Gas mass distribution")
    axes[0, 1].set(xlim=(0.02, 3.1), yscale="log", xlabel="3-D radius [kpc]", ylabel="Gas density [Msun/pc3]", title="Gas density profile")
    axes[0, 2].set(xlim=(0, 3.1), xlabel="3-D radius [kpc]", ylabel="Mean gas vphi [km/s]", title="Gas rotation")
    axes[1, 0].set(xlabel="Isolation time [Gyr]", ylabel="Gas mass [Msun]", title="Central gas retention")
    axes[1, 1].set(xlabel="Isolation time [Gyr]", ylabel="Scale [kpc]", title="Gas size and thickness")
    axes[1, 2].set(xlabel="Isolation time [Gyr]", ylabel="SFR [Msun/yr]", title="Star formation history")
    for axis in axes.flat:
        axis.grid(alpha=0.2)
        axis.legend(frameon=False, fontsize=7)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", action="append", required=True, help="NAME=PATH")
    parser.add_argument("--snapshots", default="0,15,50,150,300,386,600")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args()
    models = []
    for item in args.model:
        name, path = item.split("=", 1)
        models.append((name, path))
    snapshots = [int(value) for value in args.snapshots.split(",")]
    args.output.mkdir(parents=True, exist_ok=True)
    if args.plot_only:
        summaries = read_csv(args.output / "isolation_2047_4001_summary.csv")
        profiles = read_csv(args.output / "isolation_2047_4001_profiles.csv")
    else:
        tasks = [(name, path, snapshot) for name, path in models for snapshot in snapshots]
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            results = list(pool.map(analyze_snapshot, tasks))
        summaries = [result["summary"] for result in results]
        profiles = [row for result in results for row in result["profiles"]]
        summaries.sort(key=lambda row: (row["model"], row["snapshot"]))
        profiles.sort(key=lambda row: (row["model"], row["snapshot"], row["radial_bin_index"]))
        write_csv(args.output / "isolation_2047_4001_summary.csv", summaries)
        write_csv(args.output / "isolation_2047_4001_profiles.csv", profiles)
        with (args.output / "isolation_2047_4001_metadata.json").open("w", encoding="utf-8") as handle:
            json.dump({"models": models, "snapshots": snapshots, "workers": args.workers}, handle, indent=2)
    make_figure(summaries, profiles, args.output / "isolation_2047_4001_diagnosis.png")
    make_gas_figure(summaries, profiles, args.output / "isolation_2047_4001_gas_diagnosis.png")


if __name__ == "__main__":
    main()
