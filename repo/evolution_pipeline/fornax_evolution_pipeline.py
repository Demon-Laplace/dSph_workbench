#!/usr/bin/env python3
"""Configuration-driven physical-evolution analysis for GIZMO dwarf runs.

The extractor reuses the dSph_workbench snapshot context and observational
kinematics definitions.  The resulting CSV is self-contained for all derived
time-series calculations and plotting; the ``plot`` subcommand never reopens a
simulation snapshot.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
import tempfile
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

HERE = Path(__file__).resolve().parent
REPO_ROOT = Path(os.environ.get("FORNAX_REPO_ROOT", str(HERE.parent)))
sys.path.insert(0, str(REPO_ROOT))

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".mplconfig"))

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, gaussian_filter1d
from scipy.signal import savgol_filter

from basefunc import Analysis
from snapshot_context import prepare_snapshot_context
from snapshot_metrics import (
    detrended_dispersion_in_aperture,
    old_dwarf_star_local_mask,
    old_star_projected_kinematics,
    compute_snapshot_summary,
)


G_KPC_KMS2_PER_MSUN = 4.30091e-6
KPC_CM = 3.0856775814913673e21
MSUN_G = 1.98847e33
PROTON_MASS_G = 1.67262192369e-24
BOLTZMANN_ERG_PER_K = 1.380649e-16
KPC_PER_KMS_TO_GYR = 0.9777922216807892
NHI_PER_MSUN_PC2 = 1.248e20


def _nested(config: Mapping[str, Any], path: str, default: Any = None) -> Any:
    value: Any = config
    for part in path.split("."):
        if not isinstance(value, Mapping) or part not in value:
            return default
        value = value[part]
    return value


def _required(config: Mapping[str, Any], path: str) -> Any:
    value = _nested(config, path, None)
    if value is None:
        raise ValueError(f"Missing required configuration value: {path}")
    return value


def load_config(path: Path) -> tuple[dict[str, Any], Path]:
    path = path.resolve()
    with path.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    if not isinstance(config, dict):
        raise ValueError("The configuration root must be a JSON object")
    return config, path


def resolve_paths(config: Mapping[str, Any], config_path: Path) -> dict[str, Path]:
    run_dir = Path(str(_required(config, "paths.run_dir"))).expanduser()
    if not run_dir.is_absolute():
        run_dir = (config_path.parent / run_dir).resolve()
    else:
        run_dir = run_dir.resolve()

    snapshot_dir = Path(str(_nested(config, "paths.snapshot_dir", "output")))
    if not snapshot_dir.is_absolute():
        snapshot_dir = run_dir / snapshot_dir

    output_dir = Path(str(_nested(config, "paths.output_dir", "evolution_analysis")))
    if not output_dir.is_absolute():
        output_dir = run_dir / output_dir

    csv_name = str(_nested(config, "paths.timeseries_csv", "fornax_evolution_timeseries.csv"))
    figure_stem = str(_nested(config, "paths.figure_stem", "fornax_physical_evolution"))
    return {
        "run_dir": run_dir,
        "snapshot_dir": snapshot_dir.resolve(),
        "output_dir": output_dir.resolve(),
        "csv": (output_dir / csv_name).resolve(),
        "figure_stem": (output_dir / figure_stem).resolve(),
        "metadata": (output_dir / "fornax_evolution_metadata.json").resolve(),
    }


def analysis_config_hash(config: Mapping[str, Any]) -> str:
    """Hash choices that require rereading particle data."""
    relevant = {
        key: config.get(key)
        for key in (
            "snapshots",
            "particle_types",
            "dwarf_selection",
            "projection",
            "stellar",
            "gas",
            "hi",
            "cgm",
            "enclosed_mass",
        )
    }
    payload = json.dumps(relevant, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def derivation_config_hash(config: Mapping[str, Any]) -> str:
    """Hash table-only choices that can be changed without reopening snapshots."""
    relevant = {
        key: config.get(key)
        for key in ("smoothing", "comparison_epoch", "pericentre_detection")
    }
    payload = json.dumps(relevant, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def discover_snapshots(snapshot_dir: Path, config: Mapping[str, Any]) -> list[tuple[int, Path]]:
    glob_pattern = str(_nested(config, "snapshots.glob", "snapshot_*.hdf5"))
    number_regex = re.compile(
        str(_nested(config, "snapshots.number_regex", r"snapshot_(\d+)\.hdf5$"))
    )
    start = int(_nested(config, "snapshots.start", 0))
    stop_raw = _nested(config, "snapshots.stop", None)
    stop = None if stop_raw is None else int(stop_raw)
    step = int(_nested(config, "snapshots.step", 1))
    if step <= 0:
        raise ValueError("snapshots.step must be positive")

    discovered: list[tuple[int, Path]] = []
    for path in sorted(snapshot_dir.glob(glob_pattern)):
        match = number_regex.search(path.name)
        if match is None:
            continue
        number = int(match.group(1))
        if number < start or (stop is not None and number > stop):
            continue
        if (number - start) % step != 0:
            continue
        discovered.append((number, path.resolve()))
    if not discovered:
        raise FileNotFoundError(
            f"No snapshots matched {glob_pattern!r} in {snapshot_dir} "
            f"for range start={start}, stop={stop}, step={step}"
        )
    return discovered


def weighted_mean_vectors(vectors: np.ndarray, weights: np.ndarray) -> np.ndarray:
    vectors = np.asarray(vectors, dtype=float)
    weights = np.asarray(weights, dtype=float)
    good = np.all(np.isfinite(vectors), axis=1) & np.isfinite(weights) & (weights > 0.0)
    if np.count_nonzero(good) == 0:
        return np.full(3, np.nan)
    return np.average(vectors[good], axis=0, weights=weights[good])


def rotation_basis(inclination_deg: float, azimuth_deg: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    inc = np.radians(float(inclination_deg))
    azi = np.radians(float(azimuth_deg))
    los = np.array([np.sin(inc) * np.cos(azi), np.sin(inc) * np.sin(azi), np.cos(inc)])
    x_axis = np.array([-np.sin(azi), np.cos(azi), 0.0])
    y_axis = np.cross(los, x_axis)
    x_axis /= np.linalg.norm(x_axis)
    y_axis /= np.linalg.norm(y_axis)
    los /= np.linalg.norm(los)
    return x_axis, y_axis, los


def projected_old_star_observables(
    snapshot: Mapping[str, Any],
    summary: Mapping[str, Any],
    projection: Mapping[str, Any],
) -> dict[str, float]:
    mode = str(projection.get("mode", "native")).lower()
    if mode == "native":
        return {
            "re_major_kpc": float(summary["rhalf"]),
            "re_circular_kpc": float(summary["rhalf_circular"]),
            "sigma_los_kms": float(summary["sigma_re_circular"]),
            "axis_ratio": float(1.0 - summary["eps"]),
            "pa_rad": float(summary["pa"]),
            "shape_center_x_kpc": float(summary["shape_center_x_kpc"]),
            "shape_center_y_kpc": float(summary["shape_center_y_kpc"]),
            "sigma_n_old_stars": float(summary["sigma_re_noldstar"]),
            "sigma_gradient_kms_per_kpc": float(summary["sigma_gradient_kms_per_kpc"]),
        }
    if mode != "euler":
        raise ValueError("projection.mode must be 'native' or 'euler'")

    df = snapshot["df"]
    dwarf_mask = np.asarray(snapshot["total_dw_star_mask"], dtype=bool)
    old_local = old_dwarf_star_local_mask(snapshot)
    old_global_indices = np.flatnonzero(dwarf_mask)[old_local]
    center = np.array([snapshot["dw_xc"], snapshot["dw_yc"], snapshot["dw_zc"]], dtype=float)
    positions = df.loc[old_global_indices, ["x", "y", "z"]].to_numpy(dtype=float) - center
    velocities = df.loc[old_global_indices, ["vx", "vy", "vz"]].to_numpy(dtype=float)
    masses = df.loc[old_global_indices, "m"].to_numpy(dtype=float)

    inclination = float(_required(projection, "inclination_deg"))
    azimuth = float(_required(projection, "azimuth_deg"))
    x_axis, y_axis, los = rotation_basis(inclination, azimuth)
    x_kpc = positions @ x_axis
    y_kpc = positions @ y_axis
    vlos = velocities @ los

    shape = Analysis.calculate_projected_shape(x_kpc, y_kpc, mass=masses, n_neighbors=30)
    eps = float(shape["eps"]) if np.isfinite(shape["eps"]) else 0.0
    pa = float(shape["pa"]) if np.isfinite(shape["pa"]) else 0.0
    center_x = float(shape["center_x"]) if np.isfinite(shape["center_x"]) else 0.0
    center_y = float(shape["center_y"]) if np.isfinite(shape["center_y"]) else 0.0
    re_major = Analysis.half_light_radius(
        x_kpc, y_kpc, masses, ep=eps, pa=pa, center_x=center_x, center_y=center_y
    )
    re_circular = Analysis.half_light_radius(
        x_kpc, y_kpc, masses, ep=0.0, pa=0.0, center_x=center_x, center_y=center_y
    )
    dispersion = detrended_dispersion_in_aperture(
        x_kpc,
        y_kpc,
        vlos,
        re_circular,
        center_x_kpc=center_x,
        center_y_kpc=center_y,
        circular=True,
    )
    return {
        "re_major_kpc": float(re_major),
        "re_circular_kpc": float(re_circular),
        "sigma_los_kms": float(dispersion["sigma"]),
        "axis_ratio": float(1.0 - eps),
        "pa_rad": pa,
        "shape_center_x_kpc": center_x,
        "shape_center_y_kpc": center_y,
        "sigma_n_old_stars": float(dispersion["nstar"]),
        "sigma_gradient_kms_per_kpc": float(dispersion["gradient"]["grad_amp"]),
    }


def read_gas_smoothing_lengths(snapshot_path: Path, gas_count: int, fallback_kpc: float) -> np.ndarray:
    with h5py.File(snapshot_path, "r") as handle:
        if "PartType0" not in handle:
            return np.empty(0, dtype=float)
        group = handle["PartType0"]
        if "SmoothingLength" in group:
            values = np.asarray(group["SmoothingLength"], dtype=float)
        else:
            values = np.full(gas_count, float(fallback_kpc), dtype=float)
    if values.size != gas_count:
        raise RuntimeError(
            f"Gas smoothing-length count {values.size} does not match loaded gas count {gas_count}"
        )
    return values


def adaptive_sph_map(
    x_deg: np.ndarray,
    y_deg: np.ndarray,
    weights_msun: np.ndarray,
    smoothing_length_kpc: np.ndarray,
    distance_kpc: float,
    half_width_deg: float,
    npix: int,
    minimum_sigma_pixels: float,
) -> np.ndarray:
    limits = [[-half_width_deg, half_width_deg], [-half_width_deg, half_width_deg]]
    inside = (
        np.isfinite(x_deg)
        & np.isfinite(y_deg)
        & np.isfinite(weights_msun)
        & np.isfinite(smoothing_length_kpc)
        & (weights_msun > 0.0)
        & (np.abs(x_deg) <= half_width_deg)
        & (np.abs(y_deg) <= half_width_deg)
    )
    if not np.any(inside):
        return np.zeros((npix, npix), dtype=float)
    x = x_deg[inside]
    y = y_deg[inside]
    weights = weights_msun[inside]
    h = np.clip(smoothing_length_kpc[inside], 1.0e-4, None)
    quantile_edges = np.unique(np.quantile(h, np.linspace(0.0, 1.0, 9)))
    if quantile_edges.size < 2:
        quantile_edges = np.array([h.min(), np.nextafter(h.max(), np.inf)])
    pixel_deg = 2.0 * half_width_deg / npix
    projected_mass = np.zeros((npix, npix), dtype=float)
    for index in range(quantile_edges.size - 1):
        lower, upper = quantile_edges[index : index + 2]
        if index == quantile_edges.size - 2:
            group = (h >= lower) & (h <= upper)
        else:
            group = (h >= lower) & (h < upper)
        if not np.any(group):
            continue
        image, _, _ = np.histogram2d(
            x[group], y[group], bins=npix, range=limits, weights=weights[group]
        )
        sigma_deg = np.rad2deg(np.median(h[group]) / distance_kpc) / 2.0
        sigma_pixels = np.clip(sigma_deg / pixel_deg, minimum_sigma_pixels, npix / 3.0)
        projected_mass += gaussian_filter(image, sigma_pixels, mode="constant")
    return projected_mass


def contour_hi_mass(
    snapshot: Mapping[str, Any],
    snapshot_path: Path,
    hi_config: Mapping[str, Any],
) -> float:
    df = snapshot["df"]
    gas_global = np.flatnonzero(df["tp"].to_numpy(dtype=int) == 0)
    hsml = read_gas_smoothing_lengths(
        snapshot_path,
        gas_global.size,
        float(hi_config.get("fallback_smoothing_length_kpc", 0.05)),
    )
    hsml_global = np.full(len(df), np.nan, dtype=float)
    hsml_global[gas_global] = hsml

    cold_mask = np.asarray(snapshot["dw_cold_gas_mask"], dtype=bool)
    if np.count_nonzero(cold_mask) == 0:
        return 0.0
    mass = df.loc[cold_mask, "m"].to_numpy(dtype=float)
    neutral = df.loc[cold_mask, "nh"].to_numpy(dtype=float)
    smooth = hsml_global[cold_mask]
    distance_kpc = float(snapshot["d_mean"])
    x_deg = np.degrees(np.asarray(snapshot["cold_gas_x_kpc"], dtype=float) / distance_kpc)
    y_deg = np.degrees(np.asarray(snapshot["cold_gas_y_kpc"], dtype=float) / distance_kpc)
    half_width = float(hi_config.get("field_half_width_deg", 2.1))
    npix = int(hi_config.get("map_pixels", 520))
    smoothing_pixels = float(hi_config.get("minimum_smoothing_pixels", 2.85))
    projected_mass = adaptive_sph_map(
        x_deg,
        y_deg,
        mass * neutral,
        smooth,
        distance_kpc,
        half_width,
        npix,
        smoothing_pixels,
    )
    pixel_deg = 2.0 * half_width / npix
    pixel_kpc = np.deg2rad(pixel_deg) * distance_kpc
    pixel_area_pc2 = (pixel_kpc * 1000.0) ** 2
    surface_density = projected_mass / pixel_area_pc2
    threshold = float(hi_config.get("contour_threshold_nhi_cm2", 5.0e18)) / NHI_PER_MSUN_PC2
    selected = np.isfinite(surface_density) & (surface_density >= threshold)
    return float(np.sum(surface_density[selected]) * pixel_area_pc2)


def local_cgm_measurement(
    df: pd.DataFrame,
    dwarf_center: np.ndarray,
    dwarf_velocity: np.ndarray,
    config: Mapping[str, Any],
) -> dict[str, float]:
    particle_type = int(config.get("particle_type", 0))
    gas = df["tp"].to_numpy(dtype=int) == particle_type
    positions = df[["x", "y", "z"]].to_numpy(dtype=float)
    velocities = df[["vx", "vy", "vz"]].to_numpy(dtype=float)
    masses = df["m"].to_numpy(dtype=float)
    temperature = df["temp"].to_numpy(dtype=float)
    neutral = df["nh"].to_numpy(dtype=float)
    radius = np.linalg.norm(positions - dwarf_center, axis=1)

    inner = float(config.get("exclusion_radius_kpc", 30.0))
    outer = float(config.get("search_radius_kpc", 60.0))
    hot_min = float(config.get("temperature_min_k", 2.0e4))
    neutral_max = float(config.get("neutral_fraction_max", 1.0))
    candidates = (
        gas
        & np.isfinite(radius)
        & (radius >= inner)
        & (radius <= outer)
        & np.isfinite(temperature)
        & (temperature >= hot_min)
        & np.isfinite(neutral)
        & (neutral <= neutral_max)
        & np.isfinite(masses)
        & (masses > 0.0)
    )
    indices = np.flatnonzero(candidates)
    method = str(config.get("method", "knn_shell")).lower()
    if method == "knn_shell" and indices.size:
        neighbours = int(config.get("neighbours", 128))
        order = np.argsort(radius[indices])
        indices = indices[order[: min(neighbours, indices.size)]]
    elif method != "fixed_shell":
        raise ValueError("cgm.method must be 'knn_shell' or 'fixed_shell'")

    minimum_particles = int(config.get("minimum_particles", 32))
    if indices.size < minimum_particles:
        return {
            "cgm_particle_count": float(indices.size),
            "cgm_effective_outer_radius_kpc": np.nan,
            "cgm_density_msun_kpc3": np.nan,
            "cgm_density_g_cm3": np.nan,
            "cgm_hydrogen_number_density_cm3": np.nan,
            "cgm_total_number_density_cm3": np.nan,
            "cgm_velocity_x_kms": np.nan,
            "cgm_velocity_y_kms": np.nan,
            "cgm_velocity_z_kms": np.nan,
            "v_rel_cgm_x_kms": np.nan,
            "v_rel_cgm_y_kms": np.nan,
            "v_rel_cgm_z_kms": np.nan,
            "v_rel_cgm_kms": np.nan,
            "ram_pressure_dyn_cm2": np.nan,
            "ram_pressure_over_kb_k_cm3": np.nan,
        }

    effective_outer = float(np.max(radius[indices])) if method == "knn_shell" else outer
    volume_kpc3 = 4.0 * np.pi / 3.0 * (effective_outer**3 - inner**3)
    density_msun_kpc3 = float(np.sum(masses[indices]) / volume_kpc3)
    density_g_cm3 = density_msun_kpc3 * MSUN_G / (KPC_CM**3)
    hydrogen_fraction = float(config.get("hydrogen_mass_fraction", 0.76))
    mean_molecular_weight = float(config.get("mean_molecular_weight", 0.61))
    nh_cm3 = hydrogen_fraction * density_g_cm3 / PROTON_MASS_G
    ntotal_cm3 = density_g_cm3 / (mean_molecular_weight * PROTON_MASS_G)
    cgm_velocity = weighted_mean_vectors(velocities[indices], masses[indices])
    relative = dwarf_velocity - cgm_velocity
    speed = float(np.linalg.norm(relative))
    pressure = density_g_cm3 * (speed * 1.0e5) ** 2
    return {
        "cgm_particle_count": float(indices.size),
        "cgm_effective_outer_radius_kpc": effective_outer,
        "cgm_density_msun_kpc3": density_msun_kpc3,
        "cgm_density_g_cm3": density_g_cm3,
        "cgm_hydrogen_number_density_cm3": nh_cm3,
        "cgm_total_number_density_cm3": ntotal_cm3,
        "cgm_velocity_x_kms": float(cgm_velocity[0]),
        "cgm_velocity_y_kms": float(cgm_velocity[1]),
        "cgm_velocity_z_kms": float(cgm_velocity[2]),
        "v_rel_cgm_x_kms": float(relative[0]),
        "v_rel_cgm_y_kms": float(relative[1]),
        "v_rel_cgm_z_kms": float(relative[2]),
        "v_rel_cgm_kms": speed,
        "ram_pressure_dyn_cm2": pressure,
        "ram_pressure_over_kb_k_cm3": pressure / BOLTZMANN_ERG_PER_K,
    }


def aperture_label(radius_kpc: float) -> str:
    text = f"{float(radius_kpc):g}".replace("-", "m").replace(".", "p")
    return f"{text}kpc"


def gas_fraction(stellar_mass: float, gas_mass: float) -> float:
    denominator = stellar_mass + gas_mass
    return float(gas_mass / denominator) if denominator > 0.0 else np.nan


def process_snapshot(
    number: int,
    snapshot_path: Path,
    config: Mapping[str, Any],
    config_hash: str,
) -> dict[str, Any]:
    selection = dict(_nested(config, "dwarf_selection", {}))
    gas_config = dict(_nested(config, "gas", {}))
    core_radius = float(_required(selection, "core_radius_kpc"))
    gas_aperture = float(gas_config.get("aperture_kpc", 20.0))
    hot_split = float(gas_config.get("temperature_split_k", 2.0e4))
    snapshot = prepare_snapshot_context(
        folder_path=str(snapshot_path.parent),
        snapshot_num=number,
        core_radius=core_radius,
        r_exclude=float(selection.get("mw_exclusion_radius_kpc", 5.0)),
        dwarf_radius_factor=float(selection.get("dwarf_radius_factor", 3.0)),
        k_density=int(selection.get("density_neighbours", 16)),
        dwarf_gas_radius=gas_aperture,
        gas_temperature_split=hot_split,
        include_mw_gas=True,
        mw_gas_radius=float(_nested(config, "cgm.mw_gas_radius_kpc", 500.0)),
        include_dark_matter=True,
        include_star_birth=True,
    )
    summary = compute_snapshot_summary(snapshot, number)
    projection = projected_old_star_observables(
        snapshot, summary, dict(_nested(config, "projection", {"mode": "native"}))
    )
    df = snapshot["df"]
    center = np.array([snapshot["dw_xc"], snapshot["dw_yc"], snapshot["dw_zc"]], dtype=float)
    positions = df[["x", "y", "z"]].to_numpy(dtype=float)
    velocities = df[["vx", "vy", "vz"]].to_numpy(dtype=float)
    masses = df["m"].to_numpy(dtype=float)
    particle_types = df["tp"].to_numpy(dtype=int)
    distance_from_dwarf = np.linalg.norm(positions - center, axis=1)

    star_types = np.asarray(_nested(config, "particle_types.stars", [2, 3, 4]), dtype=int)
    gas_type = int(_nested(config, "particle_types.gas", 0))
    star_mask = np.asarray(snapshot["total_dw_star_mask"], dtype=bool) & np.isin(particle_types, star_types)
    stellar_aperture = float(_nested(config, "stellar.mass_aperture_kpc", 20.0))
    star_mass_mask = star_mask & (distance_from_dwarf <= stellar_aperture)
    stellar_mass = float(np.sum(masses[star_mass_mask]))

    configured_com_aperture = _nested(config, "stellar.com_velocity_aperture_kpc", None)
    if configured_com_aperture is None:
        com_aperture = float(
            _nested(config, "stellar.com_velocity_aperture_re_multiple", 8.0)
        ) * float(projection["re_major_kpc"])
    else:
        com_aperture = float(configured_com_aperture)
    com_mask = star_mask & (distance_from_dwarf <= com_aperture)
    dwarf_velocity = weighted_mean_vectors(velocities[com_mask], masses[com_mask])

    dwarf_gas_mask = np.asarray(snapshot["dw_gas_mask"], dtype=bool) & (particle_types == gas_type)
    total_gas_mass = float(np.sum(masses[dwarf_gas_mask]))
    neutral = df["nh"].to_numpy(dtype=float)
    temperature = df["temp"].to_numpy(dtype=float)
    hi_mask = (
        dwarf_gas_mask
        & np.isfinite(neutral)
        & (neutral >= float(_nested(config, "hi.neutral_fraction_min", 0.0)))
        & np.isfinite(temperature)
        & (temperature < float(_nested(config, "hi.temperature_max_k", 2.0e4)))
    )
    hi_particle_mass = float(np.sum(masses[hi_mask] * neutral[hi_mask]))
    hi_config = dict(_nested(config, "hi", {}))
    hi_contour_mass = (
        contour_hi_mass(snapshot, snapshot_path, hi_config)
        if bool(hi_config.get("calculate_contour_mass", True))
        else np.nan
    )
    hi_definition = str(hi_config.get("mass_definition", "contour")).lower()
    if hi_definition == "contour":
        hi_mass = hi_contour_mass
    elif hi_definition == "particle":
        hi_mass = hi_particle_mass
    else:
        raise ValueError("hi.mass_definition must be 'contour' or 'particle'")

    cgm = local_cgm_measurement(
        df,
        center,
        dwarf_velocity,
        dict(_nested(config, "cgm", {})),
    )
    rgc = float(np.linalg.norm(center))
    mw_mass = float(summary["mw_mass_r"])
    tidal = G_KPC_KMS2_PER_MSUN * mw_mass / (rgc**3) if rgc > 0.0 else np.nan
    tidal_gyr2 = tidal * (1.0 / KPC_PER_KMS_TO_GYR) ** 2

    row: dict[str, Any] = {
        "analysis_config_sha256": config_hash,
        "snapshot": number,
        "snapshot_file": snapshot_path.name,
        "time_gyr": float(snapshot["tsnap"]),
        "dwarf_center_x_kpc": float(center[0]),
        "dwarf_center_y_kpc": float(center[1]),
        "dwarf_center_z_kpc": float(center[2]),
        "dwarf_velocity_com_x_kms": float(dwarf_velocity[0]),
        "dwarf_velocity_com_y_kms": float(dwarf_velocity[1]),
        "dwarf_velocity_com_z_kms": float(dwarf_velocity[2]),
        "distance_galactocentric_kpc": rgc,
        "distance_heliocentric_kpc": float(summary["distance"]),
        "mw_enclosed_mass_msun": mw_mass,
        "tidal_proxy_kms2_kpc2": tidal,
        "tidal_proxy_gyr2": tidal_gyr2,
        "stellar_mass_msun": stellar_mass,
        "stellar_mass_old_observational_aperture_msun": float(summary["star_mass"]),
        "gas_mass_msun": total_gas_mass,
        "hi_mass_msun": float(hi_mass),
        "hi_mass_particle_msun": hi_particle_mass,
        "hi_mass_contour_msun": float(hi_contour_mass),
        "re_major_kpc": projection["re_major_kpc"],
        "re_circular_kpc": projection["re_circular_kpc"],
        "sigma_los_kms": projection["sigma_los_kms"],
        "sigma_n_old_stars": projection["sigma_n_old_stars"],
        "sigma_gradient_kms_per_kpc": projection["sigma_gradient_kms_per_kpc"],
        "axis_ratio": projection["axis_ratio"],
        "position_angle_rad": projection["pa_rad"],
        "position_angle_deg": float(np.degrees(projection["pa_rad"])),
        "shape_center_x_kpc": projection["shape_center_x_kpc"],
        "shape_center_y_kpc": projection["shape_center_y_kpc"],
        "projection_mode": str(_nested(config, "projection.mode", "native")),
        "projection_inclination_deg": _nested(config, "projection.inclination_deg", np.nan),
        "projection_azimuth_deg": _nested(config, "projection.azimuth_deg", np.nan),
        "stellar_particle_count": int(np.count_nonzero(star_mass_mask)),
        "gas_particle_count": int(np.count_nonzero(dwarf_gas_mask)),
        "hi_particle_count": int(np.count_nonzero(hi_mask)),
    }
    row.update(cgm)

    radii = [float(value) for value in _nested(config, "enclosed_mass.radii_kpc", [0.5, 1.0, 2.0])]
    for radius_kpc in radii:
        label = aperture_label(radius_kpc)
        inside = distance_from_dwarf <= radius_kpc
        star_enclosed = float(np.sum(masses[inside & star_mask]))
        gas_enclosed = float(np.sum(masses[inside & (particle_types == gas_type)]))
        row[f"stellar_mass_3d_lt_{label}_msun"] = star_enclosed
        row[f"gas_mass_3d_lt_{label}_msun"] = gas_enclosed
        row[f"gas_fraction_3d_lt_{label}"] = gas_fraction(star_enclosed, gas_enclosed)

    if bool(_nested(config, "enclosed_mass.include_effective_radius", True)):
        re = float(projection["re_major_kpc"])
        inside = distance_from_dwarf <= re
        star_enclosed = float(np.sum(masses[inside & star_mask]))
        gas_enclosed = float(np.sum(masses[inside & (particle_types == gas_type)]))
        row["stellar_mass_3d_lt_re_msun"] = star_enclosed
        row["gas_mass_3d_lt_re_msun"] = gas_enclosed
        row["gas_fraction_3d_lt_re"] = gas_fraction(star_enclosed, gas_enclosed)

    snapshot["simulation"].df = None
    return row


def atomic_write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(
            "w", encoding="utf-8", newline="", dir=path.parent, prefix=f".{path.name}.", delete=False
        ) as handle:
            tmp_path = Path(handle.name)
            frame.to_csv(handle, index=False, float_format="%.10g")
        os.replace(tmp_path, path)
    finally:
        if tmp_path is not None and tmp_path.exists():
            tmp_path.unlink()


def smooth_series(values: np.ndarray, config: Mapping[str, Any]) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.size < 2:
        return values.copy()
    finite = np.isfinite(values)
    if not np.any(finite):
        return np.full_like(values, np.nan)
    filled = pd.Series(values).interpolate(limit_direction="both").to_numpy(dtype=float)
    method = str(config.get("method", "savgol")).lower()
    if method == "none":
        smoothed = filled
    elif method == "savgol":
        window = int(_required(config, "window_snapshots"))
        polyorder = int(config.get("polyorder", 2))
        if window % 2 == 0:
            raise ValueError("smoothing.window_snapshots must be odd for Savitzky-Golay")
        if window <= polyorder:
            raise ValueError("smoothing.window_snapshots must exceed smoothing.polyorder")
        if values.size < window:
            # Checkpoint files remain usable while a resumable run has not yet
            # accumulated the configured number of snapshots.  The requested
            # filter is applied automatically once enough rows exist.
            smoothed = filled
        else:
            smoothed = savgol_filter(filled, window_length=window, polyorder=polyorder, mode="interp")
    elif method == "gaussian":
        sigma = float(_required(config, "gaussian_sigma_snapshots"))
        if sigma <= 0.0:
            raise ValueError("smoothing.gaussian_sigma_snapshots must be positive")
        smoothed = gaussian_filter1d(filled, sigma=sigma, mode="nearest")
    elif method in {"rolling_mean", "rolling_median"}:
        window = int(_required(config, "window_snapshots"))
        if window <= 0:
            raise ValueError("smoothing.window_snapshots must be positive")
        series = pd.Series(filled)
        rolling = series.rolling(window=window, center=True, min_periods=1)
        smoothed = (
            rolling.mean().to_numpy(dtype=float)
            if method == "rolling_mean"
            else rolling.median().to_numpy(dtype=float)
        )
    else:
        raise ValueError(
            "smoothing.method must be one of none, savgol, gaussian, rolling_mean, rolling_median"
        )
    smoothed[~finite] = np.nan
    return smoothed


def actual_crossing_index(
    distance: np.ndarray,
    target: float,
    tolerance: float,
    branch: str,
) -> Optional[int]:
    distance = np.asarray(distance, dtype=float)
    good_indices = np.flatnonzero(np.isfinite(distance))
    crossings: list[int] = []
    for left, right in zip(good_indices[:-1], good_indices[1:]):
        a = distance[left] - target
        b = distance[right] - target
        if a == 0.0:
            crossings.append(int(left))
        elif b == 0.0 or a * b < 0.0:
            crossings.append(int(left if abs(a) <= abs(b) else right))
    if not crossings:
        return None
    branch_lower = branch.lower()
    index = crossings[0] if branch_lower in {"first", "first_crossing", "inbound"} else crossings[-1]
    return index if abs(distance[index] - target) <= tolerance else None


def reached_pericentre_index(radius: np.ndarray, config: Mapping[str, Any]) -> Optional[int]:
    if not bool(config.get("enabled", True)):
        return None
    radius = np.asarray(radius, dtype=float)
    finite = np.isfinite(radius)
    if np.count_nonzero(finite) < 5:
        return None
    minimum_post_points = int(config.get("minimum_post_points", 3))
    minimum_rise = float(config.get("minimum_rise_kpc", 0.5))
    last_allowed = radius.size - minimum_post_points
    candidates: list[int] = []
    for index in range(1, last_allowed):
        if not (np.isfinite(radius[index - 1]) and np.isfinite(radius[index]) and np.isfinite(radius[index + 1])):
            continue
        if radius[index] <= radius[index - 1] and radius[index] < radius[index + 1]:
            post = radius[index + 1 : index + 1 + minimum_post_points]
            if post.size == minimum_post_points and np.all(np.isfinite(post)):
                if np.nanmax(post) - radius[index] >= minimum_rise:
                    candidates.append(index)
    return min(candidates, key=lambda item: radius[item]) if candidates else None


def add_derived_columns(frame: pd.DataFrame, config: Mapping[str, Any]) -> pd.DataFrame:
    frame = frame.sort_values(["time_gyr", "snapshot"]).drop_duplicates("snapshot", keep="last").reset_index(drop=True)
    time = frame["time_gyr"].to_numpy(dtype=float)
    if time.size > 1 and np.any(np.diff(time) <= 0.0):
        raise ValueError("Simulation times must be strictly increasing after sorting")
    gas = frame["gas_mass_msun"].to_numpy(dtype=float)
    smoothed = smooth_series(gas, dict(_nested(config, "smoothing", {})))
    frame["gas_mass_smoothed_msun"] = smoothed
    if time.size > 1:
        raw_derivative = np.gradient(gas, time)
        smooth_derivative = np.gradient(smoothed, time)
    else:
        raw_derivative = np.full_like(gas, np.nan)
        smooth_derivative = np.full_like(gas, np.nan)
    frame["dgas_dt_raw_msun_per_gyr"] = raw_derivative
    frame["dgas_dt_smoothed_msun_per_gyr"] = smooth_derivative
    minimum_rate = float(_nested(config, "smoothing.minimum_abs_rate_msun_per_gyr", 0.0))
    raw_valid = np.isfinite(raw_derivative) & (np.abs(raw_derivative) > minimum_rate)
    smooth_valid = np.isfinite(smooth_derivative) & (np.abs(smooth_derivative) > minimum_rate)
    frame["tau_gas_raw_gyr"] = np.divide(
        gas,
        np.abs(raw_derivative),
        out=np.full_like(gas, np.nan),
        where=raw_valid,
    )
    frame["tau_gas_smoothed_gyr"] = np.divide(
        smoothed,
        np.abs(smooth_derivative),
        out=np.full_like(gas, np.nan),
        where=smooth_valid,
    )
    re = frame["re_major_kpc"].to_numpy(dtype=float)
    sigma = frame["sigma_los_kms"].to_numpy(dtype=float)
    tdyn = np.divide(
        KPC_PER_KMS_TO_GYR * re,
        sigma,
        out=np.full_like(re, np.nan),
        where=np.isfinite(re) & np.isfinite(sigma) & (sigma > 0.0),
    )
    frame["stellar_dynamical_time_gyr"] = tdyn
    frame["tau_gas_over_tdyn"] = np.divide(
        frame["tau_gas_smoothed_gyr"].to_numpy(dtype=float),
        tdyn,
        out=np.full_like(tdyn, np.nan),
        where=np.isfinite(tdyn) & (tdyn > 0.0),
    )

    frame["is_interaction_start"] = False
    frame["is_comparison_epoch"] = False
    frame["is_pericentre"] = False
    if len(frame):
        configured_start = _nested(config, "comparison_epoch.interaction_start_time_gyr", None)
        start_index = 0 if configured_start is None else int(np.nanargmin(np.abs(time - float(configured_start))))
        frame.loc[start_index, "is_interaction_start"] = True

    comparison = dict(_nested(config, "comparison_epoch", {}))
    method = str(comparison.get("method", "heliocentric_distance")).lower()
    comparison_index: Optional[int] = None
    if method == "heliocentric_distance":
        comparison_index = actual_crossing_index(
            frame["distance_heliocentric_kpc"].to_numpy(dtype=float),
            float(_required(comparison, "target_heliocentric_distance_kpc")),
            float(comparison.get("tolerance_kpc", 0.25)),
            str(comparison.get("branch", "first_crossing")),
        )
    elif method == "time":
        comparison_index = int(np.nanargmin(np.abs(time - float(_required(comparison, "time_gyr")))))
    elif method != "none":
        raise ValueError("comparison_epoch.method must be heliocentric_distance, time, or none")
    if comparison_index is not None:
        frame.loc[comparison_index, "is_comparison_epoch"] = True

    pericentre_index = reached_pericentre_index(
        frame["distance_galactocentric_kpc"].to_numpy(dtype=float),
        dict(_nested(config, "pericentre_detection", {})),
    )
    if pericentre_index is not None:
        frame.loc[pericentre_index, "is_pericentre"] = True
    frame["derivation_config_sha256"] = derivation_config_hash(config)
    return frame


def save_metadata(
    frame: pd.DataFrame,
    config: Mapping[str, Any],
    config_path: Path,
    paths: Mapping[str, Path],
) -> None:
    def marker(column: str) -> Optional[dict[str, float]]:
        selected = frame.loc[frame[column].astype(bool)]
        if selected.empty:
            return None
        row = selected.iloc[0]
        return {
            "snapshot": int(row["snapshot"]),
            "time_gyr": float(row["time_gyr"]),
            "distance_galactocentric_kpc": float(row["distance_galactocentric_kpc"]),
            "distance_heliocentric_kpc": float(row["distance_heliocentric_kpc"]),
        }

    metadata = {
        "schema_version": 1,
        "analysis_config_sha256": analysis_config_hash(config),
        "derivation_config_sha256": derivation_config_hash(config),
        "config_file": str(config_path),
        "run_dir": str(paths["run_dir"]),
        "timeseries_csv": str(paths["csv"]),
        "row_count": int(len(frame)),
        "snapshot_min": int(frame["snapshot"].min()) if len(frame) else None,
        "snapshot_max": int(frame["snapshot"].max()) if len(frame) else None,
        "definitions": {
            "centre": "dSph_workbench shrinking 3D stellar centre after standard MW centring and dwarf classification",
            "sigma_los": "old stars; circular half-light aperture; planar LOS velocity gradient fitted and removed",
            "re_major": "projected old-star semi-major-axis half-light radius",
            "enclosed_masses": "three-dimensional spherical apertures centred on the adopted dwarf centre",
            "hi_particle": "sum(mass * neutral fraction) below the configured temperature threshold",
            "hi_contour": "adaptive projected H I map integrated above the configured fixed N_HI contour",
            "tidal_proxy": "G M_MW(<R_GC) / R_GC^3",
            "gas_timescale": "absolute M_gas / (dM_gas/dt) after the configured smoothing",
            "stellar_dynamical_time": "0.9777922217 Gyr * R_e[kpc] / sigma_los[km/s]",
        },
        "events": {
            "interaction_start": marker("is_interaction_start"),
            "comparison_epoch": marker("is_comparison_epoch"),
            "pericentre": marker("is_pericentre"),
        },
        "config": config,
    }
    paths["metadata"].write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def extract(config: Mapping[str, Any], config_path: Path, overwrite: bool = False) -> pd.DataFrame:
    paths = resolve_paths(config, config_path)
    paths["output_dir"].mkdir(parents=True, exist_ok=True)
    discovered = discover_snapshots(paths["snapshot_dir"], config)
    config_hash = analysis_config_hash(config)
    hash_mismatch = False
    if paths["csv"].exists():
        existing = pd.read_csv(paths["csv"])
        if "analysis_config_sha256" in existing.columns and len(existing):
            hashes = set(existing["analysis_config_sha256"].dropna().astype(str))
            hash_mismatch = bool(hashes and hashes != {config_hash})
            if hash_mismatch:
                if not overwrite:
                    raise RuntimeError(
                        "Existing CSV was produced with a different analysis configuration. "
                        "Use --overwrite to rebuild it."
                    )
                # Rows made with different scientific choices must never be
                # mixed into one time series, even when the new range is only
                # a subset of the old range.
                existing = pd.DataFrame()
    else:
        existing = pd.DataFrame()
    existing_numbers = set(existing.get("snapshot", pd.Series(dtype=int)).astype(int)) if not overwrite else set()
    pending = [(number, path) for number, path in discovered if number not in existing_numbers]
    print(
        f"[evolution] discovered={len(discovered)} existing={len(existing_numbers)} pending={len(pending)}",
        flush=True,
    )
    frame = existing.copy()
    checkpoint_every = int(_nested(config, "processing.checkpoint_every", 1))
    for count, (number, path) in enumerate(pending, start=1):
        print(f"[evolution] snapshot {number}: {path.name}", flush=True)
        row = process_snapshot(number, path, config, config_hash)
        if overwrite and not frame.empty and "snapshot" in frame:
            frame = frame.loc[frame["snapshot"].astype(int) != number]
        frame = pd.concat([frame, pd.DataFrame([row])], ignore_index=True, sort=False)
        if count % checkpoint_every == 0 or count == len(pending):
            checkpoint = add_derived_columns(frame, config)
            atomic_write_csv(checkpoint, paths["csv"])
    if not pending:
        frame = add_derived_columns(frame, config)
        atomic_write_csv(frame, paths["csv"])
    else:
        frame = pd.read_csv(paths["csv"])
    save_metadata(frame, config, config_path, paths)
    return frame


def _normalise(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    finite = np.flatnonzero(np.isfinite(values) & (values > 0.0))
    if finite.size == 0:
        return np.full_like(values, np.nan)
    return values / values[finite[0]]


def draw_event_lines(axes: Sequence[plt.Axes], frame: pd.DataFrame) -> tuple[list[Any], list[str]]:
    events = [
        ("is_interaction_start", "MW interaction", "#303030", (0, (1.2, 2.0))),
        ("is_comparison_epoch", "comparison epoch", "#777777", (0, (4.0, 2.8))),
        ("is_pericentre", "pericentre", "#9b4b3f", (0, (2.0, 2.0))),
    ]
    handles: list[Any] = []
    labels: list[str] = []
    for column, label, color, linestyle in events:
        if column not in frame:
            continue
        selected = frame.loc[frame[column].astype(bool)]
        if selected.empty:
            continue
        time = float(selected.iloc[0]["time_gyr"])
        for axis in axes:
            axis.axvline(time, color=color, lw=0.85, ls=linestyle, zorder=1)
        handles.append(axes[0].plot([], [], color=color, lw=0.9, ls=linestyle)[0])
        labels.append(label)
    return handles, labels


def plot_timeseries(config: Mapping[str, Any], config_path: Path) -> list[Path]:
    paths = resolve_paths(config, config_path)
    if not paths["csv"].exists():
        raise FileNotFoundError(paths["csv"])
    frame = pd.read_csv(paths["csv"]).sort_values("time_gyr")
    time = frame["time_gyr"].to_numpy(dtype=float)

    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 8.0,
            "axes.labelsize": 8.2,
            "axes.titlesize": 8.4,
            "legend.fontsize": 7.0,
            "xtick.labelsize": 7.2,
            "ytick.labelsize": 7.2,
            "axes.linewidth": 0.75,
        }
    )
    blue = "#2f6f9f"
    teal = "#2a9d8f"
    amber = "#d08c32"
    purple = "#7b3294"
    red = "#b6534c"
    fig, axes = plt.subplots(4, 1, figsize=(6.9, 7.9), sharex=True)

    ax = axes[0]
    orbit_line = ax.plot(time, frame["distance_galactocentric_kpc"], color=blue, lw=1.55, label=r"$R_{\rm GC}$")[0]
    ax.set_ylabel(r"$R_{\rm GC}$ (kpc)")
    ax2 = ax.twinx()
    pressure = frame["ram_pressure_dyn_cm2"].to_numpy(dtype=float)
    pressure_line = ax2.plot(time, pressure, color=red, lw=1.35, label=r"$P_{\rm ram}$")[0]
    if np.any(np.isfinite(pressure) & (pressure > 0.0)):
        ax2.set_yscale("log")
    ax2.set_ylabel(r"$P_{\rm ram}$ (dyn cm$^{-2}$)")
    ax.legend(handles=[orbit_line, pressure_line], loc="best", frameon=False, ncol=2)
    ax.set_title("Environment and orbit", loc="left", fontweight="semibold")

    ax = axes[1]
    gas_line = ax.plot(time, _normalise(frame["gas_mass_msun"]), color=teal, lw=1.55, label=r"$M_{\rm gas}/M_{\rm gas,0}$")[0]
    hi_line = ax.plot(time, _normalise(frame["hi_mass_msun"]), color=amber, lw=1.45, label=r"$M_{\rm H\,I}/M_{\rm H\,I,0}$")[0]
    ax.set_ylabel("Retained fraction")
    ax.set_ylim(bottom=0.0)
    ax.legend(handles=[gas_line, hi_line], loc="best", frameon=False, ncol=2)
    ax.set_title("Gas evolution", loc="left", fontweight="semibold")

    ax = axes[2]
    re_line = ax.plot(time, frame["re_major_kpc"], color=blue, lw=1.55, label=r"$R_e$")[0]
    ax.set_ylabel(r"$R_e$ (kpc)")
    ax2 = ax.twinx()
    fraction_line = ax2.plot(
        time,
        frame["gas_fraction_3d_lt_re"],
        color=purple,
        lw=1.35,
        label=r"$f_{\rm gas}(<R_e)$",
    )[0]
    ax2.set_ylabel(r"$f_{\rm gas}(<R_e)$")
    ax2.set_ylim(0.0, 1.0)
    ax.legend(handles=[re_line, fraction_line], loc="best", frameon=False, ncol=2)
    ax.set_title("Stellar structural response", loc="left", fontweight="semibold")

    ax = axes[3]
    ax.plot(time, frame["sigma_los_kms"], color=blue, lw=1.55, label=r"$\sigma_{\rm los}$")
    ax.set_ylabel(r"$\sigma_{\rm los}$ (km s$^{-1}$)")
    ax.set_xlabel("Simulation time (Gyr)")
    ax.legend(loc="best", frameon=False)
    ax.set_title("Stellar kinematics", loc="left", fontweight="semibold")

    event_handles, event_labels = draw_event_lines(axes, frame)
    axes[0].legend(
        [orbit_line, pressure_line] + event_handles,
        [r"$R_{\rm GC}$", r"$P_{\rm ram}$"] + event_labels,
        loc="best",
        frameon=False,
        ncol=2,
    )
    for axis in axes:
        axis.grid(color="#d8d8d8", lw=0.45, alpha=0.55)
        axis.tick_params(direction="in", top=True, right=False)
    fig.subplots_adjust(left=0.105, right=0.885, bottom=0.075, top=0.975, hspace=0.25)

    formats = [str(value).lower() for value in _nested(config, "plot.formats", ["pdf", "png"])]
    dpi = int(_nested(config, "plot.dpi", 350))
    outputs = []
    for extension in formats:
        output = paths["figure_stem"].with_suffix(f".{extension}")
        fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.03)
        outputs.append(output)
        print(output, flush=True)
    plt.close(fig)
    return outputs


def derive_only(config: Mapping[str, Any], config_path: Path) -> pd.DataFrame:
    paths = resolve_paths(config, config_path)
    frame = pd.read_csv(paths["csv"])
    frame = add_derived_columns(frame, config)
    atomic_write_csv(frame, paths["csv"])
    save_metadata(frame, config, config_path, paths)
    return frame


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("run", "extract", "derive", "plot"):
        subparser = subparsers.add_parser(command)
        subparser.add_argument("--config", type=Path, required=True)
        if command in {"run", "extract"}:
            subparser.add_argument(
                "--overwrite",
                action="store_true",
                help="Reprocess selected snapshots and permit a changed analysis configuration",
            )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    config, config_path = load_config(args.config)
    if args.command in {"run", "extract"}:
        extract(config, config_path, overwrite=bool(args.overwrite))
    elif args.command == "derive":
        derive_only(config, config_path)
    if args.command in {"run", "plot"}:
        plot_timeseries(config, config_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
