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
        for key in (
            "smoothing",
            "diagnostics",
            "comparison_epoch",
            "pericentre_detection",
        )
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


def projected_elliptical_aperture_mask(
    x_kpc: np.ndarray,
    y_kpc: np.ndarray,
    re_major_kpc: float,
    re_multiple: float,
    axis_ratio: float,
    pa_rad: float,
    center_x_kpc: float,
    center_y_kpc: float,
) -> np.ndarray:
    """Select a projected ellipse expressed as a multiple of old-star Re."""
    if not np.isfinite(re_major_kpc) or re_major_kpc <= 0.0:
        return np.zeros_like(np.asarray(x_kpc, dtype=float), dtype=bool)
    return Analysis.get_elliptical_radial_mask(
        x_kpc,
        y_kpc,
        float(re_multiple) * float(re_major_kpc),
        ep=1.0 - float(axis_ratio),
        pa=float(pa_rad),
        center_x=float(center_x_kpc),
        center_y=float(center_y_kpc),
    )


def projected_dwarf_star_coordinates(
    snapshot: Mapping[str, Any],
    projection_config: Mapping[str, Any],
    star_types: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return global indices and projected coordinates for all dwarf stars."""
    df = snapshot["df"]
    dwarf_mask = np.asarray(snapshot["total_dw_star_mask"], dtype=bool)
    dwarf_indices = np.flatnonzero(dwarf_mask)
    dwarf_types = df.loc[dwarf_indices, "tp"].to_numpy(dtype=int)
    star_local = np.isin(dwarf_types, np.asarray(star_types, dtype=int))
    star_indices = dwarf_indices[star_local]

    mode = str(projection_config.get("mode", "native")).lower()
    if mode == "native":
        x_kpc = np.asarray(snapshot["x_kpc"], dtype=float)[star_local]
        y_kpc = np.asarray(snapshot["y_kpc"], dtype=float)[star_local]
    elif mode == "euler":
        center = np.array(
            [snapshot["dw_xc"], snapshot["dw_yc"], snapshot["dw_zc"]], dtype=float
        )
        positions = df.loc[star_indices, ["x", "y", "z"]].to_numpy(dtype=float) - center
        x_axis, y_axis, _ = rotation_basis(
            float(_required(projection_config, "inclination_deg")),
            float(_required(projection_config, "azimuth_deg")),
        )
        x_kpc = positions @ x_axis
        y_kpc = positions @ y_axis
    else:
        raise ValueError("projection.mode must be 'native' or 'euler'")
    return star_indices, x_kpc, y_kpc


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


def prepare_analysis_snapshot(
    number: int,
    snapshot_path: Path,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    selection = dict(_nested(config, "dwarf_selection", {}))
    gas_config = dict(_nested(config, "gas", {}))
    return prepare_snapshot_context(
        folder_path=str(snapshot_path.parent),
        snapshot_num=number,
        core_radius=float(_required(selection, "core_radius_kpc")),
        r_exclude=float(selection.get("mw_exclusion_radius_kpc", 5.0)),
        dwarf_radius_factor=float(selection.get("dwarf_radius_factor", 3.0)),
        k_density=int(selection.get("density_neighbours", 16)),
        dwarf_gas_radius=float(gas_config.get("aperture_kpc", 20.0)),
        gas_temperature_split=float(gas_config.get("temperature_split_k", 2.0e4)),
        include_mw_gas=True,
        mw_gas_radius=float(_nested(config, "cgm.mw_gas_radius_kpc", 500.0)),
        include_dark_matter=True,
        include_star_birth=True,
    )


def read_gas_particle_ids(snapshot_path: Path, gas_count: int) -> np.ndarray:
    with h5py.File(snapshot_path, "r") as handle:
        if "PartType0" not in handle or "ParticleIDs" not in handle["PartType0"]:
            raise RuntimeError(f"PartType0/ParticleIDs is required in {snapshot_path}")
        values = np.asarray(handle["PartType0"]["ParticleIDs"], dtype=np.uint64)
    if values.size != gas_count:
        raise RuntimeError(
            f"Gas ParticleID count {values.size} does not match loaded gas count {gas_count}"
        )
    return values


def find_snapshot_path(
    snapshot_dir: Path,
    config: Mapping[str, Any],
    number: int,
) -> Path:
    glob_pattern = str(_nested(config, "snapshots.glob", "snapshot_*.hdf5"))
    number_regex = re.compile(
        str(_nested(config, "snapshots.number_regex", r"snapshot_(\d+)\.hdf5$"))
    )
    for path in sorted(snapshot_dir.glob(glob_pattern)):
        match = number_regex.search(path.name)
        if match is not None and int(match.group(1)) == int(number):
            return path.resolve()
    raise FileNotFoundError(f"Reference snapshot {number} was not found in {snapshot_dir}")


def initial_dwarf_gas_membership(
    snapshot_dir: Path,
    config: Mapping[str, Any],
) -> tuple[np.ndarray, dict[str, Any]]:
    exclusion = dict(_nested(config, "cgm.dwarf_gas_exclusion", {}))
    method = str(exclusion.get("method", "none")).lower()
    if method == "none":
        return np.empty(0, dtype=np.uint64), {
            "method": "none",
            "reference_snapshot": None,
            "reference_radius_kpc": None,
            "particle_count": 0,
            "particle_ids_sha256": None,
        }
    if method != "initial_particle_ids":
        raise ValueError(
            "cgm.dwarf_gas_exclusion.method must be 'none' or 'initial_particle_ids'"
        )

    reference_number = int(exclusion.get("reference_snapshot", 0))
    reference_path = find_snapshot_path(snapshot_dir, config, reference_number)
    snapshot = prepare_analysis_snapshot(reference_number, reference_path, config)
    df = snapshot["df"]
    gas_type = int(_nested(config, "particle_types.gas", 0))
    particle_types = df["tp"].to_numpy(dtype=int)
    gas_global = np.flatnonzero(particle_types == gas_type)
    gas_ids = read_gas_particle_ids(reference_path, gas_global.size)
    ids_global = np.zeros(len(df), dtype=np.uint64)
    ids_global[gas_global] = gas_ids

    center = np.array(
        [snapshot["dw_xc"], snapshot["dw_yc"], snapshot["dw_zc"]], dtype=float
    )
    positions = df[["x", "y", "z"]].to_numpy(dtype=float)
    radius = np.linalg.norm(positions - center, axis=1)
    reference_radius = float(
        exclusion.get(
            "reference_radius_kpc",
            _nested(config, "gas.aperture_kpc", 20.0),
        )
    )
    membership_mask = (
        (particle_types == gas_type)
        & np.isfinite(radius)
        & (radius <= reference_radius)
    )
    membership_ids = np.unique(ids_global[membership_mask])
    membership_ids.sort()
    fingerprint = hashlib.sha256(
        np.ascontiguousarray(membership_ids).tobytes()
    ).hexdigest()
    snapshot["simulation"].df = None
    return membership_ids, {
        "method": method,
        "reference_snapshot": reference_number,
        "reference_radius_kpc": reference_radius,
        "particle_count": int(membership_ids.size),
        "particle_ids_sha256": fingerprint,
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
    excluded_dwarf_gas_mask: Optional[np.ndarray] = None,
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
    shell_candidates = (
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
    if excluded_dwarf_gas_mask is None:
        excluded_dwarf_gas_mask = np.zeros(len(df), dtype=bool)
    else:
        excluded_dwarf_gas_mask = np.asarray(excluded_dwarf_gas_mask, dtype=bool)
        if excluded_dwarf_gas_mask.size != len(df):
            raise ValueError("excluded_dwarf_gas_mask must have one value per dataframe row")
    excluded_in_shell = int(np.count_nonzero(shell_candidates & excluded_dwarf_gas_mask))
    candidates = shell_candidates & ~excluded_dwarf_gas_mask
    indices = np.flatnonzero(candidates)
    method = str(config.get("method", "knn_shell")).lower()
    if method not in {"knn_shell", "fixed_shell"}:
        raise ValueError("cgm.method must be 'knn_shell' or 'fixed_shell'")
    if method == "knn_shell" and indices.size:
        neighbours = int(config.get("neighbours", 128))
        order = np.argsort(radius[indices])
        indices = indices[order[: min(neighbours, indices.size)]]

    minimum_particles = int(config.get("minimum_particles", 32))
    if indices.size < minimum_particles:
        return {
            "cgm_particle_count": float(indices.size),
            "cgm_excluded_dwarf_gas_particle_count": float(excluded_in_shell),
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
        "cgm_excluded_dwarf_gas_particle_count": float(excluded_in_shell),
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
    excluded_cgm_gas_ids: Optional[np.ndarray] = None,
    cgm_membership_info: Optional[Mapping[str, Any]] = None,
) -> dict[str, Any]:
    snapshot = prepare_analysis_snapshot(number, snapshot_path, config)
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
    gas_global = np.flatnonzero(particle_types == gas_type)
    gas_ids = read_gas_particle_ids(snapshot_path, gas_global.size)
    excluded_global = np.zeros(len(df), dtype=bool)
    if excluded_cgm_gas_ids is not None and len(excluded_cgm_gas_ids):
        excluded_global[gas_global] = np.isin(
            gas_ids,
            np.asarray(excluded_cgm_gas_ids, dtype=np.uint64),
            assume_unique=False,
        )
    star_mask = np.asarray(snapshot["total_dw_star_mask"], dtype=bool) & np.isin(particle_types, star_types)
    projection_config = dict(_nested(config, "projection", {"mode": "native"}))
    star_indices, star_x_kpc, star_y_kpc = projected_dwarf_star_coordinates(
        snapshot, projection_config, star_types
    )
    mass_aperture_re_multiple = float(
        _nested(config, "stellar.mass_aperture_re_multiple", 8.0)
    )
    sensitivity_re_multiple = float(
        _nested(config, "stellar.sensitivity_aperture_re_multiple", 10.0)
    )
    aperture_local = projected_elliptical_aperture_mask(
        star_x_kpc,
        star_y_kpc,
        projection["re_major_kpc"],
        mass_aperture_re_multiple,
        projection["axis_ratio"],
        projection["pa_rad"],
        projection["shape_center_x_kpc"],
        projection["shape_center_y_kpc"],
    )
    sensitivity_local = projected_elliptical_aperture_mask(
        star_x_kpc,
        star_y_kpc,
        projection["re_major_kpc"],
        sensitivity_re_multiple,
        projection["axis_ratio"],
        projection["pa_rad"],
        projection["shape_center_x_kpc"],
        projection["shape_center_y_kpc"],
    )
    star_mass_mask = np.zeros(len(df), dtype=bool)
    star_mass_mask[star_indices[aperture_local]] = True
    star_sensitivity_mask = np.zeros(len(df), dtype=bool)
    star_sensitivity_mask[star_indices[sensitivity_local]] = True

    if "birth" in df.columns:
        old_star_mask = star_mask & (df["birth"].to_numpy(dtype=float) == 0.0)
    else:
        new_star_type = int(_nested(config, "particle_types.new_stars", 4))
        old_star_mask = star_mask & (particle_types != new_star_type)
    new_star_mask = star_mask & ~old_star_mask

    stellar_mass = float(np.sum(masses[star_mass_mask]))
    stellar_mass_old = float(np.sum(masses[star_mass_mask & old_star_mask]))
    stellar_mass_new = float(np.sum(masses[star_mass_mask & new_star_mask]))
    stellar_mass_sensitivity = float(np.sum(masses[star_sensitivity_mask]))

    preselection_radius = float(
        _nested(
            config,
            "stellar.preselection_radius_kpc",
            _nested(config, "stellar.mass_aperture_kpc", 20.0),
        )
    )
    star_preselection_mask = star_mask & (distance_from_dwarf <= preselection_radius)
    stellar_mass_preselection = float(np.sum(masses[star_preselection_mask]))

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
        excluded_dwarf_gas_mask=excluded_global,
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
        "stellar_mass_old_tracer_aperture_msun": stellar_mass_old,
        "stellar_mass_new_aperture_msun": stellar_mass_new,
        "stellar_mass_sensitivity_aperture_msun": stellar_mass_sensitivity,
        "stellar_mass_preselection_msun": stellar_mass_preselection,
        "stellar_mass_aperture_re_multiple": mass_aperture_re_multiple,
        "stellar_mass_sensitivity_re_multiple": sensitivity_re_multiple,
        "stellar_mass_aperture_major_kpc": (
            mass_aperture_re_multiple * float(projection["re_major_kpc"])
        ),
        "stellar_mass_preselection_radius_kpc": preselection_radius,
        # Deprecated compatibility alias.  It is explicitly an old-star
        # tracer mass, never the total gravitating stellar mass.
        "stellar_mass_old_observational_aperture_msun": stellar_mass_old,
        "gas_mass_msun": total_gas_mass,
        "hi_mass_msun": float(hi_mass),
        "hi_fraction_total_stars": gas_fraction(stellar_mass, float(hi_mass)),
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
        "stellar_particle_count_old_tracer": int(
            np.count_nonzero(star_mass_mask & old_star_mask)
        ),
        "stellar_particle_count_new": int(np.count_nonzero(star_mass_mask & new_star_mask)),
        "stellar_particle_count_preselection": int(np.count_nonzero(star_preselection_mask)),
        "gas_particle_count": int(np.count_nonzero(dwarf_gas_mask)),
        "hi_particle_count": int(np.count_nonzero(hi_mask)),
        "cgm_dwarf_gas_membership_count": int(
            (cgm_membership_info or {}).get("particle_count", 0)
        ),
        "cgm_dwarf_gas_membership_sha256": (
            (cgm_membership_info or {}).get("particle_ids_sha256") or ""
        ),
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
        smoothed = filled.copy()
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


def add_central_gas_timescale_columns(
    frame: pd.DataFrame,
    time: np.ndarray,
    tdyn: np.ndarray,
    config: Mapping[str, Any],
) -> pd.DataFrame:
    """Add central enclosed-gas timescales without reopening snapshots.

    The raw enclosed masses are produced by :func:`process_snapshot` using the
    adopted three-dimensional dwarf centre.  Fixed physical apertures and the
    evolving projected semi-major-axis effective radius are intentionally kept
    as separate series because the latter changes the enclosing volume.
    """
    default_smoothing = dict(_nested(config, "smoothing", {}))
    central_smoothing = dict(
        _nested(config, "smoothing.central_gas", default_smoothing)
    )
    minimum_rate = float(
        central_smoothing.get(
            "minimum_abs_rate_msun_per_gyr",
            _nested(config, "smoothing.minimum_abs_rate_msun_per_gyr", 0.0),
        )
    )
    aperture_specs: list[tuple[str, Optional[float]]] = []
    for radius_kpc in _nested(config, "enclosed_mass.radii_kpc", [0.5, 1.0, 2.0]):
        radius = float(radius_kpc)
        if np.isclose(radius, 0.5) or np.isclose(radius, 1.0):
            aperture_specs.append((aperture_label(radius), radius))
    if bool(_nested(config, "enclosed_mass.include_effective_radius", True)):
        aperture_specs.append(("re", None))

    for label, fixed_radius_kpc in aperture_specs:
        gas_column = f"gas_mass_3d_lt_{label}_msun"
        stellar_column = f"stellar_mass_3d_lt_{label}_msun"
        if gas_column not in frame or stellar_column not in frame:
            continue
        gas = frame[gas_column].to_numpy(dtype=float)
        stellar = frame[stellar_column].to_numpy(dtype=float)
        smoothed_unclipped = smooth_series(gas, central_smoothing)
        smoothed = np.maximum(smoothed_unclipped, 0.0)
        raw_derivative = (
            np.gradient(gas, time) if time.size > 1 else np.full_like(gas, np.nan)
        )
        smooth_derivative = (
            np.gradient(smoothed, time)
            if time.size > 1
            else np.full_like(gas, np.nan)
        )
        raw_valid = (
            np.isfinite(gas)
            & (gas > 0.0)
            & np.isfinite(raw_derivative)
            & (np.abs(raw_derivative) > minimum_rate)
        )
        smooth_valid = (
            np.isfinite(smoothed)
            & (smoothed > 0.0)
            & np.isfinite(smooth_derivative)
            & (np.abs(smooth_derivative) > minimum_rate)
        )
        tau_raw = np.divide(
            gas,
            np.abs(raw_derivative),
            out=np.full_like(gas, np.nan),
            where=raw_valid,
        )
        tau_smoothed = np.divide(
            smoothed,
            np.abs(smooth_derivative),
            out=np.full_like(gas, np.nan),
            where=smooth_valid,
        )
        frame[f"gas_mass_3d_lt_{label}_smoothed_unclipped_msun"] = smoothed_unclipped
        frame[f"gas_mass_3d_lt_{label}_smoothed_msun"] = smoothed
        frame[f"dgas_3d_lt_{label}_dt_raw_msun_per_gyr"] = raw_derivative
        frame[f"dgas_3d_lt_{label}_dt_smoothed_msun_per_gyr"] = smooth_derivative
        frame[f"tau_gas_3d_lt_{label}_raw_gyr"] = tau_raw
        frame[f"tau_gas_3d_lt_{label}_smoothed_gyr"] = tau_smoothed
        frame[f"tau_gas_3d_lt_{label}_over_tdyn"] = np.divide(
            tau_smoothed,
            tdyn,
            out=np.full_like(tdyn, np.nan),
            where=np.isfinite(tdyn) & (tdyn > 0.0),
        )

        if fixed_radius_kpc is not None:
            factor = G_KPC_KMS2_PER_MSUN / fixed_radius_kpc
            gas_nonnegative = np.where(np.isfinite(gas) & (gas >= 0.0), gas, np.nan)
            stellar_nonnegative = np.where(
                np.isfinite(stellar) & (stellar >= 0.0), stellar, np.nan
            )
            frame[f"vcirc_gas_3d_lt_{label}_kms"] = np.sqrt(factor * gas_nonnegative)
            frame[f"vcirc_stellar_3d_lt_{label}_kms"] = np.sqrt(
                factor * stellar_nonnegative
            )
            frame[f"vcirc_baryon_3d_lt_{label}_kms"] = np.sqrt(
                factor * (gas_nonnegative + stellar_nonnegative)
            )
    return frame


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


def duration_in_value_interval(
    time: np.ndarray,
    values: np.ndarray,
    lower: float,
    upper: float,
) -> float:
    """Integrate time within [lower, upper] using linear interpolation."""
    time = np.asarray(time, dtype=float)
    values = np.asarray(values, dtype=float)
    if lower > upper:
        raise ValueError("interval lower bound must not exceed upper bound")
    if time.size != values.size:
        raise ValueError("time and values must have the same length")
    duration = 0.0
    for t0, t1, y0, y1 in zip(time[:-1], time[1:], values[:-1], values[1:]):
        if not np.all(np.isfinite([t0, t1, y0, y1])) or t1 <= t0:
            continue
        breaks = [0.0, 1.0]
        if y1 != y0:
            for boundary in (lower, upper):
                fraction = (boundary - y0) / (y1 - y0)
                if 0.0 < fraction < 1.0:
                    breaks.append(float(fraction))
        breaks = sorted(set(breaks))
        for left, right in zip(breaks[:-1], breaks[1:]):
            midpoint = 0.5 * (left + right)
            value = y0 + midpoint * (y1 - y0)
            if lower <= value <= upper:
                duration += (right - left) * (t1 - t0)
    return float(duration)


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
    frame = add_central_gas_timescale_columns(frame, time, tdyn, config)
    pressure = frame["ram_pressure_dyn_cm2"].to_numpy(dtype=float)
    pressure_smoothing = dict(
        _nested(config, "smoothing.ram_pressure", {"method": "none"})
    )
    frame["ram_pressure_smoothed_dyn_cm2"] = smooth_series(
        pressure, pressure_smoothing
    )
    sigma_interval = dict(
        _nested(
            config,
            "diagnostics.sigma_los_interval",
            {"lower_kms": 9.0, "upper_kms": 11.0},
        )
    )
    sigma_lower = float(sigma_interval.get("lower_kms", 9.0))
    sigma_upper = float(sigma_interval.get("upper_kms", 11.0))
    if sigma_lower > sigma_upper:
        raise ValueError("diagnostics.sigma_los_interval lower_kms must not exceed upper_kms")
    frame["sigma_in_fornax_like_interval"] = (
        np.isfinite(sigma) & (sigma >= sigma_lower) & (sigma <= sigma_upper)
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

    sigma_interval = dict(
        _nested(
            config,
            "diagnostics.sigma_los_interval",
            {"lower_kms": 9.0, "upper_kms": 11.0},
        )
    )
    sigma_lower = float(sigma_interval.get("lower_kms", 9.0))
    sigma_upper = float(sigma_interval.get("upper_kms", 11.0))
    sigma_duration = duration_in_value_interval(
        frame["time_gyr"].to_numpy(dtype=float),
        frame["sigma_los_kms"].to_numpy(dtype=float),
        sigma_lower,
        sigma_upper,
    )
    membership_count = (
        int(frame["cgm_dwarf_gas_membership_count"].iloc[0])
        if len(frame) and "cgm_dwarf_gas_membership_count" in frame
        else 0
    )
    membership_hash = (
        str(frame["cgm_dwarf_gas_membership_sha256"].iloc[0])
        if len(frame) and "cgm_dwarf_gas_membership_sha256" in frame
        else None
    )

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
            "hi_plot": "raw configured H I mass; no temporal smoothing",
            "local_cgm": "configured hot-gas shell after excluding all gas particles tagged as dwarf members in the reference snapshot",
            "ram_pressure_plot": "configured light smoothing of the raw rho_CGM * v_rel^2 series",
            "tidal_proxy": "G M_MW(<R_GC) / R_GC^3",
            "gas_timescale": "absolute M_gas / (dM_gas/dt) after the configured smoothing",
            "central_gas_timescale": "the same smoothed absolute mass-loss timescale applied separately to 3D gas masses within 0.5 kpc, 1.0 kpc, and the evolving R_e(t) aperture",
            "central_apertures": "0.5 and 1.0 kpc are fixed 3D spheres; R_e(t) is an evolving 3D sphere whose radius is the projected old-star semi-major-axis effective radius",
            "baryonic_circular_velocity": "sqrt(G M(<r) / r) for gas, stars, and their sum at the fixed 0.5 and 1.0 kpc apertures",
            "stellar_dynamical_time": "0.9777922217 Gyr * R_e[kpc] / sigma_los[km/s]",
        },
        "events": {
            "interaction_start": marker("is_interaction_start"),
            "comparison_epoch": marker("is_comparison_epoch"),
            "pericentre": marker("is_pericentre"),
        },
        "diagnostics": {
            "sigma_los_interval_kms": [sigma_lower, sigma_upper],
            "duration_in_sigma_los_interval_gyr": sigma_duration,
            "duration_method": "piecewise-linear interpolation between consecutive snapshots",
            "cgm_dwarf_gas_exclusion_particle_count": membership_count,
            "cgm_dwarf_gas_exclusion_particle_ids_sha256": membership_hash,
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
    excluded_cgm_gas_ids = np.empty(0, dtype=np.uint64)
    cgm_membership_info: dict[str, Any] = {}
    if pending:
        excluded_cgm_gas_ids, cgm_membership_info = initial_dwarf_gas_membership(
            paths["snapshot_dir"], config
        )
        if cgm_membership_info["method"] != "none":
            print(
                "[evolution] CGM exclusion: "
                f"{cgm_membership_info['particle_count']} initial dwarf-gas ParticleIDs "
                f"from snapshot {cgm_membership_info['reference_snapshot']}",
                flush=True,
            )
    frame = existing.copy()
    checkpoint_every = int(_nested(config, "processing.checkpoint_every", 1))
    for count, (number, path) in enumerate(pending, start=1):
        print(f"[evolution] snapshot {number}: {path.name}", flush=True)
        row = process_snapshot(
            number,
            path,
            config,
            config_hash,
            excluded_cgm_gas_ids=excluded_cgm_gas_ids,
            cgm_membership_info=cgm_membership_info,
        )
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


def draw_event_lines(axes: Sequence[plt.Axes], frame: pd.DataFrame) -> Optional[float]:
    events = [
        ("is_comparison_epoch", "comparison epoch", "#777777", (0, (4.0, 2.8))),
    ]
    comparison_time: Optional[float] = None
    for column, label, color, linestyle in events:
        if column not in frame:
            continue
        selected = frame.loc[frame[column].astype(bool)]
        if selected.empty:
            continue
        time = float(selected.iloc[0]["time_gyr"])
        for axis in axes:
            axis.axvline(time, color=color, lw=0.85, ls=linestyle, zorder=1)
        if column == "is_comparison_epoch":
            comparison_time = time
    return comparison_time


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
    blue = "#3b78a8"
    teal = "#4a9b8e"
    cyan = "#35b9c5"
    purple = "#80679b"
    red = "#c35b42"
    fig, axes = plt.subplots(4, 1, figsize=(6.9, 6.85), sharex=True)

    ax = axes[0]
    orbit_line = ax.plot(
        time,
        frame["distance_galactocentric_kpc"].to_numpy(dtype=float),
        color=blue,
        lw=1.65,
        label=r"$R_{\rm GC}$",
    )[0]
    ax.set_ylabel(r"$R_{\rm GC}$ (kpc)")
    ax2 = ax.twinx()
    pressure = frame["ram_pressure_smoothed_dyn_cm2"].to_numpy(dtype=float)
    pressure_line = ax2.plot(time, pressure, color=red, lw=1.45, label=r"$P_{\rm ram}$")[0]
    if np.any(np.isfinite(pressure) & (pressure > 0.0)):
        ax2.set_yscale("log")
    ax2.set_ylabel(r"$P_{\rm ram}$ (dyn cm$^{-2}$)")
    ax.legend(handles=[orbit_line, pressure_line], loc="lower left", frameon=False, ncol=2)
    ax.set_title("(a) Environment", loc="left", fontweight="semibold")

    ax = axes[1]
    gas_line = ax.plot(time, _normalise(frame["gas_mass_msun"]), color=teal, lw=1.65, label=r"$M_{\rm gas}/M_{\rm gas,0}$")[0]
    hi_line = ax.plot(time, _normalise(frame["hi_mass_msun"]), color=cyan, lw=1.55, label=r"$M_{\rm H\,I}/M_{\rm H\,I,0}$")[0]
    ax.set_ylabel("Normalized mass")
    ax.set_ylim(bottom=0.0)
    ax.legend(handles=[gas_line, hi_line], loc="lower left", frameon=False, ncol=2)
    ax.set_title("(b) Gas evolution", loc="left", fontweight="semibold")

    ax = axes[2]
    re_line = ax.plot(
        time,
        frame["re_major_kpc"].to_numpy(dtype=float),
        color=blue,
        lw=1.65,
        label=r"$R_e$",
    )[0]
    ax.set_ylabel(r"$R_e$ (kpc)")
    ax2 = ax.twinx()
    fraction_line = ax2.plot(
        time,
        frame["gas_fraction_3d_lt_re"].to_numpy(dtype=float),
        color=purple,
        lw=1.45,
        label=r"$f_{\rm gas}(<R_e)$",
    )[0]
    ax2.set_ylabel(r"$f_{\rm gas}(<R_e)$")
    ax2.set_ylim(0.0, 1.0)
    ax.legend(handles=[re_line, fraction_line], loc="upper left", frameon=False, ncol=2)
    ax.set_title("(c) Stellar structure", loc="left", fontweight="semibold")

    ax = axes[3]
    ax.plot(
        time,
        frame["sigma_los_kms"].to_numpy(dtype=float),
        color=blue,
        lw=1.65,
        label=r"$\sigma_{\rm los}$",
    )
    ax.set_ylabel(r"$\sigma_{\rm los}$ (km s$^{-1}$)")
    ax.set_xlabel("Simulation time (Gyr)")
    ax.legend(loc="upper right", frameon=False)
    ax.set_title("(d) Stellar kinematics", loc="left", fontweight="semibold")

    comparison_time = draw_event_lines(axes, frame)
    if comparison_time is not None:
        axes[0].text(
            comparison_time,
            0.96,
            "comparison epoch",
            transform=axes[0].get_xaxis_transform(),
            rotation=90,
            ha="right",
            va="top",
            color="#666666",
            fontsize=6.8,
        )
    for axis in axes:
        axis.tick_params(direction="in", top=True, right=False)
        axis.minorticks_on()
        axis.label_outer()
    fig.subplots_adjust(left=0.105, right=0.885, bottom=0.080, top=0.975, hspace=0.16)

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
