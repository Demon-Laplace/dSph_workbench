#!/usr/bin/env python3
"""Convert Fornax RGB-candidate counts to an equivalent V-band map.

This workflow treats the input candidate catalogue as star counts, calibrates
the background-subtracted density empirically against Muñoz et al. V-band
surface-brightness anchors, and writes reproducible maps, profiles, and
diagnostics.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from matplotlib.path import Path as MplPath
from scipy.ndimage import convolve1d
from scipy.spatial import ConvexHull


CODE_VERSION = "v0.2.3"

DEFAULT_INPUT = Path("fixtures_realistic/local_pc/Fornax/Sample1_candidates.csv")
DEFAULT_RADIAL_INPUT = Path("fixtures_realistic/local_pc/Fornax/Sample1_density.csv")
DEFAULT_OUTPUT_DIR = Path("sandbox_runs/yang_fornax_surface_brightness")
DEFAULT_PREPARED_OBSERVATION_PROFILE_NAME = "Fornax_surface_brightness_profile.csv"

FORNAX_CENTER_RA = "02h39m50.9s"
FORNAX_CENTER_DEC = "-34d30m54s"
FORNAX_ELLIPTICITY = 0.317
FORNAX_Q = 1.0 - FORNAX_ELLIPTICITY
FORNAX_PA_DEG = 47.3

MUNOZ_CENTER_MU_V = 23.59
MUNOZ_CENTER_RADIUS_KPC = 0.05
MUNOZ_RHALF_MU_V = 24.77
MUNOZ_RHALF_RADIUS_KPC = 0.791
MUNOZ_MU_V_ERR = 0.16
DEFAULT_DISTANCE_KPC = 139.6
DEFAULT_PROFILE_SELECTION_MAX_RADIUS_KPC = 6.0
DEFAULT_PROFILE_PLOT_MAX_RADIUS_KPC = 5.1
DEFAULT_MUNOZ_ANCHOR_HALF_WIDTH_KPC = 0.12
REFERENCE_RADIUS_DEG = 0.3
YANG_S1_REFERENCE_DENSITY = 11.901
YANG_S1_BACKGROUND_DENSITY = 0.421e-3
YANG_S1_1D_BACKGROUND_UNCERTAINTY = 0.024e-3
YANG_S1_1D_LIMIT_MAG_ARCSEC2 = 36.59
YANG_STYLE_DETECTION_SIGMA = 3.0

BACKGROUND_R_INNER_DEG = 6.0
BACKGROUND_R_OUTER_DEG = 10.0
SMOOTH_FWHM_DEG = 0.25
CENTRAL_PLOT_HALF_WIDTH_DEG = 3.0


@dataclass(frozen=True)
class FornaxGeometry:
    center: SkyCoord
    ellipticity: float = FORNAX_ELLIPTICITY
    q: float = FORNAX_Q
    pa_deg: float = FORNAX_PA_DEG

    @property
    def pa_rad(self) -> float:
        return float(np.deg2rad(self.pa_deg))


@dataclass
class GridMap:
    x_edges: np.ndarray
    y_edges: np.ndarray
    x_centers: np.ndarray
    y_centers: np.ndarray
    xx: np.ndarray
    yy: np.ndarray
    counts: np.ndarray
    area_arcmin2: np.ndarray
    footprint_fraction: np.ndarray
    r_ell: np.ndarray

    @property
    def pixel_area_arcmin2(self) -> float:
        return float(np.diff(self.x_edges).mean() * np.diff(self.y_edges).mean() * 3600.0)

    @property
    def pixel_size_deg(self) -> float:
        return float(np.diff(self.x_edges).mean())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert Fornax RGB-candidate counts into a Muñoz-calibrated "
            "equivalent V-band surface-brightness map."
        )
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Candidate catalogue CSV.")
    parser.add_argument(
        "--radial-input",
        type=Path,
        default=DEFAULT_RADIAL_INPUT,
        help="Optional existing 1D density profile used only for validation.",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--prepared-observation-output",
        type=Path,
        default=None,
        help=(
            "Preprocessed 1D observation CSV written for PlotFig.py. "
            "Defaults to the parent directory of the directory containing PlotFig.py."
        ),
    )
    parser.add_argument(
        "--pixel-size-deg",
        type=float,
        default=0.05,
        help="Regular map pixel size in degrees. Default samples the 0.25 deg FWHM kernel with 5 pixels.",
    )
    parser.add_argument(
        "--coverage-threshold",
        type=float,
        default=0.35,
        help="Minimum smoothed effective area fraction for reliable map cells.",
    )
    parser.add_argument(
        "--footprint-subpixels",
        type=int,
        default=5,
        help="Subpixel samples per axis for convex-hull footprint area fractions.",
    )
    parser.add_argument(
        "--distance-kpc",
        type=float,
        default=DEFAULT_DISTANCE_KPC,
        help="Distance used only for deg-to-kpc radial-profile conversion.",
    )
    parser.add_argument(
        "--profile-selection-max-kpc",
        "--profile-max-kpc",
        dest="profile_selection_max_kpc",
        type=float,
        default=DEFAULT_PROFILE_SELECTION_MAX_RADIUS_KPC,
        help="Maximum semimajor-axis radius saved in the radial profile.",
    )
    parser.add_argument(
        "--profile-plot-max-kpc",
        type=float,
        default=DEFAULT_PROFILE_PLOT_MAX_RADIUS_KPC,
        help="Maximum semimajor-axis radius shown in the radial-profile figure.",
    )
    parser.add_argument(
        "--munoz-anchor-half-width-kpc",
        type=float,
        default=DEFAULT_MUNOZ_ANCHOR_HALF_WIDTH_KPC,
        help="Half-width around the Muñoz half-light anchor used for map-density averaging.",
    )
    return parser.parse_args()


def default_prepared_observation_output() -> Path:
    plotfig_path = Path(__file__).resolve().with_name("PlotFig.py")
    return plotfig_path.parent.parent / DEFAULT_PREPARED_OBSERVATION_PROFILE_NAME


def load_and_validate_candidates(path: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    if not path.exists():
        raise FileNotFoundError(path)

    df = pd.read_csv(path)
    required = {"ra", "dec"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Candidate file is missing required coordinate columns: {missing}")

    finite = np.isfinite(df["ra"].to_numpy(dtype=float)) & np.isfinite(df["dec"].to_numpy(dtype=float))
    duplicate_source_count = None
    if "source_id" in df.columns:
        duplicate_source_count = int(df["source_id"].duplicated().sum())

    diagnostics = {
        "input_path": str(path),
        "input_kind": "individual_star_catalogue",
        "identified_columns": {
            "ra_deg": "ra",
            "dec_deg": "dec",
            "count_per_row": 1,
            "optional_projected_x_deg": "x" if "x" in df.columns else None,
            "optional_projected_y_deg": "y" if "y" in df.columns else None,
            "optional_elliptical_radius_deg": "rmajor" if "rmajor" in df.columns else None,
        },
        "coordinate_units": "ra/dec in degrees; projected x/y recomputed in tangent-plane degrees",
        "row_count": int(len(df)),
        "finite_ra_dec_rows": int(finite.sum()),
        "missing_ra_dec_rows": int((~finite).sum()),
        "duplicate_source_id_rows": duplicate_source_count,
    }
    return df.loc[finite].copy(), diagnostics


def project_to_tangent_plane(df: pd.DataFrame, geometry: FornaxGeometry) -> pd.DataFrame:
    coords = SkyCoord(df["ra"].to_numpy(dtype=float) * u.deg, df["dec"].to_numpy(dtype=float) * u.deg)
    offset_frame = SkyOffsetFrame(origin=geometry.center)
    offsets = coords.transform_to(offset_frame)

    out = df.copy()
    out["x_east_deg"] = offsets.lon.to(u.deg).value
    out["y_north_deg"] = offsets.lat.to(u.deg).value
    out["r_ell_deg"] = elliptical_radius(out["x_east_deg"], out["y_north_deg"], geometry)
    return out


def elliptical_components(x_deg: np.ndarray, y_deg: np.ndarray, geometry: FornaxGeometry) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x_deg, dtype=float)
    y = np.asarray(y_deg, dtype=float)
    pa = geometry.pa_rad
    x_major = x * np.sin(pa) + y * np.cos(pa)
    y_minor = x * np.cos(pa) - y * np.sin(pa)
    return x_major, y_minor


def elliptical_radius(x_deg: np.ndarray, y_deg: np.ndarray, geometry: FornaxGeometry) -> np.ndarray:
    x_major, y_minor = elliptical_components(x_deg, y_deg, geometry)
    return np.sqrt(x_major**2 + (y_minor / geometry.q) ** 2)


def deg_to_kpc(angle_deg: np.ndarray | float, distance_kpc: float) -> np.ndarray | float:
    return np.radians(angle_deg) * distance_kpc


def kpc_to_deg(radius_kpc: float, distance_kpc: float) -> float:
    return float(np.degrees(radius_kpc / distance_kpc))


def verify_pa_convention(geometry: FornaxGeometry) -> dict[str, Any]:
    test_geometry = FornaxGeometry(center=geometry.center, ellipticity=geometry.ellipticity, q=geometry.q, pa_deg=0.0)
    north_major, north_minor = elliptical_components(np.array([0.0]), np.array([1.0]), test_geometry)
    east_major, east_minor = elliptical_components(np.array([1.0]), np.array([0.0]), test_geometry)
    return {
        "pa_zero_north_major_component": float(north_major[0]),
        "pa_zero_north_minor_component": float(north_minor[0]),
        "pa_zero_east_major_component": float(east_major[0]),
        "pa_zero_east_minor_component": float(east_minor[0]),
        "pa_zero_major_axis_points_north": bool(np.isclose(north_major[0], 1.0) and np.isclose(north_minor[0], 0.0)),
    }


def build_regular_grid(
    stars: pd.DataFrame,
    geometry: FornaxGeometry,
    pixel_size_deg: float,
    footprint_subpixels: int,
) -> GridMap:
    x = stars["x_east_deg"].to_numpy(dtype=float)
    y = stars["y_north_deg"].to_numpy(dtype=float)
    pad = pixel_size_deg
    x_edges = make_edges(x.min() - pad, x.max() + pad, pixel_size_deg)
    y_edges = make_edges(y.min() - pad, y.max() + pad, pixel_size_deg)

    counts, _, _ = np.histogram2d(x, y, bins=[x_edges, y_edges])
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    xx, yy = np.meshgrid(x_centers, y_centers, indexing="ij")

    footprint_fraction = convex_hull_coverage_fraction(x, y, x_edges, y_edges, footprint_subpixels)
    pixel_area_arcmin2 = pixel_size_deg * pixel_size_deg * 3600.0
    area_arcmin2 = footprint_fraction * pixel_area_arcmin2
    r_ell = elliptical_radius(xx, yy, geometry)

    return GridMap(
        x_edges=x_edges,
        y_edges=y_edges,
        x_centers=x_centers,
        y_centers=y_centers,
        xx=xx,
        yy=yy,
        counts=counts,
        area_arcmin2=area_arcmin2,
        footprint_fraction=footprint_fraction,
        r_ell=r_ell,
    )


def make_edges(vmin: float, vmax: float, step: float) -> np.ndarray:
    lo = np.floor(vmin / step) * step
    hi = np.ceil(vmax / step) * step
    n = int(np.ceil((hi - lo) / step))
    return lo + np.arange(n + 1) * step


def convex_hull_coverage_fraction(
    x: np.ndarray,
    y: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    subpixels: int,
) -> np.ndarray:
    if subpixels <= 0:
        raise ValueError("footprint_subpixels must be positive")

    hull = ConvexHull(np.column_stack([x, y]))
    hull_path = MplPath(np.column_stack([x, y])[hull.vertices])

    nx = len(x_edges) - 1
    ny = len(y_edges) - 1
    coverage = np.zeros((nx, ny), dtype=float)
    offsets = (np.arange(subpixels) + 0.5) / subpixels
    dx = np.diff(x_edges)
    dy = np.diff(y_edges)
    sample_x_by_cell = x_edges[:-1, None] + offsets[None, :] * dx[:, None]

    for j in range(ny):
        sample_y = y_edges[j] + offsets * dy[j]
        sx = np.repeat(sample_x_by_cell, subpixels, axis=1).reshape(-1)
        sy = np.tile(sample_y, nx * subpixels)
        inside = hull_path.contains_points(np.column_stack([sx, sy]))
        coverage[:, j] = inside.reshape(nx, subpixels, subpixels).mean(axis=(1, 2))

    return coverage


def density_from_counts(counts: np.ndarray, area_arcmin2: np.ndarray) -> np.ndarray:
    density = np.full(counts.shape, np.nan, dtype=float)
    good = area_arcmin2 > 0
    density[good] = counts[good] / area_arcmin2[good]
    return density


def sigma_clip_mask(values: np.ndarray, sigma: float = 3.0, max_iter: int = 5) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    keep = np.isfinite(values)
    for _ in range(max_iter):
        current = values[keep]
        if current.size < 3:
            break
        center = np.median(current)
        scatter = 1.4826 * np.median(np.abs(current - center))
        if not np.isfinite(scatter) or scatter <= 0:
            scatter = np.std(current)
        if not np.isfinite(scatter) or scatter <= 0:
            break
        new_keep = np.isfinite(values) & (np.abs(values - center) <= sigma * scatter)
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep
    return keep


def weighted_mean_and_scatter(values: np.ndarray, weights: np.ndarray) -> tuple[float, float, float, float]:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if good.sum() == 0:
        return np.nan, np.nan, np.nan, np.nan
    v = values[good]
    w = weights[good]
    mean = float(np.sum(w * v) / np.sum(w))
    var = float(np.sum(w * (v - mean) ** 2) / np.sum(w))
    n_eff = float(np.sum(w) ** 2 / np.sum(w**2))
    unc_mean = float(np.sqrt(var / n_eff)) if n_eff > 0 else np.nan
    return mean, float(np.sqrt(var)), unc_mean, n_eff


def estimate_background(grid: GridMap, raw_density: np.ndarray) -> dict[str, Any]:
    bg_mask = (
        (grid.r_ell >= BACKGROUND_R_INNER_DEG)
        & (grid.r_ell < BACKGROUND_R_OUTER_DEG)
        & (grid.area_arcmin2 > 0)
    )
    counts_bg = float(np.nansum(grid.counts[bg_mask]))
    area_bg = float(np.nansum(grid.area_arcmin2[bg_mask]))
    if area_bg <= 0:
        raise ValueError("No effective area in requested Yang background annulus 6 < r_ell < 10 deg.")

    bg1 = counts_bg / area_bg
    bg1_unc = np.sqrt(counts_bg) / area_bg if counts_bg > 0 else 0.0

    density_bg = raw_density[bg_mask]
    area_cells = grid.area_arcmin2[bg_mask]
    finite = np.isfinite(density_bg) & (area_cells > 0)
    bg2, bg2_rms, bg2_unc_mean, n_eff = weighted_mean_and_scatter(density_bg[finite], area_cells[finite])

    agree_sigma = np.sqrt(bg1_unc**2 + bg2_unc_mean**2) if np.isfinite(bg2_unc_mean) else bg1_unc
    methods_agree = bool(np.isfinite(bg2) and (abs(bg1 - bg2) <= max(agree_sigma, np.finfo(float).eps)))
    adopted = float(bg1)
    adopted_unc_mean = float(np.nanmax([bg1_unc, bg2_unc_mean]))

    return {
        "background_mask": bg_mask,
        "method1_density": float(bg1),
        "method1_uncertainty": float(bg1_unc),
        "method1_counts": counts_bg,
        "method1_area_arcmin2": area_bg,
        "method2_density": float(bg2),
        "method2_uncertainty_on_mean": float(bg2_unc_mean),
        "method2_spatial_rms": float(bg2_rms),
        "method2_effective_cell_count": float(n_eff),
        "method2_used_cell_count": int(finite.sum()),
        "method2_total_cell_count": int(bg_mask.sum()),
        "method2_note": (
            "area-weighted all-cell estimate; no sigma clipping is applied because "
            "the raw background pixels are sparse and zero-inflated"
        ),
        "methods_agree_within_mean_uncertainties": methods_agree,
        "adopted_density": adopted,
        "adopted_uncertainty_on_mean": adopted_unc_mean,
        "adopted_note": "count-based Method 1 mean; Method 2 used as a spatial-fluctuation diagnostic",
    }


def gaussian_kernel1d(sigma_pix: float, truncate: float = 4.0) -> np.ndarray:
    if sigma_pix <= 0:
        return np.array([1.0], dtype=float)
    radius = int(truncate * sigma_pix + 0.5)
    x = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-0.5 * (x / sigma_pix) ** 2)
    kernel /= kernel.sum()
    return kernel


def convolve_separable(image: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    tmp = convolve1d(image, kernel, axis=0, mode="constant", cval=0.0)
    return convolve1d(tmp, kernel, axis=1, mode="constant", cval=0.0)


def smooth_counts_area_variance(grid: GridMap) -> dict[str, np.ndarray | float]:
    sigma_deg = SMOOTH_FWHM_DEG / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    sigma_pix = sigma_deg / grid.pixel_size_deg
    kernel = gaussian_kernel1d(sigma_pix)

    counts_smooth = convolve_separable(grid.counts, kernel)
    area_smooth = convolve_separable(grid.area_arcmin2, kernel)
    variance_counts_smooth = convolve_separable(grid.counts, kernel**2)

    density_smooth = density_from_counts(counts_smooth, area_smooth)
    variance_density_smooth = np.full_like(density_smooth, np.nan, dtype=float)
    good = area_smooth > 0
    variance_density_smooth[good] = variance_counts_smooth[good] / area_smooth[good] ** 2

    coverage_fraction_smooth = area_smooth / grid.pixel_area_arcmin2

    return {
        "sigma_deg": float(sigma_deg),
        "sigma_pix": float(sigma_pix),
        "kernel_size": int(kernel.size),
        "kernel_sum": float(kernel.sum()),
        "kernel_squared_sum_2d": float(np.sum(kernel**2) ** 2),
        "counts_smooth": counts_smooth,
        "area_smooth": area_smooth,
        "density_smooth": density_smooth,
        "variance_density_smooth": variance_density_smooth,
        "coverage_fraction_smooth": coverage_fraction_smooth,
    }


def calibrate_surface_brightness(
    grid: GridMap,
    smooth: dict[str, np.ndarray | float],
    bg: dict[str, Any],
    coverage_threshold: float,
    distance_kpc: float,
    anchor_half_width_kpc: float,
) -> dict[str, Any]:
    density_smooth = np.asarray(smooth["density_smooth"], dtype=float)
    variance_density = np.asarray(smooth["variance_density_smooth"], dtype=float)
    area_smooth = np.asarray(smooth["area_smooth"], dtype=float)
    coverage_fraction_smooth = np.asarray(smooth["coverage_fraction_smooth"], dtype=float)

    net_density = density_smooth - float(bg["adopted_density"])
    reliable_coverage = coverage_fraction_smooth >= coverage_threshold
    r_ell_kpc = deg_to_kpc(grid.r_ell, distance_kpc)

    yang_ref_mask = (grid.r_ell < REFERENCE_RADIUS_DEG) & reliable_coverage & np.isfinite(net_density)
    if not np.any(yang_ref_mask):
        raise ValueError("No reliable cells for r_ell < 0.3 deg validation aperture.")
    yang_ref = measure_density_anchor(
        net_density,
        variance_density,
        area_smooth,
        yang_ref_mask,
        bg,
    )

    center_aperture_kpc = max(
        MUNOZ_CENTER_RADIUS_KPC,
        float(deg_to_kpc(grid.pixel_size_deg, distance_kpc)),
    )
    center_mask = (
        (r_ell_kpc <= center_aperture_kpc)
        & reliable_coverage
        & np.isfinite(net_density)
    )
    rhalf_half_width = max(float(anchor_half_width_kpc), float(deg_to_kpc(grid.pixel_size_deg, distance_kpc)))
    rhalf_mask = (
        (np.abs(r_ell_kpc - MUNOZ_RHALF_RADIUS_KPC) <= rhalf_half_width)
        & reliable_coverage
        & np.isfinite(net_density)
    )
    if not np.any(center_mask):
        raise ValueError("No reliable cells in the Muñoz central calibration aperture.")
    if not np.any(rhalf_mask):
        raise ValueError("No reliable cells near the Muñoz half-light calibration radius.")

    anchor_measurements = {
        "center": measure_density_anchor(
            net_density,
            variance_density,
            area_smooth,
            center_mask,
            bg,
            target_mu=MUNOZ_CENTER_MU_V,
            target_radius_kpc=MUNOZ_CENTER_RADIUS_KPC,
            aperture=f"r_ell_kpc <= {center_aperture_kpc:.6g} kpc",
        ),
        "rhalf": measure_density_anchor(
            net_density,
            variance_density,
            area_smooth,
            rhalf_mask,
            bg,
            target_mu=MUNOZ_RHALF_MU_V,
            target_radius_kpc=MUNOZ_RHALF_RADIUS_KPC,
            aperture=f"|r_ell_kpc - rhalf| <= {rhalf_half_width:.6g} kpc",
        ),
    }
    calibration_candidates = build_munoz_calibration_candidates(anchor_measurements)
    valid_candidates = [
        candidate for candidate in calibration_candidates
        if np.isfinite(candidate["rms_residual_mag"])
    ]
    if not valid_candidates:
        raise ValueError("No valid Muñoz calibration candidate could be evaluated.")
    selected_calibration = min(valid_candidates, key=lambda candidate: candidate["rms_residual_mag"])

    density_unc_mean_bg = np.sqrt(variance_density + float(bg["adopted_uncertainty_on_mean"]) ** 2)
    bg_fluctuation = estimate_smoothed_background_fluctuation(grid, net_density, reliable_coverage)
    density_unc_detection = np.sqrt(variance_density + bg_fluctuation["spatial_rms"] ** 2)

    mu_v = np.full_like(net_density, np.nan, dtype=float)
    positive = net_density > 0
    mu_v[positive] = density_to_mu(net_density[positive], selected_calibration["zero_point_offset"])

    mu_err_stat = np.full_like(net_density, np.nan, dtype=float)
    good_mu = positive & np.isfinite(density_unc_mean_bg)
    mu_err_stat[good_mu] = (2.5 / np.log(10.0)) * density_unc_mean_bg[good_mu] / net_density[good_mu]
    mu_err_total = np.sqrt(
        mu_err_stat**2
        + MUNOZ_MU_V_ERR**2
        + selected_calibration["zero_point_density_uncertainty_mag"] ** 2
    )

    snr_detection = np.full_like(net_density, np.nan, dtype=float)
    good_snr = np.isfinite(density_unc_detection) & (density_unc_detection > 0)
    snr_detection[good_snr] = net_density[good_snr] / density_unc_detection[good_snr]

    density_3sigma = 3.0 * bg_fluctuation["spatial_rms"]
    mu_3sigma = np.nan
    if density_3sigma > 0:
        mu_3sigma = float(density_to_mu(density_3sigma, selected_calibration["zero_point_offset"]))

    valid_science = positive & reliable_coverage & (snr_detection >= 3.0)
    faintest_reliable_mu = float(np.nanmax(mu_v[valid_science])) if np.any(valid_science) else np.nan
    anchor_predictions = evaluate_anchor_predictions(anchor_measurements, selected_calibration["zero_point_offset"])

    return {
        "net_density": net_density,
        "density_unc_mean_bg": density_unc_mean_bg,
        "density_unc_detection": density_unc_detection,
        "mu_v": mu_v,
        "mu_err_stat": mu_err_stat,
        "mu_err_total": mu_err_total,
        "snr_detection": snr_detection,
        "reliable_coverage": reliable_coverage,
        "mask_non_positive": ~positive,
        "mask_low_snr_3": ~(snr_detection >= 3.0),
        "valid_science": valid_science,
        "yang_sigma_ref_0p3": yang_ref["density"],
        "yang_sigma_ref_0p3_uncertainty": yang_ref["density_uncertainty"],
        "yang_sigma_ref_0p3_mu_if_24p77": 24.77,
        "anchor_measurements": anchor_measurements,
        "calibration_candidates": calibration_candidates,
        "selected_calibration": selected_calibration,
        "selected_anchor_predictions": anchor_predictions,
        "distance_kpc": float(distance_kpc),
        "center_anchor_aperture_kpc": float(center_aperture_kpc),
        "rhalf_anchor_half_width_kpc": float(rhalf_half_width),
        "background_fluctuation_smoothed": bg_fluctuation,
        "density_3sigma_threshold": float(density_3sigma),
        "mu_3sigma_threshold": mu_3sigma,
        "faintest_reliable_mu": faintest_reliable_mu,
    }


def measure_density_anchor(
    net_density: np.ndarray,
    variance_density: np.ndarray,
    area_smooth: np.ndarray,
    mask: np.ndarray,
    bg: dict[str, Any],
    target_mu: float | None = None,
    target_radius_kpc: float | None = None,
    aperture: str | None = None,
) -> dict[str, Any]:
    good = mask & np.isfinite(net_density) & np.isfinite(area_smooth) & (area_smooth > 0)
    if not np.any(good):
        return {
            "density": np.nan,
            "density_uncertainty": np.nan,
            "density_mu_uncertainty_mag": np.nan,
            "area_arcmin2": 0.0,
            "cell_count": 0,
            "target_mu": target_mu,
            "target_radius_kpc": target_radius_kpc,
            "aperture": aperture,
        }
    area = area_smooth[good]
    density = float(np.sum(net_density[good] * area) / np.sum(area))
    variance = variance_density[good] + float(bg["adopted_uncertainty_on_mean"]) ** 2
    density_unc = float(np.sqrt(np.nansum((area**2) * variance)) / np.sum(area))
    density_mu_unc = (
        float((2.5 / np.log(10.0)) * density_unc / density)
        if np.isfinite(density) and density > 0
        else np.nan
    )
    return {
        "density": density,
        "density_uncertainty": density_unc,
        "density_mu_uncertainty_mag": density_mu_unc,
        "area_arcmin2": float(np.sum(area)),
        "cell_count": int(np.sum(good)),
        "target_mu": target_mu,
        "target_radius_kpc": target_radius_kpc,
        "aperture": aperture,
    }


def density_to_mu(density: np.ndarray | float, zero_point_offset: float) -> np.ndarray | float:
    return -2.5 * np.log10(density) + zero_point_offset


def zero_point_from_anchor(anchor: dict[str, Any]) -> float:
    density = float(anchor["density"])
    target_mu = float(anchor["target_mu"])
    if not np.isfinite(density) or density <= 0:
        return np.nan
    return float(target_mu + 2.5 * np.log10(density))


def build_munoz_calibration_candidates(anchor_measurements: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    center_offset = zero_point_from_anchor(anchor_measurements["center"])
    rhalf_offset = zero_point_from_anchor(anchor_measurements["rhalf"])
    offsets = {
        "munoz_center_anchor": center_offset,
        "munoz_rhalf_anchor": rhalf_offset,
    }
    if np.isfinite(center_offset) and np.isfinite(rhalf_offset):
        offsets["munoz_two_point_best"] = float(0.5 * (center_offset + rhalf_offset))

    candidates = []
    for name, offset in offsets.items():
        predictions = evaluate_anchor_predictions(anchor_measurements, offset)
        residuals = [
            value["residual_mag"]
            for value in predictions.values()
            if np.isfinite(value["residual_mag"])
        ]
        density_mu_uncertainties = [
            float(anchor["density_mu_uncertainty_mag"])
            for anchor in anchor_measurements.values()
            if np.isfinite(anchor["density_mu_uncertainty_mag"])
        ]
        candidates.append(
            {
                "name": name,
                "zero_point_offset": float(offset),
                "rms_residual_mag": (
                    float(np.sqrt(np.mean(np.square(residuals))))
                    if residuals
                    else np.nan
                ),
                "max_abs_residual_mag": (
                    float(np.max(np.abs(residuals)))
                    if residuals
                    else np.nan
                ),
                "zero_point_density_uncertainty_mag": (
                    float(np.sqrt(np.mean(np.square(density_mu_uncertainties))))
                    if density_mu_uncertainties
                    else np.nan
                ),
                "anchor_predictions": predictions,
            }
        )
    return candidates


def evaluate_anchor_predictions(
    anchor_measurements: dict[str, dict[str, Any]],
    zero_point_offset: float,
) -> dict[str, dict[str, float]]:
    predictions = {}
    for name, anchor in anchor_measurements.items():
        density = float(anchor["density"])
        target_mu = anchor["target_mu"]
        mu_pred = (
            float(density_to_mu(density, zero_point_offset))
            if np.isfinite(density) and density > 0 and np.isfinite(zero_point_offset)
            else np.nan
        )
        residual = (
            float(mu_pred - float(target_mu))
            if target_mu is not None and np.isfinite(mu_pred)
            else np.nan
        )
        predictions[name] = {
            "predicted_mu": mu_pred,
            "target_mu": float(target_mu) if target_mu is not None else np.nan,
            "residual_mag": residual,
        }
    return predictions


def estimate_smoothed_background_fluctuation(
    grid: GridMap,
    net_density: np.ndarray,
    reliable_coverage: np.ndarray,
) -> dict[str, Any]:
    bg_mask = (
        (grid.r_ell >= BACKGROUND_R_INNER_DEG)
        & (grid.r_ell < BACKGROUND_R_OUTER_DEG)
        & reliable_coverage
        & np.isfinite(net_density)
    )
    if not np.any(bg_mask):
        return {"mean": np.nan, "spatial_rms": np.nan, "uncertainty_on_mean": np.nan, "cell_count": 0}

    clipped_values = net_density[bg_mask]
    mean = float(np.mean(clipped_values))
    rms = float(np.std(clipped_values, ddof=1)) if clipped_values.size > 1 else np.nan
    unc_mean = float(rms / np.sqrt(clipped_values.size)) if clipped_values.size > 1 else np.nan
    return {
        "mean": mean,
        "spatial_rms": rms,
        "uncertainty_on_mean": unc_mean,
        "cell_count": int(clipped_values.size),
        "total_background_cell_count": int(bg_mask.sum()),
        "note": "RMS of all reliable smoothed background cells; no clipping is applied to sparse counts.",
    }


def estimate_yang_style_1d_detection_limit(
    bg: dict[str, Any],
    zero_point_offset: float,
) -> dict[str, Any]:
    method_uncertainties = {
        "method1_poisson_counts_in_6_10deg": float(bg.get("method1_uncertainty", np.nan)),
        "method2_weighted_mean_in_6_10deg": float(bg.get("method2_uncertainty_on_mean", np.nan)),
    }
    valid_methods = {
        name: value
        for name, value in method_uncertainties.items()
        if np.isfinite(value) and value >= 0
    }
    if not valid_methods:
        return {
            "background_density_stars_arcmin2": float(bg.get("adopted_density", np.nan)),
            "density_1sigma_background_stars_arcmin2": np.nan,
            "density_3sigma_threshold_stars_arcmin2": np.nan,
            "mu_3sigma_threshold_mag_arcsec2": np.nan,
            "sigma_level": YANG_STYLE_DETECTION_SIGMA,
            "adopted_uncertainty_method": None,
            "method_uncertainties_stars_arcmin2": method_uncertainties,
            "note": "No finite 1D background uncertainty was available.",
        }

    adopted_method = max(valid_methods, key=valid_methods.get)
    density_1sigma = valid_methods[adopted_method]
    density_threshold = YANG_STYLE_DETECTION_SIGMA * density_1sigma
    mu_threshold = (
        float(density_to_mu(density_threshold, zero_point_offset))
        if density_threshold > 0 and np.isfinite(zero_point_offset)
        else np.nan
    )
    return {
        "background_density_stars_arcmin2": float(bg.get("adopted_density", np.nan)),
        "density_1sigma_background_stars_arcmin2": float(density_1sigma),
        "density_3sigma_threshold_stars_arcmin2": float(density_threshold),
        "mu_3sigma_threshold_mag_arcsec2": mu_threshold,
        "sigma_level": YANG_STYLE_DETECTION_SIGMA,
        "adopted_uncertainty_method": adopted_method,
        "method_uncertainties_stars_arcmin2": method_uncertainties,
        "yang_s1_reference_density_1sigma_stars_arcmin2": YANG_S1_1D_BACKGROUND_UNCERTAINTY,
        "yang_s1_reference_mu_3sigma_mag_arcsec2": YANG_S1_1D_LIMIT_MAG_ARCSEC2,
        "note": (
            "Yang-style 1D limit: use the larger 6-10 deg background mean uncertainty "
            "from Method 1 and Method 2, then convert 3 sigma above the mean background "
            "to equivalent V-band surface brightness."
        ),
    }


def area_weighted_mean(values: np.ndarray, area: np.ndarray) -> float:
    good = np.isfinite(values) & np.isfinite(area) & (area > 0)
    if not np.any(good):
        return np.nan
    return float(np.sum(values[good] * area[good]) / np.sum(area[good]))


def make_map_table(
    grid: GridMap,
    raw_density: np.ndarray,
    smooth: dict[str, np.ndarray | float],
    sb: dict[str, Any],
    distance_kpc: float,
) -> pd.DataFrame:
    x_major, y_minor = elliptical_components(grid.xx, grid.yy, FornaxGeometry(center=SkyCoord(FORNAX_CENTER_RA, FORNAX_CENTER_DEC)))
    table = pd.DataFrame(
        {
            "x_east_deg": grid.xx.ravel(),
            "y_north_deg": grid.yy.ravel(),
            "x_major_deg": x_major.ravel(),
            "y_minor_deg": y_minor.ravel(),
            "r_ell_deg": grid.r_ell.ravel(),
            "r_ell_kpc": np.asarray(deg_to_kpc(grid.r_ell, distance_kpc)).ravel(),
            "counts_raw": grid.counts.ravel(),
            "area_arcmin2": grid.area_arcmin2.ravel(),
            "footprint_fraction": grid.footprint_fraction.ravel(),
            "density_raw_stars_arcmin2": raw_density.ravel(),
            "counts_smooth": np.asarray(smooth["counts_smooth"]).ravel(),
            "area_smooth_arcmin2": np.asarray(smooth["area_smooth"]).ravel(),
            "coverage_fraction_smooth": np.asarray(smooth["coverage_fraction_smooth"]).ravel(),
            "density_smooth_stars_arcmin2": np.asarray(smooth["density_smooth"]).ravel(),
            "density_net_stars_arcmin2": np.asarray(sb["net_density"]).ravel(),
            "density_unc_stars_arcmin2": np.asarray(sb["density_unc_mean_bg"]).ravel(),
            "snr_detection": np.asarray(sb["snr_detection"]).ravel(),
            "mu_v_mag_arcsec2": np.asarray(sb["mu_v"]).ravel(),
            "mu_v_err_stat": np.asarray(sb["mu_err_stat"]).ravel(),
            "mu_v_err_total": np.asarray(sb["mu_err_total"]).ravel(),
            "mask_non_positive": np.asarray(sb["mask_non_positive"]).ravel(),
            "mask_low_snr_3": np.asarray(sb["mask_low_snr_3"]).ravel(),
            "mask_low_coverage": (~np.asarray(sb["reliable_coverage"])).ravel(),
            "valid_science": np.asarray(sb["valid_science"]).ravel(),
        }
    )
    net_density = np.asarray(sb["net_density"], dtype=float)
    positive = net_density > 0
    for candidate in sb["calibration_candidates"]:
        mu_candidate = np.full_like(net_density, np.nan, dtype=float)
        mu_candidate[positive] = density_to_mu(
            net_density[positive],
            candidate["zero_point_offset"],
        )
        table[f"mu_v_{candidate['name']}"] = mu_candidate.ravel()
    return table


def make_radial_profile(
    grid: GridMap,
    bg: dict[str, Any],
    zero_point_offset: float,
    distance_kpc: float,
    profile_max_kpc: float,
    calibration_mu_uncertainty: float,
    min_counts: int = 9,
) -> pd.DataFrame:
    r_kpc_all = np.asarray(deg_to_kpc(grid.r_ell, distance_kpc), dtype=float)
    finite = (
        np.isfinite(grid.r_ell)
        & np.isfinite(r_kpc_all)
        & (grid.area_arcmin2 > 0)
        & (r_kpc_all <= profile_max_kpc)
    )
    order = np.argsort(r_kpc_all[finite])
    r_deg = grid.r_ell[finite][order]
    r_kpc = r_kpc_all[finite][order]
    counts = grid.counts[finite][order]
    area = grid.area_arcmin2[finite][order]

    rows: list[dict[str, float]] = []
    start = 0
    while start < len(r_kpc):
        count_sum = 0.0
        area_sum = 0.0
        end = start
        while end < len(r_kpc) and (count_sum < min_counts or end == start):
            count_sum += counts[end]
            area_sum += area[end]
            end += 1
        if area_sum > 0:
            r1_deg = float(r_deg[start])
            r2_deg = float(r_deg[end - 1])
            r1_kpc = float(r_kpc[start])
            r2_kpc = float(r_kpc[end - 1])
            density_raw = count_sum / area_sum
            density_net = density_raw - float(bg["adopted_density"])
            density_unc = np.sqrt(count_sum) / area_sum if count_sum > 0 else 0.0
            density_net_unc = np.sqrt(density_unc**2 + float(bg["adopted_uncertainty_on_mean"]) ** 2)
            mu = np.nan
            mu_err = np.nan
            if density_net > 0:
                mu = density_to_mu(density_net, zero_point_offset)
                mu_err = np.sqrt(
                    ((2.5 / np.log(10.0)) * density_net_unc / density_net) ** 2
                    + MUNOZ_MU_V_ERR**2
                    + calibration_mu_uncertainty**2
                )
            rows.append(
                {
                    "r_inner_deg": r1_deg,
                    "r_outer_deg": r2_deg,
                    "r_mid_deg": 0.5 * (r1_deg + r2_deg),
                    "r_inner_kpc": r1_kpc,
                    "r_outer_kpc": r2_kpc,
                    "r_mid_kpc": 0.5 * (r1_kpc + r2_kpc),
                    "counts": count_sum,
                    "area_arcmin2": area_sum,
                    "density_raw_stars_arcmin2": density_raw,
                    "density_net_stars_arcmin2": density_net,
                    "density_net_unc_stars_arcmin2": density_net_unc,
                    "mu_v_mag_arcsec2": mu,
                    "mu_v_err_total": mu_err,
                }
            )
        start = end

    return pd.DataFrame(rows)


def annotate_profile_with_yang_style_1d_limit(
    profile: pd.DataFrame,
    yang_1d_detection: dict[str, Any],
) -> pd.DataFrame:
    out = profile.copy()
    sigma_density = float(yang_1d_detection.get("density_1sigma_background_stars_arcmin2", np.nan))
    threshold_density = float(yang_1d_detection.get("density_3sigma_threshold_stars_arcmin2", np.nan))

    if np.isfinite(sigma_density) and sigma_density > 0:
        out["snr_yang_style_1d_background"] = out["density_net_stars_arcmin2"] / sigma_density
    else:
        out["snr_yang_style_1d_background"] = np.nan

    out["above_yang_style_1d_3sigma"] = (
        np.isfinite(out["density_net_stars_arcmin2"])
        & np.isfinite(threshold_density)
        & (out["density_net_stars_arcmin2"] >= threshold_density)
    )
    return out


def summarize_yang_style_1d_profile_detection(
    profile: pd.DataFrame,
    plot_max_kpc: float,
) -> dict[str, Any]:
    base = np.isfinite(profile["mu_v_mag_arcsec2"]) & profile["above_yang_style_1d_3sigma"]
    within_plot = base & (profile["r_mid_kpc"] <= plot_max_kpc)
    if not np.any(base):
        return {
            "detected_bin_count": 0,
            "detected_bin_count_within_plot_radius": 0,
            "faintest_detected_mu_mag_arcsec2": np.nan,
            "faintest_detected_mu_within_plot_radius_mag_arcsec2": np.nan,
            "max_detected_r_mid_kpc": np.nan,
            "max_detected_r_mid_within_plot_radius_kpc": np.nan,
        }

    summary = {
        "detected_bin_count": int(base.sum()),
        "detected_bin_count_within_plot_radius": int(within_plot.sum()),
        "faintest_detected_mu_mag_arcsec2": float(np.nanmax(profile.loc[base, "mu_v_mag_arcsec2"])),
        "max_detected_r_mid_kpc": float(np.nanmax(profile.loc[base, "r_mid_kpc"])),
    }
    if np.any(within_plot):
        summary.update(
            {
                "faintest_detected_mu_within_plot_radius_mag_arcsec2": float(
                    np.nanmax(profile.loc[within_plot, "mu_v_mag_arcsec2"])
                ),
                "max_detected_r_mid_within_plot_radius_kpc": float(
                    np.nanmax(profile.loc[within_plot, "r_mid_kpc"])
                ),
            }
        )
    else:
        summary.update(
            {
                "faintest_detected_mu_within_plot_radius_mag_arcsec2": np.nan,
                "max_detected_r_mid_within_plot_radius_kpc": np.nan,
            }
        )
    return summary


def compare_existing_radial_profile(path: Path | None) -> dict[str, Any]:
    if path is None or not path.exists():
        return {"available": False, "path": None}
    dprof = pd.read_csv(path)
    required = {"xr", "area", "counts", "density"}
    if not required.issubset(dprof.columns):
        return {"available": False, "path": str(path), "reason": "missing required columns"}

    out: dict[str, Any] = {"available": True, "path": str(path), "area_unit_inferred": "deg^2"}
    for label, lo, hi in [
        ("central_r_ell_lt_0p3", 0.0, REFERENCE_RADIUS_DEG),
        ("background_6_10", BACKGROUND_R_INNER_DEG, BACKGROUND_R_OUTER_DEG),
    ]:
        mask = (dprof["xr"] >= lo) & (dprof["xr"] < hi)
        counts = float(dprof.loc[mask, "counts"].sum())
        area_deg2 = float(dprof.loc[mask, "area"].sum())
        density = counts / (area_deg2 * 3600.0) if area_deg2 > 0 else np.nan
        out[label] = {
            "bin_count": int(mask.sum()),
            "counts": counts,
            "area_deg2": area_deg2,
            "density_stars_arcmin2": float(density),
        }
    return out


def make_preprocessed_observation_profile(
    path: Path | None,
    zero_point_offset: float,
    calibration_name: str,
    calibration_mu_uncertainty: float,
    distance_kpc: float,
) -> tuple[pd.DataFrame | None, dict[str, Any]]:
    if path is None or not path.exists():
        return None, {"available": False, "path": None if path is None else str(path)}

    dprof = pd.read_csv(path)
    required = {"xr", "xr1", "xr2", "density", "density_error", "area", "counts"}
    missing = sorted(required - set(dprof.columns))
    if missing:
        return None, {"available": False, "path": str(path), "reason": f"missing columns: {missing}"}

    out = dprof.copy()
    area_arcmin2 = out["area"].to_numpy(dtype=float) * 3600.0
    counts = out["counts"].to_numpy(dtype=float)
    density_raw = np.divide(
        counts,
        area_arcmin2,
        out=np.full_like(counts, np.nan, dtype=float),
        where=area_arcmin2 > 0,
    )
    density_error = out["density_error"].to_numpy(dtype=float)

    bg_mask = (
        (out["xr"].to_numpy(dtype=float) >= BACKGROUND_R_INNER_DEG)
        & (out["xr"].to_numpy(dtype=float) < BACKGROUND_R_OUTER_DEG)
        & np.isfinite(area_arcmin2)
        & (area_arcmin2 > 0)
    )
    bg_counts = float(np.nansum(counts[bg_mask]))
    bg_area = float(np.nansum(area_arcmin2[bg_mask]))
    if bg_area <= 0:
        raise ValueError("No radial-profile area in 6 < r_ell < 10 deg for prepared observation CSV.")
    background_density = bg_counts / bg_area
    background_uncertainty = np.sqrt(bg_counts) / bg_area if bg_counts > 0 else 0.0

    density_net = density_raw - background_density
    density_net_unc = np.sqrt(density_error**2 + background_uncertainty**2)
    mu_v = np.full_like(density_net, np.nan, dtype=float)
    mu_err = np.full_like(density_net, np.nan, dtype=float)
    positive = density_net > 0
    mu_v[positive] = density_to_mu(density_net[positive], zero_point_offset)
    mu_err[positive] = np.sqrt(
        ((2.5 / np.log(10.0)) * density_net_unc[positive] / density_net[positive]) ** 2
        + MUNOZ_MU_V_ERR**2
        + calibration_mu_uncertainty**2
    )

    density_3sigma = YANG_STYLE_DETECTION_SIGMA * background_uncertainty
    mu_3sigma = (
        float(density_to_mu(density_3sigma, zero_point_offset))
        if density_3sigma > 0
        else np.nan
    )

    out["r_mid_deg"] = out["xr"].to_numpy(dtype=float)
    out["r_inner_deg"] = out["xr1"].to_numpy(dtype=float)
    out["r_outer_deg"] = out["xr2"].to_numpy(dtype=float)
    out["r_mid_kpc"] = deg_to_kpc(out["r_mid_deg"].to_numpy(dtype=float), distance_kpc)
    out["r_inner_kpc"] = deg_to_kpc(out["r_inner_deg"].to_numpy(dtype=float), distance_kpc)
    out["r_outer_kpc"] = deg_to_kpc(out["r_outer_deg"].to_numpy(dtype=float), distance_kpc)
    out["area_arcmin2"] = area_arcmin2
    out["density_raw_stars_arcmin2"] = density_raw
    out["background_density_stars_arcmin2"] = background_density
    out["background_uncertainty_1sigma_stars_arcmin2"] = background_uncertainty
    out["background_bin_count"] = int(np.sum(bg_mask))
    out["density_net_stars_arcmin2"] = density_net
    out["density_net_unc_stars_arcmin2"] = density_net_unc
    out["mu_v_mag_arcsec2"] = mu_v
    out["mu_v_err_total"] = mu_err
    out["snr_yang_style_1d_background"] = density_net / background_uncertainty if background_uncertainty > 0 else np.nan
    out["above_yang_style_1d_3sigma"] = density_net >= density_3sigma
    out["mu_yang_style_1d_3sigma_limit_mag_arcsec2"] = mu_3sigma
    out["zero_point_offset"] = zero_point_offset
    out["calibration_name"] = calibration_name
    out["source_profile_path"] = str(path)
    out["code_version"] = CODE_VERSION

    preferred_order = [
        "r_inner_kpc",
        "r_outer_kpc",
        "r_mid_kpc",
        "r_inner_deg",
        "r_outer_deg",
        "r_mid_deg",
        "xr",
        "xr1",
        "xr2",
        "counts",
        "area",
        "area_arcmin2",
        "density_raw_stars_arcmin2",
        "density_error",
        "background_density_stars_arcmin2",
        "background_uncertainty_1sigma_stars_arcmin2",
        "density_net_stars_arcmin2",
        "density_net_unc_stars_arcmin2",
        "mu_v_mag_arcsec2",
        "mu_v_err_total",
        "snr_yang_style_1d_background",
        "above_yang_style_1d_3sigma",
        "mu_yang_style_1d_3sigma_limit_mag_arcsec2",
        "zero_point_offset",
        "calibration_name",
        "source_profile_path",
        "code_version",
    ]
    remaining = [column for column in out.columns if column not in preferred_order]
    out = out[[column for column in preferred_order if column in out.columns] + remaining]

    diagnostics = {
        "available": True,
        "path": str(path),
        "row_count": int(len(out)),
        "background_density_stars_arcmin2": float(background_density),
        "background_uncertainty_1sigma_stars_arcmin2": float(background_uncertainty),
        "background_bin_count": int(np.sum(bg_mask)),
        "mu_yang_style_1d_3sigma_limit_mag_arcsec2": mu_3sigma,
        "positive_mu_rows": int(np.sum(np.isfinite(mu_v))),
    }
    return out, diagnostics


def contour_levels(mu_v: np.ndarray, valid: np.ndarray, mu_limit: float) -> np.ndarray:
    values = mu_v[valid & np.isfinite(mu_v)]
    if values.size == 0 or not np.isfinite(mu_limit):
        return np.array([], dtype=float)
    bright = np.nanmin(values)
    faint = min(np.nanmax(values), mu_limit)
    start = np.ceil(bright * 2.0) / 2.0
    stop = np.floor(faint * 2.0) / 2.0
    if stop < start:
        return np.array([], dtype=float)
    return np.arange(start, stop + 0.25, 0.5)


def ellipse_points(radius_deg: float, geometry: FornaxGeometry, n: int = 300) -> tuple[np.ndarray, np.ndarray]:
    theta = np.linspace(0.0, 2.0 * np.pi, n)
    x_major = radius_deg * np.cos(theta)
    y_minor = geometry.q * radius_deg * np.sin(theta)
    pa = geometry.pa_rad
    x = x_major * np.sin(pa) + y_minor * np.cos(pa)
    y = x_major * np.cos(pa) - y_minor * np.sin(pa)
    return x, y


def plot_surface_brightness_map(
    grid: GridMap,
    sb: dict[str, Any],
    geometry: FornaxGeometry,
    output_png: Path,
    output_pdf: Path,
) -> list[float]:
    mu_v = np.asarray(sb["mu_v"], dtype=float)
    valid = np.asarray(sb["valid_science"], dtype=bool)
    mu_plot = np.where(valid, mu_v, np.nan)
    central = (
        (grid.xx >= -CENTRAL_PLOT_HALF_WIDTH_DEG)
        & (grid.xx <= CENTRAL_PLOT_HALF_WIDTH_DEG)
        & (grid.yy >= -CENTRAL_PLOT_HALF_WIDTH_DEG)
        & (grid.yy <= CENTRAL_PLOT_HALF_WIDTH_DEG)
    )
    levels = contour_levels(mu_v, valid & central, float(sb["mu_3sigma_threshold"]))

    cmap = plt.get_cmap("magma_r").copy()
    cmap.set_bad("0.82")
    fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    image = ax.imshow(
        mu_plot.T,
        origin="lower",
        extent=[grid.x_edges[0], grid.x_edges[-1], grid.y_edges[0], grid.y_edges[-1]],
        cmap=cmap,
        vmin=np.nanpercentile(mu_plot[central & np.isfinite(mu_plot)], 2) if np.any(central & np.isfinite(mu_plot)) else None,
        vmax=float(sb["mu_3sigma_threshold"]) if np.isfinite(sb["mu_3sigma_threshold"]) else None,
        interpolation="nearest",
        aspect="equal",
    )
    cbar = fig.colorbar(image, ax=ax, pad=0.02)
    cbar.set_label(r"$\mu_V$ (mag arcsec$^{-2}$)")
    cbar.ax.invert_yaxis()

    if levels.size > 0:
        contours = ax.contour(grid.xx, grid.yy, mu_plot, levels=levels, colors="black", linewidths=0.8)
        ax.clabel(contours, fmt="%.1f", fontsize=8)

    for radius, color in [(0.3, "white"), (1.3, "cyan"), (2.1, "lime")]:
        ex, ey = ellipse_points(radius, geometry)
        ax.plot(ex, ey, color=color, lw=1.0, ls="--", alpha=0.9)

    ax.scatter([0.0], [0.0], marker="+", s=80, c="red", lw=1.6, label="Fornax centre")
    ax.annotate("N", xy=(2.55, -2.55), xytext=(2.55, -2.9), arrowprops={"arrowstyle": "->", "lw": 1.2}, ha="center")
    ax.annotate("E", xy=(2.9, -2.55), xytext=(2.55, -2.55), arrowprops={"arrowstyle": "->", "lw": 1.2}, va="center")
    ax.set_xlim(-CENTRAL_PLOT_HALF_WIDTH_DEG, CENTRAL_PLOT_HALF_WIDTH_DEG)
    ax.set_ylim(-CENTRAL_PLOT_HALF_WIDTH_DEG, CENTRAL_PLOT_HALF_WIDTH_DEG)
    ax.set_xlabel("East offset (deg)")
    ax.set_ylabel("North offset (deg)")
    ax.set_title("Fornax equivalent V-band surface brightness")
    ax.legend(loc="upper right", fontsize=8)
    ax.text(
        0.02,
        0.02,
        "gray: masked or <3 sigma",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 2},
    )
    fig.savefig(output_png, dpi=220)
    fig.savefig(output_pdf)
    plt.close(fig)
    return [float(level) for level in levels]


def plot_radial_profile(
    profile: pd.DataFrame,
    bg: dict[str, Any],
    sb: dict[str, Any],
    output_png: Path,
    plot_max_kpc: float,
) -> None:
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(7.0, 7.5), sharex=True, constrained_layout=True)
    display_profile = profile[profile["r_mid_kpc"] <= plot_max_kpc].copy()

    ax0.errorbar(
        display_profile["r_mid_kpc"],
        display_profile["density_raw_stars_arcmin2"],
        yerr=np.sqrt(display_profile["counts"]) / display_profile["area_arcmin2"],
        fmt="o",
        ms=3,
        lw=0.8,
        label="raw density",
    )
    ax0.errorbar(
        display_profile["r_mid_kpc"],
        display_profile["density_net_stars_arcmin2"],
        yerr=display_profile["density_net_unc_stars_arcmin2"],
        fmt="s",
        ms=3,
        lw=0.8,
        label="background-subtracted",
    )
    bg_density = float(bg["adopted_density"])
    yang_1d = sb.get("yang_style_1d_detection", {})
    one_d_sigma = float(yang_1d.get("density_1sigma_background_stars_arcmin2", np.nan))
    ax0.axhline(bg_density, color="gray", ls="--", lw=1.2, label="adopted background")
    if np.isfinite(one_d_sigma):
        ax0.axhline(bg_density + one_d_sigma, color="gray", ls=":", lw=1.0, label="Yang-style 1D 1 sigma")
        ax0.axhline(
            bg_density + YANG_STYLE_DETECTION_SIGMA * one_d_sigma,
            color="gray",
            ls="-.",
            lw=1.0,
            label="Yang-style 1D 3 sigma",
        )
    ax0.set_yscale("symlog", linthresh=1e-4)
    ax0.set_ylabel(r"$\Sigma$ (stars arcmin$^{-2}$)")
    ax0.legend(fontsize=8)

    mu_limit = float(yang_1d.get("mu_3sigma_threshold_mag_arcsec2", sb["mu_3sigma_threshold"]))
    if "above_yang_style_1d_3sigma" in display_profile.columns:
        good = np.isfinite(display_profile["mu_v_mag_arcsec2"]) & display_profile["above_yang_style_1d_3sigma"]
    else:
        density_limit = float(yang_1d.get("density_3sigma_threshold_stars_arcmin2", sb["density_3sigma_threshold"]))
        good = (
            np.isfinite(display_profile["mu_v_mag_arcsec2"])
            & (display_profile["density_net_stars_arcmin2"] >= density_limit)
        )
    ax1.errorbar(
        display_profile.loc[good, "r_mid_kpc"],
        display_profile.loc[good, "mu_v_mag_arcsec2"],
        yerr=display_profile.loc[good, "mu_v_err_total"],
        fmt="o",
        ms=3,
        lw=0.8,
        color="black",
    )
    if np.isfinite(mu_limit):
        ax1.axhline(mu_limit, color="gray", ls="-.", lw=1.0, label="Yang-style 1D 3 sigma limit")
    ax1.scatter(
        [MUNOZ_CENTER_RADIUS_KPC, MUNOZ_RHALF_RADIUS_KPC],
        [MUNOZ_CENTER_MU_V, MUNOZ_RHALF_MU_V],
        facecolors="none",
        edgecolors="orange",
        s=55,
        lw=1.4,
        label="Muñoz anchors",
        zorder=5,
    )
    ax1.invert_yaxis()
    if np.any(good):
        bright_limit = np.floor(float(np.nanmin(display_profile.loc[good, "mu_v_mag_arcsec2"])) - 0.5)
        faint_limit = np.ceil(mu_limit + 0.5) if np.isfinite(mu_limit) else np.ceil(
            float(np.nanmax(display_profile.loc[good, "mu_v_mag_arcsec2"])) + 0.5
        )
        ax1.set_ylim(faint_limit, bright_limit)
    ax1.set_xlim(0.0, plot_max_kpc)
    ax1.set_xlabel(r"Elliptical semimajor radius $r_{\rm ell}$ (kpc)")
    ax1.set_ylabel(r"$\mu_V$ (mag arcsec$^{-2}$)")
    ax1.legend(fontsize=8)
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


def calibration_comparison_table(candidates: list[dict[str, Any]]) -> pd.DataFrame:
    rows = []
    for candidate in candidates:
        row = {
            "name": candidate["name"],
            "zero_point_offset": candidate["zero_point_offset"],
            "rms_residual_mag": candidate["rms_residual_mag"],
            "max_abs_residual_mag": candidate["max_abs_residual_mag"],
            "zero_point_density_uncertainty_mag": candidate["zero_point_density_uncertainty_mag"],
        }
        for anchor_name, prediction in candidate["anchor_predictions"].items():
            row[f"{anchor_name}_target_mu"] = prediction["target_mu"]
            row[f"{anchor_name}_predicted_mu"] = prediction["predicted_mu"]
            row[f"{anchor_name}_residual_mag"] = prediction["residual_mag"]
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_profile_annulus(
    profile: pd.DataFrame,
    r_inner_kpc: float,
    r_outer_kpc: float,
) -> dict[str, Any]:
    overlap = (profile["r_outer_kpc"] > r_inner_kpc) & (profile["r_inner_kpc"] < r_outer_kpc)
    if not np.any(overlap):
        return {
            "r_inner_kpc": float(r_inner_kpc),
            "r_outer_kpc": float(r_outer_kpc),
            "bin_count": 0,
            "counts": 0.0,
            "area_arcmin2": 0.0,
            "density_raw_stars_arcmin2": np.nan,
            "density_net_stars_arcmin2": np.nan,
        }

    subset = profile.loc[overlap]
    counts = float(subset["counts"].sum())
    area = float(subset["area_arcmin2"].sum())
    raw_density = counts / area if area > 0 else np.nan
    bg_weighted_net_density = (
        float(np.sum(subset["density_net_stars_arcmin2"] * subset["area_arcmin2"]) / area)
        if area > 0
        else np.nan
    )
    return {
        "r_inner_kpc": float(r_inner_kpc),
        "r_outer_kpc": float(r_outer_kpc),
        "bin_count": int(overlap.sum()),
        "counts": counts,
        "area_arcmin2": area,
        "density_raw_stars_arcmin2": float(raw_density),
        "density_net_stars_arcmin2": bg_weighted_net_density,
    }


def write_diagnostics(path: Path, diagnostics: dict[str, Any]) -> None:
    path.write_text(json.dumps(to_jsonable(diagnostics), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def to_jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): to_jsonable(v) for k, v in value.items() if not isinstance(v, np.ndarray)}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    return value


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    geometry = FornaxGeometry(center=SkyCoord(FORNAX_CENTER_RA, FORNAX_CENTER_DEC))
    stars_raw, input_diag = load_and_validate_candidates(args.input)
    stars = project_to_tangent_plane(stars_raw, geometry)
    grid = build_regular_grid(stars, geometry, args.pixel_size_deg, args.footprint_subpixels)
    raw_density = density_from_counts(grid.counts, grid.area_arcmin2)
    bg = estimate_background(grid, raw_density)
    smooth = smooth_counts_area_variance(grid)
    sb = calibrate_surface_brightness(
        grid,
        smooth,
        bg,
        args.coverage_threshold,
        args.distance_kpc,
        args.munoz_anchor_half_width_kpc,
    )
    zero_point_offset = float(sb["selected_calibration"]["zero_point_offset"])
    yang_1d_detection = estimate_yang_style_1d_detection_limit(bg, zero_point_offset)
    sb["yang_style_1d_detection"] = yang_1d_detection
    profile = make_radial_profile(
        grid,
        bg,
        zero_point_offset,
        args.distance_kpc,
        args.profile_selection_max_kpc,
        float(sb["selected_calibration"]["zero_point_density_uncertainty_mag"]),
    )
    profile = annotate_profile_with_yang_style_1d_limit(profile, yang_1d_detection)
    prepared_observation_profile, prepared_observation_diag = make_preprocessed_observation_profile(
        args.radial_input,
        zero_point_offset,
        str(sb["selected_calibration"]["name"]),
        float(sb["selected_calibration"]["zero_point_density_uncertainty_mag"]),
        args.distance_kpc,
    )

    prefix = f"yang_fornax_sb_{CODE_VERSION}"
    map_csv = args.output_dir / f"{prefix}_processed_map.csv.gz"
    profile_csv = args.output_dir / f"{prefix}_radial_profile.csv"
    calibration_csv = args.output_dir / f"{prefix}_calibration_comparison.csv"
    map_png = args.output_dir / f"{prefix}_contour_map.png"
    map_pdf = args.output_dir / f"{prefix}_contour_map.pdf"
    profile_png = args.output_dir / f"{prefix}_radial_profile.png"
    diagnostics_json = args.output_dir / f"{prefix}_diagnostics.json"
    prepared_observation_csv = (
        args.prepared_observation_output
        if args.prepared_observation_output is not None
        else default_prepared_observation_output()
    )

    map_table = make_map_table(grid, raw_density, smooth, sb, args.distance_kpc)
    map_table.to_csv(map_csv, index=False)
    profile.to_csv(profile_csv, index=False)
    if prepared_observation_profile is not None:
        prepared_observation_csv.parent.mkdir(parents=True, exist_ok=True)
        prepared_observation_profile.to_csv(prepared_observation_csv, index=False)
        prepared_observation_diag["output_csv"] = str(prepared_observation_csv)
    calibration_comparison_table(sb["calibration_candidates"]).to_csv(calibration_csv, index=False)
    contour_level_values = plot_surface_brightness_map(grid, sb, geometry, map_png, map_pdf)
    plot_radial_profile(profile, bg, sb, profile_png, args.profile_plot_max_kpc)

    existing_profile_validation = compare_existing_radial_profile(args.radial_input)
    central_profile_density = np.nan
    background_profile_density = np.nan
    if existing_profile_validation.get("available"):
        central_profile_density = existing_profile_validation["central_r_ell_lt_0p3"]["density_stars_arcmin2"]
        background_profile_density = existing_profile_validation["background_6_10"]["density_stars_arcmin2"]
    reliable = np.asarray(sb["valid_science"], dtype=bool)
    non_positive = np.asarray(sb["mask_non_positive"], dtype=bool)
    low_snr = np.asarray(sb["mask_low_snr_3"], dtype=bool)
    low_coverage = ~np.asarray(sb["reliable_coverage"], dtype=bool)
    total_cells = int(grid.counts.size)
    outer_profile_contamination = summarize_profile_annulus(
        profile,
        args.profile_plot_max_kpc,
        args.profile_selection_max_kpc,
    )
    yang_1d_profile_summary = summarize_yang_style_1d_profile_detection(
        profile,
        args.profile_plot_max_kpc,
    )

    diagnostics = {
        "code_version": CODE_VERSION,
        "input": input_diag,
        "outputs": {
            "script": str(Path(__file__).resolve()),
            "processed_map_csv_gz": str(map_csv),
            "radial_profile_csv": str(profile_csv),
            "calibration_comparison_csv": str(calibration_csv),
            "contour_map_png": str(map_png),
            "contour_map_pdf": str(map_pdf),
            "radial_profile_png": str(profile_png),
            "diagnostics_json": str(diagnostics_json),
            "prepared_observation_profile_csv": (
                str(prepared_observation_csv)
                if prepared_observation_profile is not None
                else None
            ),
        },
        "geometry": {
            "center_ra": FORNAX_CENTER_RA,
            "center_dec": FORNAX_CENTER_DEC,
            "center_ra_deg": float(geometry.center.ra.deg),
            "center_dec_deg": float(geometry.center.dec.deg),
            "ellipticity": geometry.ellipticity,
            "axis_ratio_q": geometry.q,
            "pa_deg_east_of_north": geometry.pa_deg,
            "pa_convention_test": verify_pa_convention(geometry),
        },
        "grid": {
            "nx": int(len(grid.x_centers)),
            "ny": int(len(grid.y_centers)),
            "pixel_size_deg": grid.pixel_size_deg,
            "pixel_area_arcmin2": grid.pixel_area_arcmin2,
            "x_range_deg": [float(grid.x_edges[0]), float(grid.x_edges[-1])],
            "y_range_deg": [float(grid.y_edges[0]), float(grid.y_edges[-1])],
            "footprint_method": "convex hull of candidate tangent-plane positions with subpixel area sampling",
            "footprint_subpixels_per_axis": int(args.footprint_subpixels),
            "effective_coverage_area_deg2": float(np.nansum(grid.area_arcmin2) / 3600.0),
            "cells_with_positive_effective_area": int(np.sum(grid.area_arcmin2 > 0)),
            "distance_kpc_for_profile": float(args.distance_kpc),
            "profile_selection_max_radius_kpc": float(args.profile_selection_max_kpc),
            "profile_plot_max_radius_kpc": float(args.profile_plot_max_kpc),
            "outer_profile_contamination": outer_profile_contamination,
        },
        "prepared_observation_profile": prepared_observation_diag,
        "background": {k: v for k, v in bg.items() if k != "background_mask"},
        "smoothing": {
            "fwhm_deg": SMOOTH_FWHM_DEG,
            "sigma_deg": smooth["sigma_deg"],
            "sigma_pix": smooth["sigma_pix"],
            "kernel_size": smooth["kernel_size"],
            "kernel_sum": smooth["kernel_sum"],
            "kernel_squared_sum_2d": smooth["kernel_squared_sum_2d"],
            "coverage_threshold": float(args.coverage_threshold),
        },
        "calibration": {
            "munoz_center_mu_v_mag_arcsec2": MUNOZ_CENTER_MU_V,
            "munoz_center_radius_kpc": MUNOZ_CENTER_RADIUS_KPC,
            "munoz_rhalf_mu_v_mag_arcsec2": MUNOZ_RHALF_MU_V,
            "munoz_rhalf_radius_kpc": MUNOZ_RHALF_RADIUS_KPC,
            "munoz_mu_v_external_uncertainty_mag": MUNOZ_MU_V_ERR,
            "rhalf_anchor_half_width_kpc": sb["rhalf_anchor_half_width_kpc"],
            "anchor_measurements": sb["anchor_measurements"],
            "calibration_candidates": sb["calibration_candidates"],
            "selected_calibration": sb["selected_calibration"],
            "selected_anchor_predictions": sb["selected_anchor_predictions"],
            "yang_reference_radius_deg_for_validation": REFERENCE_RADIUS_DEG,
            "yang_sigma_ref_0p3_stars_arcmin2": sb["yang_sigma_ref_0p3"],
            "yang_sigma_ref_0p3_uncertainty_stars_arcmin2": sb["yang_sigma_ref_0p3_uncertainty"],
            "yang_s1_reference_density_stars_arcmin2": YANG_S1_REFERENCE_DENSITY,
            "fraction_of_yang_s1_reference_density": float(sb["yang_sigma_ref_0p3"] / YANG_S1_REFERENCE_DENSITY),
        },
        "detection": {
            "two_dimensional_smoothed_map": {
                "background_fluctuation_smoothed": sb["background_fluctuation_smoothed"],
                "density_3sigma_threshold_stars_arcmin2": sb["density_3sigma_threshold"],
                "mu_3sigma_threshold_mag_arcsec2": sb["mu_3sigma_threshold"],
                "faintest_reliable_mu_mag_arcsec2": sb["faintest_reliable_mu"],
                "contour_levels_mag_arcsec2": contour_level_values,
                "note": "Used only for the 2D smoothed contour map.",
            },
            "yang_style_1d_radial_profile": yang_1d_detection,
            "yang_style_1d_profile_summary": yang_1d_profile_summary,
            "background_fluctuation_smoothed": sb["background_fluctuation_smoothed"],
            "density_3sigma_threshold_stars_arcmin2": sb["density_3sigma_threshold"],
            "mu_3sigma_threshold_mag_arcsec2": sb["mu_3sigma_threshold"],
            "faintest_reliable_mu_mag_arcsec2": sb["faintest_reliable_mu"],
            "contour_levels_mag_arcsec2": contour_level_values,
        },
        "masking": {
            "total_cells": total_cells,
            "non_positive_net_density_cells": int(np.sum(non_positive)),
            "non_positive_net_density_fraction": float(np.sum(non_positive) / total_cells),
            "low_snr_3_cells": int(np.sum(low_snr)),
            "low_snr_3_fraction": float(np.sum(low_snr) / total_cells),
            "low_coverage_cells": int(np.sum(low_coverage)),
            "low_coverage_fraction": float(np.sum(low_coverage) / total_cells),
            "valid_science_cells": int(np.sum(reliable)),
            "valid_science_fraction": float(np.sum(reliable) / total_cells),
        },
        "validation_checks": {
            "density_units": "stars per square arcminute",
            "central_calibration_area_weighted": True,
            "background_subtracted_in_linear_density_space": True,
            "gaussian_smoothing_before_log_conversion": True,
            "radial_profile_uses_yang_style_1d_limit": True,
            "contour_map_uses_2d_smoothed_spatial_rms_limit": True,
            "no_log_for_non_positive_density": bool(np.all(~np.isfinite(np.asarray(sb["mu_v"])[np.asarray(sb["net_density"]) <= 0]))),
            "contours_above_3sigma_limit": bool(
                len(contour_level_values) == 0
                or max(contour_level_values) <= float(sb["mu_3sigma_threshold"]) + 1e-9
            ),
            "existing_radial_profile_validation": existing_profile_validation,
        },
        "published_value_comparison": {
            "yang_s1_reference_density_stars_arcmin2": YANG_S1_REFERENCE_DENSITY,
            "yang_style_smoothed_r_lt_0p3_density_stars_arcmin2": sb["yang_sigma_ref_0p3"],
            "yang_style_smoothed_r_lt_0p3_fraction_of_yang": float(sb["yang_sigma_ref_0p3"] / YANG_S1_REFERENCE_DENSITY),
            "existing_unsmoothed_profile_reference_density_stars_arcmin2": float(central_profile_density),
            "existing_unsmoothed_profile_fraction_of_yang": (
                float(central_profile_density / YANG_S1_REFERENCE_DENSITY)
                if np.isfinite(central_profile_density)
                else np.nan
            ),
            "yang_s1_background_density_stars_arcmin2": YANG_S1_BACKGROUND_DENSITY,
            "yang_s1_1d_background_uncertainty_stars_arcmin2": YANG_S1_1D_BACKGROUND_UNCERTAINTY,
            "yang_s1_1d_limit_mag_arcsec2": YANG_S1_1D_LIMIT_MAG_ARCSEC2,
            "measured_background_density_stars_arcmin2": bg["adopted_density"],
            "measured_yang_style_1d_background_uncertainty_stars_arcmin2": yang_1d_detection[
                "density_1sigma_background_stars_arcmin2"
            ],
            "measured_yang_style_1d_limit_mag_arcsec2": yang_1d_detection[
                "mu_3sigma_threshold_mag_arcsec2"
            ],
            "measured_minus_yang_s1_1d_limit_mag": (
                float(yang_1d_detection["mu_3sigma_threshold_mag_arcsec2"] - YANG_S1_1D_LIMIT_MAG_ARCSEC2)
                if np.isfinite(yang_1d_detection["mu_3sigma_threshold_mag_arcsec2"])
                else np.nan
            ),
            "existing_profile_background_density_stars_arcmin2": float(background_profile_density),
            "interpretation": (
                "The background agrees with the published S1 value. The final map is no "
                "longer normalized to the Yang r_ell<0.3 deg value; it is normalized by "
                "selecting the Muñoz calibration candidate with the smallest RMS residual "
                "against the central and half-light-radius Muñoz surface-brightness anchors."
            ),
        },
        "assumptions_and_limitations": [
            "The candidate file is an individual-star catalogue with no explicit survey mask or cell-area column.",
            "The 2D effective area is approximated by the convex hull of candidate positions; this can overestimate coverage if the true footprint has holes or concave edges.",
            "The published Yang S1 density and background are validation targets only and are not forced.",
            "The final map is calibrated from the smoothed background-subtracted density using Muñoz anchors.",
            "The 2D contour-map limit and Yang-style 1D radial-profile limit are reported separately.",
            "Surface-brightness values are not distance-modulus shifted; the 139.6 kpc distance is used only for angular-to-kpc radii.",
        ],
    }
    write_diagnostics(diagnostics_json, diagnostics)

    print(json.dumps(to_jsonable({
        "background_density": bg["adopted_density"],
        "background_uncertainty_on_mean": bg["adopted_uncertainty_on_mean"],
        "selected_calibration": sb["selected_calibration"]["name"],
        "selected_calibration_rms_mag": sb["selected_calibration"]["rms_residual_mag"],
        "yang_sigma_ref_0p3": sb["yang_sigma_ref_0p3"],
        "mu_2d_smoothed_3sigma_threshold": sb["mu_3sigma_threshold"],
        "mu_yang_style_1d_3sigma_threshold": yang_1d_detection["mu_3sigma_threshold_mag_arcsec2"],
        "yang_s1_reference_1d_limit": YANG_S1_1D_LIMIT_MAG_ARCSEC2,
        "faintest_reliable_mu": sb["faintest_reliable_mu"],
        "faintest_yang_style_1d_profile_mu_within_plot_radius": yang_1d_profile_summary[
            "faintest_detected_mu_within_plot_radius_mag_arcsec2"
        ],
        "outputs": diagnostics["outputs"],
    }), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
