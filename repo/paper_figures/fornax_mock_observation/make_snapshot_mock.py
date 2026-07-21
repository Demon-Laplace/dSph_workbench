"""Render one GIZMO snapshot as a paper-ready mock stellar observation.

The particle selection follows the existing Fornax analysis pipeline:

1. centre the simulation on the Milky Way with a shrinking sphere;
2. locate the densest stellar system outside the Milky Way;
3. select the dwarf within three Fornax core radii;
4. transform to ICRS and project onto a local R.A./Dec. tangent plane;
5. select gas within 20 kpc of the dwarf stellar centre.

The display layer is deliberately different from the diagnostic PlotFig.py
panels: stellar particles are converted to a PSF-smoothed light map, while
total gas and cold neutral gas are shown as blue emission and cyan contours.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEPS = HERE / ".deps"
if DEPS.exists():
    sys.path.insert(0, str(DEPS))

# Keep rendering headless and all caches inside this self-contained demo.
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".mplconfig"))

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import (
    CartesianDifferential,
    CartesianRepresentation,
    Galactocentric,
    Galactic,
    ICRS,
    SkyCoord,
)
from astropy.utils import iers
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree


BOX_CENTER_FRACTION = 0.5
FORNAX_CORE_RADIUS_KPC = 5.4
DWARF_RADIUS_FACTOR = 3.0
MW_EXCLUSION_RADIUS_KPC = 5.0
K_DENSITY = 16
MASS_TO_LIGHT = 2.6
SOLAR_ABS_MAG_V = 4.83
SUN_X_KPC = 8.122
SUN_Z_KPC = 0.0208
FIELD_HALF_WIDTH_DEG = 2.1
DWARF_GAS_RADIUS_KPC = 20.0
COLD_GAS_TEMPERATURE_K = 2.0e4
HI_MASS_FIELD_BUFFER = 1.10
ADOPTED_DISTANCE_KPC = 139.6
SNAPSHOT_CADENCE_GYR = 0.01

iers.conf.auto_download = False


def read_particles(snapshot_path: Path) -> tuple[dict[int, dict[str, np.ndarray]], float, float]:
    """Read the stellar particle fields required for imaging and proper motion."""
    particles: dict[int, dict[str, np.ndarray]] = {}
    with h5py.File(snapshot_path, "r") as handle:
        header = handle["Header"].attrs
        box_size = float(header["BoxSize"])
        time = float(header["Time"])
        offset = BOX_CENTER_FRACTION * box_size

        for ptype in (2, 3, 4):
            group_name = f"PartType{ptype}"
            if group_name not in handle:
                continue
            group = handle[group_name]
            particles[ptype] = {
                "position": np.asarray(group["Coordinates"], dtype=np.float64) - offset,
                "mass": np.asarray(group["Masses"], dtype=np.float64) * 1.0e10,
                "velocity": np.asarray(group["Velocities"], dtype=np.float64),
            }

    if not particles:
        raise RuntimeError("No stellar particle groups (PartType2/3/4) were found.")
    return particles, box_size, time


def shrinking_sphere_mw_center(
    particles: dict[int, dict[str, np.ndarray]],
    r_init: float = 300.0,
    r_min: float = 5.0,
    shrink: float = 0.8,
) -> np.ndarray:
    """Match the server pipeline's mass-weighted Milky Way centring."""
    selected = [particles[p] for p in (2, 4) if p in particles]
    positions = np.concatenate([item["position"] for item in selected])
    masses = np.concatenate([item["mass"] for item in selected])
    center = np.average(positions, axis=0, weights=masses)

    radius = r_init
    while radius > r_min:
        distances = np.linalg.norm(positions - center, axis=1)
        inside = distances < radius
        if np.count_nonzero(inside) < 100:
            break
        center = np.average(positions[inside], axis=0, weights=masses[inside])
        radius *= shrink
    return center


def concatenate_stars(
    particles: dict[int, dict[str, np.ndarray]], mw_center: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    positions = []
    masses = []
    particle_types = []
    velocities = []
    for ptype in sorted(particles):
        count = particles[ptype]["position"].shape[0]
        positions.append(particles[ptype]["position"] - mw_center)
        masses.append(particles[ptype]["mass"])
        particle_types.append(np.full(count, ptype, dtype=np.int8))
        velocities.append(particles[ptype]["velocity"])
    return (
        np.concatenate(positions),
        np.concatenate(masses),
        np.concatenate(particle_types),
        np.concatenate(velocities),
    )


def snapshot_time_gyr(snapshot_path: Path, header_time_gyr: float) -> tuple[int | None, float]:
    """Use the run's 0.01 Gyr snapshot cadence, falling back to Header/Time."""
    match = re.search(r"snapshot_(\d+)", snapshot_path.stem)
    if match is None:
        return None, float(header_time_gyr)
    snapshot_number = int(match.group(1))
    return snapshot_number, SNAPSHOT_CADENCE_GYR * snapshot_number


def mean_proper_motion(
    positions_kpc: np.ndarray, velocities_kms: np.ndarray
) -> tuple[float, float]:
    """Return the simulated mean ICRS proper motion, as in the original PlotFig.py."""
    position = CartesianRepresentation(
        x=positions_kpc[:, 0] * u.kpc,
        y=positions_kpc[:, 1] * u.kpc,
        z=positions_kpc[:, 2] * u.kpc,
    )
    velocity = CartesianDifferential(
        d_x=velocities_kms[:, 0] * u.km / u.s,
        d_y=velocities_kms[:, 1] * u.km / u.s,
        d_z=velocities_kms[:, 2] * u.km / u.s,
    )
    coordinate = SkyCoord(position.with_differentials(velocity), frame="galactocentric")
    icrs = coordinate.transform_to(ICRS())
    return (
        float(np.mean(icrs.pm_ra_cosdec.to_value(u.mas / u.yr))),
        float(np.mean(icrs.pm_dec.to_value(u.mas / u.yr))),
    )


def find_dwarf_mask(
    positions: np.ndarray,
    exclusion_radius_kpc: float = MW_EXCLUSION_RADIUS_KPC,
    dwarf_radius_kpc: float = DWARF_RADIUS_FACTOR * FORNAX_CORE_RADIUS_KPC,
    k_density: int = K_DENSITY,
) -> tuple[np.ndarray, np.ndarray]:
    """Locate the outer stellar density peak using the existing kNN idea."""
    galactocentric_radius = np.linalg.norm(positions, axis=1)
    outer_indices = np.flatnonzero(galactocentric_radius > exclusion_radius_kpc)
    if outer_indices.size <= k_density:
        raise RuntimeError("Too few stellar particles outside the Milky Way exclusion radius.")

    tree = cKDTree(positions)
    # Asking only for the kth neighbour avoids allocating an N x k matrix.
    kth_distance = tree.query(
        positions[outer_indices], k=[k_density + 1], workers=-1
    )[0][:, 0]
    peak_index = outer_indices[np.argmin(kth_distance)]
    peak = positions[peak_index]
    dwarf_mask = np.linalg.norm(positions - peak, axis=1) < dwarf_radius_kpc
    return dwarf_mask, peak


def gas_temperature(internal_energy: np.ndarray, electron_abundance: np.ndarray) -> np.ndarray:
    """Use the same internal-energy-to-temperature conversion as zim.py."""
    boltzmann = 1.3806e-16
    proton_mass = 1.6726e-24
    hydrogen_fraction = 0.76
    gamma = 5.0 / 3.0
    velocity_unit_squared = 1.0e10
    mean_weight = (
        4.0
        / (3.0 * hydrogen_fraction + 1.0 + 4.0 * hydrogen_fraction * electron_abundance)
        * proton_mass
    )
    return mean_weight / boltzmann * (gamma - 1.0) * internal_energy * velocity_unit_squared


def read_dwarf_gas(
    snapshot_path: Path,
    box_size: float,
    mw_center: np.ndarray,
    dwarf_center: np.ndarray,
    radius_kpc: float = DWARF_GAS_RADIUS_KPC,
    chunk_size: int = 400_000,
) -> dict[str, np.ndarray]:
    """Read gas in chunks and retain only particles near the dwarf."""
    position_chunks = []
    mass_chunks = []
    neutral_chunks = []
    temperature_chunks = []
    smoothing_length_chunks = []

    with h5py.File(snapshot_path, "r") as handle:
        if "PartType0" not in handle:
            return {
                "position": np.empty((0, 3)),
                "mass": np.empty(0),
                "neutral_fraction": np.empty(0),
                "temperature": np.empty(0),
                "smoothing_length": np.empty(0),
            }
        group = handle["PartType0"]
        count = int(group["Coordinates"].shape[0])
        offset = BOX_CENTER_FRACTION * box_size

        for start in range(0, count, chunk_size):
            stop = min(start + chunk_size, count)
            positions = (
                np.asarray(group["Coordinates"][start:stop], dtype=np.float64)
                - offset
                - mw_center
            )
            nearby = np.sum((positions - dwarf_center) ** 2, axis=1) < radius_kpc**2
            if not np.any(nearby):
                continue

            masses = np.asarray(group["Masses"][start:stop], dtype=np.float64)[nearby] * 1.0e10
            if "NeutralHydrogenAbundance" in group:
                neutral = np.asarray(
                    group["NeutralHydrogenAbundance"][start:stop], dtype=np.float64
                )[nearby]
            else:
                neutral = np.zeros(np.count_nonzero(nearby), dtype=np.float64)
            internal_energy = np.asarray(
                group["InternalEnergy"][start:stop], dtype=np.float64
            )[nearby]
            if "ElectronAbundance" in group:
                electron_abundance = np.asarray(
                    group["ElectronAbundance"][start:stop], dtype=np.float64
                )[nearby]
            else:
                electron_abundance = np.zeros(np.count_nonzero(nearby), dtype=np.float64)

            position_chunks.append(positions[nearby])
            mass_chunks.append(masses)
            neutral_chunks.append(neutral)
            temperature_chunks.append(gas_temperature(internal_energy, electron_abundance))
            if "SmoothingLength" in group:
                smoothing_length_chunks.append(
                    np.asarray(group["SmoothingLength"][start:stop], dtype=np.float64)[nearby]
                )
            else:
                smoothing_length_chunks.append(
                    np.full(np.count_nonzero(nearby), 0.05, dtype=np.float64)
                )

    if not position_chunks:
        return {
            "position": np.empty((0, 3)),
            "mass": np.empty(0),
            "neutral_fraction": np.empty(0),
            "temperature": np.empty(0),
            "smoothing_length": np.empty(0),
        }
    return {
        "position": np.concatenate(position_chunks),
        "mass": np.concatenate(mass_chunks),
        "neutral_fraction": np.concatenate(neutral_chunks),
        "temperature": np.concatenate(temperature_chunks),
        "smoothing_length": np.concatenate(smoothing_length_chunks),
    }


def positions_to_icrs(
    positions: np.ndarray, chunk_size: int = 100_000
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Match the server pipeline's Galactocentric-to-ICRS transformation."""
    frame = Galactocentric()
    x_sun = frame.galcen_distance.to(u.kpc).value
    z_sun = frame.z_sun.to(u.kpc).value
    heliocentric = positions - np.array([x_sun, 0.0, z_sun])
    distance = np.linalg.norm(heliocentric, axis=1)

    ra_chunks = []
    dec_chunks = []
    for start in range(0, positions.shape[0], chunk_size):
        stop = min(start + chunk_size, positions.shape[0])
        coords = Galactocentric(
            x=heliocentric[start:stop, 0] * u.kpc,
            y=heliocentric[start:stop, 1] * u.kpc,
            z=heliocentric[start:stop, 2] * u.kpc,
        ).transform_to(Galactic())
        icrs = coords.transform_to(ICRS())
        ra_chunks.append(icrs.ra.deg)
        dec_chunks.append(icrs.dec.deg)

    return np.concatenate(ra_chunks), np.concatenate(dec_chunks), distance, heliocentric


def find_sky_center(
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
    shrink_rate: float = 0.95,
    threshold_deg: float = 0.1,
) -> tuple[float, float]:
    """Match Analysis.find_center_2d(..., units='degree')."""
    histogram, ra_edges, dec_edges = np.histogram2d(
        ra_deg,
        dec_deg,
        bins=[360, 180],
        range=[[0.0, 360.0], [-90.0, 90.0]],
    )
    i, j = np.unravel_index(np.argmax(histogram), histogram.shape)
    center_ra = float(np.mean(ra_edges[i : i + 2]))
    center_dec = float(np.mean(dec_edges[j : j + 2]))

    range_ra = FORNAX_CORE_RADIUS_KPC / 3.0
    range_dec = FORNAX_CORE_RADIUS_KPC / 3.0
    while range_ra > threshold_deg and range_dec > threshold_deg:
        inside = (
            (ra_deg >= center_ra - range_ra)
            & (ra_deg <= center_ra + range_ra)
            & (dec_deg >= center_dec - range_dec)
            & (dec_deg <= center_dec + range_dec)
        )
        if not np.any(inside):
            break
        center_ra = float(np.median(ra_deg[inside]))
        center_dec = float(np.median(dec_deg[inside]))
        range_ra *= shrink_rate
        range_dec *= shrink_rate
    return center_ra, center_dec


def sky_offsets(
    ra_deg: np.ndarray, dec_deg: np.ndarray, center: tuple[float, float]
) -> tuple[np.ndarray, np.ndarray]:
    """Return local angular offsets in degrees, matching rotate_to_sky()."""
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    center_ra, center_dec = np.deg2rad(center)
    delta_ra = np.cos(dec) * np.sin(ra - center_ra)
    delta_dec = (
        np.sin(dec) * np.cos(center_dec)
        - np.cos(dec) * np.sin(center_dec) * np.cos(ra - center_ra)
    )
    return np.rad2deg(delta_ra), np.rad2deg(delta_dec)


def project_to_sky(
    positions: np.ndarray,
    center: tuple[float, float] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, tuple[float, float]]:
    ra_deg, dec_deg, distance, heliocentric = positions_to_icrs(positions)
    if center is None:
        center = find_sky_center(ra_deg, dec_deg)
    x_deg, y_deg = sky_offsets(ra_deg, dec_deg, center)
    return x_deg, y_deg, distance, heliocentric, center


def stellar_colormap() -> LinearSegmentedColormap:
    """Black-to-ivory palette resembling deep resolved-stellar imaging."""
    return LinearSegmentedColormap.from_list(
        "fornax_starlight",
        [
            (0.00, "#000000"),
            (0.10, "#060814"),
            (0.28, "#171a3a"),
            (0.48, "#49375b"),
            (0.68, "#9b6969"),
            (0.84, "#e3b78e"),
            (1.00, "#fff7dc"),
        ],
        N=256,
    )


def gas_colormap() -> LinearSegmentedColormap:
    """Deep blue-to-cyan palette for projected gas surface density."""
    return LinearSegmentedColormap.from_list(
        "fornax_gas",
        [
            (0.00, "#000000"),
            (0.18, "#031227"),
            (0.42, "#073d66"),
            (0.70, "#087f9d"),
            (1.00, "#67e8ef"),
        ],
        N=256,
    )


def make_light_map(
    x_deg: np.ndarray,
    y_deg: np.ndarray,
    masses: np.ndarray,
    distance_kpc: float,
    half_width_deg: float,
    npix: int,
    smoothing_sigma_pixels: float,
    mass_to_light: float = MASS_TO_LIGHT,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return a PSF-smoothed V-band surface-brightness map in mag/arcsec^2."""
    limits = [[-half_width_deg, half_width_deg], [-half_width_deg, half_width_deg]]
    mass_image, x_edges, y_edges = np.histogram2d(
        x_deg, y_deg, bins=npix, range=limits, weights=masses
    )

    fine = gaussian_filter(mass_image, smoothing_sigma_pixels, mode="constant")
    diffuse = gaussian_filter(mass_image, 3.5 * smoothing_sigma_pixels, mode="constant")
    smoothed_mass = (0.78 * fine) + (0.22 * diffuse)

    pixel_deg = 2.0 * half_width_deg / npix
    pixel_area_arcsec2 = (pixel_deg * 3600.0) ** 2

    surface_brightness = np.full_like(smoothed_mass, np.inf, dtype=np.float64)
    positive = smoothed_mass > 0.0
    luminosity = smoothed_mass[positive] / mass_to_light
    absolute_magnitude = SOLAR_ABS_MAG_V - 2.5 * np.log10(luminosity)
    distance_modulus = 5.0 * np.log10(distance_kpc * 1000.0) - 5.0
    apparent_magnitude = absolute_magnitude + distance_modulus
    surface_brightness[positive] = apparent_magnitude + 2.5 * np.log10(pixel_area_arcsec2)
    return surface_brightness, x_edges, y_edges, smoothed_mass


def adaptive_sph_map(
    x_deg: np.ndarray,
    y_deg: np.ndarray,
    weights: np.ndarray,
    smoothing_length_kpc: np.ndarray,
    distance_kpc: float,
    half_width_deg: float,
    npix: int,
    minimum_sigma_pixels: float,
) -> np.ndarray:
    """Approximate adaptive SPH projection by smoothing h-binned mass maps."""
    limits = [[-half_width_deg, half_width_deg], [-half_width_deg, half_width_deg]]
    inside = (
        np.isfinite(x_deg)
        & np.isfinite(y_deg)
        & np.isfinite(weights)
        & np.isfinite(smoothing_length_kpc)
        & (weights > 0.0)
        & (np.abs(x_deg) <= half_width_deg)
        & (np.abs(y_deg) <= half_width_deg)
    )
    if not np.any(inside):
        return np.zeros((npix, npix), dtype=np.float64)

    x = x_deg[inside]
    y = y_deg[inside]
    w = weights[inside]
    h = np.clip(smoothing_length_kpc[inside], 1.0e-4, None)
    quantile_edges = np.unique(np.quantile(h, np.linspace(0.0, 1.0, 9)))
    if quantile_edges.size < 2:
        quantile_edges = np.array([h.min(), np.nextafter(h.max(), np.inf)])

    pixel_deg = 2.0 * half_width_deg / npix
    projected_mass = np.zeros((npix, npix), dtype=np.float64)
    for index in range(quantile_edges.size - 1):
        lower = quantile_edges[index]
        upper = quantile_edges[index + 1]
        if index == quantile_edges.size - 2:
            group = (h >= lower) & (h <= upper)
        else:
            group = (h >= lower) & (h < upper)
        if not np.any(group):
            continue
        mass_image, _, _ = np.histogram2d(
            x[group], y[group], bins=npix, range=limits, weights=w[group]
        )
        # GIZMO's h is a kernel scale; sigma=h/2 gives a smooth projected proxy.
        sigma_deg = np.rad2deg(np.median(h[group]) / distance_kpc) / 2.0
        sigma_pixels = np.clip(
            sigma_deg / pixel_deg, minimum_sigma_pixels, npix / 3.0
        )
        projected_mass += gaussian_filter(mass_image, sigma_pixels, mode="constant")
    return projected_mass


def make_hi_map(
    x_deg: np.ndarray,
    y_deg: np.ndarray,
    masses: np.ndarray,
    neutral_fraction: np.ndarray,
    temperature: np.ndarray,
    smoothing_length_kpc: np.ndarray,
    distance_kpc: float,
    half_width_deg: float,
    npix: int,
    smoothing_sigma_pixels: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the NH-weighted cold H I surface density in Msun/pc^2."""
    cold = temperature < COLD_GAS_TEMPERATURE_K
    neutral_mass = adaptive_sph_map(
        x_deg[cold],
        y_deg[cold],
        masses[cold] * neutral_fraction[cold],
        smoothing_length_kpc[cold],
        distance_kpc,
        half_width_deg,
        npix,
        0.75 * smoothing_sigma_pixels,
    )

    pixel_deg = 2.0 * half_width_deg / npix
    pixel_kpc = np.deg2rad(pixel_deg) * distance_kpc
    pixel_area_pc2 = (pixel_kpc * 1000.0) ** 2
    return neutral_mass / pixel_area_pc2, cold


def hi_to_display(surface_density: np.ndarray) -> np.ndarray:
    """Compress H I dynamic range without applying a density cut."""
    peak = float(np.nanmax(surface_density)) if surface_density.size else 0.0
    if not np.isfinite(peak) or peak <= 0.0:
        return np.zeros_like(surface_density)
    display = np.log1p(250.0 * np.clip(surface_density / peak, 0.0, None)) / np.log1p(250.0)
    return display


def magnitude_to_display(
    surface_brightness: np.ndarray,
    bright_limit: float,
    faint_limit: float,
) -> np.ndarray:
    """Map magnitudes to [0, 1], preserving faint tidal structure."""
    brightness = (faint_limit - surface_brightness) / (faint_limit - bright_limit)
    # Empty pixels have mu=+inf and therefore brightness=-inf: keep them black.
    # The opposite limit (formally mu=-inf) is saturated at the bright end.
    brightness = np.nan_to_num(brightness, nan=0.0, posinf=1.0, neginf=0.0)
    brightness = np.clip(brightness, 0.0, 1.0)
    return np.arcsinh(4.5 * brightness) / np.arcsinh(4.5)


def style_axis(axis: plt.Axes, half_width_deg: float) -> None:
    axis.set_facecolor("black")
    for spine in axis.spines.values():
        spine.set_color("#747474")
        spine.set_linewidth(0.55)
    tick_limit = int(np.floor(half_width_deg))
    ticks = np.arange(-tick_limit, tick_limit + 1, 1)
    axis.set_xticks(ticks)
    axis.set_yticks(ticks)
    axis.tick_params(colors="#d4d4d4", width=0.55, length=3.0, labelsize=6.8, direction="in")
    axis.set_xlabel("RA", color="#e5e5e5", fontsize=7.4)
    axis.set_ylabel("Dec", color="#e5e5e5", fontsize=7.4)
    # Astronomical images conventionally show R.A. increasing to the left.
    axis.set_xlim(half_width_deg, -half_width_deg)
    axis.set_ylim(-half_width_deg, half_width_deg)
    axis.set_aspect("equal")


def mass_to_tex(mass_msun: float) -> str:
    """Format a mass compactly for an in-panel math label."""
    if not np.isfinite(mass_msun) or mass_msun <= 0.0:
        return "0"
    exponent = int(np.floor(np.log10(mass_msun)))
    mantissa = mass_msun / (10.0**exponent)
    return rf"{mantissa:.2f}\times10^{{{exponent}}}"


def draw_proper_motion_arrow(
    axis: plt.Axes,
    pmra_masyr: float,
    pmdec_masyr: float,
    arrow_length_deg: float = 0.72,
) -> None:
    """Draw the projected motion direction; +RA points left on the reversed axis."""
    norm = float(np.hypot(pmra_masyr, pmdec_masyr))
    if not np.isfinite(norm) or norm == 0.0:
        return
    dx = pmra_masyr / norm * arrow_length_deg
    dy = pmdec_masyr / norm * arrow_length_deg
    axis.arrow(
        0.0,
        0.0,
        dx,
        dy,
        width=0.009,
        head_width=0.105,
        head_length=0.13,
        length_includes_head=True,
        facecolor="#f1f3f2",
        edgecolor="#17191b",
        linewidth=0.35,
        alpha=0.92,
        zorder=18,
    )


def save_mock_panel(
    display_image: np.ndarray,
    hi_display: np.ndarray,
    hi_surface_density: np.ndarray,
    half_width_deg: float,
    time_gyr: float,
    stellar_mass_msun: float,
    hi_mass_msun: float,
    distance_kpc: float,
    pmra_masyr: float,
    pmdec_masyr: float,
    output_png: Path,
    output_pdf: Path,
    dpi: int,
) -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 8.0,
            "savefig.facecolor": "black",
            "figure.facecolor": "black",
        }
    )
    fig, axis = plt.subplots(figsize=(3.45, 3.45))
    extent = [-half_width_deg, half_width_deg, -half_width_deg, half_width_deg]
    hi_alpha = 0.46 * np.power(np.clip(hi_display.T, 0.0, 1.0), 0.72)
    axis.imshow(
        hi_display.T,
        origin="lower",
        extent=extent,
        interpolation="bilinear",
        cmap=gas_colormap(),
        vmin=0.0,
        vmax=1.0,
        alpha=hi_alpha,
        rasterized=True,
    )

    stellar_rgba = stellar_colormap()(np.clip(display_image.T, 0.0, 1.0))
    stellar_rgba[..., 3] = np.power(np.clip(display_image.T, 0.0, 1.0), 0.62)
    axis.imshow(
        stellar_rgba,
        origin="lower",
        extent=extent,
        interpolation="bilinear",
        rasterized=True,
    )
    hi_peak = float(np.nanmax(hi_surface_density))
    if np.isfinite(hi_peak) and hi_peak > 0.0:
        axis.contour(
            hi_surface_density.T,
            levels=hi_peak * np.array([0.02, 0.10, 0.35, 0.70]),
            origin="lower",
            extent=extent,
            colors=["#249fb4", "#42cbd3", "#72ebeb", "#d3ffff"],
            linewidths=[0.42, 0.52, 0.64, 0.78],
            alpha=0.92,
        )
    style_axis(axis, half_width_deg)
    draw_proper_motion_arrow(axis, pmra_masyr, pmdec_masyr)
    axis.text(
        0.035,
        0.965,
        rf"$t={time_gyr:.2f}\,\mathrm{{Gyr}}$",
        transform=axis.transAxes,
        color="#f4f1ea",
        fontsize=7.7,
        ha="left",
        va="top",
    )
    axis.text(
        0.035,
        0.910,
        rf"$M_\star={mass_to_tex(stellar_mass_msun)}\,M_\odot$",
        transform=axis.transAxes,
        color="#f0dfcf",
        fontsize=6.3,
        ha="left",
        va="top",
    )
    axis.text(
        0.035,
        0.865,
        rf"$M_{{\mathrm{{H\,I}}}}={mass_to_tex(hi_mass_msun)}\,M_\odot$",
        transform=axis.transAxes,
        color="#94f3f3",
        fontsize=6.0,
        ha="left",
        va="top",
    )
    axis.text(
        0.035,
        0.820,
        rf"$D_\odot={distance_kpc:.2f}\,\mathrm{{kpc}}$",
        transform=axis.transAxes,
        color="#d7dadd",
        fontsize=5.9,
        ha="left",
        va="top",
    )
    legend = axis.legend(
        handles=[
            Line2D([0], [0], color="#a7ffff", lw=0.9, label="HI cloud"),
        ],
        loc="upper right",
        frameon=False,
        fontsize=6.1,
        handlelength=1.45,
        handletextpad=0.45,
        borderpad=0.15,
        labelcolor="#e8eeee",
    )
    legend.set_zorder(20)
    fig.subplots_adjust(left=0.15, right=0.985, bottom=0.14, top=0.985)
    fig.savefig(output_png, dpi=dpi, facecolor="black", bbox_inches="tight", pad_inches=0.01)
    fig.savefig(output_pdf, dpi=dpi, facecolor="black", bbox_inches="tight", pad_inches=0.01)
    plt.close(fig)


def save_diagnostic(
    surface_brightness: np.ndarray,
    hi_surface_density: np.ndarray,
    half_width_deg: float,
    bright_limit: float,
    faint_limit: float,
    output_path: Path,
    dpi: int,
) -> None:
    masked = np.ma.masked_greater(surface_brightness, faint_limit)
    cmap = stellar_colormap().reversed()
    cmap.set_bad("black")
    fig, axis = plt.subplots(figsize=(4.35, 3.75), facecolor="black")
    axis.set_facecolor("black")
    extent = [-half_width_deg, half_width_deg, -half_width_deg, half_width_deg]
    image = axis.imshow(
        masked.T,
        origin="lower",
        extent=extent,
        cmap=cmap,
        vmin=bright_limit,
        vmax=faint_limit,
        interpolation="bilinear",
    )
    hi_peak = float(np.nanmax(hi_surface_density))
    if np.isfinite(hi_peak) and hi_peak > 0.0:
        axis.contour(
            hi_surface_density.T,
            levels=hi_peak * np.array([0.02, 0.10, 0.35, 0.70]),
            origin="lower",
            extent=extent,
            colors="#b7ffff",
            linewidths=[0.42, 0.52, 0.64, 0.78],
        )
    axis.set_xlabel("RA", color="#dedede")
    axis.set_ylabel("Dec", color="#dedede")
    axis.tick_params(colors="#cfcfcf", width=0.6, labelsize=7)
    for spine in axis.spines.values():
        spine.set_color("#777777")
        spine.set_linewidth(0.6)
    axis.set_xlim(half_width_deg, -half_width_deg)
    axis.set_ylim(-half_width_deg, half_width_deg)
    axis.set_aspect("equal")
    colorbar = fig.colorbar(image, ax=axis, fraction=0.046, pad=0.035)
    colorbar.set_label(r"$\mu_V\ [\mathrm{mag\,arcsec^{-2}}]$", color="#dedede")
    colorbar.ax.tick_params(colors="#cfcfcf", labelsize=7)
    colorbar.ax.invert_yaxis()
    fig.tight_layout(pad=0.7)
    fig.savefig(output_path, dpi=dpi, facecolor="black", bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path, help="Path to a GIZMO HDF5 snapshot")
    parser.add_argument("--outdir", type=Path, default=HERE / "output")
    parser.add_argument("--npix", type=int, default=720)
    parser.add_argument("--dpi", type=int, default=450)
    parser.add_argument("--field-half-deg", type=float, default=FIELD_HALF_WIDTH_DEG)
    parser.add_argument("--smooth-pixels", type=float, default=2.0)
    parser.add_argument("--gas-smooth-pixels", type=float, default=5.0)
    parser.add_argument("--bright-limit", type=float, default=22.0)
    parser.add_argument("--faint-limit", type=float, default=34.0)
    parser.add_argument(
        "--distance-kpc",
        type=float,
        default=ADOPTED_DISTANCE_KPC,
        help="Adopted heliocentric distance for angular and photometric scaling",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    snapshot_path = args.snapshot.resolve()
    output_dir = args.outdir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    pdf_dir = output_dir / "pdf"
    pdf_dir.mkdir(parents=True, exist_ok=True)

    particles, box_size, header_time_gyr = read_particles(snapshot_path)
    snapshot_number, time_gyr = snapshot_time_gyr(snapshot_path, header_time_gyr)
    mw_center = shrinking_sphere_mw_center(particles)
    positions, masses, particle_types, velocities = concatenate_stars(particles, mw_center)
    dwarf_mask, dwarf_peak = find_dwarf_mask(positions)
    dwarf_center = np.average(positions[dwarf_mask], axis=0, weights=masses[dwarf_mask])
    gas = read_dwarf_gas(
        snapshot_path=snapshot_path,
        box_size=box_size,
        mw_center=mw_center,
        dwarf_center=dwarf_center,
    )

    # Match the analysis pipeline's old-star definition: pre-existing types 2/3.
    old_dwarf = dwarf_mask & (particle_types != 4)
    selected_positions = positions[old_dwarf]
    selected_masses = masses[old_dwarf]
    selected_velocities = velocities[old_dwarf]
    if selected_positions.shape[0] == 0:
        raise RuntimeError("The dwarf selection contains no old stellar particles.")

    x_deg, y_deg, stellar_distance, stellar_heliocentric, sky_center = project_to_sky(
        selected_positions
    )
    transformed_distance_kpc = float(np.average(stellar_distance, weights=selected_masses))
    distance_kpc = float(args.distance_kpc)
    angular_rescale = transformed_distance_kpc / distance_kpc
    x_deg *= angular_rescale
    y_deg *= angular_rescale
    heliocentric_center = np.average(
        stellar_heliocentric, axis=0, weights=selected_masses
    )
    pmra_masyr, pmdec_masyr = mean_proper_motion(
        selected_positions, selected_velocities
    )
    half_width_kpc = float(np.tan(np.deg2rad(args.field_half_deg)) * distance_kpc)
    surface_brightness, _, _, smoothed_mass = make_light_map(
        x_deg=x_deg,
        y_deg=y_deg,
        masses=selected_masses,
        distance_kpc=distance_kpc,
        half_width_deg=args.field_half_deg,
        npix=args.npix,
        smoothing_sigma_pixels=args.smooth_pixels,
    )
    display_image = magnitude_to_display(
        surface_brightness, args.bright_limit, args.faint_limit
    )

    if gas["position"].shape[0] > 0:
        gas_x_deg, gas_y_deg, _, _, _ = project_to_sky(
            gas["position"], center=sky_center
        )
        gas_x_deg *= angular_rescale
        gas_y_deg *= angular_rescale
        hi_surface_density, cold_gas = make_hi_map(
            x_deg=gas_x_deg,
            y_deg=gas_y_deg,
            masses=gas["mass"],
            neutral_fraction=gas["neutral_fraction"],
            temperature=gas["temperature"],
            smoothing_length_kpc=gas["smoothing_length"],
            distance_kpc=distance_kpc,
            half_width_deg=args.field_half_deg,
            npix=args.npix,
            smoothing_sigma_pixels=args.gas_smooth_pixels,
        )
    else:
        gas_x_deg = np.empty(0)
        gas_y_deg = np.empty(0)
        cold_gas = np.zeros(0, dtype=bool)
        hi_surface_density = np.zeros_like(surface_brightness)
    hi_display = hi_to_display(hi_surface_density)
    stellar_mass_msun = float(selected_masses.sum())
    hi_mass_half_width_deg = HI_MASS_FIELD_BUFFER * args.field_half_deg
    hi_mass_aperture = (
        cold_gas
        & (np.abs(gas_x_deg) <= hi_mass_half_width_deg)
        & (np.abs(gas_y_deg) <= hi_mass_half_width_deg)
    )
    hi_mass_msun = float(
        np.sum(gas["mass"][hi_mass_aperture] * gas["neutral_fraction"][hi_mass_aperture])
    )

    stem = snapshot_path.stem
    panel_png = output_dir / f"{stem}_mock_panel.png"
    panel_pdf = pdf_dir / f"{stem}_mock_panel.pdf"
    diagnostic_png = output_dir / f"{stem}_surface_brightness_diagnostic.png"
    metadata_json = output_dir / f"{stem}_mock_metadata.json"

    save_mock_panel(
        display_image,
        hi_display,
        hi_surface_density,
        args.field_half_deg,
        time_gyr,
        stellar_mass_msun,
        hi_mass_msun,
        distance_kpc,
        pmra_masyr,
        pmdec_masyr,
        panel_png,
        panel_pdf,
        args.dpi,
    )
    save_diagnostic(
        surface_brightness,
        hi_surface_density,
        args.field_half_deg,
        args.bright_limit,
        args.faint_limit,
        diagnostic_png,
        args.dpi,
    )

    finite_mu = surface_brightness[np.isfinite(surface_brightness)]
    metadata = {
        "snapshot": str(snapshot_path),
        "snapshot_number": snapshot_number,
        "snapshot_time_gyr": time_gyr,
        "header_time_gyr": header_time_gyr,
        "box_size_kpc": box_size,
        "mw_center_raw_kpc": mw_center.tolist(),
        "dwarf_density_peak_galactocentric_kpc": dwarf_peak.tolist(),
        "dwarf_stellar_center_galactocentric_kpc": dwarf_center.tolist(),
        "dwarf_stellar_center_galactocentric_radius_kpc": float(
            np.linalg.norm(dwarf_center)
        ),
        "dwarf_heliocentric_center_kpc": heliocentric_center.tolist(),
        "transformed_particle_distance_kpc": transformed_distance_kpc,
        "adopted_heliocentric_distance_kpc": distance_kpc,
        "distance_reference": "Li et al. (2021), distance modulus 20.72 mag",
        "angular_rescale_factor": angular_rescale,
        "mean_proper_motion_mas_per_yr": {
            "pmra_cosdec": pmra_masyr,
            "pmdec": pmdec_masyr,
        },
        "sky_center_icrs_deg": {"ra": sky_center[0], "dec": sky_center[1]},
        "selection_radius_kpc": DWARF_RADIUS_FACTOR * FORNAX_CORE_RADIUS_KPC,
        "selected_old_star_particles": int(np.count_nonzero(old_dwarf)),
        "selected_old_stellar_mass_msun": stellar_mass_msun,
        "field_half_width_deg": args.field_half_deg,
        "field_half_width_kpc": half_width_kpc,
        "ra_axis_increases_to_left": True,
        "gas": {
            "rendered_component": "NH-weighted cold H I only",
            "selection_radius_kpc": DWARF_GAS_RADIUS_KPC,
            "particle_count": int(gas["mass"].size),
            "total_mass_msun": float(gas["mass"].sum()),
            "cold_temperature_threshold_k": COLD_GAS_TEMPERATURE_K,
            "cold_particle_count": int(np.count_nonzero(cold_gas)),
            "cold_mass_msun": float(gas["mass"][cold_gas].sum()),
            "cold_neutral_mass_msun": hi_mass_msun,
            "mass_aperture": {
                "shape": "projected square",
                "field_buffer_factor": HI_MASS_FIELD_BUFFER,
                "half_width_deg": hi_mass_half_width_deg,
                "particle_count": int(np.count_nonzero(hi_mass_aperture)),
            },
        },
        "surface_brightness_percentiles_mag_arcsec2": {
            "p01": float(np.percentile(finite_mu, 1)),
            "p10": float(np.percentile(finite_mu, 10)),
            "p50": float(np.percentile(finite_mu, 50)),
            "p90": float(np.percentile(finite_mu, 90)),
            "p99": float(np.percentile(finite_mu, 99)),
        },
        "render": {
            "npix": args.npix,
            "dpi": args.dpi,
            "smoothing_sigma_pixels": args.smooth_pixels,
            "gas_smoothing_sigma_pixels": args.gas_smooth_pixels,
            "hi_density_cut_applied": False,
            "mass_to_light_v": MASS_TO_LIGHT,
            "bright_limit_mag_arcsec2": args.bright_limit,
            "faint_limit_mag_arcsec2": args.faint_limit,
        },
        "outputs": {
            "paper_panel_png": str(panel_png),
            "paper_panel_pdf": str(panel_pdf),
            "diagnostic_png": str(diagnostic_png),
        },
    }
    metadata_json.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
