"""Build a paper-ready multi-snapshot evolution figure for the Fornax run."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import make_snapshot_mock as mock


PANEL_COUNT = 6
DEFAULT_LOOKBACK_STEP_GYR = 0.5


def load_elinfo(path: Path | None) -> dict[int, dict[str, str]]:
    """Index elinfo rows by numsp while ignoring the two comment lines."""
    if path is None:
        return {}
    with path.open(encoding="utf-8") as handle:
        rows = csv.DictReader(line for line in handle if not line.startswith("#"))
        return {int(float(row["numsp"])): row for row in rows}


def finite_float(row: dict[str, str] | None, key: str) -> float | None:
    if row is None or key not in row:
        return None
    try:
        value = float(row[key])
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def present_reference_from_distance(
    elinfo: dict[int, dict[str, str]], target_distance_kpc: float
) -> dict[str, float | int]:
    """Find the elinfo time whose heliocentric distance best matches the target."""
    candidates: list[tuple[float, float, float, int]] = []
    for snapshot_number, row in elinfo.items():
        distance_kpc = finite_float(row, "distance")
        age_gyr = finite_float(row, "age")
        if distance_kpc is None or age_gyr is None:
            continue
        candidates.append(
            (
                abs(distance_kpc - target_distance_kpc),
                distance_kpc,
                age_gyr,
                snapshot_number,
            )
        )
    if not candidates:
        raise ValueError("elinfo contains no finite heliocentric distance/age rows")

    _, distance_kpc, age_gyr, snapshot_number = min(candidates, key=lambda item: item[0])
    return {
        "snapshot": snapshot_number,
        "distance_kpc": distance_kpc,
        "time_gyr": age_gyr,
    }


def select_default_snapshots(
    elinfo: dict[int, dict[str, str]],
    present_time_gyr: float,
    step_gyr: float,
    panel_count: int = PANEL_COUNT,
) -> list[int]:
    """Select a recent regular lookback grid, using the initial frame for a remainder."""
    if step_gyr <= 0.0:
        raise ValueError("--lookback-step-gyr must be positive")

    timed_rows = []
    for snapshot_number, row in elinfo.items():
        age_gyr = finite_float(row, "age")
        if age_gyr is not None and age_gyr <= present_time_gyr + 1.0e-8:
            timed_rows.append((snapshot_number, age_gyr))
    if len(timed_rows) < panel_count:
        raise ValueError(f"elinfo contains fewer than {panel_count} usable snapshots")

    earliest_time_gyr = min(age for _, age in timed_rows)
    available_lookback_gyr = present_time_gyr - earliest_time_gyr
    full_span_gyr = step_gyr * (panel_count - 1)

    if available_lookback_gyr >= full_span_gyr - 1.0e-8:
        target_lookbacks = [
            step_gyr * index for index in range(panel_count - 1, -1, -1)
        ]
    else:
        regular = [
            step_gyr * index
            for index in range(int(np.floor(available_lookback_gyr / step_gyr)) + 1)
        ]
        if not np.isclose(regular[-1], available_lookback_gyr, atol=1.0e-8):
            regular.append(available_lookback_gyr)
        target_lookbacks = sorted(regular, reverse=True)
        if len(target_lookbacks) != panel_count:
            raise ValueError(
                f"A {step_gyr:.3g} Gyr grid spans only {len(target_lookbacks)} panels; "
                "pass --snapshots explicitly for this shorter run."
            )

    selected = []
    for lookback_gyr in target_lookbacks:
        target_time_gyr = present_time_gyr - lookback_gyr
        snapshot_number, _ = min(
            timed_rows, key=lambda item: abs(item[1] - target_time_gyr)
        )
        selected.append(snapshot_number)
    if len(set(selected)) != panel_count:
        raise ValueError("The lookback grid mapped to duplicate snapshots")
    return selected


def locate_snapshot(model_dir: Path, number: int) -> Path:
    candidates = (
        model_dir / "output" / f"snapshot_{number:03d}.hdf5",
        model_dir / f"snapshot_{number:03d}.hdf5",
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"snapshot_{number:03d}.hdf5 not found below {model_dir}")


def process_snapshot(
    snapshot_path: Path,
    snapshot_number: int,
    elinfo_row: dict[str, str] | None,
    args: argparse.Namespace,
) -> dict[str, object]:
    particles, box_size, header_time_gyr = mock.read_particles(snapshot_path)
    mw_center = mock.shrinking_sphere_mw_center(particles)
    positions, masses, particle_types, velocities = mock.concatenate_stars(
        particles, mw_center
    )
    dwarf_mask, dwarf_peak = mock.find_dwarf_mask(positions)
    dwarf_center = np.average(positions[dwarf_mask], axis=0, weights=masses[dwarf_mask])
    gas = mock.read_dwarf_gas(
        snapshot_path=snapshot_path,
        box_size=box_size,
        mw_center=mw_center,
        dwarf_center=dwarf_center,
    )

    old_dwarf = dwarf_mask & (particle_types != 4)
    selected_positions = positions[old_dwarf]
    selected_masses = masses[old_dwarf]
    selected_velocities = velocities[old_dwarf]
    if selected_positions.shape[0] == 0:
        raise RuntimeError(f"No old dwarf stars in {snapshot_path.name}")

    x_deg, y_deg, particle_distance, _, sky_center = mock.project_to_sky(
        selected_positions
    )
    transformed_distance_kpc = float(
        np.average(particle_distance, weights=selected_masses)
    )
    angular_rescale = transformed_distance_kpc / args.distance_kpc
    x_deg *= angular_rescale
    y_deg *= angular_rescale

    surface_brightness, _, _, _ = mock.make_light_map(
        x_deg=x_deg,
        y_deg=y_deg,
        masses=selected_masses,
        distance_kpc=args.distance_kpc,
        half_width_deg=args.field_half_deg,
        npix=args.npix,
        smoothing_sigma_pixels=args.smooth_pixels,
    )
    stellar_display = mock.magnitude_to_display(
        surface_brightness, args.bright_limit, args.faint_limit
    )

    if gas["position"].shape[0] > 0:
        gas_x_deg, gas_y_deg, _, _, _ = mock.project_to_sky(
            gas["position"], center=sky_center
        )
        gas_x_deg *= angular_rescale
        gas_y_deg *= angular_rescale
        hi_surface_density, cold_gas = mock.make_hi_map(
            x_deg=gas_x_deg,
            y_deg=gas_y_deg,
            masses=gas["mass"],
            neutral_fraction=gas["neutral_fraction"],
            temperature=gas["temperature"],
            smoothing_length_kpc=gas["smoothing_length"],
            distance_kpc=args.distance_kpc,
            half_width_deg=args.field_half_deg,
            npix=args.npix,
            smoothing_sigma_pixels=args.gas_smooth_pixels,
        )
    else:
        gas_x_deg = np.empty(0)
        gas_y_deg = np.empty(0)
        cold_gas = np.zeros(0, dtype=bool)
        hi_surface_density = np.zeros_like(surface_brightness)

    hi_display = mock.hi_to_display(hi_surface_density)
    hi_aperture_half_deg = mock.HI_MASS_FIELD_BUFFER * args.field_half_deg
    hi_aperture = (
        cold_gas
        & (np.abs(gas_x_deg) <= hi_aperture_half_deg)
        & (np.abs(gas_y_deg) <= hi_aperture_half_deg)
    )
    hi_mass_msun = float(
        np.sum(gas["mass"][hi_aperture] * gas["neutral_fraction"][hi_aperture])
    )

    pmra = finite_float(elinfo_row, "pmra")
    pmdec = finite_float(elinfo_row, "pmdec")
    pm_source = "elinfo"
    if pmra is None or pmdec is None:
        pmra, pmdec = mock.mean_proper_motion(selected_positions, selected_velocities)
        pm_source = "stellar particles"

    time_gyr = finite_float(elinfo_row, "age")
    if time_gyr is None:
        _, time_gyr = mock.snapshot_time_gyr(snapshot_path, header_time_gyr)

    return {
        "snapshot": snapshot_number,
        "snapshot_path": str(snapshot_path),
        "time_gyr": float(time_gyr),
        "header_time_gyr": float(header_time_gyr),
        "stellar_display": stellar_display,
        "hi_display": hi_display,
        "hi_surface_density": hi_surface_density,
        "stellar_mass_msun": float(selected_masses.sum()),
        "hi_mass_msun": hi_mass_msun,
        "pmra_masyr": float(pmra),
        "pmdec_masyr": float(pmdec),
        "pm_source": pm_source,
        "dwarf_center_galactocentric_kpc": dwarf_center.tolist(),
        "dwarf_radius_galactocentric_kpc": float(np.linalg.norm(dwarf_center)),
        "elinfo_distance_gal_kpc": finite_float(elinfo_row, "distance_gal"),
        "elinfo_heliocentric_distance_kpc": finite_float(elinfo_row, "distance"),
        "elinfo_cold_gas_mass_msun": finite_float(elinfo_row, "coldgas_mass"),
        "transformed_particle_distance_kpc": transformed_distance_kpc,
        "adopted_heliocentric_distance_kpc": float(args.distance_kpc),
        "dwarf_density_peak_galactocentric_kpc": dwarf_peak.tolist(),
        "old_star_particle_count": int(np.count_nonzero(old_dwarf)),
        "gas_particle_count_within_20_kpc": int(gas["mass"].size),
    }


def draw_panel(
    axis: plt.Axes,
    panel: dict[str, object],
    args: argparse.Namespace,
    panel_index: int,
) -> None:
    extent = [
        -args.field_half_deg,
        args.field_half_deg,
        -args.field_half_deg,
        args.field_half_deg,
    ]
    hi_display = np.asarray(panel["hi_display"])
    stellar_display = np.asarray(panel["stellar_display"])
    hi_surface_density = np.asarray(panel["hi_surface_density"])

    hi_alpha = 0.46 * np.power(np.clip(hi_display.T, 0.0, 1.0), 0.72)
    axis.imshow(
        hi_display.T,
        origin="lower",
        extent=extent,
        interpolation="bilinear",
        cmap=mock.gas_colormap(),
        vmin=0.0,
        vmax=1.0,
        alpha=hi_alpha,
        rasterized=True,
    )
    stellar_rgba = mock.stellar_colormap()(np.clip(stellar_display.T, 0.0, 1.0))
    stellar_rgba[..., 3] = np.power(np.clip(stellar_display.T, 0.0, 1.0), 0.62)
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
            linewidths=[0.38, 0.48, 0.58, 0.70],
            alpha=0.92,
        )
    mock.style_axis(axis, args.field_half_deg)
    mock.draw_proper_motion_arrow(
        axis,
        float(panel["pmra_masyr"]),
        float(panel["pmdec_masyr"]),
        arrow_length_deg=0.62,
    )

    row, column = divmod(panel_index, 3)
    if row == 0:
        axis.set_xlabel("")
        axis.tick_params(labelbottom=False)
    if column != 0:
        axis.set_ylabel("")
        axis.tick_params(labelleft=False)

    axis.text(
        0.035,
        0.955,
        rf"$M_\star={mock.mass_to_tex(float(panel['stellar_mass_msun']))}\,M_\odot$",
        transform=axis.transAxes,
        color="#f0dfcf",
        fontsize=5.35,
        ha="left",
        va="top",
    )
    axis.text(
        0.035,
        0.900,
        rf"$M_{{\mathrm{{H\,I}}}}={mock.mass_to_tex(float(panel['hi_mass_msun']))}\,M_\odot$",
        transform=axis.transAxes,
        color="#94f3f3",
        fontsize=5.15,
        ha="left",
        va="top",
    )
    panel_distance = panel.get("elinfo_heliocentric_distance_kpc")
    if panel_distance is None:
        panel_distance = panel["adopted_heliocentric_distance_kpc"]
    axis.text(
        0.035,
        0.035,
        rf"$t_{{\mathrm{{lookback}}}}={float(panel['lookback_time_gyr']):.2f}\,\mathrm{{Gyr}}$"
        "\n"
        rf"$D={float(panel_distance):.1f}\,\mathrm{{kpc}}$",
        transform=axis.transAxes,
        color="#d7dadd",
        fontsize=5.25,
        ha="left",
        va="bottom",
        linespacing=1.22,
    )
    axis.text(
        0.965,
        0.955,
        f"({chr(ord('a') + panel_index)})",
        transform=axis.transAxes,
        color="#e8e8e8",
        fontsize=6.2,
        ha="right",
        va="top",
    )

    if panel_index == 0:
        axis.legend(
            handles=[Line2D([0], [0], color="#a7ffff", lw=0.9, label="HI cloud")],
            loc="lower right",
            frameon=False,
            fontsize=5.35,
            handlelength=1.25,
            handletextpad=0.38,
            labelcolor="#e8eeee",
        )


def public_metadata(panel: dict[str, object]) -> dict[str, object]:
    return {
        key: value
        for key, value in panel.items()
        if key not in {"stellar_display", "hi_display", "hi_surface_density"}
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model_dir", type=Path)
    parser.add_argument("--elinfo", type=Path)
    parser.add_argument(
        "--snapshots",
        nargs="+",
        type=int,
        help="Explicit snapshot numbers; defaults to an elinfo-derived lookback grid",
    )
    parser.add_argument("--outdir", type=Path, default=Path("output"))
    parser.add_argument("--distance-kpc", type=float, default=mock.ADOPTED_DISTANCE_KPC)
    parser.add_argument(
        "--present-time-gyr",
        type=float,
        help="Optional override for lookback-time zero; defaults to the elinfo distance match",
    )
    parser.add_argument(
        "--lookback-step-gyr",
        type=float,
        default=DEFAULT_LOOKBACK_STEP_GYR,
        help="Spacing of the automatically selected lookback-time grid",
    )
    parser.add_argument("--field-half-deg", type=float, default=mock.FIELD_HALF_WIDTH_DEG)
    parser.add_argument("--npix", type=int, default=520)
    parser.add_argument("--dpi", type=int, default=450)
    parser.add_argument("--smooth-pixels", type=float, default=1.7)
    parser.add_argument("--gas-smooth-pixels", type=float, default=3.8)
    parser.add_argument("--bright-limit", type=float, default=22.0)
    parser.add_argument("--faint-limit", type=float, default=34.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    model_dir = args.model_dir.resolve()
    output_dir = args.outdir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    elinfo_path = args.elinfo
    if elinfo_path is None:
        candidates = sorted(model_dir.glob("elinfo*.csv"))
        elinfo_path = candidates[0] if candidates else None
    elinfo = load_elinfo(elinfo_path.resolve() if elinfo_path else None)

    present_reference = None
    if args.present_time_gyr is None:
        present_reference = present_reference_from_distance(elinfo, args.distance_kpc)
        present_time_gyr = float(present_reference["time_gyr"])
    else:
        present_time_gyr = float(args.present_time_gyr)

    snapshot_numbers = (
        list(args.snapshots)
        if args.snapshots is not None
        else select_default_snapshots(
            elinfo,
            present_time_gyr,
            args.lookback_step_gyr,
        )
    )
    if len(snapshot_numbers) != PANEL_COUNT:
        raise ValueError(
            f"The paper layout currently requires exactly {PANEL_COUNT} snapshots."
        )

    panels = []
    for number in snapshot_numbers:
        snapshot_path = locate_snapshot(model_dir, number)
        print(f"Processing {snapshot_path.name}", flush=True)
        panels.append(process_snapshot(snapshot_path, number, elinfo.get(number), args))

    for panel in panels:
        panel["lookback_time_gyr"] = present_time_gyr - float(panel["time_gyr"])

    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 7.0,
            "savefig.facecolor": "black",
            "figure.facecolor": "black",
        }
    )
    fig, axes = plt.subplots(2, 3, figsize=(7.10, 4.78), facecolor="black")
    for index, (axis, panel) in enumerate(zip(axes.flat, panels)):
        draw_panel(axis, panel, args, index)
    fig.subplots_adjust(left=0.073, right=0.995, bottom=0.090, top=0.995, wspace=0.035, hspace=0.035)

    stem = f"{model_dir.name.lower()}_evolution_" + "_".join(
        f"{number:03d}" for number in snapshot_numbers
    )
    png_path = output_dir / f"{stem}.png"
    metadata_path = output_dir / f"{stem}_metadata.json"
    fig.savefig(png_path, dpi=args.dpi, facecolor="black", bbox_inches="tight", pad_inches=0.01)
    plt.close(fig)

    metadata = {
        "model_dir": str(model_dir),
        "elinfo": str(elinfo_path.resolve()) if elinfo_path else None,
        "snapshots": snapshot_numbers,
        "snapshot_cadence_gyr": mock.SNAPSHOT_CADENCE_GYR,
        "lookback_step_gyr": args.lookback_step_gyr,
        "present_time_gyr": present_time_gyr,
        "lookback_time_reference": {
            "method": "closest elinfo heliocentric distance"
            if present_reference is not None
            else "explicit present time",
            "target_distance_kpc": args.distance_kpc,
            "matched_snapshot": present_reference["snapshot"]
            if present_reference is not None
            else None,
            "matched_distance_kpc": present_reference["distance_kpc"]
            if present_reference is not None
            else None,
        },
        "adopted_heliocentric_distance_kpc": args.distance_kpc,
        "distance_reference": f"Li et al. (2021), dm=20.72 mag ({args.distance_kpc:g} kpc)",
        "field_half_width_deg": args.field_half_deg,
        "ra_axis_increases_to_left": True,
        "rendered_gas": "NH-weighted cold HI, no density cut",
        "gas_selection_radius_kpc": mock.DWARF_GAS_RADIUS_KPC,
        "panels": [public_metadata(panel) for panel in panels],
        "outputs": {"png": str(png_path)},
    }
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == "__main__":
    main()
