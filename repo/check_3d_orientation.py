import argparse
import json
import os
import re
import warnings
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig-dsph")
warnings.filterwarnings(
    "ignore",
    message=r"mergecube\.par not found under .*; falling back to legacy cube_path=.*",
    category=RuntimeWarning,
)

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objs as go

from basefunc import Analysis
from snapshot_context import prepare_snapshot_context
from variable import cd_threshold, fornax_core_radius


CODE_VERSION = "v0.3.1"
DEFAULT_INITIAL_PLACEMENT_ANGLES_DEG = (83.0, 45.0, 0.0)

SCRIPT_DIR = Path(__file__).resolve().parent
WORKSPACE_ROOT = SCRIPT_DIR.parent


def parse_snapshot_path(snapshot_path):
    path = Path(snapshot_path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Snapshot file not found: {path}")
    match = re.fullmatch(r"snapshot_(\d+)\.hdf5", path.name)
    if match is None:
        raise ValueError(f"Expected snapshot filename like snapshot_211.hdf5, got: {path.name}")
    return path, path.parent, int(match.group(1))


def default_output_dir():
    sandbox_root = WORKSPACE_ROOT / "sandbox_runs"
    if sandbox_root.exists():
        return (sandbox_root / "3dcheck").resolve()
    return (Path.cwd().resolve() / "3dcheck_output").resolve()


def make_output_paths(snapshot_path, output_dir=None):
    snapshot_path, _, snapshot_num = parse_snapshot_path(snapshot_path)
    output_dir = default_output_dir() if output_dir is None else Path(output_dir).expanduser().resolve()
    stem = f"3dcheck_{snapshot_path.parent.name}_{snapshot_num:03d}_{CODE_VERSION}"
    return {
        "output_dir": output_dir,
        "html": output_dir / f"{stem}.html",
        "gas_png": output_dir / f"{stem}_gas_contour.png",
        "yz_png": output_dir / f"{stem}_yz_orientation.png",
        "summary_json": output_dir / f"{stem}_summary.json",
    }


def load_snapshot_context_from_path(
    snapshot_path,
    core_radius=fornax_core_radius,
    r_exclude=5.0,
    dwarf_radius_factor=3.0,
    k_density=16,
    dwarf_gas_radius=30.0,
    gas_temperature_split=20000.0,
):
    snapshot_path, folder_path, snapshot_num = parse_snapshot_path(snapshot_path)
    context = prepare_snapshot_context(
        folder_path=f"{folder_path}{os.sep}",
        snapshot_num=snapshot_num,
        core_radius=core_radius,
        r_exclude=r_exclude,
        dwarf_radius_factor=dwarf_radius_factor,
        k_density=k_density,
        dwarf_gas_radius=dwarf_gas_radius,
        gas_temperature_split=gas_temperature_split,
    )
    context["snapshot_path"] = snapshot_path
    context["snapshot_num"] = snapshot_num
    context["model_dir"] = folder_path
    context["code_version"] = CODE_VERSION
    return context


def _sample_indices(size, max_points=None, seed=0):
    if max_points is None or max_points <= 0 or size <= max_points:
        return np.arange(size)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(size, size=max_points, replace=False))


def heliocentric_stellar_center(context):
    star_coords = context["star_coords"]
    star_mask = np.asarray(context["total_dw_star_mask"], dtype=bool)
    star_mass = context["df"].loc[star_mask, "m"].to_numpy(dtype=float)
    return np.asarray(
        Analysis.find_center_3d(
            star_coords["xh"],
            star_coords["yh"],
            star_coords["zh"],
            mass=star_mass,
        ),
        dtype=float,
    )


def galactocentric_stellar_center(context):
    return np.asarray(
        [context["dw_xc"], context["dw_yc"], context["dw_zc"]],
        dtype=float,
    )


def _safe_unit_vector(vector):
    vector = np.asarray(vector, dtype=float)
    norm = np.linalg.norm(vector)
    if not np.isfinite(norm) or norm <= 0:
        raise ValueError(f"Cannot normalize vector with norm={norm}: {vector}")
    return vector / norm, norm


def _rotation_matrix_from_vectors(source, target):
    source_unit, _ = _safe_unit_vector(source)
    target_unit, _ = _safe_unit_vector(target)
    cross = np.cross(source_unit, target_unit)
    dot = float(np.clip(np.dot(source_unit, target_unit), -1.0, 1.0))

    if np.isclose(dot, 1.0):
        return np.eye(3)

    if np.isclose(dot, -1.0):
        basis = np.eye(3)
        axis = basis[np.argmin(np.abs(source_unit))]
        axis = axis - np.dot(axis, source_unit) * source_unit
        axis, _ = _safe_unit_vector(axis)
        return 2.0 * np.outer(axis, axis) - np.eye(3)

    skew = np.array(
        [
            [0.0, -cross[2], cross[1]],
            [cross[2], 0.0, -cross[0]],
            [-cross[1], cross[0], 0.0],
        ]
    )
    return np.eye(3) + skew + skew @ skew * ((1.0 - dot) / np.dot(cross, cross))


def rotation_x(angle):
    return np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, np.cos(angle), -np.sin(angle)],
            [0.0, np.sin(angle), np.cos(angle)],
        ]
    )


def rotation_y(angle):
    return np.array(
        [
            [np.cos(angle), 0.0, np.sin(angle)],
            [0.0, 1.0, 0.0],
            [-np.sin(angle), 0.0, np.cos(angle)],
        ]
    )


def rotation_z(angle):
    return np.array(
        [
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )


def rotation_matrix_xyz(angles):
    ax, ay, az = np.asarray(angles, dtype=float)
    return rotation_z(az) @ rotation_y(ay) @ rotation_x(ax)


def euler_xyz_from_rotation_matrix(matrix):
    """Return angles for R = Rz(az) @ Ry(ay) @ Rx(ax)."""
    matrix = np.asarray(matrix, dtype=float)
    sy = -matrix[2, 0]
    ay = np.arcsin(np.clip(sy, -1.0, 1.0))
    cy = np.cos(ay)

    if abs(cy) > 1e-12:
        ax = np.arctan2(matrix[2, 1], matrix[2, 2])
        az = np.arctan2(matrix[1, 0], matrix[0, 0])
    else:
        ax = 0.0
        az = np.arctan2(-matrix[0, 1], matrix[1, 1])

    return np.asarray([ax, ay, az], dtype=float)


def wrap_degrees_180(angles_deg):
    return (np.asarray(angles_deg, dtype=float) + 180.0) % 360.0 - 180.0


def rotate_coords(x, y, z, angles, center=None):
    coords = np.vstack([np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float)])
    if center is not None:
        center = np.asarray(center, dtype=float).reshape(3, 1)
        coords = coords - center
    rotated = rotation_matrix_xyz(angles) @ coords
    if center is not None:
        rotated = rotated + center
    return rotated[0], rotated[1], rotated[2]


def apply_ic_rotation_correction(x, y, z, correction, center=None):
    matrix = np.asarray(correction["correction_matrix"], dtype=float)
    coords = np.vstack([np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float)])
    if center is not None:
        center = np.asarray(center, dtype=float).reshape(3, 1)
        coords = coords - center
    rotated = matrix @ coords
    if center is not None:
        rotated = rotated + center
    return rotated[0], rotated[1], rotated[2]


def infer_ic_rotation_correction(center_xyz):
    """Infer the legacy correction angles that align a center vector with +Z.

    The legacy notebook returned (0, ay, az), then verified the correction with
    Ry(-ay) @ Rz(-az) @ n.  This function keeps both representations explicit.
    """
    center = np.asarray(center_xyz, dtype=float)
    unit_vector, radius = _safe_unit_vector(center)

    alpha_z = np.arctan2(unit_vector[1], unit_vector[0])
    r_xy = np.hypot(unit_vector[0], unit_vector[1])
    alpha_y = np.arctan2(r_xy, unit_vector[2])
    alpha_x = 0.0

    legacy_angles = np.asarray([alpha_x, alpha_y, alpha_z], dtype=float)
    correction_matrix = rotation_y(-alpha_y) @ rotation_z(-alpha_z)
    corrected_unit_vector = correction_matrix @ unit_vector

    return {
        "center_xyz": center,
        "radius": radius,
        "unit_vector": unit_vector,
        "legacy_angles_rad": legacy_angles,
        "legacy_angles_deg": np.degrees(legacy_angles),
        "correction_sequence_rad": {
            "first_Rz": -alpha_z,
            "then_Ry": -alpha_y,
        },
        "correction_sequence_deg": {
            "first_Rz": float(np.degrees(-alpha_z)),
            "then_Ry": float(np.degrees(-alpha_y)),
        },
        "correction_matrix": correction_matrix,
        "corrected_unit_vector": corrected_unit_vector,
    }


def orientation_report(context):
    return {
        "heliocentric": infer_ic_rotation_correction(heliocentric_stellar_center(context)),
        "galactocentric": infer_ic_rotation_correction(galactocentric_stellar_center(context)),
    }


def _relative_positions(context, mask):
    df = context["df"]
    mask = np.asarray(mask, dtype=bool)
    center = galactocentric_stellar_center(context)
    return np.column_stack(
        [
            df.loc[mask, "x"].to_numpy(dtype=float) - center[0],
            df.loc[mask, "y"].to_numpy(dtype=float) - center[1],
            df.loc[mask, "z"].to_numpy(dtype=float) - center[2],
        ]
    )


def stellar_major_axis_3d(context, radius=20.0):
    star_mask = np.asarray(context["total_dw_star_mask"], dtype=bool)
    coords = _relative_positions(context, star_mask)
    mass = context["df"].loc[star_mask, "m"].to_numpy(dtype=float)

    r = np.linalg.norm(coords, axis=1)
    finite = np.isfinite(coords).all(axis=1) & np.isfinite(mass) & (mass > 0)
    if radius is not None:
        finite &= r <= radius
    if finite.sum() < 3:
        raise ValueError("Need at least three finite stellar particles to estimate a 3D major axis")

    coords = coords[finite]
    mass = mass[finite]
    weighted_mean = np.average(coords, axis=0, weights=mass)
    coords = coords - weighted_mean
    tensor = (coords * mass[:, None]).T @ coords / mass.sum()

    eigenvalues, eigenvectors = np.linalg.eigh(tensor)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    major_axis = eigenvectors[:, 0]
    major_axis, _ = _safe_unit_vector(major_axis)

    return {
        "major_axis": major_axis,
        "eigenvalues": eigenvalues,
        "eigenvectors": eigenvectors,
        "radius": radius,
        "particle_count": int(finite.sum()),
    }


def line_of_sight_unit(context):
    los, _ = _safe_unit_vector(heliocentric_stellar_center(context))
    return los


def long_axis_orientation_correction(context, radius=20.0):
    axis_info = stellar_major_axis_3d(context, radius=radius)
    major_axis = axis_info["major_axis"]
    los = line_of_sight_unit(context)

    signed_los = los if np.dot(major_axis, los) >= 0 else -los
    correction_matrix = _rotation_matrix_from_vectors(major_axis, signed_los)
    corrected_major_axis = correction_matrix @ major_axis
    angle_rad = np.arccos(np.clip(abs(np.dot(major_axis, los)), -1.0, 1.0))

    return {
        "axis_radius": radius,
        "major_axis": major_axis,
        "line_of_sight_sun_to_dwarf": los,
        "signed_target_los": signed_los,
        "angle_rad": float(angle_rad),
        "angle_deg": float(np.degrees(angle_rad)),
        "correction_matrix": correction_matrix,
        "corrected_major_axis": corrected_major_axis,
        "axis_particle_count": axis_info["particle_count"],
        "axis_eigenvalues": axis_info["eigenvalues"],
    }


def placement_angle_update(correction, initial_angles_deg=DEFAULT_INITIAL_PLACEMENT_ANGLES_DEG):
    initial_angles_deg = np.asarray(initial_angles_deg, dtype=float)
    initial_matrix = rotation_matrix_xyz(np.radians(initial_angles_deg))
    correction_matrix = np.asarray(correction["correction_matrix"], dtype=float)

    # The correction is an orientation-only rotation in the simulation frame, so it
    # left-multiplies the existing placement matrix under this explicit convention.
    corrected_matrix = correction_matrix @ initial_matrix
    corrected_angles_rad = euler_xyz_from_rotation_matrix(corrected_matrix)
    corrected_angles_deg = np.degrees(corrected_angles_rad)

    legacy_initial_matrix = rotation_matrix_xyz(np.radians(-initial_angles_deg))
    legacy_corrected_matrix = correction_matrix @ legacy_initial_matrix
    legacy_corrected_active_angles_deg = np.degrees(euler_xyz_from_rotation_matrix(legacy_corrected_matrix))
    legacy_corrected_angles_deg = -legacy_corrected_active_angles_deg

    return {
        "convention": "active extrinsic XYZ: R = Rz(az) @ Ry(ay) @ Rx(ax), angles_deg=(ax, ay, az)",
        "composition": "corrected_matrix = long_axis_correction_matrix @ initial_matrix",
        "initial_angles_deg": initial_angles_deg,
        "corrected_angles_deg": corrected_angles_deg,
        "corrected_angles_wrapped_deg": wrap_degrees_180(corrected_angles_deg),
        "initial_matrix": initial_matrix,
        "corrected_matrix": corrected_matrix,
        "legacy_rot3d_coordinate_system": {
            "convention": (
                "coordinate-system rotation like mysubs.rot_3d/zim.rot_3d called in X,Y,Z order; "
                "effective active matrix is Rz(-az) @ Ry(-ay) @ Rx(-ax)"
            ),
            "corrected_angles_deg": legacy_corrected_angles_deg,
            "corrected_angles_wrapped_deg": wrap_degrees_180(legacy_corrected_angles_deg),
            "corrected_matrix": legacy_corrected_matrix,
        },
    }


def _sample_coords(coords, max_points, seed):
    indices = _sample_indices(coords.shape[0], max_points=max_points, seed=seed)
    return coords[indices]


def _plot_yz_components(ax, components, radius, title, major_axis, los_axis):
    colors = {
        "stars": "gray",
        "cold_gas": "aquamarine",
        "hot_gas": "orange",
    }
    sizes = {
        "stars": 2,
        "cold_gas": 4,
        "hot_gas": 3,
    }
    labels = {
        "stars": "Stars",
        "cold_gas": "HI Gas",
        "hot_gas": "HII Gas",
    }

    for name, coords in components.items():
        if coords.size == 0:
            continue
        ax.scatter(coords[:, 1], coords[:, 2], s=sizes[name], c=colors[name], lw=0, alpha=0.75, label=labels[name])

    axis_len = 0.75 * radius
    ax.plot(
        [-major_axis[1] * axis_len, major_axis[1] * axis_len],
        [-major_axis[2] * axis_len, major_axis[2] * axis_len],
        color="dodgerblue",
        lw=2.0,
        label="3D stellar major axis",
    )
    ax.plot(
        [-los_axis[1] * axis_len, los_axis[1] * axis_len],
        [-los_axis[2] * axis_len, los_axis[2] * axis_len],
        color="black",
        lw=2.0,
        ls="--",
        label="Sun LOS axis",
    )
    ax.set_title(title)
    ax.set_xlabel("Y - Yc (kpc)")
    ax.set_ylabel("Z - Zc (kpc)")
    ax.set_xlim(-radius, radius)
    ax.set_ylim(-radius, radius)
    ax.set_aspect("equal", adjustable="box")
    ax.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax.grid(True, alpha=0.25)


def render_yz_orientation_figure(
    context,
    correction=None,
    radius=20.0,
    max_star_points=200000,
    max_gas_points=100000,
    seed=0,
):
    if correction is None:
        correction = long_axis_orientation_correction(context, radius=radius)
    matrix = np.asarray(correction["correction_matrix"], dtype=float)

    component_masks = {
        "stars": np.asarray(context["total_dw_star_mask"], dtype=bool),
        "cold_gas": np.asarray(context["dw_cold_gas_mask"], dtype=bool),
        "hot_gas": np.asarray(context["dw_hot_gas_mask"], dtype=bool),
    }
    max_points = {
        "stars": max_star_points,
        "cold_gas": max_gas_points,
        "hot_gas": max_gas_points,
    }

    current_components = {}
    rotated_components = {}
    for idx, (name, mask) in enumerate(component_masks.items()):
        coords = _relative_positions(context, mask)
        if coords.size == 0:
            current_components[name] = coords
            rotated_components[name] = coords
            continue
        r = np.linalg.norm(coords, axis=1)
        in_radius = np.isfinite(coords).all(axis=1) & (r <= radius)
        coords = _sample_coords(coords[in_radius], max_points[name], seed + idx)
        current_components[name] = coords
        rotated_components[name] = (matrix @ coords.T).T

    major_axis = np.asarray(correction["major_axis"], dtype=float)
    corrected_major_axis = np.asarray(correction["corrected_major_axis"], dtype=float)
    los_axis = np.asarray(correction["line_of_sight_sun_to_dwarf"], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.6), dpi=140, sharex=True, sharey=True)
    _plot_yz_components(
        axes[0],
        current_components,
        radius=radius,
        title="Current orientation",
        major_axis=major_axis,
        los_axis=los_axis,
    )
    _plot_yz_components(
        axes[1],
        rotated_components,
        radius=radius,
        title=f"After long-axis correction ({correction['angle_deg']:.2f} deg)",
        major_axis=corrected_major_axis,
        los_axis=los_axis,
    )
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, frameon=False)
    fig.suptitle(
        f"{context['model_dir'].name} snapshot {context['snapshot_num']:03d}: centered Y-Z view",
        y=0.98,
    )
    fig.tight_layout(rect=[0, 0.08, 1, 0.94])
    return fig, correction


def save_yz_orientation_plot(context, output_path, correction=None, radius=20.0, max_star_points=200000, max_gas_points=100000, seed=0):
    fig, correction = render_yz_orientation_figure(
        context,
        correction=correction,
        radius=radius,
        max_star_points=max_star_points,
        max_gas_points=max_gas_points,
        seed=seed,
    )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)
    return output_path, correction


def build_3d_figure(
    context,
    max_star_points=200000,
    max_cold_gas_points=100000,
    cold_gas_radius=5.0,
    seed=0,
    velocity_arrow_scale=1.0e4,
    sun_arrow_length=5.0,
):
    star_coords = context["star_coords"]
    cold_gas_coords = context["cold_gas_coords"]
    center = heliocentric_stellar_center(context)

    star_pos = np.column_stack([star_coords["xh"], star_coords["yh"], star_coords["zh"]])
    star_idx = _sample_indices(star_pos.shape[0], max_points=max_star_points, seed=seed)

    sun_unit = -center / np.linalg.norm(center)
    star_depth = -(star_pos[star_idx] - center) @ sun_unit

    fig = go.Figure()
    fig.add_trace(
        go.Scatter3d(
            x=star_pos[star_idx, 0],
            y=star_pos[star_idx, 1],
            z=star_pos[star_idx, 2],
            mode="markers",
            marker={
                "size": 2,
                "color": star_depth,
                "colorscale": "Bluered",
                "opacity": 0.85,
                "colorbar": {"title": "LOS depth", "x": 0.05},
            },
            name="Dwarf stars",
        )
    )

    mean_velocity = np.array(
        [
            np.mean(star_coords["vxh"]),
            np.mean(star_coords["vyh"]),
            np.mean(star_coords["vzh"]),
        ],
        dtype=float,
    )
    legacy_km_s_to_kpc_display = 1.0227e-6
    velocity_delta = mean_velocity * legacy_km_s_to_kpc_display * velocity_arrow_scale
    velocity_end = center + velocity_delta
    fig.add_trace(
        go.Scatter3d(
            x=[center[0], velocity_end[0]],
            y=[center[1], velocity_end[1]],
            z=[center[2], velocity_end[2]],
            mode="lines+markers",
            line={"color": "palegreen", "width": 8},
            marker={"size": 5, "color": "palegreen"},
            name="Mean heliocentric velocity",
        )
    )

    sun_end = center + sun_unit * sun_arrow_length
    fig.add_trace(
        go.Scatter3d(
            x=[center[0], sun_end[0]],
            y=[center[1], sun_end[1]],
            z=[center[2], sun_end[2]],
            mode="lines+markers",
            line={"color": "black", "width": 8, "dash": "dash"},
            marker={"size": 4, "color": "black"},
            name="Direction to Sun",
        )
    )

    cold_pos = np.column_stack([cold_gas_coords["xh"], cold_gas_coords["yh"], cold_gas_coords["zh"]])
    cold_radius = np.linalg.norm(cold_pos - center, axis=1) if cold_pos.size else np.array([])
    cold_near = np.flatnonzero(cold_radius < cold_gas_radius)
    cold_idx_local = _sample_indices(cold_near.size, max_points=max_cold_gas_points, seed=seed + 1)
    cold_idx = cold_near[cold_idx_local]
    if cold_idx.size > 0:
        fig.add_trace(
            go.Scatter3d(
                x=cold_pos[cold_idx, 0],
                y=cold_pos[cold_idx, 1],
                z=cold_pos[cold_idx, 2],
                mode="markers",
                marker={"size": 2, "color": "cyan", "opacity": 0.9},
                name=f"Cold gas < {cold_gas_radius:g} kpc",
            )
        )

    fig.update_layout(
        title=(
            f"{context['model_dir'].name} snapshot {context['snapshot_num']:03d} "
            f"(T={context['tsnap']:.3f} Gyr)"
        ),
        scene={
            "xaxis_title": "Xh (kpc)",
            "yaxis_title": "Yh (kpc)",
            "zaxis_title": "Zh (kpc)",
            "aspectmode": "data",
        },
        margin={"l": 0, "r": 0, "b": 0, "t": 35},
    )

    info = {
        "star_particles_total": int(star_pos.shape[0]),
        "star_particles_plotted": int(star_idx.size),
        "cold_gas_particles_total": int(cold_pos.shape[0]),
        "cold_gas_particles_within_radius": int(cold_near.size),
        "cold_gas_particles_plotted": int(cold_idx.size),
        "heliocentric_center": center,
        "mean_heliocentric_velocity": mean_velocity,
    }
    return fig, info


def gas_contour_payload(context, box=fornax_core_radius, nbins=50, smooth_sigma=2.0, unit="deg"):
    gas_mask = np.asarray(context["dw_cold_gas_mask"], dtype=bool)
    cold_gas_mass_nh = (
        context["df"].loc[gas_mask, "m"].to_numpy(dtype=float)
        * context["df"].loc[gas_mask, "nh"].to_numpy(dtype=float)
    )
    mass_high, x_grid, y_grid, cd_smooth = Analysis.est_gas_contour(
        context["cold_gas_x_kpc"],
        context["cold_gas_y_kpc"],
        cold_gas_mass_nh,
        box=box,
        nbins=nbins,
        cd_threshold=cd_threshold,
        smooth_sigma=smooth_sigma,
        unit=unit,
        d_mean=context["d_mean"],
    )
    return {
        "mass_above_threshold": float(mass_high),
        "x_grid": x_grid,
        "y_grid": y_grid,
        "cd_smooth": cd_smooth,
        "unit": unit,
        "cd_threshold": float(cd_threshold),
    }


def save_gas_contour_plot(context, output_path, gas_payload=None):
    payload = gas_contour_payload(context) if gas_payload is None else gas_payload
    x_grid = payload["x_grid"]
    y_grid = payload["y_grid"]
    cd_smooth = payload["cd_smooth"]
    has_finite_density = cd_smooth.size > 0 and np.isfinite(cd_smooth).any()
    density_max = np.nanmax(cd_smooth) if has_finite_density else np.nan

    fig, ax = plt.subplots(figsize=(5.5, 5.0), dpi=140)
    if has_finite_density and density_max > payload["cd_threshold"]:
        levels = np.linspace(payload["cd_threshold"], density_max, 5)
        ax.contour(x_grid, y_grid, cd_smooth.T, levels=levels, colors="powderblue")
    unit = payload["unit"]
    ax.scatter(
        context["x_kpc"] if unit == "kpc" else np.degrees(context["x_kpc"] / context["d_mean"]),
        context["y_kpc"] if unit == "kpc" else np.degrees(context["y_kpc"] / context["d_mean"]),
        s=1,
        c="black",
        alpha=0.12,
    )
    ax.set_title("Cold gas column-density contour")
    ax.set_xlabel(f"X ({unit})")
    ax.set_ylabel(f"Y ({unit})")
    ax.set_xlim(2, -2) if unit == "deg" else ax.set_xlim(5, -5)
    ax.set_ylim(-2, 2) if unit == "deg" else ax.set_ylim(-5, 5)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)
    return output_path


def _json_ready(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    return value


def write_summary_json(summary, output_path):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(_json_ready(summary), indent=2), encoding="utf-8")
    return output_path


def run_check(
    snapshot_path,
    output_dir=None,
    save_html=True,
    save_gas_png=True,
    save_yz_png=True,
    initial_placement_angles_deg=DEFAULT_INITIAL_PLACEMENT_ANGLES_DEG,
    max_star_points=200000,
    max_cold_gas_points=100000,
    seed=0,
):
    context = load_snapshot_context_from_path(snapshot_path)
    paths = make_output_paths(snapshot_path, output_dir=output_dir)
    paths["output_dir"].mkdir(parents=True, exist_ok=True)

    figure, scene_info = build_3d_figure(
        context,
        max_star_points=max_star_points,
        max_cold_gas_points=max_cold_gas_points,
        seed=seed,
    )
    gas_payload = gas_contour_payload(context)
    orientation = orientation_report(context)
    long_axis_correction = long_axis_orientation_correction(context, radius=20.0)
    placement_update = placement_angle_update(
        long_axis_correction,
        initial_angles_deg=initial_placement_angles_deg,
    )

    if save_html:
        figure.write_html(paths["html"], include_plotlyjs="cdn")
    else:
        paths["html"] = None

    if save_gas_png:
        save_gas_contour_plot(context, paths["gas_png"], gas_payload=gas_payload)
    else:
        paths["gas_png"] = None

    if save_yz_png:
        save_yz_orientation_plot(
            context,
            paths["yz_png"],
            correction=long_axis_correction,
            radius=20.0,
            max_star_points=max_star_points,
            max_gas_points=max_cold_gas_points,
            seed=seed,
        )
    else:
        paths["yz_png"] = None

    summary = {
        "code_version": CODE_VERSION,
        "snapshot_path": context["snapshot_path"],
        "model_dir": context["model_dir"],
        "snapshot_num": context["snapshot_num"],
        "tsnap": context["tsnap"],
        "scene": scene_info,
        "gas_contour": {
            "mass_above_threshold": gas_payload["mass_above_threshold"],
            "unit": gas_payload["unit"],
            "cd_threshold": gas_payload["cd_threshold"],
        },
        "orientation": orientation,
        "long_axis_correction": long_axis_correction,
        "placement_angle_update": placement_update,
        "outputs": paths,
    }
    write_summary_json(summary, paths["summary_json"])
    return context, figure, summary


def parse_args():
    parser = argparse.ArgumentParser(
        description="3D visual inspection and IC rotation-correction check for one snapshot."
    )
    parser.add_argument("snapshot_path", help="Path to snapshot_###.hdf5")
    parser.add_argument(
        "--output-dir",
        help="Directory for generated files. Defaults to ../sandbox_runs/3dcheck.",
    )
    parser.add_argument("--no-html", action="store_true", help="Do not write the Plotly 3D HTML file.")
    parser.add_argument("--no-gas-png", action="store_true", help="Do not write the gas contour PNG.")
    parser.add_argument("--no-yz-png", action="store_true", help="Do not write the centered Y-Z orientation PNG.")
    parser.add_argument(
        "--initial-placement-angles",
        nargs=3,
        type=float,
        default=DEFAULT_INITIAL_PLACEMENT_ANGLES_DEG,
        metavar=("AX", "AY", "AZ"),
        help="Original dwarf placement angles in degrees. Default: 83 45 0.",
    )
    parser.add_argument("--max-plot-stars", type=int, default=200000, help="Maximum stars drawn in the 3D HTML.")
    parser.add_argument("--max-plot-gas", type=int, default=100000, help="Maximum cold gas points drawn in the 3D HTML.")
    parser.add_argument("--seed", type=int, default=0, help="Deterministic sampling seed.")
    return parser.parse_args()


def main():
    args = parse_args()
    _, _, summary = run_check(
        args.snapshot_path,
        output_dir=args.output_dir,
        save_html=not args.no_html,
        save_gas_png=not args.no_gas_png,
        save_yz_png=not args.no_yz_png,
        initial_placement_angles_deg=args.initial_placement_angles,
        max_star_points=args.max_plot_stars,
        max_cold_gas_points=args.max_plot_gas,
        seed=args.seed,
    )

    heliocentric = summary["orientation"]["heliocentric"]
    print(f"Snapshot: {summary['snapshot_path']}")
    print(f"Time: {summary['tsnap']:.3f} Gyr")
    print(f"3D stars plotted: {summary['scene']['star_particles_plotted']} / {summary['scene']['star_particles_total']}")
    print(
        "Heliocentric legacy angles deg (ax, ay, az): "
        f"{np.asarray(heliocentric['legacy_angles_deg'])}"
    )
    print(
        "Heliocentric correction sequence deg: "
        f"Rz({heliocentric['correction_sequence_deg']['first_Rz']:.6f}), "
        f"then Ry({heliocentric['correction_sequence_deg']['then_Ry']:.6f})"
    )
    print(f"Corrected unit vector: {np.asarray(heliocentric['corrected_unit_vector'])}")
    print(f"Long-axis LOS angle: {summary['long_axis_correction']['angle_deg']:.6f} deg")
    placement = summary["placement_angle_update"]
    print(
        "Initial placement angles deg (ax, ay, az): "
        f"{np.asarray(placement['initial_angles_deg'])}"
    )
    print(
        "Corrected placement angles deg (ax, ay, az): "
        f"{np.asarray(placement['corrected_angles_wrapped_deg'])}"
    )
    print(
        "Corrected placement angles deg for legacy rot_3d convention (ax, ay, az): "
        f"{np.asarray(placement['legacy_rot3d_coordinate_system']['corrected_angles_wrapped_deg'])}"
    )
    print(f"Placement angle convention: {placement['convention']}")
    print(f"Summary JSON: {summary['outputs']['summary_json']}")
    if summary["outputs"]["html"] is not None:
        print(f"3D HTML: {summary['outputs']['html']}")
    if summary["outputs"]["gas_png"] is not None:
        print(f"Gas contour PNG: {summary['outputs']['gas_png']}")
    if summary["outputs"]["yz_png"] is not None:
        print(f"Y-Z orientation PNG: {summary['outputs']['yz_png']}")


if __name__ == "__main__":
    main()
