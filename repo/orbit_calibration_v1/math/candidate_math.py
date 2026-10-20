"""Numerical helpers for explicitly physical orbit/orientation *candidates*.

No observer, camera, snapshot, IC or production parameter file is changed here.
The empirical C_E/k_L mapping is a local calibration, not a new evolved orbit.
Position and velocity units must be consistent within each calculation. Merger
placement helpers accept the pre-unit-conversion native code units.
"""
from __future__ import annotations

import itertools
import math
import warnings
import numpy as np
from scipy.spatial.transform import Rotation

LEGACY_RAD_PER_DEG = 0.017453


def _vector(value, name):
    value = np.asarray(value, dtype=float)
    if value.shape != (3,) or not np.isfinite(value).all():
        raise ValueError(name + ' must be a finite 3-vector')
    return value


def unit(value):
    value = _vector(value, 'vector')
    length = float(np.linalg.norm(value))
    if length <= 0:
        raise ValueError('zero vector has no direction')
    return value / length


def decompose_state(position, velocity):
    r, v = _vector(position, 'position'), _vector(velocity, 'velocity')
    rhat = unit(r)
    vr = float(v @ rhat)
    tangential = v - vr * rhat
    vt = float(np.linalg.norm(tangential))
    if vt <= 1e-12 * max(1.0, float(np.linalg.norm(v))):
        raise ValueError('radial state has no well-conditioned orbit plane')
    that = tangential / vt
    lhat = unit(np.cross(rhat, that))
    return dict(radius=float(np.linalg.norm(r)), speed=float(np.linalg.norm(v)),
                vr=vr, vt=vt, rhat=rhat, that=that, lhat=lhat,
                angular_momentum=float(np.linalg.norm(np.cross(r, v))))


def empirical_coefficients(r0, v0, rf, vf):
    initial, final = decompose_state(r0, v0), decompose_state(rf, vf)
    return dict(C_E=final['speed'] ** 2 - initial['speed'] ** 2,
                k_L=final['angular_momentum'] / initial['angular_momentum'],
                initial=initial, final=final)


def solve_initial_velocity(r0, v0, target_r, target_v, C_E, k_L, radial_sign=None):
    """Keep seed radial/tangential directions; solve empirical final E/L targets.

    radial_sign defaults to the seed branch (zero defaults inward). Pass -1 for
    an explicit inward branch. No pericentre condition is imposed or certified.
    """
    seed, target = decompose_state(r0, v0), decompose_state(target_r, target_v)
    if not np.isfinite([C_E, k_L]).all() or k_L <= 0:
        raise ValueError('C_E must be finite and k_L positive')
    vt0 = target['angular_momentum'] / (float(k_L) * seed['radius'])
    speed0_squared = target['speed'] ** 2 - float(C_E)
    vr0_squared = speed0_squared - vt0 ** 2
    if speed0_squared < 0 or vr0_squared < 0:
        raise ValueError('empirical target has no real initial radial velocity')
    sign = (-1 if seed['vr'] <= 0 else 1) if radial_sign is None else radial_sign
    if sign not in (-1, 1):
        raise ValueError('radial_sign must be -1 or +1')
    vr0 = sign * math.sqrt(vr0_squared)
    velocity = vr0 * seed['rhat'] + vt0 * seed['that']
    return dict(velocity=velocity, vr0=vr0, vt0=vt0,
                speed0=math.sqrt(speed0_squared), seed=seed, target=target,
                C_E=float(C_E), k_L=float(k_L), pericentre_constrained=False,
                reconstructed_target_speed_squared=float(velocity @ velocity) + C_E,
                reconstructed_target_angular_momentum=k_L * seed['radius'] * vt0)


def orbit_frame(position, velocity):
    state = decompose_state(position, velocity)
    return np.column_stack([state[key] for key in ['rhat', 'that', 'lhat']])


def orbit_plane_rotation(current_r, current_v, target_r, target_v):
    """Proper SO(3) map of radial/tangential/normal directions only.

    Different speed, radial speed and radius remain separate constraints. This
    is a physical registration candidate; a non-spherical live MW is not
    guaranteed dynamically invariant under this rotation.
    """
    source, target = orbit_frame(current_r, current_v), orbit_frame(target_r, target_v)
    result = target @ source.T
    if np.linalg.det(result) < 0.999999999 or np.max(abs(result.T @ result - np.eye(3))) > 1e-10:
        raise ValueError('orbit frame map is not a proper rotation')
    return result


def axis_angle_rotation(axis, angle_deg):
    axis = unit(axis)
    angle = np.deg2rad(float(angle_deg))
    x, y, z = axis
    cross = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]], dtype=float)
    return np.eye(3) + np.sin(angle) * cross + (1 - np.cos(angle)) * (cross @ cross)


def axis_pa_deg(vector, east, north):
    return float(np.rad2deg(np.arctan2(np.dot(vector, east), np.dot(vector, north))) % 180.0)


def pa_difference_deg(target, current):
    return float((float(target) - float(current) + 90.0) % 180.0 - 90.0)


def pa_roll_rotation(los, delta_pa_deg):
    """Increase north->east PA by delta, with E cross N = outward LOS.

    Positive standard right-handed rotation around LOS decreases that PA;
    therefore the Rodrigues angle is minus the requested PA change.
    """
    return axis_angle_rotation(los, -float(delta_pa_deg))


def legacy_axis_matrices(angles_deg, float32_trig=True):
    """Exact angle parsing/trig construction used by current ic_merge.c."""
    angles = _vector(angles_deg, 'angles_deg')
    if float32_trig:
        angles = angles.astype(np.float32)
    matrices = []
    for axis in range(3):
        angle = float(angles[axis]) * LEGACY_RAD_PER_DEG
        s, c = math.sin(angle), math.cos(angle)
        dtype = np.float32 if float32_trig else float
        s, c = dtype(s), dtype(c)
        matrix = np.zeros((3, 3), dtype=dtype)
        i, j = (axis + 1) % 3, (axis + 2) % 3
        matrix[axis, axis] = 1
        matrix[i, i] = matrix[j, j] = c
        matrix[i, j], matrix[j, i] = s, -s
        matrices.append(matrix)
    return matrices


def legacy_matrix(angles_deg, float32_trig=True):
    """Linear matrix Q/G, before per-operation and intermediate ASCII rounding."""
    rx, ry, rz = [x.astype(float) for x in legacy_axis_matrices(angles_deg, float32_trig)]
    return rz @ ry @ rx


def _nearest_proper_rotation(matrix):
    if np.linalg.det(matrix) <= 0 or np.max(abs(matrix.T @ matrix - np.eye(3))) > 2e-5:
        raise ValueError('expected a proper rotation with only small float-rounding error')
    u, _, vt = np.linalg.svd(matrix)
    sign = np.eye(3)
    sign[-1, -1] = np.linalg.det(u @ vt)
    return u @ sign @ vt


def legacy_angles_from_matrix(matrix, reference_angles=None):
    """Invert actual sign/order/.017453, then optimize nearby float32 choices.

    A conjugation by a float-rounded Q is slightly non-orthogonal. The target
    is first projected to SO(3), with its projection error reported. Periodicity
    is 2*pi/.017453, not exactly 360 parameter degrees.
    """
    matrix = np.asarray(matrix, dtype=float)
    if matrix.shape != (3, 3) or not np.isfinite(matrix).all():
        raise ValueError('matrix must be a finite 3x3 array')
    proper = _nearest_proper_rotation(matrix)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        euler = Rotation.from_matrix(proper).as_euler('xyz')
    alternatives = [euler, np.array([euler[0] + np.pi, np.pi - euler[1], euler[2] + np.pi])]
    reference = np.zeros(3) if reference_angles is None else _vector(reference_angles, 'reference_angles')
    period = 2 * np.pi / LEGACY_RAD_PER_DEG
    candidates = []
    for radians in alternatives:
        parameters = -radians / LEGACY_RAD_PER_DEG
        parameters += np.round((reference - parameters) / period) * period
        rounded = parameters.astype(np.float32)
        choices = [[np.nextafter(x, np.float32(-np.inf)), x, np.nextafter(x, np.float32(np.inf))] for x in rounded]
        for triplet in itertools.product(*choices):
            angles = np.array(triplet, dtype=np.float32)
            realized = legacy_matrix(angles)
            error = float(np.linalg.norm(realized - matrix))
            distance = float(np.linalg.norm(angles.astype(float) - reference))
            candidates.append((error, distance, angles, realized))
    # Only use the close Euler branch unless it materially worsens the matrix.
    near = min(candidates, key=lambda item: item[1])
    close_branch = [item for item in candidates if item[1] <= near[1] + 1e-3]
    error, distance, angles, realized = min(close_branch, key=lambda item: item[0])
    return dict(angles_deg=angles.astype(float), matrix=realized,
                max_abs_matrix_error=float(np.max(abs(realized - matrix))),
                frobenius_error=error, reference_distance_deg=distance,
                target_SO3_projection_error=float(np.max(abs(proper - matrix))),
                realized_orthogonality_error=float(np.max(abs(realized.T @ realized - np.eye(3)))))


def corrected_internal_angles(old_angles, cube_angles, internal_rotation, orbit_rotation=None):
    q, g = legacy_matrix(cube_angles), legacy_matrix(old_angles)
    a = np.asarray(internal_rotation, dtype=float)
    p = np.eye(3) if orbit_rotation is None else np.asarray(orbit_rotation, dtype=float)
    for name, rotation in [('internal_rotation', a), ('orbit_rotation', p)]:
        if rotation.shape != (3, 3) or not np.isfinite(rotation).all():
            raise ValueError(name + ' must be a finite 3x3 proper rotation')
        _nearest_proper_rotation(rotation)
    target = np.linalg.solve(q, a @ p @ q @ g)
    result = legacy_angles_from_matrix(target, reference_angles=old_angles)
    result.update(target_matrix=target, cube_matrix=q, original_galaxy_matrix=g,
                  native_body_rotation=a @ p,
                  native_matrix_error=float(np.max(abs(q @ result['matrix'] - a @ p @ q @ g))))
    return result


def initial_state_from_merge(galaxy_parameters, cube_parameters, mw_center, mw_bulk, dwarf_center, dwarf_bulk):
    """Relative MW/dwarf tracer state in final native code units, linear limit."""
    galaxy, cube = np.asarray(galaxy_parameters, dtype=float), np.asarray(cube_parameters, dtype=float)
    if galaxy.shape != (12,) or cube.shape != (12,):
        raise ValueError('merger parameter vectors must each contain 12 numbers')
    q, g = legacy_matrix(cube[9:12]), legacy_matrix(galaxy[9:12])
    r = q @ (g @ _vector(dwarf_center, 'dwarf_center') + galaxy[:3] - _vector(mw_center, 'mw_center'))
    v = q @ (g @ _vector(dwarf_bulk, 'dwarf_bulk') + galaxy[3:6] - _vector(mw_bulk, 'mw_bulk'))
    return r, v


def placement_for_target(target_r, target_v, new_angles, cube_angles,
                         mw_center, mw_bulk, dwarf_center, dwarf_bulk):
    """Solve mergegalaxy first 6 values with all four input offsets included.

    Common cube translation/boost cancels in MW-relative vectors. Target vectors
    must already be converted back from kpc/km/s to native code units.
    """
    q, g = legacy_matrix(cube_angles), legacy_matrix(new_angles)
    r = np.linalg.solve(q, _vector(target_r, 'target_r')) + _vector(mw_center, 'mw_center') - g @ _vector(dwarf_center, 'dwarf_center')
    v = np.linalg.solve(q, _vector(target_v, 'target_v')) + _vector(mw_bulk, 'mw_bulk') - g @ _vector(dwarf_bulk, 'dwarf_bulk')
    return np.concatenate([r, v])


def legacy_transform_float32(values, parameters):
    """Replay current C transform's operation order, with no FMA or fast-math."""
    values = np.asarray(values, dtype=np.float32).copy()
    if values.shape[-1] < 6:
        raise ValueError('last axis must contain at least position and velocity')
    parameters = np.asarray(parameters, dtype=np.float32)
    if parameters.shape != (12,):
        raise ValueError('parameters must have 12 values')
    matrices = legacy_axis_matrices(parameters[9:12])
    for offset in (0, 3):
        vector = values[..., offset:offset + 3].copy()
        for matrix in matrices:
            output = np.zeros_like(vector)
            for i in range(3):
                for j in range(3):
                    product = np.multiply(vector[..., j], matrix[i, j], dtype=np.float32)
                    output[..., i] = np.add(output[..., i], product, dtype=np.float32)
            vector = output
        values[..., offset:offset + 3] = vector
    values[..., :6] = np.add(values[..., :6], parameters[:6], dtype=np.float32)
    return values


def round_ascii_7g(values):
    values = np.asarray(values, dtype=np.float32)
    result = np.array([np.float32(format(float(x), '.7g')) for x in values.ravel()], dtype=np.float32)
    return result.reshape(values.shape)


def emulate_one_pass_particle(values8, role, galaxy_parameters, cube_parameters, final_ascii_round=True):
    """Return current C merger's particle values, excluding IDs/type ordering.

    role=0 background cube, 1 MW first galaxy, 2 dwarf second galaxy. The
    intermediate %.7g applies to all eight values even for an unrotated MW.
    """
    result = np.asarray(values8, dtype=np.float32).copy()
    if result.shape[-1] != 8 or role not in (0, 1, 2):
        raise ValueError('expected eight float fields and role 0, 1 or 2')
    if role == 2:
        result = legacy_transform_float32(result, galaxy_parameters)
    if role != 0:
        result = legacy_transform_float32(round_ascii_7g(result), cube_parameters)
    return round_ascii_7g(result) if final_ascii_round else result
