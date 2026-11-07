"""Static, COM-preserving SO(3) morphology search; never writes ICs or submits jobs.

The caller supplies a verified native catalogue, a fixed canonical E/N/LOS
basis, and explicit reference data. Search scores are engineering distances,
not observational likelihoods. A snapshot rotation is not an initial pose fit.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation
from scipy.stats import qmc


def proper_rotation(matrix):
    matrix = np.asarray(matrix, dtype=float)
    if (matrix.shape != (3, 3) or not np.isfinite(matrix).all()
            or not np.allclose(matrix.T @ matrix, np.eye(3), atol=2e-10, rtol=0)
            or abs(np.linalg.det(matrix) - 1) > 2e-10):
        raise ValueError('Expected a finite proper SO(3) rotation')
    return matrix


def sobol_rotations(power, seed):
    """Scrambled Sobol points mapped to Haar-uniform unit quaternions (xyzw)."""
    u = qmc.Sobol(d=3, scramble=True, seed=int(seed)).random_base2(int(power))
    a, b = np.sqrt(1 - u[:, 0]), np.sqrt(u[:, 0])
    quaternion = np.column_stack((a * np.sin(2 * np.pi * u[:, 1]),
                                 a * np.cos(2 * np.pi * u[:, 1]),
                                 b * np.sin(2 * np.pi * u[:, 2]),
                                 b * np.cos(2 * np.pi * u[:, 2])))
    return Rotation.from_quat(quaternion).as_matrix()


def project_tan(offset, rotation, basis, distance_kpc, *, with_distance=False):
    """Perspective TAN projection with each star's rotated line-of-sight depth.

    basis has native east, north, LOS as rows. Coordinates are plane degrees,
    not spherical arc angles. Float32 offset storage is retained; small basis
    matrices and the perspective denominator are evaluated with checked precision.
    """
    offset = np.asarray(offset)
    if offset.dtype.kind != 'f':
        offset = offset.astype(np.float32)
    rotation = proper_rotation(rotation)
    basis = proper_rotation(basis)
    distance = float(distance_kpc)
    if (offset.ndim != 2 or offset.shape[1] != 3 or not np.isfinite(offset).all()
            or not np.isfinite(distance) or distance <= 0):
        raise ValueError('Finite Nx3 offsets and positive COM distance required')
    directions = (basis @ rotation).astype(offset.dtype)
    xyz = offset @ directions.T
    depth = xyz[:, 2].astype(np.float64) + distance
    if np.any(depth <= 0):
        raise ValueError('Source member behind TAN plane; cannot silently drop it')
    xy = np.rad2deg(xyz[:, :2] / depth[:, None])
    if with_distance:
        return xy, np.sqrt(np.sum(xyz.astype(np.float64) ** 2, axis=1)
                           + distance * distance + 2 * distance * xyz[:, 2])
    return xy


def cic_image(xy, weight, half_width_deg, npix):
    """Cloud-in-cell image, with fixed sky pixels and no image recentering."""
    xy, weight = np.asarray(xy), np.asarray(weight)
    if (xy.ndim != 2 or xy.shape[1] != 2 or weight.shape != (len(xy),)
            or not np.isfinite(xy).all() or not np.isfinite(weight).all()
            or np.any(weight < 0) or half_width_deg <= 0 or int(npix) < 2):
        raise ValueError('Invalid image inputs')
    npix = int(npix)
    pixel = (xy + half_width_deg) * (npix / (2 * half_width_deg)) - .5
    lower = np.floor(pixel).astype(np.int64)
    fraction = pixel - lower
    image = np.zeros(npix * npix, dtype=float)
    for i, j in ((0, 0), (0, 1), (1, 0), (1, 1)):
        x, y = lower[:, 0] + i, lower[:, 1] + j
        use = (x >= 0) & (x < npix) & (y >= 0) & (y < npix)
        w = weight * (fraction[:, 0] if i else 1 - fraction[:, 0])
        w = w * (fraction[:, 1] if j else 1 - fraction[:, 1])
        image += np.bincount(x[use] * npix + y[use], weights=w[use], minlength=npix * npix)
    return image.reshape(npix, npix)


class DensityObjective:
    """Multi-window, smoothed mass-fraction Hellinger distances.

    Normalize by total source mass before applying any field window. Thus lost
    or redistributed light is penalized; no window is independently renormalized.
    Absolute mass/brightness and canonical shape require separate final checks.
    """
    def __init__(self, offset, mass, basis, distance_kpc, reference_xy,
                 reference_mass, windows, npix=48, sigma_pixels=1.2):
        self.offset = np.asarray(offset)
        mass, reference_mass = np.asarray(mass), np.asarray(reference_mass)
        if (mass.shape != (len(self.offset),) or np.any(mass <= 0)
                or not np.isfinite(mass).all() or reference_mass.ndim != 1
                or np.any(reference_mass <= 0) or not np.isfinite(reference_mass).all()
                or not np.isfinite(sigma_pixels) or sigma_pixels <= 0):
            raise ValueError('Positive finite source/reference masses and smoothing required')
        self.weight = mass / mass.sum(dtype=np.float64)
        self.basis = proper_rotation(basis)
        self.distance = float(distance_kpc)
        self.windows = [float(w) for w in windows]
        if not self.windows or any(not np.isfinite(w) or w <= 0 for w in self.windows):
            raise ValueError('Explicit positive image windows required')
        self.npix, self.sigma = int(npix), float(sigma_pixels)
        ref_weight = reference_mass / reference_mass.sum(dtype=np.float64)
        self.reference = [self.image(reference_xy, ref_weight, w) for w in self.windows]
        self.reference_sqrt = [np.sqrt(image) for image in self.reference]
        self.reference_enclosed = [image.sum() for image in self.reference]
        if min(self.reference_enclosed) <= 0:
            raise ValueError('Reference has no light in an evaluation window')

    def image(self, xy, weight, window):
        return gaussian_filter(cic_image(xy, weight, window, self.npix),
                               self.sigma, mode='constant')

    def components(self, rotation):
        xy = project_tan(self.offset, rotation, self.basis, self.distance)
        return np.array([np.sum((np.sqrt(self.image(xy, self.weight, w)) - ref) ** 2) / norm
                         for w, ref, norm in zip(self.windows, self.reference_sqrt,
                                                self.reference_enclosed)])

    def __call__(self, rotation):
        return float(np.mean(self.components(rotation)))


def rotation_separation_deg(left, right):
    return float(np.rad2deg(Rotation.from_matrix(proper_rotation(left)
                            @ proper_rotation(right).T).magnitude()))


def distinct_candidates(matrices, scores, count, minimum_separation_deg):
    selected = []
    for index in np.argsort(scores, kind='stable'):
        if not np.isfinite(scores[index]):
            continue
        if all(rotation_separation_deg(matrices[index], matrices[old]) >= minimum_separation_deg
               for old in selected):
            selected.append(int(index))
            if len(selected) >= count:
                break
    return selected


def global_search(objective, *, power, seed, extra_rotations=(), progress=None):
    matrices = np.concatenate((np.eye(3)[None], sobol_rotations(power, seed),
                               np.asarray(extra_rotations).reshape(-1, 3, 3)))
    scores = np.empty(len(matrices))
    for index, matrix in enumerate(matrices):
        try:
            scores[index] = objective(matrix)
        except ValueError:
            scores[index] = np.inf
        if progress and (index % 512 == 0 or index == len(matrices) - 1):
            progress(dict(stage='global', completed=index + 1, total=len(matrices),
                          best_score=float(np.min(scores[:index + 1]))))
    if not np.isfinite(scores).any():
        raise ValueError('No projectable rotations in the search')
    return matrices, scores


def refine(objective, matrices, *, radius_deg, max_evaluations, xtol_rad=1e-4,
           ftol=1e-7, progress=None):
    radius = np.deg2rad(float(radius_deg))
    if radius <= 0 or not np.isfinite(radius):
        raise ValueError('Positive local rotation-vector bounds required')
    results = []
    for index, matrix in enumerate(matrices):
        matrix = proper_rotation(matrix)
        initial_score = objective(matrix)
        def function(delta):
            try:
                return objective(Rotation.from_rotvec(delta).as_matrix() @ matrix)
            except ValueError:
                return np.inf
        fit = minimize(function, np.zeros(3), method='Powell',
                       bounds=[(-radius, radius)] * 3,
                       options=dict(maxfev=int(max_evaluations), xtol=float(xtol_rad), ftol=float(ftol)))
        # An exhausted or unsuccessful local search cannot make a seed worse.
        improved = np.isfinite(fit.fun) and fit.fun < initial_score
        final = Rotation.from_rotvec(fit.x).as_matrix() @ matrix if improved else matrix
        result = dict(matrix=final, score=objective(final), seed_score=initial_score,
                      solver_success=bool(fit.success), solver_message=str(fit.message),
                      evaluations=int(fit.nfev), local_delta_rotvec_rad=fit.x if improved else np.zeros(3))
        results.append(result)
        if progress:
            progress(dict(stage='refine', completed=index + 1, total=len(matrices),
                          score=result['score'], solver_success=result['solver_success']))
    return results
