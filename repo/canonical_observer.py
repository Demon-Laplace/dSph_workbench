"""Versioned mock-observation measurements with an explicit physical observer.

All functions are pure: they neither write files nor alter supplied arrays.
Positions, velocities and the observer use one caller-specified native frame.
The canonical per-star LOS is calculated in that frame before any camera.
No target position, proper motion, physical model rotation or observer fitting
is implemented.  The Walker implementation is limited to simulated catalogues
with membership probabilities one and velocity measurement errors zero.
"""
from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from typing import Any

import numpy as np


METHOD_VERSION = "observing_v1"
ASTROPY_NATIVE_DIRECT_CANONICAL_ID = "astropy_native_direct_v1"
GALACTIC_CENTER_ICRS_DEG = (266.4051, -28.936175)
GALACTOCENTRIC_ROLL_DEG = 0.0
ROTATION_ATOL = 1.0e-10
WALKER_VARIANCE_RELATIVE_TOL = 64.0 * np.finfo(np.float64).eps
PM_KMS_PER_MASYR_KPC = 4.740470463533349


def _vector(value: Any, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if result.shape != (3,) or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be a finite vector of shape (3,)")
    return result


def _xyz(value: Any, name: str, n: int | None = None) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if result.ndim != 2 or result.shape[1] != 3 or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be a finite array of shape (N, 3)")
    if n is not None and len(result) != n:
        raise ValueError(f"{name} must have {n} rows")
    return result


def _ids(value: Any) -> np.ndarray:
    result = np.asarray(value)
    if result.ndim != 1 or result.dtype != np.dtype("uint64"):
        raise ValueError("particle_id must be a one-dimensional uint64 array")
    if len(np.unique(result)) != len(result):
        raise ValueError("particle_id must be unique")
    return result


def _scalar(value: Any, name: str, *, positive: bool = False) -> float:
    if np.ndim(value) != 0:
        raise ValueError(f"{name} must be a finite scalar")
    result = float(value)
    if not np.isfinite(result) or (positive and result <= 0):
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    return result


def _series(value: Any, name: str, n: int) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if result.shape != (n,) or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be a finite array of shape ({n},)")
    return result


def _mask(value: Any, n: int, name: str) -> np.ndarray:
    result = np.asarray(value)
    if result.shape != (n,) or result.dtype != np.dtype("bool"):
        raise ValueError(f"{name} must be a boolean array of shape ({n},)")
    return result


def _proper_rotation(value: Any, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if result.shape != (3, 3) or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be a finite 3 by 3 matrix")
    if (not np.allclose(result.T @ result, np.eye(3), atol=ROTATION_ATOL, rtol=0)
            or abs(float(np.linalg.det(result)) - 1.0) > ROTATION_ATOL):
        raise ValueError(f"{name} must be orthogonal with determinant +1")
    return result


def _icrs_basis(unit_vectors: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    ra = np.arctan2(unit_vectors[..., 1], unit_vectors[..., 0])
    dec = np.arctan2(unit_vectors[..., 2], np.hypot(unit_vectors[..., 0], unit_vectors[..., 1]))
    east = np.stack((-np.sin(ra), np.cos(ra), np.zeros_like(ra)), axis=-1)
    north = np.stack((-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)), axis=-1)
    return ra, dec, east, north


def native_observables(
    position_raw: Any, velocity_raw: Any, particle_id: Any, *,
    host_center_raw: Any, host_bulk_raw: Any,
    length_to_kpc: float, velocity_to_kms: float,
    sun_position_kpc: Any, sun_velocity_kms: Any,
) -> dict[str, Any]:
    """Form the physical heliocentric six-vector and LOS without a camera.

    ``host_center_raw`` and ``host_bulk_raw`` have the same raw units as the
    corresponding particle arrays.  The Sun is relative to the MW, in true
    kpc and km/s in the same native axes.  Both unit factors are mandatory.
    Thus h=(p-H)*length_to_kpc-S; w=(v-B)*velocity_to_kms; u=w-Vs.
    ``vlos_native_kms=sum(u*(h/D), axis=1)`` is the canonical scalar to retain
    with the returned uint64 IDs.  Solar reflex and perspective are retained.
    """
    ids = _ids(particle_id)
    p = _xyz(position_raw, "position_raw", len(ids))
    v = _xyz(velocity_raw, "velocity_raw", len(ids))
    H = _vector(host_center_raw, "host_center_raw")
    B = _vector(host_bulk_raw, "host_bulk_raw")
    S = _vector(sun_position_kpc, "sun_position_kpc")
    Vs = _vector(sun_velocity_kms, "sun_velocity_kms")
    length = _scalar(length_to_kpc, "length_to_kpc", positive=True)
    speed = _scalar(velocity_to_kms, "velocity_to_kms", positive=True)
    r = (p - H) * length
    w = (v - B) * speed
    h = r - S
    u = w - Vs
    distance = np.sqrt(np.sum(h * h, axis=1))
    if np.any(distance <= 0):
        raise ValueError("a particle coincides with the observer; its LOS is undefined")
    ray = h / distance[:, None]
    los = np.sum(u * ray, axis=1)
    return dict(
        method_version=METHOD_VERSION, particle_id=ids.copy(),
        host_relative_position_kpc=r, grf_velocity_native_kms=w,
        h_native_kpc=h, u_native_kms=u, distance_kpc=distance,
        ray_native=ray, vlos_native_kms=los, vlos_source="native_direct",
        observer=dict(host_center_raw=H.copy(), host_bulk_raw=B.copy(),
                      length_to_kpc=length, velocity_to_kms=speed,
                      sun_position_kpc=S.copy(), sun_velocity_kms=Vs.copy()),
    )


def attach_frozen_vlos(
    native_result: Mapping[str, Any], frozen_particle_id: Any,
    frozen_vlos_kms: Any, *, atol: float = 1.0e-10,
) -> dict[str, Any]:
    """Verify complete ID coverage, then keep the saved LOS as primary.

    IDs must be unique uint64 on both sides, with identical complete sets.
    The newly calculated native scalar is used only for the absolute-error
    comparison; the returned ``vlos_native_kms`` copies the frozen floating
    array in native-result ID order, preserving its stored dtype and values.
    The recomputed scalar and its difference are separately named diagnostics.
    All returned state is copied, so neither caller input is mutated/aliased.
    An existing saved LOS is never silently replaced by its recomputation.
    """
    ids = _ids(native_result["particle_id"])
    frozen_ids = _ids(frozen_particle_id)
    if len(frozen_ids) != len(ids):
        raise ValueError("frozen IDs must completely cover the native catalogue, without extras")
    frozen = np.asarray(frozen_vlos_kms)
    if (frozen.shape != (len(ids),) or frozen.dtype.kind != "f"
            or not np.all(np.isfinite(frozen))):
        raise ValueError("frozen_vlos_kms must be a finite floating array with one value per ID")
    recomputed = _series(native_result["vlos_native_kms"], "native vlos_native_kms", len(ids))
    tolerance = _scalar(atol, "atol")
    if tolerance < 0:
        raise ValueError("atol must be nonnegative")
    order = np.argsort(frozen_ids)
    locations = np.searchsorted(frozen_ids[order], ids)
    if np.any(locations >= len(ids)):
        raise ValueError("frozen IDs do not match the native catalogue")
    joined = order[locations]
    if not np.array_equal(frozen_ids[joined], ids):
        raise ValueError("frozen IDs do not match the native catalogue")
    primary = frozen[joined].copy()
    difference = recomputed - primary
    maximum = float(np.max(abs(difference))) if len(ids) else 0.0
    if maximum > tolerance:
        raise ValueError(f"frozen/native LOS mismatch: max absolute difference {maximum} km/s exceeds {tolerance}")
    result = deepcopy(dict(native_result))
    result.update(vlos_native_kms=primary, vlos_source="frozen_id_matched",
                  vlos_native_recomputed_kms=recomputed.copy(),
                  vlos_recomputed_minus_frozen_kms=difference.copy(),
                  frozen_vlos_max_abs_difference_kms=maximum,
                  frozen_vlos_comparison_atol_kms=tolerance,
                  frozen_vlos_complete_ID_coverage=True)
    return result


def physical_sky_bridge(
    *, disc_north_native: Any, sun_position_kpc: Any, sun_velocity_kms: Any,
) -> dict[str, Any]:
    """Build a proper native-heliocentric-to-ICRS rotation from the disc/Sun.

    The caller supplies the already signed north pole.  Its length is
    normalized, but its sign is never selected or flipped using the target.
    a is the Sun's in-disc radial unit vector; b=k cross a.  The intermediate
    Galactocentric basis is T=stack(-a,-b,k).  The full Astropy transform uses
    |S|, z_sun=S dot k and galcen_v_sun=T Vs, preserving the actual Sun.
    The fixed conventional GC ICRS direction and roll zero are recorded;
    there is no observed Fornax target or proper-motion registration.
    Only this function requires Astropy; no runtime paths are inserted.
    """
    from astropy import units as units
    from astropy.coordinates import CartesianDifferential, CartesianRepresentation, Galactocentric, ICRS, SkyCoord

    k = _vector(disc_north_native, "disc_north_native")
    knorm = float(np.linalg.norm(k))
    if knorm == 0:
        raise ValueError("disc_north_native must be nonzero and already signed")
    k = k / knorm
    S = _vector(sun_position_kpc, "sun_position_kpc")
    Vs = _vector(sun_velocity_kms, "sun_velocity_kms")
    height = float(S @ k)
    in_plane = S - height * k
    radius = float(np.linalg.norm(in_plane))
    if radius <= np.finfo(float).eps * max(float(np.linalg.norm(S)), 1.0):
        raise ValueError("Sun-GC in-disc direction is undefined for a Sun on the disc pole")
    a = in_plane / radius
    b = np.cross(k, a)
    T = _proper_rotation(np.vstack((-a, -b, k)), "native_to_disc_GC")
    Sd = T @ S
    Vd = T @ Vs
    attributes = dict(
        galcen_distance=float(np.linalg.norm(S)) * units.kpc,
        z_sun=height * units.kpc,
        galcen_v_sun=CartesianDifferential(Vd * units.km / units.s),
        galcen_coord=SkyCoord(ra=GALACTIC_CENTER_ICRS_DEG[0] * units.deg,
                             dec=GALACTIC_CENTER_ICRS_DEG[1] * units.deg, frame="icrs"),
        roll=GALACTOCENTRIC_ROLL_DEG * units.deg,
    )

    def transform_gc_points(points: np.ndarray) -> np.ndarray:
        coords = points @ T.T
        frame = Galactocentric(CartesianRepresentation(coords.T * units.kpc), **attributes)
        return frame.transform_to(ICRS()).cartesian.xyz.to_value(units.kpc).T

    zero = transform_gc_points(S[None, :])[0]
    if not np.allclose(zero, 0.0, atol=ROTATION_ATOL, rtol=0):
        raise ValueError("explicit Galactocentric frame does not place the supplied Sun at the origin")
    M = _proper_rotation((transform_gc_points(S[None, :] + np.eye(3)) - zero).T,
                         "physical_to_icrs")
    return dict(
        method_version=METHOD_VERSION, physical_to_icrs=M.copy(),
        native_to_disc_GC=T.copy(), disc_north_native=k.copy(),
        disc_north_sign_policy="Caller supplied; never inferred from the dwarf or observed PM.",
        sun_position_kpc=S.copy(), sun_velocity_kms=Vs.copy(),
        sun_position_disc_kpc=Sd, sun_velocity_disc_kms=Vd,
        galcen_distance_kpc=float(np.linalg.norm(S)), z_sun_kpc=height,
        sun_in_plane_radius_kpc=radius,
        galcen_icrs_deg=tuple(GALACTIC_CENTER_ICRS_DEG), roll_deg=GALACTOCENTRIC_ROLL_DEG,
        physical_to_icrs_determinant=float(np.linalg.det(M)),
        bridge_orthogonality_max_error=float(np.max(abs(M.T @ M - np.eye(3)))),
    )


def project_sky(
    h_native_kpc: Any, particle_id: Any, *, center_h_native_kpc: Any,
    physical_to_icrs: Any, aperture_kpc: float = 5.0,
) -> dict[str, Any]:
    """Project the physical rays; no velocity or canonical LOS is modified.

    The caller selects the centre (normally the full pre-existing-star mass
    centre) independently.  ``east_deg/north_deg`` resolve true spherical
    rho along the centre's physical E/N; they are the Walker coordinates.
    ``delta_RA_deg`` is the raw wrapped RA difference, with no cos(dec).
    Physical x/y are gnomonic: D0*(ray dot E/N)/(ray dot n0).  Therefore
    Rproj=D0*tan(rho), and the circular mask is Rproj<aperture_kpc in front
    of the camera, not a 3-D spherical cut.  Per-star ICRS E/N are returned
    for optional proper-motion output from the same physical M and u.
    """
    ids = _ids(particle_id)
    h = _xyz(h_native_kpc, "h_native_kpc", len(ids))
    center = _vector(center_h_native_kpc, "center_h_native_kpc")
    M = _proper_rotation(physical_to_icrs, "physical_to_icrs")
    aperture = _scalar(aperture_kpc, "aperture_kpc", positive=True)
    D0 = float(np.linalg.norm(center))
    if D0 == 0:
        raise ValueError("center_h_native_kpc must not coincide with the observer")
    distance = np.sqrt(np.sum(h * h, axis=1))
    if np.any(distance <= 0):
        raise ValueError("a particle coincides with the observer")
    rays = h / distance[:, None]
    n0 = center / D0
    ra0, dec0, eeq, neq = _icrs_basis(M @ n0)
    east = M.T @ eeq
    north = M.T @ neq
    R = _proper_rotation(np.vstack((east, north, n0)), "pointing_camera")
    ex = rays @ east
    ny = rays @ north
    denominator = rays @ n0
    transverse = np.hypot(ex, ny)
    rho = np.arctan2(transverse, denominator)
    factor = np.divide(np.degrees(rho), transverse, out=np.zeros_like(transverse), where=transverse > 0)
    xk = np.divide(D0 * ex, denominator, out=np.full_like(ex, np.inf), where=denominator != 0)
    yk = np.divide(D0 * ny, denominator, out=np.full_like(ny, np.inf), where=denominator != 0)
    radius = np.hypot(xk, yk)
    eqrays = rays @ M.T
    ra, dec, per_east, per_north = _icrs_basis(eqrays)
    ra_deg = np.degrees(ra) % 360.0
    dec_deg = np.degrees(dec)
    ra0_deg = float(np.degrees(ra0) % 360.0)
    dec0_deg = float(np.degrees(dec0))
    return dict(
        method_version=METHOD_VERSION, particle_id=ids.copy(), distance_kpc=distance,
        sky_source="explicit_native_physical_bridge",
        center_h_native_kpc=center.copy(), center_distance_kpc=D0,
        center_RA_deg=ra0_deg, center_Dec_deg=dec0_deg,
        physical_to_icrs=M.copy(), pointing_camera=R.copy(),
        east_native=east, north_native=north, LOS_native=n0,
        ray_icrs=eqrays, east_icrs=per_east, north_icrs=per_north,
        RA_deg=ra_deg, Dec_deg=dec_deg,
        delta_RA_deg=(ra_deg - ra0_deg + 180.0) % 360.0 - 180.0,
        delta_Dec_deg=dec_deg - dec0_deg,
        east_deg=ex * factor, north_deg=ny * factor, rho_deg=np.degrees(rho),
        east_kpc=xk, north_kpc=yk, projected_radius_kpc=radius,
        front_mask=denominator > 0,
        fit_R5_mask=(denominator > 0) & (radius < aperture),
        aperture_kpc=aperture,
    )


def project_catalog_sky(
    ra_deg: Any, dec_deg: Any, distance_kpc: Any, particle_id: Any, *,
    center_ra_deg: float, center_dec_deg: float, center_distance_kpc: float,
    physical_to_icrs: Any | None = None, aperture_kpc: float = 5.0,
) -> dict[str, Any]:
    """Use existing source-catalogue ICRS coordinates without redefining sky.

    RA/Dec/D and the centre are caller-supplied source data, not fitting
    targets.  The supplied RA/Dec values are copied unchanged.  There is no
    Sun or disc reconstruction, coordinate registration, or observed-PM roll.
    The same spherical Walker components and gnomonic radius as project_sky
    are derived directly from these ICRS rays.  Raw wrapped DeltaRA does not
    contain a cos(dec) factor.

    Optional physical_to_icrs must be an independently verified source-frame
    bridge; it is never inferred from the catalogue.  It only expresses
    native GRF vectors and the camera in the already supplied sky.  If it is
    absent, native E/N/LOS and the bridge/camera are None.  The native-GRF
    arrow routine then refuses to guess those missing physical axes.
    """
    ids = _ids(particle_id)
    n = len(ids)
    ra_input = _series(ra_deg, "ra_deg", n)
    dec_input = _series(dec_deg, "dec_deg", n)
    distance = _series(distance_kpc, "distance_kpc", n)
    if np.any(abs(dec_input) > 90.0) or np.any(distance <= 0):
        raise ValueError("catalogue declinations must lie in [-90,90] and distances must be positive")
    ra0 = _scalar(center_ra_deg, "center_ra_deg")
    dec0 = _scalar(center_dec_deg, "center_dec_deg")
    if abs(dec0) > 90.0:
        raise ValueError("center_dec_deg must lie in [-90,90]")
    D0 = _scalar(center_distance_kpc, "center_distance_kpc", positive=True)
    aperture = _scalar(aperture_kpc, "aperture_kpc", positive=True)
    ra, dec = np.deg2rad(ra_input), np.deg2rad(dec_input)
    rays = np.stack((np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)), axis=1)
    star_east = np.stack((-np.sin(ra), np.cos(ra), np.zeros_like(ra)), axis=1)
    star_north = np.stack((-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)), axis=1)
    r0, d0 = np.deg2rad([ra0, dec0])
    n0eq = np.array([np.cos(d0) * np.cos(r0), np.cos(d0) * np.sin(r0), np.sin(d0)])
    e0eq = np.array([-np.sin(r0), np.cos(r0), 0.0])
    north0eq = np.array([-np.sin(d0) * np.cos(r0), -np.sin(d0) * np.sin(r0), np.cos(d0)])
    ex, ny, denominator = rays @ e0eq, rays @ north0eq, rays @ n0eq
    transverse = np.hypot(ex, ny)
    rho = np.arctan2(transverse, denominator)
    factor = np.divide(np.degrees(rho), transverse, out=np.zeros_like(transverse), where=transverse > 0)
    xk = np.divide(D0 * ex, denominator, out=np.full_like(ex, np.inf), where=denominator != 0)
    yk = np.divide(D0 * ny, denominator, out=np.full_like(ny, np.inf), where=denominator != 0)
    radius = np.hypot(xk, yk)
    M = None if physical_to_icrs is None else _proper_rotation(physical_to_icrs, "physical_to_icrs")
    east = None if M is None else M.T @ e0eq
    north = None if M is None else M.T @ north0eq
    n0native = None if M is None else M.T @ n0eq
    camera = None if M is None else _proper_rotation(np.vstack((east, north, n0native)), "pointing_camera")
    return dict(
        method_version=METHOD_VERSION, sky_source="source_catalog_ICRS",
        particle_id=ids.copy(), distance_kpc=distance.copy(),
        center_h_native_kpc=None if M is None else D0 * n0native,
        center_distance_kpc=D0, center_RA_deg=ra0, center_Dec_deg=dec0,
        physical_to_icrs=None if M is None else M.copy(), pointing_camera=camera,
        east_native=east, north_native=north, LOS_native=n0native,
        center_ray_icrs=n0eq, center_east_icrs=e0eq, center_north_icrs=north0eq,
        ray_icrs=rays, east_icrs=star_east, north_icrs=star_north,
        RA_deg=ra_input.copy(), Dec_deg=dec_input.copy(),
        delta_RA_deg=(ra_input - ra0 + 180.0) % 360.0 - 180.0,
        delta_Dec_deg=dec_input - dec0,
        east_deg=ex * factor, north_deg=ny * factor, rho_deg=np.degrees(rho),
        east_kpc=xk, north_kpc=yk, projected_radius_kpc=radius,
        front_mask=denominator > 0,
        fit_R5_mask=(denominator > 0) & (radius < aperture), aperture_kpc=aperture,
    )


def walker61_zero_error(
    projected: Mapping[str, Any], vlos_kms: Any, *, fit_mask: Any,
) -> dict[str, Any]:
    """Walker directional fits for P_i=1 and measurement error zero only.

    For theta0=0,3,...,180 degrees north toward east, fit
    v0+k*(rho_deg*cos(theta-theta0)).  The equivalent coordinate is formed
    from the spherical E/N components, never raw DeltaRA/DeltaDec.  Select
    the largest absolute directional k, not the norm of a 2-D plane fit.
    sigma_mle is sqrt(mean(residual**2)), with divisor N and no ddof
    correction.  This is not the observational membership/error MLE.
    A direction with variance no larger than 64*machine-epsilon times the
    full centred 2-D spatial variance is numerically singular and invalid.
    Exact 0/90/180 trigonometric values prevent artificial endpoint leverage.
    """
    ids = _ids(projected["particle_id"])
    n = len(ids)
    los = _series(vlos_kms, "vlos_kms", n)
    use = _mask(fit_mask, n, "fit_mask")
    if int(use.sum()) < 3:
        raise ValueError("Walker fit requires at least three selected stars")
    x = _series(projected["east_deg"], "east_deg", n)[use]
    y = _series(projected["north_deg"], "north_deg", n)[use]
    v = los[use]
    angles = np.arange(0.0, 181.0, 3.0)
    k = np.full(61, np.nan)
    v0 = np.full(61, np.nan)
    sigma = np.full(61, np.nan)
    valid = np.zeros(61, dtype=bool)
    leverage = np.zeros(61)
    dx = x - np.mean(x)
    dy = y - np.mean(y)
    spatial_variance_sum = float(np.dot(dx, dx) + np.dot(dy, dy))
    leverage_floor = WALKER_VARIANCE_RELATIVE_TOL * spatial_variance_sum
    sine = np.sin(np.deg2rad(angles))
    cosine = np.cos(np.deg2rad(angles))
    sine[0] = sine[-1] = 0.0
    cosine[0], cosine[-1] = 1.0, -1.0
    sine[30], cosine[30] = 1.0, 0.0
    vmean = float(np.mean(v))
    dv = v - vmean
    for j in range(61):
        coordinate = x * sine[j] + y * cosine[j]
        mean = float(np.mean(coordinate))
        dt = coordinate - mean
        variance_sum = float(np.dot(dt, dt))
        leverage[j] = variance_sum
        if variance_sum <= leverage_floor:
            continue
        k[j] = float(np.dot(dt, dv) / variance_sum)
        v0[j] = vmean - k[j] * mean
        residual = v - (v0[j] + k[j] * coordinate)
        sigma[j] = float(np.sqrt(np.mean(residual * residual)))
        valid[j] = True
    if not np.any(valid):
        raise ValueError("all Walker directions have zero spatial leverage")
    best = int(np.argmax(np.where(valid, abs(k), -np.inf)))
    return dict(
        method_version=METHOD_VERSION, membership_probability=1.0,
        velocity_measurement_error_kms=0.0, N_fit=int(use.sum()),
        fit_particle_id=ids[use].copy(), theta0_scan_deg=angles,
        k_scan_kms_per_deg=k, v0_scan_kms=v0, sigma_mle_scan_kms=sigma,
        valid_scan=valid, invalid_scan=~valid,
        scan_leverage_sum_deg2=leverage,
        spatial_variance_sum_deg2=spatial_variance_sum,
        numerical_leverage_floor_deg2=leverage_floor,
        variance_relative_tolerance=WALKER_VARIANCE_RELATIVE_TOL,
        invalid_direction_definition="Centred directional variance at or below relative numerical leverage floor.",
        theta_index=best, theta0_best_deg=float(angles[best]),
        k_abs=float(abs(k[best])), k_signed=float(k[best]),
        v0_kms=float(v0[best]), sigma_mle_kms=float(sigma[best]),
        sigma_raw_mle_kms=float(np.std(v, ddof=0)),
        dispersion_definition="MLE residual RMS, divisor N; no sample ddof correction.",
    )


def velocity_grid40(
    projected: Mapping[str, Any], vlos_kms: Any, mass_msun: Any,
    grf_velocity_native_kms: Any, *, display_v0_kms: float,
    include_grf_arrow: bool = True,
) -> dict[str, Any]:
    """40x40 occupied-cell display and independent same-aperture GRF arrow.

    The square display is [-5,+5] kpc on each gnomonic axis.  Every occupied
    cell is retained, including stars outside the circular fitting aperture.
    Only one supplied constant is subtracted for display.  Counts/means are
    indexed [east_bin,north_bin]; a plot uses mean.T with the returned XY
    angular vertex mesh, east left and north up, and fixed +/-12 km/s colors.
    The arrow uses projected.fit_R5_mask only, mass-mean native w=v-B,
    projected on the physical centre E/N.  No Sun velocity is subtracted.
    With include_grf_arrow=False, no native bridge, mass or GRF velocity is
    required: mass_msun and grf_velocity_native_kms may be None.  Sky vertices
    then come directly from the source ICRS centre.  All arrow values are None
    and its supported/direction flags are False; no bridge is guessed.
    """
    ids = _ids(projected["particle_id"])
    n = len(ids)
    los = _series(vlos_kms, "vlos_kms", n)
    if not isinstance(include_grf_arrow, (bool, np.bool_)):
        raise ValueError("include_grf_arrow must be boolean")
    x = np.asarray(projected["east_kpc"], dtype=np.float64)
    y = np.asarray(projected["north_kpc"], dtype=np.float64)
    if x.shape != (n,) or y.shape != (n,):
        raise ValueError("projected physical axes must have one value per ID")
    front = _mask(projected["front_mask"], n, "front_mask")
    zero = _scalar(display_v0_kms, "display_v0_kms")
    edges = np.linspace(-5.0, 5.0, 41)
    shown = front & np.isfinite(x) & np.isfinite(y) & (abs(x) <= 5.0) & (abs(y) <= 5.0)
    count = np.histogram2d(x[shown], y[shown], bins=(edges, edges))[0].astype(np.int64)
    total = np.histogram2d(x[shown], y[shown], bins=(edges, edges), weights=los[shown] - zero)[0]
    mean = np.divide(total, count, out=np.full((40, 40), np.nan), where=count > 0)
    D0 = _scalar(projected["center_distance_kpc"], "center_distance_kpc", positive=True)
    xx, yy = np.meshgrid(edges, edges, indexing="xy")
    ra0 = _scalar(projected["center_RA_deg"], "center_RA_deg")
    dec0 = _scalar(projected["center_Dec_deg"], "center_Dec_deg")
    meanw = motion = motion_unit = arrow_ids = None
    amplitude = 0.0
    if include_grf_arrow:
        if projected.get("physical_to_icrs") is None:
            raise ValueError("native GRF arrow requires an explicit verified source physical_to_icrs bridge; catalogue sky does not determine it")
        mass = _series(mass_msun, "mass_msun", n)
        if np.any(mass < 0):
            raise ValueError("mass_msun must be nonnegative")
        w = _xyz(grf_velocity_native_kms, "grf_velocity_native_kms", n)
        fit = _mask(projected["fit_R5_mask"], n, "fit_R5_mask")
        if not np.any(fit) or float(np.sum(mass[fit])) <= 0:
            raise ValueError("same-aperture GRF arrow requires positive selected mass")
        meanw = np.average(w[fit], axis=0, weights=mass[fit])
        east = _vector(projected["east_native"], "east_native")
        north = _vector(projected["north_native"], "north_native")
        motion = np.array([meanw @ east, meanw @ north])
        amplitude = float(np.linalg.norm(motion))
        motion_unit = motion / amplitude if amplitude > 0 else np.full(2, np.nan)
        arrow_ids = ids[fit].copy()
        M = _proper_rotation(projected["physical_to_icrs"], "physical_to_icrs")
        n0 = _vector(projected["LOS_native"], "LOS_native")
        rays = n0 + xx[..., None] * east / D0 + yy[..., None] * north / D0
        rays = rays / np.linalg.norm(rays, axis=-1)[..., None]
        ra, dec, _, _ = _icrs_basis(rays @ M.T)
    else:
        # Existing catalogue sky remains primary even without any native bridge.
        ar, dr = np.deg2rad([ra0, dec0])
        n0eq = np.array([np.cos(dr) * np.cos(ar), np.cos(dr) * np.sin(ar), np.sin(dr)])
        e0eq = np.array([-np.sin(ar), np.cos(ar), 0.0])
        north0eq = np.array([-np.sin(dr) * np.cos(ar), -np.sin(dr) * np.sin(ar), np.cos(dr)])
        rays = n0eq + xx[..., None] * e0eq / D0 + yy[..., None] * north0eq / D0
        rays = rays / np.linalg.norm(rays, axis=-1)[..., None]
        ra, dec, _, _ = _icrs_basis(rays)
    return dict(
        method_version=METHOD_VERSION, count=count, mean_HRF_minus_v0=mean,
        edges_gnomonic_kpc=edges, display_v0_kms=zero,
        display_mask=shown, display_particle_id=ids[shown].copy(),
        delta_RA_vertices_deg=(np.degrees(ra) - ra0 + 180.0) % 360.0 - 180.0,
        delta_Dec_vertices_deg=np.degrees(dec) - dec0,
        display_color_limits_kms=(-12.0, 12.0),
        GRF_massmean_native_kms=meanw, GRF_massmean_EN_kms=motion,
        GRF_motion_unit_EN=motion_unit, GRF_motion_direction_defined=amplitude > 0,
        GRF_arrow_particle_id=arrow_ids, GRF_arrow_supported=bool(include_grf_arrow),
        GRF_arrow_status="computed_from_explicit_source_bridge" if include_grf_arrow else "not_requested_no_bridge_inferred",
        arrow_label="Motion relative to the Milky Way" if include_grf_arrow else None,
        grid_array_order="east_bin,north_bin; transpose once for XY vertex mesh",
    )


def passive_camera(
    h_native_kpc: Any, u_native_kms: Any, *, camera_R: Any, physical_to_icrs: Any,
) -> dict[str, Any]:
    """Change representation only, with the synchronized output sky bridge.

    Column-vector convention: hcam=R h, ucam=R u, Mcam=M R.T.
    Native per-ID LOS remains the canonical saved value and is not replaced.
    Reflections and non-orthogonal matrices are rejected.
    """
    h = _xyz(h_native_kpc, "h_native_kpc")
    u = _xyz(u_native_kms, "u_native_kms", len(h))
    R = _proper_rotation(camera_R, "camera_R")
    M = _proper_rotation(physical_to_icrs, "physical_to_icrs")
    bridge = _proper_rotation(M @ R.T, "camera_to_icrs")
    return dict(
        method_version=METHOD_VERSION, h_camera_kpc=h @ R.T,
        u_camera_kms=u @ R.T, camera_to_icrs=bridge,
        camera_R=R.copy(), camera_determinant=float(np.linalg.det(R)),
        camera_orthogonality_max_error=float(np.max(abs(R.T @ R - np.eye(3)))),
    )


def proper_motion_components(
    h_native_kpc: Any, u_native_kms: Any, *, physical_to_icrs: Any,
) -> dict[str, Any]:
    """Actual per-point heliocentric PM using each point's own ICRS E/N/D.

    The inputs are already heliocentric physical vectors.  This function
    performs no additional Sun/bulk subtraction and no fitting or camera.
    """
    h = _xyz(h_native_kpc, "h_native_kpc")
    u = _xyz(u_native_kms, "u_native_kms", len(h))
    M = _proper_rotation(physical_to_icrs, "physical_to_icrs")
    D = np.sqrt(np.sum(h * h, axis=1))
    if np.any(D <= 0):
        raise ValueError("proper motion is undefined at the observer")
    rays = (h / D[:, None]) @ M.T
    ra, dec, east, north = _icrs_basis(rays)
    velocity = u @ M.T
    return dict(
        RA_deg=np.degrees(ra) % 360.0, Dec_deg=np.degrees(dec), distance_kpc=D,
        vlos_kms=np.sum(u * (h / D[:, None]), axis=1),
        mu_alpha_cosdec_masyr=np.sum(velocity * east, axis=1) / (PM_KMS_PER_MASYR_KPC * D),
        mu_delta_masyr=np.sum(velocity * north, axis=1) / (PM_KMS_PER_MASYR_KPC * D),
        east_icrs=east, north_icrs=north,
    )


def astropy_native_direct_observables(
    position_raw: Any, velocity_raw: Any, particle_id: Any, *,
    host_center_raw: Any, host_bulk_raw: Any,
    length_to_kpc: float, velocity_to_kms: float,
    galcen_distance_kpc: float, z_sun_kpc: float, galcen_v_sun_kms: Any,
    galcen_icrs_deg: Any, roll_deg: float, native_to_GC: Any,
) -> dict[str, Any]:
    """One complete, explicit Galactocentric 6-D to ICRS transformation.

    Raw inputs first form native r=(p-H)*L and w=(v-B)*V.  The supplied
    proper native_to_GC matrix then represents them in the declared Astropy
    Galactocentric frame, which transforms the full CartesianDifferential
    once.  The pure native LOS is separately calculated without a camera.
    Every frame parameter is mandatory.  The CLI canonical v1 contract uses
    native_to_GC=I and fixed explicit parameters; this pure function does not
    infer a Sun, flip a native axis or fit a target.
    """
    import astropy
    from astropy import units as units
    from astropy.coordinates import CartesianDifferential, CartesianRepresentation, Galactocentric, ICRS, SkyCoord

    ids = _ids(particle_id)
    p = _xyz(position_raw, "position_raw", len(ids))
    v = _xyz(velocity_raw, "velocity_raw", len(ids))
    H = _vector(host_center_raw, "host_center_raw")
    B = _vector(host_bulk_raw, "host_bulk_raw")
    length = _scalar(length_to_kpc, "length_to_kpc", positive=True)
    speed = _scalar(velocity_to_kms, "velocity_to_kms", positive=True)
    T = _proper_rotation(native_to_GC, "native_to_GC")
    distance_gc = _scalar(galcen_distance_kpc, "galcen_distance_kpc", positive=True)
    height = _scalar(z_sun_kpc, "z_sun_kpc")
    if abs(height) >= distance_gc:
        raise ValueError("abs(z_sun_kpc) must be smaller than galcen_distance_kpc")
    solar_velocity_gc = _vector(galcen_v_sun_kms, "galcen_v_sun_kms")
    gc_angles = np.asarray(galcen_icrs_deg, dtype=float)
    if gc_angles.shape != (2,) or not np.all(np.isfinite(gc_angles)) or abs(gc_angles[1]) > 90:
        raise ValueError("galcen_icrs_deg must contain a finite RA and valid Dec")
    roll = _scalar(roll_deg, "roll_deg")
    frame_attributes = dict(
        galcen_distance=distance_gc * units.kpc, z_sun=height * units.kpc,
        galcen_v_sun=CartesianDifferential(solar_velocity_gc * units.km / units.s),
        galcen_coord=SkyCoord(ra=gc_angles[0] * units.deg, dec=gc_angles[1] * units.deg, frame="icrs"),
        roll=roll * units.deg,
    )
    frame = Galactocentric(**frame_attributes)
    # Recover the observer from this exact explicit affine frame, rather
    # than confusing galcen_distance with its projected cylindrical radius.
    zero = ICRS(CartesianRepresentation(np.zeros(3) * units.kpc,
                differentials=CartesianDifferential(np.zeros(3) * units.km / units.s)))
    observer_gc = zero.transform_to(frame)
    S_gc = observer_gc.cartesian.xyz.to_value(units.kpc)
    Vs_gc = observer_gc.velocity.d_xyz.to_value(units.km / units.s)
    S, Vs = T.T @ S_gc, T.T @ Vs_gc
    native = native_observables(p, v, ids, host_center_raw=H, host_bulk_raw=B,
        length_to_kpc=length, velocity_to_kms=speed, sun_position_kpc=S, sun_velocity_kms=Vs)
    gc_position = native["host_relative_position_kpc"] @ T.T
    gc_velocity = native["grf_velocity_native_kms"] @ T.T
    representation = CartesianRepresentation(gc_position.T * units.kpc,
        differentials=CartesianDifferential(gc_velocity.T * units.km / units.s))
    full_sky = Galactocentric(representation, **frame_attributes).transform_to(ICRS())

    def position_transform(native_positions):
        q = np.asarray(native_positions) @ T.T
        return Galactocentric(CartesianRepresentation(q.T * units.kpc), **frame_attributes).transform_to(ICRS()).cartesian.xyz.to_value(units.kpc).T

    origin = position_transform(S[None, :])[0]
    M = _proper_rotation((position_transform(S[None, :] + np.eye(3)) - origin).T,
                         "physical_to_icrs")
    calculated = proper_motion_components(native["h_native_kpc"], native["u_native_kms"], physical_to_icrs=M)
    direct = dict(
        RA_deg=full_sky.ra.to_value(units.deg), Dec_deg=full_sky.dec.to_value(units.deg),
        distance_astropy_kpc=full_sky.distance.to_value(units.kpc),
        vlos_astropy_kms=full_sky.radial_velocity.to_value(units.km / units.s),
        mu_alpha_cosdec_masyr=full_sky.pm_ra_cosdec.to_value(units.mas / units.yr),
        mu_delta_masyr=full_sky.pm_dec.to_value(units.mas / units.yr),
    )
    def maximum(values):
        return float(np.max(abs(values))) if np.size(values) else 0.0
    ra_error = (calculated["RA_deg"] - direct["RA_deg"] + 180) % 360 - 180
    checks = dict(
        Sun_origin_max_error_kpc=maximum(origin),
        distance_max_error_kpc=maximum(native["distance_kpc"] - direct["distance_astropy_kpc"]),
        LOS_max_error_kms=maximum(native["vlos_native_kms"] - direct["vlos_astropy_kms"]),
        RA_max_error_deg=maximum(ra_error), Dec_max_error_deg=maximum(calculated["Dec_deg"] - direct["Dec_deg"]),
        PMRA_max_error_masyr=maximum(calculated["mu_alpha_cosdec_masyr"] - direct["mu_alpha_cosdec_masyr"]),
        PMDec_max_error_masyr=maximum(calculated["mu_delta_masyr"] - direct["mu_delta_masyr"]),
    )
    if any(value > 1e-8 for value in checks.values()):
        raise ValueError(f"full Astropy and native physical 6-D calculations disagree: {checks}")
    native.update(direct)
    is_canonical = (
        abs(distance_gc - 8.122) <= 1e-12 and abs(height - .0208) <= 1e-12
        and np.allclose(solar_velocity_gc, [12.9, 245.6, 7.78], rtol=0, atol=1e-12)
        and np.allclose(gc_angles, [266.4051, -28.936175], rtol=0, atol=1e-12)
        and abs(roll) <= 1e-12 and np.allclose(T, np.eye(3), rtol=0, atol=1e-12)
    )
    native.update(
        canonical_id=ASTROPY_NATIVE_DIRECT_CANONICAL_ID if is_canonical else None,
        physical_to_icrs=M.copy(), native_to_GC=T.copy(),
        frame_contract=dict(
            galcen_distance_kpc=distance_gc, z_sun_kpc=height,
            galcen_v_sun_kms=solar_velocity_gc.copy(), galcen_icrs_deg=gc_angles.copy(), roll_deg=roll,
            native_to_GC=T.copy(), recovered_Sun_GC_kpc=S_gc.copy(), recovered_Sun_velocity_GC_kms=Vs_gc.copy(),
            recovered_Sun_native_kpc=S.copy(), recovered_Sun_velocity_native_kms=Vs.copy(),
            astropy_version=astropy.__version__, full_GC_6D_transform_count=1,
            no_second_solar_subtraction=True,
        ),
        full_astropy_independent_check=dict(state="PASS", **checks),
    )
    return native


def aperture_motion_summaries(
    native_result: Mapping[str, Any], projected: Mapping[str, Any], mass_msun: Any, *, fit_mask: Any,
) -> dict[str, Any]:
    """Keep three distinct same-aperture PM definitions explicitly named.

    Primary display vectors share the full image centre's E/N and D0:
    same-R5 massmean GRF, solar reflex, and HRF bulk velocities.  Separately
    report PM of the R5 mass-centre six-vector at its own sightline/distance,
    and the mass-weighted mean of individual per-star PM components.
    """
    ids = _ids(native_result["particle_id"])
    if not np.array_equal(ids, _ids(projected["particle_id"])):
        raise ValueError("native and sky catalogues must already share exactly the same ID order")
    use = _mask(fit_mask, len(ids), "fit_mask")
    mass = _series(mass_msun, "mass_msun", len(ids))
    if np.any(mass < 0) or not np.any(use) or float(np.sum(mass[use])) <= 0:
        raise ValueError("motion summary requires nonnegative masses and positive selected mass")
    h = _xyz(native_result["h_native_kpc"], "h_native_kpc", len(ids))
    u = _xyz(native_result["u_native_kms"], "u_native_kms", len(ids))
    w = _xyz(native_result["grf_velocity_native_kms"], "grf_velocity_native_kms", len(ids))
    Vs = _vector(native_result["observer"]["sun_velocity_kms"], "sun_velocity_kms")
    M = _proper_rotation(projected["physical_to_icrs"], "physical_to_icrs")
    east = _vector(projected["east_native"], "east_native")
    north = _vector(projected["north_native"], "north_native")
    D0 = _scalar(projected["center_distance_kpc"], "center_distance_kpc", positive=True)
    hm = np.average(h[use], axis=0, weights=mass[use])
    um = np.average(u[use], axis=0, weights=mass[use])
    wm = np.average(w[use], axis=0, weights=mass[use])
    if not np.allclose(um, wm - Vs, atol=1e-9, rtol=0):
        raise ValueError("same-aperture HRF/GRF means do not share the declared solar velocity")
    def EN(vector):
        return np.array([vector @ east, vector @ north])
    grf, solar, hrf = EN(wm), EN(-Vs), EN(um)
    def direction(vector):
        amplitude = float(np.linalg.norm(vector))
        return vector / amplitude if amplitude > 0 else None
    if "mu_alpha_cosdec_masyr" in native_result and "mu_delta_masyr" in native_result:
        all_pm = {name: _series(native_result[name], name, len(ids))
                  for name in ("mu_alpha_cosdec_masyr", "mu_delta_masyr")}
        individual_pm_source = "Native result's explicitly supplied per-star PM columns (full Astropy in canonical mode)."
    else:
        all_pm = proper_motion_components(h, u, physical_to_icrs=M)
        individual_pm_source = "Per-star projection of supplied heliocentric h/u through the same explicit bridge."
    point = proper_motion_components(hm[None, :], um[None, :], physical_to_icrs=M)
    return dict(
        N_R5=int(use.sum()), selected_mass_msun=float(np.sum(mass[use])),
        common_display_center_bulk_motion=dict(
            definition="Same R5 massmean velocities, at the full input image-centre E/N and D0; primary display definition.",
            center_RA_deg=float(projected["center_RA_deg"]), center_Dec_deg=float(projected["center_Dec_deg"]),
            center_distance_kpc=D0, east_native=east.copy(), north_native=north.copy(),
            GRF_massmean_native_kms=wm, HRF_massmean_native_kms=um, solar_reflex_native_kms=-Vs,
            GRF_EN_kms=grf, solar_reflex_EN_kms=solar, HRF_EN_kms=hrf,
            GRF_EN_masyr=grf / (PM_KMS_PER_MASYR_KPC * D0),
            solar_reflex_EN_masyr=solar / (PM_KMS_PER_MASYR_KPC * D0),
            HRF_EN_masyr=hrf / (PM_KMS_PER_MASYR_KPC * D0),
            GRF_unit_EN=direction(grf), HRF_unit_EN=direction(hrf),
            component_sum_max_error_kms=float(np.max(abs(grf + solar - hrf))),
        ),
        R5_COM_true6D_system_PM=dict(
            definition="Proper motion of the R5 mass-centre six-vector at its own E/N and distance; not the primary display vector.",
            h_native_kpc=hm, u_native_kms=um, grf_velocity_native_kms=wm,
            RA_deg=float(point["RA_deg"][0]), Dec_deg=float(point["Dec_deg"][0]), distance_kpc=float(point["distance_kpc"][0]),
            vlos_kms=float(point["vlos_kms"][0]),
            mu_alpha_cosdec_masyr=float(point["mu_alpha_cosdec_masyr"][0]), mu_delta_masyr=float(point["mu_delta_masyr"][0]),
        ),
        massmean_individual_PM_R5=dict(
            definition="Mass-weighted mean of individual stellar PM components using each star's E/N and distance; not the COM-point PM.",
            individual_PM_source=individual_pm_source,
            mu_alpha_cosdec_masyr=float(np.average(all_pm["mu_alpha_cosdec_masyr"][use], weights=mass[use])),
            mu_delta_masyr=float(np.average(all_pm["mu_delta_masyr"][use], weights=mass[use])),
        ),
    )


# Production adapter: the scientific implementation above is a verbatim copy
# of methods/observing_v1/observation.py, pinned in the migration receipt.
CANONICAL_OBSERVER_ID = "astropy_native_direct_v1"
OBSERVER_FRAME = dict(galcen_distance_kpc=8.122, z_sun_kpc=.0208,
    galcen_v_sun_kms=[12.9,245.6,7.78],galcen_icrs_deg=[266.4051,-28.936175],
    roll_deg=0.,native_to_GC=[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
LENGTH_CODE_TO_KPC = 1.00000064815577
VELOCITY_CODE_TO_KMS = 1.0
SKY_FIELD_HALF_WIDTH_DEG = 2.1
KINEMATIC_TRACER_DEFINITION = "current_dwarf_members_old_and_new"

def observer_contract_sha256():
    import hashlib, json
    data=dict(observer_id=CANONICAL_OBSERVER_ID,frame=OBSERVER_FRAME,
        length_code_to_kpc=LENGTH_CODE_TO_KPC,velocity_code_to_kms=VELOCITY_CODE_TO_KMS,
        host_bulk="TP2/3/4 massmean within true 3D r<5kpc of shrink MW centre",
        default_motion="Heliocentric proper motion",display_half_width_deg=2.1,
        display_projection="raw DeltaRA wrapped and raw DeltaDec",walker="61direction spherical Rproj<5kpc",
        kinematic_tracer=KINEMATIC_TRACER_DEFINITION)
    return hashlib.sha256(json.dumps(data,sort_keys=True,separators=(",",":")).encode()).hexdigest()

def require_canonical_csv(path):
    """Refuse mixed/legacy CSV instead of relabelling or silently continuing it."""
    from pathlib import Path
    p=Path(path)
    if not p.exists(): return
    comments=[]
    with p.open("r",encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"): comments.append(line.strip())
            elif line.strip(): break
    if canonical_csv_comment().strip() not in comments:
        raise ValueError("Legacy or different-observer CSV: use a NEW analysis output path and rerun GetInfo; existing CSV is preserved: "+str(p))

def canonical_csv_comment():
    import os,hashlib
    from pathlib import Path
    default=Path.cwd()/'output'/'analysis'/CANONICAL_OBSERVER_ID/'cohort_ids.npz'
    cohort=Path(os.environ.get('DSPH_PREEXISTING_ID_PATH',str(default))).resolve()
    if not cohort.is_file():raise ValueError('Frozen cohort file is required to create or append canonical CSV: '+str(cohort))
    return "# canonical_observer_id="+CANONICAL_OBSERVER_ID+" observer_config_sha256="+observer_contract_sha256()+" cohort_file_sha256="+hashlib.sha256(cohort.read_bytes()).hexdigest()+"\n"

def production_observables(position_grf_kpc, velocity_grf_kms, particle_id):
    """Inputs already have H/B and units applied exactly once by center_on_MW."""
    core=astropy_native_direct_observables(position_grf_kpc,velocity_grf_kms,particle_id,
        host_center_raw=[0,0,0],host_bulk_raw=[0,0,0],length_to_kpc=1.,velocity_to_kms=1.,**OBSERVER_FRAME)
    h,u=core["h_native_kpc"],core["u_native_kms"]
    p,w=core["host_relative_position_kpc"],core["grf_velocity_native_kms"]
    result=dict(x=p[:,0],y=p[:,1],z=p[:,2],vx=w[:,0],vy=w[:,1],vz=w[:,2],
        xh=h[:,0],yh=h[:,1],zh=h[:,2],vxh=u[:,0],vyh=u[:,1],vzh=u[:,2],
        id=core["particle_id"],ra=core["RA_deg"],dec=core["Dec_deg"],rh=core["distance_kpc"],
        vlos=core["vlos_native_kms"],pmra=core["mu_alpha_cosdec_masyr"],pmdec=core["mu_delta_masyr"],
        canonical_observer_id=CANONICAL_OBSERVER_ID,observer_config_sha256=observer_contract_sha256(),_native=core)
    return result

def production_sky(coords, center_h_native_kpc):
    native=coords["_native"]
    return project_sky(native["h_native_kpc"],native["particle_id"],
        center_h_native_kpc=center_h_native_kpc,physical_to_icrs=native["physical_to_icrs"],aperture_kpc=5.)

def production_motion(coords, sky, mass_msun, selected):
    return aperture_motion_summaries(coords["_native"],sky,mass_msun,fit_mask=selected)["common_display_center_bulk_motion"]

def angle_velocity_grid(delta_ra_deg,delta_dec_deg,vlos_kms,v0_kms,npix=40):
    """Raw DeltaRA/DeltaDec display bins; footprint never defines Walker sample."""
    x,y,v=map(np.asarray,(delta_ra_deg,delta_dec_deg,vlos_kms))
    edges=np.linspace(-SKY_FIELD_HALF_WIDTH_DEG,SKY_FIELD_HALF_WIDTH_DEG,int(npix)+1)
    count,_,_=np.histogram2d(x,y,bins=(edges,edges))
    total,_,_=np.histogram2d(x,y,bins=(edges,edges),weights=v-float(v0_kms))
    field=np.divide(total,count,out=np.full_like(total,np.nan),where=count>0)
    return field,(-SKY_FIELD_HALF_WIDTH_DEG,SKY_FIELD_HALF_WIDTH_DEG,-SKY_FIELD_HALF_WIDTH_DEG,SKY_FIELD_HALF_WIDTH_DEG)

def tangent_kpc_to_raw_offsets(projected,x_kpc,y_kpc):
    """Exact inverse gnomonic mapping for contour/shape overlays on raw sky axes."""
    a,d=np.deg2rad([projected['center_RA_deg'],projected['center_Dec_deg']])
    n0=np.array([np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)])
    e=np.array([-np.sin(a),np.cos(a),0.]);north=np.array([-np.sin(d)*np.cos(a),-np.sin(d)*np.sin(a),np.cos(d)])
    x,y=np.broadcast_arrays(np.asarray(x_kpc),np.asarray(y_kpc))
    q=projected['center_distance_kpc']*n0+x[...,None]*e+y[...,None]*north
    ra=np.degrees(np.arctan2(q[...,1],q[...,0]));dec=np.degrees(np.arctan2(q[...,2],np.hypot(q[...,0],q[...,1])))
    return (ra-projected['center_RA_deg']+180)%360-180,dec-projected['center_Dec_deg']

def raw_sky_pixel_area_arcsec2(center_dec_deg,npix,half_width_deg=2.1):
    """Exact spherical area for each row of a constant raw-RA/Dec grid."""
    edges=np.radians(center_dec_deg+np.linspace(-half_width_deg,half_width_deg,npix+1))
    return np.radians(2*half_width_deg/npix)*np.diff(np.sin(edges))*(180./np.pi*3600.)**2

def frozen_cohort_mask(snapshot_particle_id,folder_path):
    """Require an explicitly selected-then-frozen tracer cohort, never infer it by density."""
    import os,hashlib
    from pathlib import Path
    default=Path(folder_path).resolve()/'analysis'/CANONICAL_OBSERVER_ID/'cohort_ids.npz'
    path=Path(os.environ.get('DSPH_PREEXISTING_ID_PATH',str(default))).resolve()
    if not path.is_file():
        raise ValueError('Canonical analysis requires a selected-then-frozen cohort_ids.npz (particle_id uint64). Set DSPH_PREEXISTING_ID_PATH or create the explicit cohort file; no dynamic-membership fallback: '+str(path))
    with np.load(str(path),allow_pickle=False) as data:
        cohort=_ids(data['particle_id'])
    if len(cohort)==0:raise ValueError('Frozen tracer cohort is empty')
    current=_ids(snapshot_particle_id)
    missing=np.setdiff1d(cohort,current,assume_unique=True)
    if len(missing):raise ValueError('Snapshot lacks frozen cohort ID(s), count='+str(len(missing)))
    mask=np.isin(current,cohort)
    return mask,dict(cohort_definition='selected_then_frozen',cohort_path=str(path),
        cohort_file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),cohort_N=len(cohort),
        sorted_cohort_id_sha256=hashlib.sha256(np.sort(cohort).astype('<u8').tobytes()).hexdigest(),
        not_claimed_to_be_all_initial_IC_stars=True)

def write_snapshot_observer_metadata(folder_path,snapshot_num,host_state,cohort,star_coords,walker,motion):
    """Small atomic sidecar; source snapshot/particle arrays are never modified."""
    import os,json,hashlib,tempfile
    from pathlib import Path
    root=Path(os.environ.get('DSPH_CANONICAL_OUTPUT_ROOT',str(Path(folder_path).resolve()/'analysis'/CANONICAL_OBSERVER_ID)))
    root.mkdir(parents=True,exist_ok=True)
    def serial(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [serial(v) for v in x]
        return x
    sources=[]
    for path in sorted(Path(folder_path).glob('snapshot_'+str(int(snapshot_num)).zfill(3)+'*.hdf5')):
        st=path.stat()
        # Source fingerprints are explicitly named: no unsupported full-snapshot SHA claim.
        sources.append(dict(path=str(path.resolve()),size_bytes=st.st_size,mtime_ns=st.st_mtime_ns,
            fingerprint_kind='path_size_mtime_only'))
    result=dict(canonical_observer_id=CANONICAL_OBSERVER_ID,observer_config_sha256=observer_contract_sha256(),
        observer_frame=OBSERVER_FRAME,host=host_state,cohort=cohort,snapshot_num=int(snapshot_num),sources=sources,
        kinematic_tracer=KINEMATIC_TRACER_DEFINITION,
        implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_selected_id_sha256=hashlib.sha256(star_coords['id'].astype('<u8').tobytes()).hexdigest(),
        walker_R5={k:v for k,v in walker.items() if k not in ('fit_particle_id',)},motion_R5=motion)
    target=root/('observer_snapshot_'+str(int(snapshot_num)).zfill(3)+'.json')
    with tempfile.NamedTemporaryFile('w',encoding='utf-8',dir=str(root),prefix='.observer-',suffix='.json',delete=False) as f:
        temp=Path(f.name);json.dump(serial(result),f,ensure_ascii=False,indent=2)
    os.replace(str(temp),str(target))
    return str(target)
