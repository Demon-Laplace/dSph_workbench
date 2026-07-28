import numpy as np
import astropy.units as u
from astropy.coordinates import CartesianDifferential, CartesianRepresentation, SkyCoord

from basefunc import Analysis
from variable import r_3d_cut, stellar_region_rhalf_multiplier


def _sum_cold_gas_mass(mass, nh, radial_mask):
    return (nh[radial_mask] * mass[radial_mask]).sum()


def _sum_hot_gas_mass(mass, nh, radial_mask):
    return (mass[radial_mask] * (1 - nh[radial_mask])).sum()


def _safe_std(values):
    values = np.asarray(values, dtype=float)
    return np.nan if values.size == 0 else np.std(values)


def _safe_mean(values):
    values = np.asarray(values, dtype=float)
    return np.nan if values.size == 0 else np.mean(values)


def old_dwarf_star_local_mask(snapshot):
    """Return a local mask selecting old dwarf stars from total_dw_star_mask arrays."""
    df = snapshot["df"]
    total_dw_star_mask = np.asarray(snapshot["total_dw_star_mask"], dtype=bool)

    if "birth" in df.columns:
        birth = df.loc[total_dw_star_mask, "birth"].to_numpy(dtype=float)
        return birth == 0

    # In the current snapshots, pre-existing collisionless dwarf stars are
    # particle types 2/3; newly formed stars are type 4 and carry birth times.
    tp = df.loc[total_dw_star_mask, "tp"].to_numpy(dtype=int)
    return tp != 4


def fit_velocity_gradient(x_kpc, y_kpc, vlos, mask=None):
    """Fit vlos = v0 + a*x + b*y and return model and residual arrays."""
    x_kpc = np.asarray(x_kpc, dtype=float)
    y_kpc = np.asarray(y_kpc, dtype=float)
    vlos = np.asarray(vlos, dtype=float)

    good = np.isfinite(x_kpc) & np.isfinite(y_kpc) & np.isfinite(vlos)
    if mask is not None:
        good &= np.asarray(mask, dtype=bool)
    if good.sum() < 3:
        return {
            "v0": np.nan,
            "a": np.nan,
            "b": np.nan,
            "grad_amp": np.nan,
            "angle_deg_from_x": np.nan,
            "v_model": np.full_like(vlos, np.nan, dtype=float),
            "v_resid": np.full_like(vlos, np.nan, dtype=float),
            "good_count": int(good.sum()),
        }

    design = np.column_stack([np.ones(good.sum()), x_kpc[good], y_kpc[good]])
    coeff, *_ = np.linalg.lstsq(design, vlos[good], rcond=None)
    v0, a, b = coeff

    v_model = np.full_like(vlos, np.nan, dtype=float)
    v_resid = np.full_like(vlos, np.nan, dtype=float)
    v_model[good] = v0 + a * x_kpc[good] + b * y_kpc[good]
    v_resid[good] = vlos[good] - v_model[good]

    return {
        "v0": float(v0),
        "a": float(a),
        "b": float(b),
        "grad_amp": float(np.hypot(a, b)),
        "angle_deg_from_x": float(np.degrees(np.arctan2(b, a))),
        "v_model": v_model,
        "v_resid": v_resid,
        "good_count": int(good.sum()),
    }


def _star_projected_kinematics(snapshot, local_mask):
    """Return projected positions and LOS velocities for a local star mask."""
    local_mask = np.asarray(local_mask, dtype=bool)
    star_coords = snapshot["star_coords"]

    x_kpc = np.asarray(snapshot["x_kpc"], dtype=float)[local_mask]
    y_kpc = np.asarray(snapshot["y_kpc"], dtype=float)[local_mask]
    vlos = Analysis.cal_vlos(
        np.asarray(star_coords["xh"], dtype=float)[local_mask],
        np.asarray(star_coords["yh"], dtype=float)[local_mask],
        np.asarray(star_coords["zh"], dtype=float)[local_mask],
        np.asarray(star_coords["vxh"], dtype=float)[local_mask],
        np.asarray(star_coords["vyh"], dtype=float)[local_mask],
        np.asarray(star_coords["vzh"], dtype=float)[local_mask],
    )
    return {
        "x_kpc": x_kpc,
        "y_kpc": y_kpc,
        "vlos": vlos,
    }


def old_star_projected_kinematics(snapshot):
    """Projected old-star arrays with raw and velocity-gradient-subtracted vlos."""
    old_local_mask = old_dwarf_star_local_mask(snapshot)
    projected = _star_projected_kinematics(snapshot, old_local_mask)
    x_kpc = projected["x_kpc"]
    y_kpc = projected["y_kpc"]
    vlos = projected["vlos"]
    gradient = fit_velocity_gradient(x_kpc, y_kpc, vlos)

    return {
        "old_local_mask": old_local_mask,
        "x_kpc": x_kpc,
        "y_kpc": y_kpc,
        "vlos": vlos,
        "vlos_detrended": gradient["v_resid"],
        "velocity_gradient": gradient,
    }


def detrended_dispersion_in_aperture(
    x_kpc,
    y_kpc,
    vlos,
    radius_kpc,
    center_x_kpc=0.0,
    center_y_kpc=0.0,
    ep=0.0,
    pa=0.0,
    circular=True,
):
    """Fit and remove a planar velocity gradient inside an aperture.

    The default aperture is circular, matching the Walker-style observational
    selection.  Set ``circular=False`` only for an explicit elliptical
    diagnostic.
    """
    x_kpc = np.asarray(x_kpc, dtype=float)
    y_kpc = np.asarray(y_kpc, dtype=float)
    vlos = np.asarray(vlos, dtype=float)
    if x_kpc.shape != y_kpc.shape or x_kpc.shape != vlos.shape:
        raise ValueError("x_kpc, y_kpc, and vlos must have matching shapes")

    if circular:
        radius = np.hypot(x_kpc - center_x_kpc, y_kpc - center_y_kpc)
    else:
        radius = Analysis.projected_elliptical_radius(
            x_kpc,
            y_kpc,
            ep=ep,
            pa=pa,
            center_x=center_x_kpc,
            center_y=center_y_kpc,
        )

    aperture_mask = (
        np.isfinite(radius)
        & np.isfinite(vlos)
        & np.isfinite(radius_kpc)
        & (radius_kpc > 0.0)
        & (radius <= radius_kpc)
    )
    gradient = fit_velocity_gradient(x_kpc, y_kpc, vlos, mask=aperture_mask)
    residual = gradient["v_resid"]
    sigma = _safe_std(residual[aperture_mask])
    return {
        "sigma": sigma,
        "nstar": int(np.count_nonzero(aperture_mask)),
        "gradient": gradient,
    }


def _mean_pm_from_arrays(x, y, z, vx, vy, vz):
    pos = CartesianRepresentation(
        x * u.kpc,
        y * u.kpc,
        z * u.kpc,
    )
    vel = CartesianDifferential(
        vx * u.km / u.s,
        vy * u.km / u.s,
        vz * u.km / u.s,
    )
    coord = SkyCoord(pos.with_differentials(vel), frame="galactocentric")
    coord_icrs = coord.transform_to("icrs")
    return np.mean(coord_icrs.pm_ra_cosdec).value, np.mean(coord_icrs.pm_dec).value


def compute_snapshot_summary(snapshot, numsp):
    df = snapshot["df"]
    tsnap = snapshot["tsnap"]
    mw_gas_mask = snapshot["mw_gas_mask"]
    dw_hot_gas_mask = snapshot["dw_hot_gas_mask"]
    dw_cold_gas_mask = snapshot["dw_cold_gas_mask"]
    total_dw_star_mask = snapshot["total_dw_star_mask"]
    total_mw_star_mask = snapshot["total_mw_star_mask"]
    dw_xc = snapshot["dw_xc"]
    dw_yc = snapshot["dw_yc"]
    dw_zc = snapshot["dw_zc"]
    star_coords = snapshot["star_coords"]
    d_mean = snapshot["d_mean"]
    d_mean_gal = snapshot["d_mean_gal"]
    x_kpc = snapshot["x_kpc"]
    y_kpc = snapshot["y_kpc"]
    hot_gas_x_kpc = snapshot["hot_gas_x_kpc"]
    hot_gas_y_kpc = snapshot["hot_gas_y_kpc"]
    cold_gas_x_kpc = snapshot["cold_gas_x_kpc"]
    cold_gas_y_kpc = snapshot["cold_gas_y_kpc"]
    rotra_dw_cold_gas = snapshot["rotra_dw_cold_gas"]
    rotdec_dw_cold_gas = snapshot["rotdec_dw_cold_gas"]
    old_kinematics = old_star_projected_kinematics(snapshot)
    old_star_local_mask = old_kinematics["old_local_mask"]

    x_all = df["x"].to_numpy()
    y_all = df["y"].to_numpy()
    z_all = df["z"].to_numpy()
    vx_all = df["vx"].to_numpy()
    vy_all = df["vy"].to_numpy()
    vz_all = df["vz"].to_numpy()
    m_all = df["m"].to_numpy()
    r_all = df["r"].to_numpy()
    nh_all = df["nh"].to_numpy()
    tp_all = df["tp"].to_numpy()

    star_x = x_all[total_dw_star_mask]
    star_y = y_all[total_dw_star_mask]
    star_z = z_all[total_dw_star_mask]
    star_vx = vx_all[total_dw_star_mask]
    star_vy = vy_all[total_dw_star_mask]
    star_vz = vz_all[total_dw_star_mask]
    star_m = m_all[total_dw_star_mask]

    obs_star_x = star_x[old_star_local_mask]
    obs_star_y = star_y[old_star_local_mask]
    obs_star_z = star_z[old_star_local_mask]
    obs_star_vx = star_vx[old_star_local_mask]
    obs_star_vy = star_vy[old_star_local_mask]
    obs_star_vz = star_vz[old_star_local_mask]
    obs_star_m = star_m[old_star_local_mask]
    obs_x_kpc = old_kinematics["x_kpc"]
    obs_y_kpc = old_kinematics["y_kpc"]
    hot_x = x_all[dw_hot_gas_mask]
    hot_y = y_all[dw_hot_gas_mask]
    hot_z = z_all[dw_hot_gas_mask]
    hot_m = m_all[dw_hot_gas_mask]
    hot_nh = nh_all[dw_hot_gas_mask]

    cold_x = x_all[dw_cold_gas_mask]
    cold_y = y_all[dw_cold_gas_mask]
    cold_z = z_all[dw_cold_gas_mask]
    cold_m = m_all[dw_cold_gas_mask]
    cold_nh = nh_all[dw_cold_gas_mask]

    mw_star_r = total_mw_star_mask & (r_all <= d_mean_gal)
    mw_gas_r = np.zeros_like(r_all, dtype=bool) if mw_gas_mask is None else (mw_gas_mask & (r_all <= d_mean_gal))
    mw_dm_r = (tp_all == 1) & (r_all <= d_mean_gal)
    mass_mw_star_r = m_all[mw_star_r].sum()
    mass_mw_gas_r = m_all[mw_gas_r].sum()
    mass_mw_dm_r = m_all[mw_dm_r].sum()
    mw_mass_r = mass_mw_star_r + mass_mw_gas_r + mass_mw_dm_r

    projected_shape = Analysis.calculate_projected_shape(
        obs_x_kpc,
        obs_y_kpc,
        mass=obs_star_m,
        n_neighbors=30,
    )
    eps = projected_shape["eps"]
    pa = projected_shape["pa"]
    shape_center_x_kpc = projected_shape["center_x"]
    shape_center_y_kpc = projected_shape["center_y"]
    if not np.isfinite(eps):
        eps = 0.0
    if not np.isfinite(pa):
        pa = 0.0
    if not np.isfinite(shape_center_x_kpc):
        shape_center_x_kpc = 0.0
    if not np.isfinite(shape_center_y_kpc):
        shape_center_y_kpc = 0.0

    r_half = Analysis.half_light_radius(
        obs_x_kpc,
        obs_y_kpc,
        obs_star_m,
        ep=eps,
        pa=pa,
        center_x=shape_center_x_kpc,
        center_y=shape_center_y_kpc,
    )
    r_half_circularized = float(Analysis.circularized_half_light_radius(r_half, eps))
    r_half_circular = Analysis.half_light_radius(
        obs_x_kpc,
        obs_y_kpc,
        obs_star_m,
        ep=0.0,
        pa=0.0,
        center_x=shape_center_x_kpc,
        center_y=shape_center_y_kpc,
    )

    r_half_3d = Analysis.calculate_half_mass_radius(
        obs_star_x,
        obs_star_y,
        obs_star_z,
        obs_star_m,
        center_x=dw_xc,
        center_y=dw_yc,
        center_z=dw_zc,
    )

    dw_star_rhalf_mask = Analysis.get_3d_radial_mask(
        obs_star_x - dw_xc,
        obs_star_y - dw_yc,
        obs_star_z - dw_zc,
        r_half_3d,
    )
    dw_hot_gas_rhalf_mask = Analysis.get_3d_radial_mask(
        hot_x - dw_xc,
        hot_y - dw_yc,
        hot_z - dw_zc,
        r_half_3d,
    )
    dw_cold_gas_rhalf_mask = Analysis.get_3d_radial_mask(
        cold_x - dw_xc,
        cold_y - dw_yc,
        cold_z - dw_zc,
        r_half_3d,
    )

    star_half_mass = obs_star_m[dw_star_rhalf_mask].sum()
    hotgas_half_mass = _sum_hot_gas_mass(hot_m, hot_nh, dw_hot_gas_rhalf_mask)
    coldgas_half_mass = _sum_cold_gas_mass(cold_m, cold_nh, dw_cold_gas_rhalf_mask)

    sigma_re_circular_result = detrended_dispersion_in_aperture(
        old_kinematics["x_kpc"],
        old_kinematics["y_kpc"],
        old_kinematics["vlos"],
        r_half_circular,
        center_x_kpc=shape_center_x_kpc,
        center_y_kpc=shape_center_y_kpc,
        circular=True,
    )
    sigma_re_elliptical_result = detrended_dispersion_in_aperture(
        old_kinematics["x_kpc"],
        old_kinematics["y_kpc"],
        old_kinematics["vlos"],
        r_half,
        center_x_kpc=shape_center_x_kpc,
        center_y_kpc=shape_center_y_kpc,
        ep=eps,
        pa=pa,
        circular=False,
    )
    sigma_fixed_500pc_result = detrended_dispersion_in_aperture(
        old_kinematics["x_kpc"],
        old_kinematics["y_kpc"],
        old_kinematics["vlos"],
        r_3d_cut,
        center_x_kpc=shape_center_x_kpc,
        center_y_kpc=shape_center_y_kpc,
        circular=True,
    )
    sigma_re_circular = sigma_re_circular_result["sigma"]
    sigma_re_elliptical = sigma_re_elliptical_result["sigma"]
    sigma_fixed_500pc = sigma_fixed_500pc_result["sigma"]
    sigma_re_noldstar = sigma_re_circular_result["nstar"]
    sigma_gradient_kms_per_kpc = sigma_re_circular_result["gradient"]["grad_amp"]

    stellar_region_radius = stellar_region_rhalf_multiplier * r_half
    dw_star_stellar_region_mask = Analysis.get_elliptical_radial_mask(
        obs_x_kpc,
        obs_y_kpc,
        stellar_region_radius,
        ep=eps,
        pa=pa,
        center_x=shape_center_x_kpc,
        center_y=shape_center_y_kpc,
    )
    dw_hot_gas_stellar_region_mask = Analysis.get_elliptical_radial_mask(
        hot_gas_x_kpc,
        hot_gas_y_kpc,
        stellar_region_radius,
        ep=eps,
        pa=pa,
        center_x=shape_center_x_kpc,
        center_y=shape_center_y_kpc,
    )
    dw_cold_gas_stellar_region_mask = Analysis.get_elliptical_radial_mask(
        cold_gas_x_kpc,
        cold_gas_y_kpc,
        stellar_region_radius,
        ep=eps,
        pa=pa,
        center_x=shape_center_x_kpc,
        center_y=shape_center_y_kpc,
    )

    v_r, v_theta, v_phi = Analysis.cartesian_to_spherical(
        obs_star_x[dw_star_stellar_region_mask],
        obs_star_y[dw_star_stellar_region_mask],
        obs_star_z[dw_star_stellar_region_mask],
        obs_star_vx[dw_star_stellar_region_mask],
        obs_star_vy[dw_star_stellar_region_mask],
        obs_star_vz[dw_star_stellar_region_mask],
    )

    dw_star_500pc_mask = Analysis.get_3d_radial_mask(
        obs_star_x - dw_xc,
        obs_star_y - dw_yc,
        obs_star_z - dw_zc,
        r_3d_cut,
    )

    cold_gas_r = np.sqrt(np.array(rotra_dw_cold_gas) ** 2 + np.array(rotdec_dw_cold_gas) ** 2)
    cold_gas_within_mask = cold_gas_r <= 2.1
    rotra_dw_cold_gas_within = np.array(rotra_dw_cold_gas)[cold_gas_within_mask]
    rotdec_dw_cold_gas_within = np.array(rotdec_dw_cold_gas)[cold_gas_within_mask]
    if len(rotra_dw_cold_gas_within) > 0:
        cold_gas_center_ra, cold_gas_center_dec = Analysis.find_center_2d(
            rotra_dw_cold_gas_within,
            rotdec_dw_cold_gas_within,
            units='offset',
        )
    else:
        cold_gas_center_ra, cold_gas_center_dec = np.nan, np.nan

    mean_pmra, mean_pmdec = _mean_pm_from_arrays(
        obs_star_x,
        obs_star_y,
        obs_star_z,
        obs_star_vx,
        obs_star_vy,
        obs_star_vz,
    )

    star_mass = obs_star_m[dw_star_stellar_region_mask].sum()
    coldgas_mass = _sum_cold_gas_mass(cold_m, cold_nh, dw_cold_gas_stellar_region_mask)
    hotgas_mass = _sum_hot_gas_mass(hot_m, hot_nh, dw_hot_gas_stellar_region_mask)
    gas_density = Analysis.GetGasdensity(
        cold_gas_x_kpc,
        cold_gas_y_kpc,
        cold_m * cold_nh,
        r_half,
        ep=eps,
        pa=pa,
        center_x=shape_center_x_kpc,
        center_y=shape_center_y_kpc,
    )
    sigma_x = _safe_std(obs_star_vx[dw_star_500pc_mask])
    sigma_y = _safe_std(obs_star_vy[dw_star_500pc_mask])
    sigma_z = _safe_std(obs_star_vz[dw_star_500pc_mask])
    sigma_xyz = np.sqrt((sigma_x**2 + sigma_y**2 + sigma_z**2) / 3)

    baryon_half_mass = star_half_mass + hotgas_half_mass + coldgas_half_mass
    if not np.isfinite(r_half_3d) or r_half_3d <= 0 or baryon_half_mass < 0:
        theoretical_sigma = np.nan
    else:
        theoretical_sigma = np.sqrt(baryon_half_mass / 930.0 / (r_half_3d * 1000.0))

    return {
        "star_mass": star_mass,
        "star_half_mass": star_half_mass,
        "hotgas_mass": hotgas_mass,
        "hotgas_half_mass": hotgas_half_mass,
        "coldgas_mass": coldgas_mass,
        "coldgas_half_mass": coldgas_half_mass,
        "mw_mass_r": mw_mass_r,
        "eps": eps,
        "pa": pa,
        "gas_density": gas_density,
        "vr": _safe_mean(v_r),
        "vtheta": _safe_mean(v_theta),
        "vphi": _safe_mean(v_phi),
        "distance_gal": d_mean_gal,
        "rhalf": r_half,
        "rhalf_circularized": r_half_circularized,
        "rhalf_circular": r_half_circular,
        "shape_center_x_kpc": shape_center_x_kpc,
        "shape_center_y_kpc": shape_center_y_kpc,
        "distance": d_mean,
        "age": tsnap,
        # ``sigma`` is the Walker-style observable: old stars inside
        # the directly measured circular half-light aperture, after fitting
        # and removing a planar LOS velocity gradient on that same sample.
        "sigma": sigma_re_circular,
        "sigma_re_circular": sigma_re_circular,
        "sigma_re_elliptical": sigma_re_elliptical,
        "sigma_fixed_500pc": sigma_fixed_500pc,
        "sigma_re_noldstar": sigma_re_noldstar,
        "sigma_gradient_kms_per_kpc": sigma_gradient_kms_per_kpc,
        "tsigma": theoretical_sigma,
        "tsigma_xyz": sigma_xyz,
        "sigma_x": sigma_x,
        "sigma_y": sigma_y,
        "sigma_z": sigma_z,
        "pmra": mean_pmra,
        "pmdec": mean_pmdec,
        "cold_gas_center_ra": cold_gas_center_ra,
        "cold_gas_center_dec": cold_gas_center_dec,
        "numsp": numsp,
    }
