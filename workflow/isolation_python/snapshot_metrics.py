import numpy as np
import astropy.units as u
from astropy.coordinates import CartesianDifferential, CartesianRepresentation, SkyCoord

from basefunc import Analysis
from variable import r_3d_cut


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

    eps, pa = Analysis.calculate_ellipticity(
        x_kpc,
        y_kpc,
        mass=star_m,
        n_neighbors=30,
    )

    rotate = np.array(
        [
            [np.cos(-pa), -np.sin(-pa)],
            [np.sin(-pa), np.cos(-pa)],
        ]
    )
    coords = np.vstack((x_kpc, y_kpc))
    rotated_coords = rotate @ coords

    x_rot, y_rot = rotated_coords[0], rotated_coords[1]
    r_half = Analysis.half_light_radius(x_rot, y_rot, star_m, ep=eps)

    r_half_3d = Analysis.calculate_half_mass_radius(
        star_x,
        star_y,
        star_z,
        star_m,
        center_x=dw_xc,
        center_y=dw_yc,
        center_z=dw_zc,
    )

    dw_star_rhalf_mask = Analysis.get_3d_radial_mask(
        star_x - dw_xc,
        star_y - dw_yc,
        star_z - dw_zc,
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

    star_half_mass = star_m[dw_star_rhalf_mask].sum()
    hotgas_half_mass = _sum_hot_gas_mass(hot_m, hot_nh, dw_hot_gas_rhalf_mask)
    coldgas_half_mass = _sum_cold_gas_mass(cold_m, cold_nh, dw_cold_gas_rhalf_mask)

    vlos_array = Analysis.cal_vlos(
        star_coords["xh"],
        star_coords["yh"],
        star_coords["zh"],
        star_coords["vxh"],
        star_coords["vyh"],
        star_coords["vzh"],
    )
    r_projected = np.sqrt(x_kpc**2 + y_kpc**2)
    sigma_los = _safe_std(vlos_array[r_projected < r_3d_cut])

    dw_star_4rhalf_mask = Analysis.get_radial_mask(x_kpc, y_kpc, 4 * r_half)
    dw_hot_gas_4rhalf_mask = Analysis.get_radial_mask(hot_gas_x_kpc, hot_gas_y_kpc, 4 * r_half)
    dw_cold_gas_4rhalf_mask = Analysis.get_radial_mask(cold_gas_x_kpc, cold_gas_y_kpc, 4 * r_half)

    v_r, v_theta, v_phi = Analysis.cartesian_to_spherical(
        star_x[dw_star_4rhalf_mask],
        star_y[dw_star_4rhalf_mask],
        star_z[dw_star_4rhalf_mask],
        star_vx[dw_star_4rhalf_mask],
        star_vy[dw_star_4rhalf_mask],
        star_vz[dw_star_4rhalf_mask],
    )

    dw_star_500pc_mask = Analysis.get_3d_radial_mask(
        star_x - dw_xc,
        star_y - dw_yc,
        star_z - dw_zc,
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

    mean_pmra, mean_pmdec = _mean_pm_from_arrays(star_x, star_y, star_z, star_vx, star_vy, star_vz)

    star_mass = star_m[dw_star_4rhalf_mask].sum()
    coldgas_mass = _sum_cold_gas_mass(cold_m, cold_nh, dw_cold_gas_4rhalf_mask)
    hotgas_mass = _sum_hot_gas_mass(hot_m, hot_nh, dw_hot_gas_4rhalf_mask)
    gas_density = Analysis.GetGasdensity(cold_gas_x_kpc, cold_gas_y_kpc, cold_m * cold_nh, r_half)
    sigma_x = _safe_std(star_vx[dw_star_500pc_mask])
    sigma_y = _safe_std(star_vy[dw_star_500pc_mask])
    sigma_z = _safe_std(star_vz[dw_star_500pc_mask])
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
        "distance": d_mean,
        "age": tsnap,
        "sigma": sigma_los,
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
