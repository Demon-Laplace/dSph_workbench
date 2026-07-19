import numpy as np

from basefunc import Analysis, GalaxySimulation


def prepare_snapshot_context(
    folder_path,
    snapshot_num,
    core_radius,
    r_exclude=5,
    dwarf_radius_factor=3.0,
    k_density=16,
    dwarf_gas_radius=30.0,
    gas_temperature_split=20000.0,
    include_mw_gas=False,
    mw_gas_radius=300.0,
    include_dark_matter=False,
    include_star_birth=False,
):
    simulation = GalaxySimulation(folder_path=folder_path, snapshot_num=snapshot_num)
    particle_types = (0, 1, 2, 3, 4) if include_dark_matter else (0, 2, 3, 4)
    dftype = 3 if include_star_birth else 2
    star_fields = ('birth',) if include_star_birth else ()
    df, tsnap = simulation.load_snapshot(
        dftype=dftype,
        gas_fields=('temp', 'nh'),
        star_fields=star_fields,
        particle_types=particle_types,
    )
    df = simulation.center_on_MW()

    mw_star_mask, dw_star_mask = simulation.classify_mw_dwarf(
        df=df,
        r_exclude=r_exclude,
        dwarf_radius=dwarf_radius_factor * core_radius,
        k_density=k_density,
    )

    dw_gas_mask = simulation.find_dwarf_gas(dw_star_mask=dw_star_mask, radius=dwarf_gas_radius)
    mw_gas_mask = None
    if include_mw_gas:
        mw_gas_mask = simulation.find_MW_gas(radius=mw_gas_radius)

    dw_hot_gas_mask = dw_gas_mask & (df.temp >= gas_temperature_split)
    dw_cold_gas_mask = dw_gas_mask & (df.temp < gas_temperature_split)
    total_dw_star_mask = dw_star_mask
    total_mw_star_mask = mw_star_mask

    star_coords = simulation.convert_coordinates_for_mask(total_dw_star_mask, df=df, with_velocity=True, include_galactic=False)
    hot_gas_coords = simulation.convert_coordinates_for_mask(dw_hot_gas_mask, df=df, with_velocity=False, include_galactic=False)
    cold_gas_coords = simulation.convert_coordinates_for_mask(dw_cold_gas_mask, df=df, with_velocity=False, include_galactic=False)

    star_df = df.loc[total_dw_star_mask]
    dw_cra, dw_cdec = Analysis.find_center_2d(
        star_coords["ra"],
        star_coords["dec"],
        units="degree",
    )
    rotra_dw_star, rotdec_dw_star = Analysis.rotate_to_sky(
        star_coords["ra"],
        star_coords["dec"],
        center=[dw_cra, dw_cdec],
    )
    rotra_dw_hot_gas, rotdec_dw_hot_gas = Analysis.rotate_to_sky(
        hot_gas_coords["ra"],
        hot_gas_coords["dec"],
        center=[dw_cra, dw_cdec],
    )
    rotra_dw_cold_gas, rotdec_dw_cold_gas = Analysis.rotate_to_sky(
        cold_gas_coords["ra"],
        cold_gas_coords["dec"],
        center=[dw_cra, dw_cdec],
    )

    d_mean = np.mean(star_coords["rh"])
    d_mean_gal = np.mean(star_df["r"])

    x_kpc = np.asarray(np.radians(rotra_dw_star) * d_mean)
    y_kpc = np.asarray(np.radians(rotdec_dw_star) * d_mean)
    hot_gas_x_kpc = np.asarray(np.radians(rotra_dw_hot_gas) * d_mean)
    hot_gas_y_kpc = np.asarray(np.radians(rotdec_dw_hot_gas) * d_mean)
    cold_gas_x_kpc = np.asarray(np.radians(rotra_dw_cold_gas) * d_mean)
    cold_gas_y_kpc = np.asarray(np.radians(rotdec_dw_cold_gas) * d_mean)

    dw_xc, dw_yc, dw_zc = Analysis.find_center_3d(
        star_df["x"],
        star_df["y"],
        star_df["z"],
        mass=star_df["m"],
    )
    cx, cy, cz = dw_xc, dw_yc, dw_zc

    return {
        "simulation": simulation,
        "df": df,
        "tsnap": tsnap,
        "mw_star_mask": mw_star_mask,
        "dw_star_mask": dw_star_mask,
        "dw_gas_mask": dw_gas_mask,
        "mw_gas_mask": mw_gas_mask,
        "dw_hot_gas_mask": dw_hot_gas_mask,
        "dw_cold_gas_mask": dw_cold_gas_mask,
        "total_dw_star_mask": total_dw_star_mask,
        "total_mw_star_mask": total_mw_star_mask,
        "dw_cra": dw_cra,
        "dw_cdec": dw_cdec,
        "dw_xc": dw_xc,
        "dw_yc": dw_yc,
        "dw_zc": dw_zc,
        "cx": cx,
        "cy": cy,
        "cz": cz,
        "star_coords": star_coords,
        "hot_gas_coords": hot_gas_coords,
        "cold_gas_coords": cold_gas_coords,
        "rotra_dw_star": rotra_dw_star,
        "rotdec_dw_star": rotdec_dw_star,
        "rotra_dw_hot_gas": rotra_dw_hot_gas,
        "rotdec_dw_hot_gas": rotdec_dw_hot_gas,
        "rotra_dw_cold_gas": rotra_dw_cold_gas,
        "rotdec_dw_cold_gas": rotdec_dw_cold_gas,
        "d_mean": d_mean,
        "d_mean_gal": d_mean_gal,
        "x_kpc": x_kpc,
        "y_kpc": y_kpc,
        "hot_gas_x_kpc": hot_gas_x_kpc,
        "hot_gas_y_kpc": hot_gas_y_kpc,
        "cold_gas_x_kpc": cold_gas_x_kpc,
        "cold_gas_y_kpc": cold_gas_y_kpc,
    }
