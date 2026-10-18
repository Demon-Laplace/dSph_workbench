import numpy as np

from basefunc import Analysis, GalaxySimulation
from canonical_observer import production_sky, production_motion, walker61_zero_error, CANONICAL_OBSERVER_ID, observer_contract_sha256, frozen_cohort_mask, write_snapshot_observer_metadata


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
        box=False,
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

    fixed_old_mask,cohort_state=frozen_cohort_mask(df['id'].to_numpy(),folder_path)
    # Frozen tracers define observed light/kinematics. Newly formed stars remain
    # a separately classified physical mass component, never relabelled as the cohort.
    if 'birth' in df:
        newly_formed=df['birth'].to_numpy()!=0
    else:
        newly_formed=df['tp'].to_numpy()==4
    dynamic_new_mask=np.asarray(dw_star_mask,dtype=bool)&newly_formed&~fixed_old_mask
    dw_star_mask=fixed_old_mask|dynamic_new_mask
    mw_star_mask=np.asarray(mw_star_mask,dtype=bool)&~dw_star_mask

    dw_gas_mask = simulation.find_dwarf_gas(dw_star_mask=dw_star_mask, radius=dwarf_gas_radius)
    mw_gas_mask = None
    if include_mw_gas:
        mw_gas_mask = simulation.find_MW_gas(radius=mw_gas_radius)

    dw_hot_gas_mask = dw_gas_mask & (df.temp >= gas_temperature_split)
    dw_cold_gas_mask = dw_gas_mask & (df.temp < gas_temperature_split)
    total_dw_star_mask = dw_star_mask
    total_mw_star_mask = mw_star_mask

    star_coords = simulation.convert_coordinates_for_mask(total_dw_star_mask, df=df, with_velocity=True, include_galactic=False)
    hot_gas_coords = simulation.convert_coordinates_for_mask(dw_hot_gas_mask, df=df, with_velocity=True, include_galactic=False)
    cold_gas_coords = simulation.convert_coordinates_for_mask(dw_cold_gas_mask, df=df, with_velocity=True, include_galactic=False)

    star_df = df.loc[total_dw_star_mask]
    # Image centre is the pre-existing stellar mass-centre; newly formed stars
    # and gas share this same physical pointing but retain their own six-vectors.
    old=fixed_old_mask[total_dw_star_mask]
    if not old.any():raise ValueError('No old stellar tracers for canonical image centre')
    star_mass=star_df['m'].to_numpy(dtype=float)
    center_h=np.average(star_coords['_native']['h_native_kpc'][old],axis=0,weights=star_mass[old])
    for coords in (star_coords,hot_gas_coords,cold_gas_coords):
        coords['_sky']=production_sky(coords,center_h)
    sky=star_coords['_sky']
    dw_cra,dw_cdec=sky['center_RA_deg'],sky['center_Dec_deg']
    rotra_dw_star,rotdec_dw_star=sky['delta_RA_deg'],sky['delta_Dec_deg']
    rotra_dw_hot_gas=hot_gas_coords['_sky']['delta_RA_deg']
    rotdec_dw_hot_gas=hot_gas_coords['_sky']['delta_Dec_deg']
    rotra_dw_cold_gas=cold_gas_coords['_sky']['delta_RA_deg']
    rotdec_dw_cold_gas=cold_gas_coords['_sky']['delta_Dec_deg']
    d_mean=sky['center_distance_kpc']
    d_mean_gal=np.linalg.norm(np.average(star_df.loc[old,['x','y','z']].to_numpy(dtype=float),axis=0,weights=star_mass[old]))
    # Physical profiles retain true gnomonic east/north, rather than treating raw dRA as east distance.
    x_kpc,y_kpc=sky['east_kpc'],sky['north_kpc']
    hot_gas_x_kpc,hot_gas_y_kpc=hot_gas_coords['_sky']['east_kpc'],hot_gas_coords['_sky']['north_kpc']
    cold_gas_x_kpc,cold_gas_y_kpc=cold_gas_coords['_sky']['east_kpc'],cold_gas_coords['_sky']['north_kpc']
    fit_R5=old & sky['fit_R5_mask']
    walker=walker61_zero_error(sky,star_coords['vlos'],fit_mask=fit_R5)
    motion=production_motion(star_coords,sky,star_mass,fit_R5)
    observer_metadata_path=write_snapshot_observer_metadata(folder_path,snapshot_num,
        simulation.observer_host_state,cohort_state,star_coords,walker,motion)

    dw_xc, dw_yc, dw_zc = Analysis.find_center_3d(
        star_df["x"],
        star_df["y"],
        star_df["z"],
        mass=star_df["m"],
    )
    cx, cy, cz = dw_xc, dw_yc, dw_zc

    return {
        "simulation": simulation,
        "canonical_observer_id": CANONICAL_OBSERVER_ID,
        "observer_config_sha256": observer_contract_sha256(),
        "observer_host_state": simulation.observer_host_state,
        "walker_R5": walker,
        "motion_R5": motion,
        "fit_R5_mask": fit_R5,
        "fixed_old_local_mask": old,
        "cohort_state": cohort_state,
        "observer_metadata_path": observer_metadata_path,
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
