import numpy as np

from basefunc import Analysis, GalaxySimulation
from canonical_observer import production_sky, production_motion, walker61_zero_error, CANONICAL_OBSERVER_ID, MEMBER_SELECTION_ID, observer_contract_sha256, write_snapshot_observer_metadata


def prepare_snapshot_context(
    folder_path,
    snapshot_num,
    core_radius,
    r_exclude=5,
    dwarf_radius_factor=3.0,
    k_density=16,
    dwarf_gas_radius=15.0,
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

    # Current spatial membership is authoritative. A source ID or formation
    # time never adds/removes a stellar row from this physical selection.
    dw_star_mask = np.asarray(dw_star_mask, dtype=bool)
    mw_star_mask = np.asarray(mw_star_mask, dtype=bool)
    particle_types_array = df['tp'].to_numpy()
    old_star_mask = dw_star_mask & np.isin(particle_types_array, [2, 3])
    star_df = df.loc[dw_star_mask]
    dw_xc, dw_yc, dw_zc = simulation.dwarf_center_grf_kpc
    dwarf_center = np.array([dw_xc, dw_yc, dw_zc], dtype=float)
    # Keep strict aperture boundaries at the precision of the measured centre;
    # only the radius work vector is double precision, not the full particle data.
    radius_squared = np.zeros(len(df), dtype=np.float64)
    for axis, center_component in zip(('x', 'y', 'z'), dwarf_center):
        offset = df[axis].to_numpy(dtype=np.float64) - center_component
        radius_squared += offset * offset
    gas_rows = particle_types_array == 0
    star_rows = np.isin(particle_types_array, [2, 3, 4])
    gas_masks = {radius: gas_rows & (radius_squared < radius ** 2) for radius in (5, 15)}
    star_masks = {radius: star_rows & (radius_squared < radius ** 2) for radius in (5, 15)}
    # A caller's former 30-kpc radius is a preselection only. Final gas,
    # temperature phases and matched-aperture fractions use at most r<15 kpc.
    preselection_radius = max(15.0, float(dwarf_gas_radius))
    dw_gas_preselection_mask = gas_rows & (radius_squared < preselection_radius ** 2)
    dw_gas_mask = gas_masks[15]
    membership_state = dict(
        member_selection_id=MEMBER_SELECTION_ID,
        stellar_definition='stellar properties and current density-centre spatial aperture',
        stellar_radius_kpc=float(dwarf_radius_factor * core_radius),
        light_tracer_particle_types=[2, 3], selected_light_tracer_N=int(old_star_mask.sum()),
        center_grf_kpc=dwarf_center.tolist(),
        center_method=simulation.dwarf_center_method,
        gas_definition='all TP0 rows in strict 3D apertures, without origin-ID filtering',
        gas_apertures_kpc=[5, 15], gas_outer_annulus_kpc=[5, 15],
        gas_preselection_radius_kpc=preselection_radius,
        gas_temperature_split_K=float(gas_temperature_split),
        velocity_reference='native MW-centred GRF after the same snapshot host bulk subtraction',
        stellar_age=getattr(simulation, 'stellar_age_metadata', []),
        catalogue_alignment='analysis_row_index; source_particle_id is optional metadata only',
        source_id_availability=getattr(simulation, 'source_id_metadata', []))
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

    # Image centre is the pre-existing stellar mass-centre; newly formed stars
    # and gas share this same physical pointing but retain their own six-vectors.
    old=old_star_mask[total_dw_star_mask]
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
    # Kinematics include every current dwarf member, including newly formed stars.
    fit_R5=sky['fit_R5_mask']
    walker=walker61_zero_error(sky,star_coords['vlos'],fit_mask=fit_R5)
    motion=production_motion(star_coords,sky,star_mass,fit_R5)
    observer_metadata_path=write_snapshot_observer_metadata(folder_path,snapshot_num,
        simulation.observer_host_state,membership_state,star_coords,walker,motion)
    cx, cy, cz = dw_xc, dw_yc, dw_zc

    return {
        "simulation": simulation,
        "canonical_observer_id": CANONICAL_OBSERVER_ID,
        "member_selection_id": MEMBER_SELECTION_ID,
        "observer_config_sha256": observer_contract_sha256(),
        "observer_host_state": simulation.observer_host_state,
        "walker_R5": walker,
        "motion_R5": motion,
        "fit_R5_mask": fit_R5,
        "old_star_local_mask": old,
        "membership_state": membership_state,
        "observer_metadata_path": observer_metadata_path,
        "df": df,
        "tsnap": tsnap,
        "mw_star_mask": mw_star_mask,
        "dw_star_mask": dw_star_mask,
        "dw_gas_mask": dw_gas_mask,
        "dw_gas_mask_3d5": gas_masks[5],
        "dw_gas_mask_3d15": gas_masks[15],
        "dw_gas_mask_3d5to15": gas_masks[15] & ~gas_masks[5],
        "dw_star_mask_3d5": star_masks[5],
        "dw_star_mask_3d15": star_masks[15],
        "dw_gas_preselection_mask": dw_gas_preselection_mask,
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
