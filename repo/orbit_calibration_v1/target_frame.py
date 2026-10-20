"""One explicit six-dimensional target conversion, without physical registration.

Li (2021) Table 1 contains heliocentric observables; its Table 2 contains
Galactocentric spherical summaries in a different axis/solar convention.
Only the former are input observables here. No active model rotation occurs.
"""
from copy import deepcopy
import numpy as np

INPUT_FRAME = 'icrs_heliocentric'
NATIVE_FRAME = 'mw_relative_native_astropy_gc'
FRAME = dict(galcen_distance_kpc=8.122, z_sun_kpc=.0208,
             galcen_v_sun_kms=[12.9,245.6,7.78],
             galcen_icrs_deg=[266.4051,-28.936175],roll_deg=0.,
             native_to_GC=[[1,0,0],[0,1,0],[0,0,1]])


def explicit_frame():
    from astropy import units as u
    from astropy.coordinates import Galactocentric,ICRS,CartesianDifferential
    f=FRAME
    return Galactocentric(galcen_distance=f['galcen_distance_kpc']*u.kpc,
        z_sun=f['z_sun_kpc']*u.kpc,
        galcen_v_sun=CartesianDifferential(np.array(f['galcen_v_sun_kms'])*u.km/u.s),
        galcen_coord=ICRS(ra=f['galcen_icrs_deg'][0]*u.deg,dec=f['galcen_icrs_deg'][1]*u.deg),
        roll=f['roll_deg']*u.deg)


def resolve_target(target):
    """Validate semantics, invert full ICRS 6D once, and retain a round trip."""
    from astropy import units as u
    from astropy.coordinates import SkyCoord,ICRS
    if target.get('input_frame') != INPUT_FRAME:
        raise ValueError('target.input_frame must be icrs_heliocentric; Li Table 2 GC velocities are not heliocentric inputs.')
    keys=['ra_deg','dec_deg','D_kpc','vlos_kms','mu_alpha_cosdec_masyr','mu_delta_masyr']
    if any(k not in target or not np.isfinite(target[k]) for k in keys):
        raise ValueError('Supply a finite, complete ICRS heliocentric six-vector.')
    if not (0<=target['ra_deg']<360 and -90<=target['dec_deg']<=90 and target['D_kpc']>0):
        raise ValueError('Invalid target position/distance.')
    if any(k in target for k in ['vtheta_kms','vphi_kms','theta_deg','phi_deg']):
        raise ValueError('A GC spherical target cannot be mixed into an ICRS observable six-vector.')
    sky=SkyCoord(ra=target['ra_deg']*u.deg,dec=target['dec_deg']*u.deg,
        distance=target['D_kpc']*u.kpc,radial_velocity=target['vlos_kms']*u.km/u.s,
        pm_ra_cosdec=target['mu_alpha_cosdec_masyr']*u.mas/u.yr,
        pm_dec=target['mu_delta_masyr']*u.mas/u.yr,frame=ICRS())
    point=sky.transform_to(explicit_frame())
    r=point.cartesian.xyz.to_value(u.kpc);v=point.velocity.d_xyz.to_value(u.km/u.s)
    derived=dict(r_grf_native_kpc=r.tolist(),v_grf_native_kms=v.tolist(),
        R_gc_kpc=float(np.linalg.norm(r)),V3D_kms=float(np.linalg.norm(v)),
        Vr_gc_kms=float(r@v/np.linalg.norm(r)),
        Vt_gc_kms=float(np.linalg.norm(np.cross(r,v))/np.linalg.norm(r)))
    for key,value in derived.items():
        if key in target and not np.allclose(target[key],value,atol=1e-8,rtol=0):
            raise ValueError('Supplied target '+key+' disagrees with its ICRS six-vector and fixed solar convention.')
    back=point.transform_to(ICRS())
    restored=dict(ra_deg=float(back.ra.deg),dec_deg=float(back.dec.deg),D_kpc=float(back.distance.to_value(u.kpc)),
        vlos_kms=float(back.radial_velocity.to_value(u.km/u.s)),
        mu_alpha_cosdec_masyr=float(back.pm_ra_cosdec.to_value(u.mas/u.yr)),
        mu_delta_masyr=float(back.pm_dec.to_value(u.mas/u.yr)))
    errors={k:restored[k]-target[k] for k in keys}
    errors['ra_deg']=(errors['ra_deg']+180)%360-180
    if max(abs(x) for x in errors.values())>1e-8:
        raise ValueError('ICRS/GC six-dimensional round trip failed.')
    return dict(input_frame=INPUT_FRAME,output_frame=NATIVE_FRAME,explicit_frame=deepcopy(FRAME),
        native_state=derived,restored_observables=restored,roundtrip_error=errors,
        physical_model_rotation_applied=False,
        note='Derived GC numbers use the current solar convention, not the Li 2021 Table 2 solar convention.')
