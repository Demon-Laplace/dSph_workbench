"""Angular display contract; independent of the 5-kpc measurement aperture."""
import numpy as np

DISPLAY_CONTRACT_ID = 'heliocentric_delta_radec_2p1_v1'

def angular_velocity_grid40(projected, vlos_kms, *, display_v0_kms):
    ids=np.asarray(projected['particle_id']);v=np.asarray(vlos_kms,dtype=float)
    x=np.asarray(projected['delta_RA_deg'],dtype=float);y=np.asarray(projected['delta_Dec_deg'],dtype=float)
    if v.shape!=ids.shape or x.shape!=ids.shape or y.shape!=ids.shape or not np.all(np.isfinite(v)):
        raise ValueError('one finite velocity and sky position per ID are required')
    edges=np.linspace(-2.1,2.1,41)
    mask=np.asarray(projected['front_mask'],dtype=bool)&(abs(x)<=2.1)&(abs(y)<=2.1)
    count=np.histogram2d(x[mask],y[mask],bins=(edges,edges))[0].astype(np.int64)
    total=np.histogram2d(x[mask],y[mask],bins=(edges,edges),weights=v[mask]-display_v0_kms)[0]
    mean=np.divide(total,count,out=np.full((40,40),np.nan),where=count>0)
    xx,yy=np.meshgrid(edges,edges,indexing='xy')
    return dict(count=count,mean_HRF_minus_v0=mean,delta_RA_vertices_deg=xx,delta_Dec_vertices_deg=yy,
        edges_delta_RA_deg=edges.copy(),edges_delta_Dec_deg=edges.copy(),display_mask=mask,
        display_particle_id=ids[mask].copy(),display_v0_kms=float(display_v0_kms),
        display_color_limits_kms=(-12.,12.),display_contract_id=DISPLAY_CONTRACT_ID,
        display_coordinate_system='raw_delta_RA_delta_Dec',field_half_width_deg=2.1,
        grid_array_order='east_bin,north_bin; transpose once for XY vertex mesh')

def center_motion_arrow(projected,unit_EN,*,length_deg=.65):
    unit=np.asarray(unit_EN,dtype=float)
    if unit.shape!=(2,) or not np.all(np.isfinite(unit)) or np.linalg.norm(unit)<=0:
        raise ValueError('a nonzero finite east/north direction is required')
    unit=unit/np.linalg.norm(unit)
    a,d=np.deg2rad([projected['center_RA_deg'],projected['center_Dec_deg']])
    n=np.array([np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)])
    e=np.array([-np.sin(a),np.cos(a),0.]);north=np.array([-np.sin(d)*np.cos(a),-np.sin(d)*np.sin(a),np.cos(d)])
    psi=np.deg2rad(length_deg);tip=n*np.cos(psi)+(unit[0]*e+unit[1]*north)*np.sin(psi)
    ra=np.degrees(np.arctan2(tip[1],tip[0]));dec=np.degrees(np.arctan2(tip[2],np.hypot(tip[0],tip[1])))
    return dict(origin_delta_RA_Dec_deg=np.array([0.,0.]),
        tip_delta_RA_Dec_deg=np.array([(ra-np.degrees(a)+180)%360-180,dec-np.degrees(d)]),
        angular_length_deg=float(length_deg),length_represents='direction_only')
