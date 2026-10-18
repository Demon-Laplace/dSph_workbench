"""Independent explicit Astropy + direct spherical reference for focused checks."""
from __future__ import annotations
import numpy as np
import astropy.units as un
from astropy.coordinates import Galactocentric, ICRS, SkyCoord, CartesianDifferential

OBSERVER_ID='astropy_native_direct_v1'
DISTANCE=8.122
Z_SUN=.0208
V_SUN=np.array([12.9,245.6,7.78])
S_SUN=np.array([-np.sqrt(DISTANCE**2-Z_SUN**2),0.,Z_SUN])
KM_PER_YEAR_PER_MASYR_KPC=4.740470463533349

def explicit_frame():
    return Galactocentric(galcen_distance=DISTANCE*un.kpc,
        z_sun=Z_SUN*un.kpc,galcen_v_sun=CartesianDifferential(V_SUN*un.km/un.s),
        galcen_coord=SkyCoord(ra=266.4051*un.deg,dec=-28.936175*un.deg,frame='icrs'),
        roll=0*un.deg)

def expected_sixd(xyz,vel):
    p=np.asarray(xyz,dtype=float); v=np.asarray(vel,dtype=float)
    sky=SkyCoord(x=p[:,0]*un.kpc,y=p[:,1]*un.kpc,z=p[:,2]*un.kpc,
        v_x=v[:,0]*un.km/un.s,v_y=v[:,1]*un.km/un.s,v_z=v[:,2]*un.km/un.s,
        frame=explicit_frame()).transform_to(ICRS())
    h=p-S_SUN; u=v-V_SUN; d=np.linalg.norm(h,axis=1)
    direct=np.sum(h*u,axis=1)/d
    return dict(ra=sky.ra.deg,dec=sky.dec.deg,rh=sky.distance.to_value(un.kpc),
                vlos=sky.radial_velocity.to_value(un.km/un.s),
                pmra=sky.pm_ra_cosdec.to_value(un.mas/un.yr),
                pmdec=sky.pm_dec.to_value(un.mas/un.yr),
                direct_vlos=direct,h=h,u=u)

def true_spherical_offsets(ra,dec,ra0,dec0,d0):
    a,d=np.deg2rad(ra),np.deg2rad(dec)
    a0,d0rad=np.deg2rad(ra0),np.deg2rad(dec0)
    rays=np.column_stack([np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)])
    n0=np.array([np.cos(d0rad)*np.cos(a0),np.cos(d0rad)*np.sin(a0),np.sin(d0rad)])
    e0=np.array([-np.sin(a0),np.cos(a0),0.])
    north0=np.cross(n0,e0)
    east=rays@e0; north=rays@north0; front=rays@n0
    sinrho=np.hypot(east,north)
    rho=np.arctan2(sinrho,front)
    factor=np.divide(np.rad2deg(rho),sinrho,out=np.ones(len(rays))*180/np.pi,where=sinrho>0)
    return dict(east_deg=east*factor,north_deg=north*factor,rho_deg=np.rad2deg(rho),
                projected_radius_kpc=d0*np.tan(rho),
                fit_R5_mask=(front>0)&(d0*np.tan(rho)<5),
                east_icrs=e0,north_icrs=north0,center_ray_icrs=n0)

def independent_walker(ra,dec,vlos,*,ra0,dec0,d0,mask=None):
    coords=true_spherical_offsets(ra,dec,ra0,dec0,d0)
    use=coords['fit_R5_mask'] if mask is None else np.asarray(mask,dtype=bool)
    v=np.asarray(vlos,dtype=float)[use]
    ee,nn=coords['east_deg'][use],coords['north_deg'][use]
    ks=[]; intercepts=[]; sigmas=[]
    for theta in np.arange(0.,181.,3.):
        angle=np.deg2rad(theta)
        sine,cosine=np.sin(angle),np.cos(angle)
        if theta==0: sine,cosine=0.,1.
        if theta==90: sine,cosine=1.,0.
        if theta==180: sine,cosine=0.,-1.
        x=ee*sine+nn*cosine
        X=np.column_stack([np.ones(len(x)),x])
        b=np.linalg.lstsq(X,v,rcond=None)[0]
        if np.ptp(x)<1e-13:
            ks.append(np.nan);intercepts.append(np.nan);sigmas.append(np.nan)
            continue
        ks.append(b[1]);intercepts.append(b[0]);sigmas.append(np.sqrt(np.mean((v-X@b)**2)))
    ks=np.array(ks); index=int(np.nanargmax(np.abs(ks)))
    return dict(k_scan=ks,v0_scan=np.array(intercepts),sigma_scan=np.array(sigmas),
                index=index,k_signed=ks[index],k_abs=abs(ks[index]),v0=intercepts[index],
                sigma=sigmas[index],sigma_raw=np.std(v,ddof=0),N=len(v))

def center_pm_from_h_and_u(h_native,u_native):
    # Obtain the rotation only from the explicit GC->ICRS transform, using differences.
    frame=explicit_frame()
    origin=SkyCoord(x=0*un.kpc,y=0*un.kpc,z=0*un.kpc,frame=frame).transform_to(ICRS()).cartesian.xyz.to_value(un.kpc)
    cols=[]
    for basis in np.eye(3):
        point=SkyCoord(x=basis[0]*un.kpc,y=basis[1]*un.kpc,z=basis[2]*un.kpc,frame=frame).transform_to(ICRS()).cartesian.xyz.to_value(un.kpc)
        cols.append(point-origin)
    M=np.column_stack(cols)
    h=M@np.asarray(h_native);u=M@np.asarray(u_native)
    distance=np.linalg.norm(h);n=h/distance
    ra=np.arctan2(n[1],n[0]);dec=np.arcsin(n[2])
    east=np.array([-np.sin(ra),np.cos(ra),0.]);north=np.array([-np.sin(dec)*np.cos(ra),-np.sin(dec)*np.sin(ra),np.cos(dec)])
    return np.array([u@east,u@north])/(KM_PER_YEAR_PER_MASYR_KPC*distance)
