"""Sky display contracts, independent of the observer and measurement aperture.

The raw RA/Dec helpers remain available for historical releases.  New renderers
use the explicitly named TAN helpers below; none of them redefines a particle's
line-of-sight velocity, physical selection, image centre, or motion direction.
"""
import numpy as np

# Kept for callers that reproduce the historical raw RA/Dec display.
DISPLAY_CONTRACT_ID = 'heliocentric_delta_radec_2p1_v1'
LEGACY_DISPLAY_CONTRACT_ID = 'heliocentric_delta_radec_2p1_v1'
TAN_DISPLAY_CONTRACT_ID = 'heliocentric_tan_xi_eta_2p1_v1'
CURRENT_DISPLAY_CONTRACT_ID = 'heliocentric_tan_xi_eta_2p1_v1'
DISPLAY_RELATIVE_DIR = 'display/heliocentric_tan_xi_eta_2p1_v1'

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


def _finite_scalar(value, name, *, positive=False):
    if np.ndim(value) != 0:
        raise ValueError(f'{name} must be a finite scalar')
    value = float(value)
    if not np.isfinite(value) or (positive and value <= 0):
        raise ValueError(f'{name} must be finite' + (' and positive' if positive else ''))
    return value


def tangent_kpc_to_display(projected, x_kpc, y_kpc):
    """Convert physical gnomonic E/N distances to WCS-TAN plane degrees.

    The distances are D0 times the ray's E/N ratios, not arc lengths.  This is
    therefore a linear change of units, without atan or a second cos(dec).
    """
    distance = _finite_scalar(projected['center_distance_kpc'],
                              'center_distance_kpc', positive=True)
    x, y = np.broadcast_arrays(np.asarray(x_kpc, dtype=float),
                               np.asarray(y_kpc, dtype=float))
    return np.rad2deg(x / distance), np.rad2deg(y / distance)


def tangent_coordinates(projected):
    """Return the existing catalogue rays in the common centre's TAN plane."""
    x = np.asarray(projected['east_kpc'], dtype=float)
    y = np.asarray(projected['north_kpc'], dtype=float)
    if x.shape != y.shape:
        raise ValueError('east_kpc and north_kpc must have the same shape')
    return tangent_kpc_to_display(projected, x, y)


def tangent_velocity_grid40(projected, vlos_kms, *, display_v0_kms):
    """Rebin frozen per-ID LOS in the +/-2.1-degree square TAN display.

    Every occupied bin has a number-weighted mean.  Returned value/count arrays
    are [xi_bin, eta_bin]; transpose once when using the XY vertex meshes.
    """
    ids = np.asarray(projected['particle_id'])
    velocity = np.asarray(vlos_kms, dtype=float)
    xi, eta = tangent_coordinates(projected)
    front = np.asarray(projected['front_mask'])
    v0 = _finite_scalar(display_v0_kms, 'display_v0_kms')
    if ids.ndim != 1 or ids.dtype != np.dtype('uint64') or len(np.unique(ids)) != len(ids):
        raise ValueError('particle_id must be a one-dimensional unique uint64 array')
    if (velocity.shape != ids.shape or xi.shape != ids.shape or eta.shape != ids.shape
            or front.shape != ids.shape or front.dtype != np.dtype('bool')
            or not np.all(np.isfinite(velocity))):
        raise ValueError('one finite frozen LOS and sky position per ID are required')
    edges = np.linspace(-2.1, 2.1, 41)
    use = (front & np.isfinite(xi) & np.isfinite(eta)
           & (np.abs(xi) <= 2.1) & (np.abs(eta) <= 2.1))
    count = np.histogram2d(xi[use], eta[use], bins=(edges, edges))[0].astype(np.int64)
    total = np.histogram2d(xi[use], eta[use], bins=(edges, edges),
                           weights=velocity[use] - v0)[0]
    mean = np.divide(total, count, out=np.full((40, 40), np.nan), where=count > 0)
    xx, yy = np.meshgrid(edges, edges, indexing='xy')
    return dict(count=count, mean_HRF_minus_v0=mean,
                xi_vertices_deg=xx, eta_vertices_deg=yy, xi_deg=xi, eta_deg=eta,
                edges_xi_deg=edges.copy(), edges_eta_deg=edges.copy(),
                display_mask=use, display_particle_id=ids[use].copy(),
                display_v0_kms=v0, display_color_limits_kms=(-12., 12.),
                display_contract_id=TAN_DISPLAY_CONTRACT_ID,
                display_coordinate_system='gnomonic_TAN_xi_eta_deg',
                field_half_width_deg=2.1, aspect='equal', north_up_east_left=True,
                grid_array_order='xi_bin,eta_bin; transpose once for XY vertex mesh')


def tangent_motion_arrow(unit_EN, *, length_deg=.65):
    """Direction-only arrow at the origin, preserving its spherical length."""
    direction = np.asarray(unit_EN, dtype=float)
    if (direction.shape != (2,) or not np.all(np.isfinite(direction))
            or np.linalg.norm(direction) <= 0):
        raise ValueError('a nonzero finite centre east/north direction is required')
    length = _finite_scalar(length_deg, 'length_deg', positive=True)
    if length >= 90:
        raise ValueError('arrow angular length must be in (0,90) degrees')
    direction = direction / np.linalg.norm(direction)
    return dict(origin_xi_eta_deg=np.array([0., 0.]),
                tip_xi_eta_deg=np.rad2deg(np.tan(np.deg2rad(length))) * direction,
                angular_length_deg=length, length_represents='direction_only',
                projection='TAN', display_contract_id=TAN_DISPLAY_CONTRACT_ID)


def raw_offsets_to_tangent(projected, delta_ra_deg, delta_dec_deg):
    """Project a historical raw-RA/Dec point without changing that point.

    Useful for existing diagnostic centres saved in CSVs.  It is not a rescale
    by cos(dec): the complete absolute ray is reconstructed.  Missing inputs
    and points outside the forward hemisphere return NaN coordinates.
    """
    ra0 = _finite_scalar(projected['center_RA_deg'], 'center_RA_deg')
    dec0 = _finite_scalar(projected['center_Dec_deg'], 'center_Dec_deg')
    if abs(dec0) > 90:
        raise ValueError('center_Dec_deg must be within [-90,90]')
    dra, ddec = np.broadcast_arrays(np.asarray(delta_ra_deg, dtype=float),
                                   np.asarray(delta_dec_deg, dtype=float))
    dec_deg = dec0 + ddec
    if np.any(np.isfinite(dec_deg) & (np.abs(dec_deg) > 90)):
        raise ValueError('raw offsets must describe a valid absolute declination')
    da, dec, d0 = np.deg2rad(dra), np.deg2rad(dec_deg), np.deg2rad(dec0)
    # The centre RA cancels in these scalar products, including across RA=0.
    denominator = np.sin(dec) * np.sin(d0) + np.cos(dec) * np.cos(d0) * np.cos(da)
    east = np.cos(dec) * np.sin(da)
    north = np.sin(dec) * np.cos(d0) - np.cos(dec) * np.sin(d0) * np.cos(da)
    valid = np.isfinite(dra) & np.isfinite(ddec) & (denominator > 0)
    xi = np.divide(east, denominator, out=np.full(dra.shape, np.nan), where=valid)
    eta = np.divide(north, denominator, out=np.full(dra.shape, np.nan), where=valid)
    return np.rad2deg(xi), np.rad2deg(eta)


def tangent_pixel_area_arcsec2(npix, half_width_deg=2.1):
    """Exact solid angle of each square-grid TAN pixel, [xi_bin, eta_bin].

    Integrates dOmega=dx*dy/(1+x*x+y*y)**(3/2) over each cell, where x/y
    are dimensionless TAN ratios.  The result is square arcseconds, not a
    projected physical kpc area or a raw RA/Dec row-area approximation.
    """
    if isinstance(npix, (bool, np.bool_)) or not isinstance(npix, (int, np.integer)) or npix <= 0:
        raise ValueError('npix must be a positive integer')
    half_width = _finite_scalar(half_width_deg, 'half_width_deg', positive=True)
    edges = np.deg2rad(np.linspace(-half_width, half_width, int(npix) + 1))
    x, y = np.meshgrid(edges, edges, indexing='ij')
    primitive = np.arctan2(x * y, np.sqrt(1. + x*x + y*y))
    steradians = np.diff(np.diff(primitive, axis=0), axis=1)
    return steradians * (180. / np.pi * 3600.)**2
