"""Construct physical, explicitly unevolved orbit/PA candidates from frozen states.

Read-only inputs: native state JSON/NPZ, fixed observer configuration and the
actual one-pass merger parameters. Outputs are new snapshot-adjacent candidate
files. No simulation, observer change or production parameter replacement occurs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'math'))
import candidate_math as cm
from .target_frame import resolve_target

TAG = 'orbit_correction_7004_7005_20261020_v1'
OBSERVER = 'astropy_native_direct_v1'
DISPLAY = 'heliocentric_tan_xi_eta_2p1_v1'
MAS_FACTOR = 4.740470463533349
API = {}


def plain(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, Path): return str(value)
    if isinstance(value, dict): return {key: plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)): return [plain(item) for item in value]
    return value


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''): h.update(chunk)
    return h.hexdigest()


def save(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(plain(value), stream, ensure_ascii=False, indent=2)
        stream.write('\n')


def vector(value):
    result = np.asarray(value, dtype=float)
    assert result.shape == (3,) and np.isfinite(result).all()
    return result


def center(state, name='old_full'):
    item = state['centers'][name]
    return vector(item['r_grf_kpc']), vector(item['v_grf_kms'])


def sky_basis(r_center, config):
    matrix = np.asarray(config['physical_to_icrs'], dtype=float)
    h = vector(r_center) - vector(config['sun_position_native_kpc'])
    distance = float(np.linalg.norm(h))
    n = h / distance
    n_icrs = matrix @ n
    ra, dec = np.arctan2(n_icrs[1], n_icrs[0]), np.arcsin(np.clip(n_icrs[2], -1, 1))
    east_icrs = np.array([-np.sin(ra), np.cos(ra), 0.0])
    north_icrs = np.cross(n_icrs, east_icrs)
    east, north = matrix.T @ east_icrs, matrix.T @ north_icrs
    assert np.max(abs(np.cross(east, north) - n)) < 2e-12
    return dict(n=n, east=east, north=north, D_kpc=distance,
                RA_deg=float(np.rad2deg(ra) % 360), Dec_deg=float(np.rad2deg(dec)))


def independent_point_observables(r, v, config):
    basis = sky_basis(r, config)
    u = vector(v) - vector(config['sun_velocity_native_kms'])
    return dict(RA_deg=basis['RA_deg'], Dec_deg=basis['Dec_deg'], D_kpc=basis['D_kpc'],
                vlos_kms=float(u @ basis['n']),
                mu_alpha_cosdec_masyr=float(u @ basis['east'] / (MAS_FACTOR * basis['D_kpc'])),
                mu_delta_masyr=float(u @ basis['north'] / (MAS_FACTOR * basis['D_kpc'])),
                R_gc_kpc=float(np.linalg.norm(r)), V3D_kms=float(np.linalg.norm(v)),
                Vr_gc_kms=float(vector(r) @ vector(v) / np.linalg.norm(r)))


def point_observables(r, v, config):
    catalogue = API['observables'](vector(r)[None, :], vector(v)[None, :], np.array([0], dtype=np.uint64))
    result = dict(RA_deg=float(catalogue['ra'][0]), Dec_deg=float(catalogue['dec'][0]),
                  D_kpc=float(catalogue['rh'][0]), vlos_kms=float(catalogue['vlos'][0]),
                  mu_alpha_cosdec_masyr=float(catalogue['pmra'][0]), mu_delta_masyr=float(catalogue['pmdec'][0]),
                  R_gc_kpc=float(np.linalg.norm(r)), V3D_kms=float(np.linalg.norm(v)),
                  Vr_gc_kms=float(vector(r) @ vector(v) / np.linalg.norm(r)))
    independent = independent_point_observables(r, v, config)
    assert max(abs(result[key] - independent[key]) for key in result) < 1e-8
    return result


def independent_project(r_center, dr, config):
    basis = sky_basis(r_center, config)
    h0 = vector(r_center) - vector(config['sun_position_native_kpc'])
    h = np.asarray(dr, dtype=float) + h0
    denominator = h @ basis['n']
    ex, ny = h @ basis['east'], h @ basis['north']
    x, y = basis['D_kpc'] * ex / denominator, basis['D_kpc'] * ny / denominator
    distance = np.linalg.norm(h, axis=1)
    rproj = basis['D_kpc'] * np.hypot(ex, ny) / denominator
    front = denominator > 0
    return dict(basis=basis, x_kpc=x, y_kpc=y,
                xi_deg=np.rad2deg(ex / denominator), eta_deg=np.rad2deg(ny / denominator),
                per_star_D_kpc=distance, Rproj_kpc=rproj, front=front,
                R5=front & (rproj < 5.0))


def project(r_center, dr, config, ids, center_velocity, dv):
    catalogue = API['observables'](np.asarray(dr, dtype=float) + vector(r_center),
                                   np.asarray(dv, dtype=float) + vector(center_velocity),
                                   np.asarray(ids, dtype=np.uint64))
    assert catalogue['canonical_observer_id'] == config['canonical_observer_id']
    assert catalogue['observer_config_sha256'] == config['observer_config_sha256']
    sky = API['sky'](catalogue, vector(r_center) - vector(config['sun_position_native_kpc']))
    xi, eta = API['tangent'](sky)
    independent = independent_project(r_center, dr, config)
    errors = dict(xi_max_error_deg=float(np.max(abs(xi - independent['xi_deg']))),
                  eta_max_error_deg=float(np.max(abs(eta - independent['eta_deg']))),
                  D_max_error_kpc=float(np.max(abs(catalogue['rh'] - independent['per_star_D_kpc']))),
                  Rproj_max_error_kpc=float(np.max(abs(sky['projected_radius_kpc'] - independent['Rproj_kpc']))))
    assert max(errors.values()) < 2e-8
    assert np.array_equal(sky['fit_R5_mask'], independent['R5'])
    native = catalogue['_native']
    assert np.max(abs(native['physical_to_icrs'] - np.asarray(config['physical_to_icrs']))) < 1e-12
    assert np.max(abs(native['observer']['sun_position_kpc'] - vector(config['sun_position_native_kpc']))) < 1e-12
    assert np.max(abs(native['observer']['sun_velocity_kms'] - vector(config['sun_velocity_native_kms']))) < 1e-12
    basis = dict(n=sky['LOS_native'], east=sky['east_native'], north=sky['north_native'],
                 D_kpc=sky['center_distance_kpc'], RA_deg=sky['center_RA_deg'], Dec_deg=sky['center_Dec_deg'])
    return dict(basis=basis, x_kpc=sky['east_kpc'], y_kpc=sky['north_kpc'], xi_deg=xi, eta_deg=eta,
                per_star_D_kpc=catalogue['rh'], Rproj_kpc=sky['projected_radius_kpc'],
                R5=sky['fit_R5_mask'], front=sky['front_mask'],
                source='canonical_observer.production_observables/production_sky + canonical_display.tangent_coordinates',
                independent_formula_QA=errors,
                no_R_LOS_QA=native['full_astropy_independent_check'])


def measure(r_center, dr, mass, config, Analysis, ids, center_velocity, dv):
    projected = project(r_center, dr, config, ids, center_velocity, dv)
    use = projected['front'] & np.isfinite(projected['x_kpc']) & np.isfinite(projected['y_kpc'])
    assert use.all(), 'Frozen old tracer behind observer: cannot silently change shape support'
    shape = Analysis.calculate_projected_shape(projected['x_kpc'], projected['y_kpc'], mass, n_neighbors=30)
    rhalf = Analysis.half_light_radius(projected['x_kpc'], projected['y_kpc'], mass,
                                     shape['eps'], pa=shape['pa'],
                                     center_x=shape['center_x'], center_y=shape['center_y'])
    result = dict(PA_deg=float((90.0 - np.degrees(shape['pa'])) % 180.0),
                  q=float(shape['q']), pa_math_rad=float(shape['pa']),
                  rhalf_major_kpc=float(rhalf), center_x_kpc=float(shape['center_x']),
                  center_y_kpc=float(shape['center_y']), N=len(mass),
                  definition='All fixed old tracers; mass times 30-neighbour density weighted covariance; canonical TAN')
    assert np.isfinite([result['PA_deg'], result['q'], result['rhalf_major_kpc']]).all()
    return result, projected


def load_frame(reference, config):
    state_path, cache_path = Path(reference['state_json']), Path(reference['cache_npz'])
    state = json.loads(state_path.read_text())
    assert state['state'] == 'PASS'
    assert state['host']['canonical_observer_id'] == config['canonical_observer_id']
    assert state['host']['observer_config_sha256'] == config['observer_config_sha256']
    assert state['cache_sha256'] == sha(cache_path)
    with np.load(cache_path, allow_pickle=False) as cache:
        old = cache['is_old'].astype(bool)
        ids = cache['id'][old]
        assert ids.dtype == np.uint64 and len(np.unique(ids)) == len(ids)
        order = np.argsort(ids)
        r, v = center(state)
        result = dict(id=ids[order], mass=cache['mass_msun'][old][order].astype(float),
                      dr=cache['dr_native_kpc'][old][order].astype(float),
                      dv=cache['dv_native_kms'][old][order].astype(float), r=r, v=v,
                      state=state, state_path=str(state_path), cache_path=str(cache_path))
    assert len(result['id']) == state['centers']['old_full']['N']
    assert (result['mass'] > 0).all()
    return result


def interpolate_frames(left, right, fraction):
    assert np.array_equal(left['id'], right['id']), 'Both endpoints must contain exactly the frozen old cohort'
    f = float(fraction)
    r, v = (1 - f) * left['r'] + f * right['r'], (1 - f) * left['v'] + f * right['v']
    mass = (1 - f) * left['mass'] + f * right['mass']
    dr = (1 - f) * (left['dr'] + left['r'] - r) + f * (right['dr'] + right['r'] - r)
    dv = (1 - f) * (left['dv'] + left['v'] - v) + f * (right['dv'] + right['v'] - v)
    # Cache quantization or actual mass change can leave tiny residual first
    # moments. Explicitly zero these relative offsets while retaining the
    # separately defined linear COM anchor; record the correction, never hide it.
    dr_bias, dv_bias = np.average(dr, axis=0, weights=mass), np.average(dv, axis=0, weights=mass)
    dr -= dr_bias
    dv -= dv_bias
    time = (1 - f) * left['state']['header_time'] + f * right['state']['header_time']
    time_unit = left['state']['unit_time_to_Gyr']
    assert abs(time_unit - right['state']['unit_time_to_Gyr']) < 1e-12
    return dict(id=left['id'], mass=mass, dr=dr, dv=dv, r=r, v=v,
                fraction=f, header_time=float(time), physical_time_Gyr=float(time * time_unit),
                interpolation_definition='Same frozen uint64 IDs; linear native position/velocity/mass; residual offset first moments explicitly re-zeroed',
                dr_recentering_kpc=dr_bias, dv_recentering_kms=dv_bias,
                endpoints=[left['state']['numsp'], right['state']['numsp']],
                endpoint_cache_sha256=[left['state']['cache_sha256'], right['state']['cache_sha256']],
                endpoint_relative_mass_change_max=float(np.max(abs(right['mass'] - left['mass']) / left['mass'])))


def radial_anchor(left, right, target_radius, origin):
    rleft, rright = left['r'] - origin, right['r'] - origin
    function = lambda f: np.linalg.norm((1 - f) * rleft + f * rright) - target_radius
    assert function(0) >= -1e-10 and function(1) <= 1e-10, 'Selected bracket is not inward across the requested radius'
    fraction = brentq(function, 0, 1, xtol=2e-14)
    result = interpolate_frames(left, right, fraction)
    assert abs(function(fraction)) < 2e-10
    return result


def compact_anchor(frame, config):
    return dict(endpoints=frame['endpoints'], fraction=frame['fraction'], r_grf_kpc=frame['r'],
                v_grf_kms=frame['v'], header_time=frame['header_time'],
                physical_time_Gyr=frame['physical_time_Gyr'],
                observables_of_old_COM=point_observables(frame['r'], frame['v'], config),
                old_tracer_N=len(frame['id']), id_sha256=hashlib.sha256(frame['id'].astype('<u8').tobytes()).hexdigest(),
                interpolation_definition=frame['interpolation_definition'],
                residual_dr_recentering_kpc=frame['dr_recentering_kpc'],
                residual_dv_recentering_kms=frame['dv_recentering_kms'],
                endpoint_relative_mass_change_max=frame['endpoint_relative_mass_change_max'],
                endpoint_cache_sha256=frame['endpoint_cache_sha256'])


def initial_biases(initial, galaxy, cube):
    host = initial['state']['host']
    length, velocity = float(host['length_to_kpc']), float(host['velocity_to_kms'])
    h, b = vector(host['host_center_input_code']), vector(host['host_bulk_input_code'])
    q, g = cm.legacy_matrix(cube[9:12]), cm.legacy_matrix(galaxy[9:12])
    cmw = np.linalg.solve(q, h - cube[:3])
    bmw = np.linalg.solve(q, b - cube[3:6])
    cd = np.linalg.solve(g, np.linalg.solve(q, initial['r'] / length + h - cube[:3]) - galaxy[:3])
    bd = np.linalg.solve(g, np.linalg.solve(q, initial['v'] / velocity + b - cube[3:6]) - galaxy[3:6])
    restored = cm.placement_for_target(initial['r'] / length, initial['v'] / velocity,
                                      galaxy[9:12], cube[9:12], cmw, bmw, cd, bd)
    assert np.max(abs(restored - galaxy[:6])) < 1e-9
    return dict(cM=cmw, bM=bmw, cD=cd, bD=bd, length_to_kpc=length, velocity_to_kms=velocity,
                original_placement_recovery_error=float(np.max(abs(restored - galaxy[:6]))),
                definition='Effective input-tracer centre/bulk offsets inferred from actual snapshot-0 old COM and canonical MW H/B; not independently measured raw-IC global COMs')


def merger_candidate(initial_r, initial_v, old_galaxy, cube, biases, body_roll, plane):
    orientation = cm.corrected_internal_angles(old_galaxy[9:12], cube[9:12], body_roll, plane)
    placement = cm.placement_for_target(initial_r / biases['length_to_kpc'], initial_v / biases['velocity_to_kms'],
                                        orientation['angles_deg'], cube[9:12],
                                        biases['cM'], biases['bM'], biases['cD'], biases['bD'])
    parameters = np.r_[placement, old_galaxy[6:9], orientation['angles_deg']]
    # This is what the C parser actually receives from the candidate text.
    parameters = np.array([float(format(float(x), '.9g')) for x in parameters], dtype=np.float32)
    actual_r, actual_v = cm.initial_state_from_merge(parameters, cube, biases['cM'], biases['bM'], biases['cD'], biases['bD'])
    dr = actual_r * biases['length_to_kpc'] - initial_r
    dv = actual_v * biases['velocity_to_kms'] - initial_v
    assert np.max(abs(dr)) < 5e-4 and np.max(abs(dv)) < 5e-4
    return dict(parameters_12=parameters.astype(float), orientation=orientation,
                initial_r_grf_kpc=initial_r, initial_v_grf_kms=initial_v,
                float32_parsed_linear_centroid_r_error_kpc=dr,
                float32_parsed_linear_centroid_v_error_kms=dv,
                preserved_first_input_angles=old_galaxy[6:9],
                note='Linear centroid closure only; does not replace full particle IC validation or live evolution')


def write_par(path, input_names, merger):
    values = ' '.join(format(float(x), '.9g') for x in merger['parameters_12'])
    text = ' '.join(input_names) + '\n' + values + '\n\n'
    text += '# CANDIDATE ONLY: not re-evolved, not installed, no pericentre certified.\n'
    text += '# Current one-pass merger: preserve associated mergecube.par unchanged.\n'
    with Path(path).open('x', encoding='utf-8') as stream: stream.write(text)


def invariants(before, after, shape_before, shape_after):
    distance_error = float(np.max(abs(before['per_star_D_kpc'] - after['per_star_D_kpc'])))
    radius_error = float(np.max(abs(before['Rproj_kpc'] - after['Rproj_kpc'])))
    q_error = abs(shape_before['q'] - shape_after['q'])
    re_error = abs(shape_before['rhalf_major_kpc'] - shape_after['rhalf_major_kpc'])
    assert np.array_equal(before['R5'], after['R5'])
    assert distance_error < 1e-9 and radius_error < 1e-8
    assert q_error < 2e-7 and re_error < 2e-6
    return dict(per_star_D_max_error_kpc=distance_error, Rproj_max_error_kpc=radius_error,
                identical_R5_particle_mask=True, R5_N=int(before['R5'].sum()),
                q_abs_error=q_error, rhalf_abs_error_kpc=re_error,
                before_projection_QA=before['independent_formula_QA'], after_projection_QA=after['independent_formula_QA'])


def check_target(config):
    from astropy import units as u
    from astropy.coordinates import SkyCoord, Galactocentric, ICRS, CartesianDifferential
    assert config['canonical_observer_id'] == OBSERVER
    matrix = np.asarray(config['physical_to_icrs'], dtype=float)
    assert np.max(abs(matrix.T @ matrix - np.eye(3))) < 1e-11 and np.linalg.det(matrix) > .999999999
    sun = vector(config['sun_position_native_kpc'])
    vs = vector(config['sun_velocity_native_kms'])
    assert np.max(abs(sun - [-np.sqrt(8.122 ** 2 - .0208 ** 2), 0, .0208])) < 1e-12
    assert np.max(abs(vs - [12.9, 245.6, 7.78])) < 1e-12
    target = config['target']
    conversion=resolve_target(target)
    frame = Galactocentric(galcen_distance=8.122 * u.kpc, z_sun=.0208 * u.kpc,
                          galcen_v_sun=CartesianDifferential(vs * u.km / u.s),
                          galcen_coord=ICRS(ra=266.4051 * u.deg, dec=-28.936175 * u.deg), roll=0 * u.deg)
    point = SkyCoord(ra=target['ra_deg'] * u.deg, dec=target['dec_deg'] * u.deg,
                     distance=target['D_kpc'] * u.kpc, radial_velocity=target['vlos_kms'] * u.km / u.s,
                     pm_ra_cosdec=target['mu_alpha_cosdec_masyr'] * u.mas / u.yr,
                     pm_dec=target['mu_delta_masyr'] * u.mas / u.yr, frame='icrs').transform_to(frame)
    r, v = vector(target['r_grf_native_kpc']), vector(target['v_grf_native_kms'])
    rerr = float(np.max(abs(point.cartesian.xyz.to_value(u.kpc) - r)))
    verr = float(np.max(abs(point.velocity.d_xyz.to_value(u.km / u.s) - v)))
    assert rerr < 1e-7 and verr < 1e-7
    obs = point_observables(r, v, config)
    for output, source in [('RA_deg', 'ra_deg'), ('Dec_deg', 'dec_deg'), ('D_kpc', 'D_kpc'),
                           ('vlos_kms', 'vlos_kms'), ('mu_alpha_cosdec_masyr', 'mu_alpha_cosdec_masyr'),
                           ('mu_delta_masyr', 'mu_delta_masyr')]:
        assert abs(obs[output] - target[source]) < 1e-7, output
    return dict(state='PASS', astropy_r_max_error_kpc=rerr, astropy_v_max_error_kms=verr,
                target_old_COM_point_observables=obs, fixed_observer=True,
                conversion_contract=conversion,
                not_legacy_140_minus37=True)


def render_figure(path, model, frames, shapes):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    edges = np.linspace(-2.1, 2.1, 121)
    grids = [np.histogram2d(item['xi_deg'], item['eta_deg'], bins=(edges, edges), weights=mass)[0]
             for item, mass in frames]
    logs = [np.where(grid > 0, np.log10(np.maximum(grid, 1e-100)), np.nan) for grid in grids]
    vmin = min(float(np.nanmin(grid)) for grid in logs)
    vmax = max(float(np.nanmax(grid)) for grid in logs)
    fig, axes = plt.subplots(1, 3, figsize=(14.8, 5.4), constrained_layout=True)
    names = ['Matched GC-radius baseline', 'Physical target 6D registration', 'Registration + internal PA=43 deg']
    for ax, (projected, _), grid, shape, title in zip(axes, frames, logs, shapes, names):
        mesh = ax.pcolormesh(edges, edges, grid.T, shading='flat', cmap='magma', vmin=vmin, vmax=vmax)
        t = np.linspace(-2.0, 2.0, 2)
        ax.plot(t * np.sin(np.deg2rad(43)), t * np.cos(np.deg2rad(43)), '--', color='cyan', lw=1.2, label='43 deg reference')
        ax.set(xlim=(2.1, -2.1), ylim=(-2.1, 2.1), xlabel=r'$\xi$ (deg)', ylabel=r'$\eta$ (deg)')
        ax.set_aspect('equal', adjustable='box')
        basis = projected['basis']
        ax.set_title(title + '\nPA=%.2f deg, q=%.3f, Re=%.3f kpc' % (shape['PA_deg'], shape['q'], shape['rhalf_major_kpc']), fontsize=10)
        ax.text(.04, .04, 'RA %.3f, Dec %.3f\nD %.3f kpc' % (basis['RA_deg'], basis['Dec_deg'], basis['D_kpc']),
                transform=ax.transAxes, color='white', fontsize=8, bbox=dict(facecolor='black', alpha=.55, edgecolor='none'))
    axes[0].legend(loc='upper right', fontsize=8)
    fig.colorbar(mesh, ax=axes, label=r'$\log_{10}(M_{\rm old}/M_\odot)$ per TAN pixel', shrink=.77)
    fig.suptitle(model + ': static candidate geometry; no re-evolved orbit\nSame fixed observer and old-star IDs; 43 deg is north toward east', fontsize=12)
    fig.savefig(path, dpi=180)
    fig.savefig(Path(path).with_suffix('.pdf'))
    plt.close(fig)


def solve(index_path, source_code=None):
    index_path = Path(index_path).resolve()
    index = json.loads(index_path.read_text())
    config, inputs = index['config'], index['server_inputs']
    if source_code:
        source = Path(source_code).resolve()
    elif config.get('analysis_source'):
        source = Path(config['analysis_source'])
        if not source.is_absolute(): source = Path('/travail/xhuang') / source
    elif (HERE.parent / 'canonical_observer.py').is_file():
        source = HERE.parent
    else:
        source = Path('/travail/xhuang') / inputs['active_code']
    sys.path.insert(0, str(source))
    from canonical_observer import production_observables, production_sky, observer_contract_sha256, CANONICAL_OBSERVER_ID
    from canonical_display import tangent_coordinates, TAN_DISPLAY_CONTRACT_ID
    assert CANONICAL_OBSERVER_ID == OBSERVER and TAN_DISPLAY_CONTRACT_ID == DISPLAY
    assert observer_contract_sha256() == config['observer_config_sha256']
    API.update(observables=production_observables, sky=production_sky, tangent=tangent_coordinates)
    target_qa = check_target(config)
    from analysis_core import Analysis
    target_r, target_v = vector(config['target']['r_grf_native_kpc']), vector(config['target']['v_grf_native_kms'])
    target_pa = float(config['target']['morphology_PA_E_of_N_deg'])
    summary = dict(state='RUNNING', canonical_observer_id=OBSERVER, display_contract_id=DISPLAY,
                   observer_config_sha256=config['observer_config_sha256'], target=config['target'],
                   target_6D_QA=target_qa, no_new_simulation=True, analysis_source=str(source),
                   projection_source_sha256={name: sha(source / name) for name in
                                             ['canonical_observer.py', 'canonical_display.py', 'analysis_core.py']}, cases=[])
    csv_rows = []
    for case in config['cases']:
        model = case['model']
        refs = index['states'][model]
        initial = load_frame(refs[str(case['initial'])], config)
        latest = load_frame(refs[str(case['latest'])], config)
        for state in [initial['state'], latest['state']]:
            for name, digest in summary['projection_source_sha256'].items():
                assert state['source_module_sha256'][name] == digest, 'Extraction/projection source mismatch: ' + name
        latest_dir = Path(case.get('output_dir', Path(latest['state_path']).parent))
        latest_dir.mkdir(parents=True, exist_ok=True)
        out, figures = latest_dir / 'candidates', latest_dir / 'figures'
        out.mkdir(exist_ok=False)
        figures.mkdir(exist_ok=False)
        gal_record, cube_record = inputs['models'][model]['mergegalaxy.par'], inputs['models'][model]['mergecube.par']
        galaxy, cube = np.array(gal_record[1], float), np.array(cube_record[1], float)
        biases = initial_biases(initial, galaxy, cube)
        selected = {}
        for name, bracket, radius, origin in [
            ('GC_radius', case['GC_radius_bracket'], float(np.linalg.norm(target_r)), np.zeros(3)),
            ('Dhel', case['crossing_bracket'], float(config['target']['D_kpc']), vector(config['sun_position_native_kpc']))]:
            left, right = [load_frame(refs[str(n)], config) for n in bracket]
            selected[name] = radial_anchor(left, right, radius, origin)
        anchor = selected['GC_radius']
        coefficient = cm.empirical_coefficients(initial['r'], initial['v'], anchor['r'], anchor['v'])
        solved = cm.solve_initial_velocity(initial['r'], initial['v'], target_r, target_v,
                                           coefficient['C_E'], coefficient['k_L'], radial_sign=-1)
        plane = cm.orbit_plane_rotation(anchor['r'], anchor['v'], target_r, target_v)
        assert np.max(abs(plane @ anchor['r'] - target_r)) < 1e-8
        baseline_shape, baseline_projection = measure(anchor['r'], anchor['dr'], anchor['mass'], config, Analysis,
                                                     anchor['id'], anchor['v'], anchor['dv'])
        registered_dr = anchor['dr'] @ plane.T
        registered_dv = anchor['dv'] @ plane.T
        registered_shape, registered_projection = measure(target_r, registered_dr, anchor['mass'], config, Analysis,
                                                         anchor['id'], target_v, registered_dv)
        delta_pa = cm.pa_difference_deg(target_pa, registered_shape['PA_deg'])
        internal = cm.pa_roll_rotation(registered_projection['basis']['n'], delta_pa)
        corrected_dr = registered_dr @ internal.T
        corrected_dv = registered_dv @ internal.T
        assert np.max(abs(np.average(corrected_dv, axis=0, weights=anchor['mass']))) < 1e-8
        corrected_shape, corrected_projection = measure(target_r, corrected_dr, anchor['mass'], config, Analysis,
                                                       anchor['id'], target_v, corrected_dv)
        assert abs(cm.pa_difference_deg(target_pa, corrected_shape['PA_deg'])) < .02
        roll_qa = invariants(registered_projection, corrected_projection, registered_shape, corrected_shape)
        new_initial_r, new_initial_v = plane @ initial['r'], plane @ solved['velocity']
        merger = merger_candidate(new_initial_r, new_initial_v, galaxy, cube, biases, internal, plane)
        joint = dict(candidate_id='joint_GC_radius_target6D_PA43', status='UNEVOLVED_STATIC_CANDIDATE',
                     model=model, observer_unchanged=True, anchor=compact_anchor(anchor, config),
                     empirical_coefficients=coefficient, initial_velocity_solution=solved,
                     target_point_6D=config['target'], orbital_rotation_native=plane,
                     internal_rotation_native=internal, requested_PA_delta_after_registration_deg=delta_pa,
                     baseline_shape=baseline_shape, registered_shape=registered_shape, corrected_shape=corrected_shape,
                     pure_internal_roll_invariants=roll_qa,
                     static_exact_target_bulk_delta_from_rotated_measured_bulk_kms=target_v - plane @ anchor['v'],
                     static_old_mean_relative_velocity_residual_kms=np.average(corrected_dv, axis=0, weights=anchor['mass']),
                     mass_profile_3D_preserved_by_rigid_registration=True,
                     q_change_registration=registered_shape['q'] - baseline_shape['q'],
                     rhalf_relative_change_registration=registered_shape['rhalf_major_kpc'] / baseline_shape['rhalf_major_kpc'] - 1,
                     merger=merger, input_offsets=biases, pericentre_measured=False,
                     limitations=['Final bulk is set exactly only in this static geometry proxy.',
                                  'Empirical C_E/k_L and rigid plane transport do not certify a new live-host trajectory.',
                                  'Changing initial body orientation may change gas interaction and later morphology.',
                                  'Only fixed old stars were interpolated; no claim of a complete new-star catalogue.',
                                  'No pericentre has been enforced or predicted.'])
        save(out / 'joint_GC_PA43.json', joint)
        write_par(out / 'mergegalaxy.joint_GC_PA43.CANDIDATE.par', gal_record[0], merger)

        d_anchor = selected['Dhel']
        d_coeff = cm.empirical_coefficients(initial['r'], initial['v'], d_anchor['r'], d_anchor['v'])
        try:
            d_solved = cm.solve_initial_velocity(initial['r'], initial['v'], target_r, target_v,
                                                 d_coeff['C_E'], d_coeff['k_L'], radial_sign=-1)
            d_feasible, d_error = True, None
        except ValueError as error:
            d_solved, d_feasible, d_error = None, False, str(error)
        sensitivity = dict(status='SENSITIVITY_ONLY_NOT_MAIN_SOLUTION', anchor=compact_anchor(d_anchor, config),
                           empirical_coefficients=d_coeff, initial_velocity_solution=d_solved,
                           real_initial_velocity_feasible=d_feasible, infeasibility_reason=d_error,
                           delta_C_E_from_GC_primary=d_coeff['C_E'] - coefficient['C_E'],
                           delta_k_L_from_GC_primary=d_coeff['k_L'] - coefficient['k_L'],
                           delta_seed_vr0_kms=(d_solved['vr0'] - solved['vr0']) if d_feasible else None,
                           delta_seed_vt0_kms=(d_solved['vt0'] - solved['vt0']) if d_feasible else None,
                           warning='Equal heliocentric D does not imply equal MW radius/potential stage; primary calibration matches Galactocentric radius.')
        save(out / 'sensitivity_Dhel_anchor.json', sensitivity)

        latest_shape, latest_projection = measure(latest['r'], latest['dr'], latest['mass'], config, Analysis,
                                                 latest['id'], latest['v'], latest['dv'])
        stored = latest['state']['shape_old_full']
        assert abs(cm.pa_difference_deg(stored['PA_deg'], latest_shape['PA_deg'])) < .02
        latest_delta = cm.pa_difference_deg(target_pa, latest_shape['PA_deg'])
        latest_a = cm.pa_roll_rotation(latest_projection['basis']['n'], latest_delta)
        latest_corrected_shape, latest_corrected_projection = measure(latest['r'], latest['dr'] @ latest_a.T,
                                                                      latest['mass'], config, Analysis,
                                                                      latest['id'], latest['v'], latest['dv'] @ latest_a.T)
        assert abs(cm.pa_difference_deg(target_pa, latest_corrected_shape['PA_deg'])) < .02
        latest_qa = invariants(latest_projection, latest_corrected_projection, latest_shape, latest_corrected_shape)
        latest_merger = merger_candidate(initial['r'], initial['v'], galaxy, cube, biases, latest_a, np.eye(3))
        local = dict(candidate_id='latest_internal_PA43_only', status='UNEVOLVED_STATIC_CANDIDATE',
                     model=model, snapshot=case['latest'], requested_PA_delta_deg=latest_delta,
                     old_shape=latest_shape, corrected_shape=latest_corrected_shape, rotation_native=latest_a,
                     latest_old_COM_r_grf_kpc=latest['r'], latest_old_COM_v_grf_kms=latest['v'],
                     COM_orbit_kept_fixed_in_static_proxy=True, observer_unchanged=True,
                     initial_r_grf_kpc=initial['r'], initial_v_grf_kms=initial['v'],
                     roll_invariants=latest_qa, merger=latest_merger,
                     caveat='Applying the inferred body rotation to snapshot-0 merger angles requires new evolution; PA43 at the old final time is not guaranteed dynamically.')
        save(out / 'latest_internal_PA43_only.json', local)
        write_par(out / 'mergegalaxy.latest_internal_PA43_only.CANDIDATE.par', gal_record[0], latest_merger)
        figure = figures / (model + '_GC_registration_PA43_static.png')
        render_figure(figure, model,
                      [(baseline_projection, anchor['mass']), (registered_projection, anchor['mass']),
                       (corrected_projection, anchor['mass'])],
                      [baseline_shape, registered_shape, corrected_shape])
        record = dict(state='PASS', model=model, latest_snapshot=case['latest'], output_dir=str(latest_dir),
                      primary=joint, Dhel_sensitivity=sensitivity, latest_internal_only=local,
                      no_measured_pericentre=True,
                      latest_old_COM_radial_velocity_kms=cm.decompose_state(latest['r'], latest['v'])['vr'],
                      figure=str(figure), figure_sha256=sha(figure),
                      source_state_paths={num: value['state_json'] for num, value in refs.items()},
                      candidate_math_sha256=sha(HERE / 'math/candidate_math.py'), solver_sha256=sha(__file__))
        save(latest_dir / 'CANDIDATES_REPORT.json', record)
        summary['cases'].append(record)
        csv_rows.append(dict(model=model, latest_snapshot=case['latest'],
                             latest_PA_deg=latest_shape['PA_deg'], latest_internal_delta_PA_deg=latest_delta,
                             GC_anchor_fraction=anchor['fraction'], GC_anchor_physical_time_Gyr=anchor['physical_time_Gyr'],
                             C_E=coefficient['C_E'], k_L=coefficient['k_L'],
                             corrected_initial_vr0_kms=solved['vr0'], corrected_initial_vt0_kms=solved['vt0'],
                             registered_PA_deg=registered_shape['PA_deg'], corrected_PA_deg=corrected_shape['PA_deg'],
                             baseline_q=baseline_shape['q'], registered_q=registered_shape['q'],
                             baseline_rhalf_kpc=baseline_shape['rhalf_major_kpc'], registered_rhalf_kpc=registered_shape['rhalf_major_kpc'],
                             pericentre_measured=False, re_evolved=False))
        print('CANDIDATE PASS', model, 'latest PA', latest_shape['PA_deg'], 'target PA', corrected_shape['PA_deg'], flush=True)
    summary['state'] = 'PASS'
    save(index_path.parent / 'ORBIT_CANDIDATES_SUMMARY.json', summary)
    with (index_path.parent / 'ORBIT_CANDIDATES_SUMMARY.csv').open('x', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(csv_rows[0]))
        writer.writeheader()
        writer.writerows(csv_rows)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', required=True, type=Path)
    parser.add_argument('--source-code', type=Path, help='Testing only: same production shape implementation')
    args = parser.parse_args()
    solve(args.index, args.source_code)
