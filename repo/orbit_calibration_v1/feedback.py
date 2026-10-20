"""Measure evolved six-vectors and propose IC updates from evolved responses.

No orbit integrator, observer fit, final-state registration, job submission or
production-file replacement is performed. Every proposal needs a new evolution.
"""
from pathlib import Path
import csv
import hashlib
import json

import numpy as np

from .cli import OBSERVER_ID, OBSERVER_SHA, SUN, SUN_VELOCITY, sha, require_close
from .math import candidate_math as cm

RESIDUAL_KEYS = ['sky_east_deg', 'sky_north_deg', 'D_kpc', 'vlos_kms',
                 'mu_alpha_cosdec_masyr', 'mu_delta_masyr']
CONTROL_KEYS = ['initial_x_kpc', 'initial_y_kpc', 'initial_z_kpc',
                'initial_vx_kms', 'initial_vy_kms', 'initial_vz_kms']
TARGET_KEYS = ['ra_deg', 'dec_deg', 'D_kpc', 'vlos_kms',
               'mu_alpha_cosdec_masyr', 'mu_delta_masyr']


def plain(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {k: plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [plain(v) for v in value]
    return value


def digest(value):
    return hashlib.sha256(json.dumps(plain(value), sort_keys=True, allow_nan=False).encode()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def six(value, name, positive=False):
    value = np.asarray(value, dtype=float)
    if value.shape != (6,) or not np.isfinite(value).all() or (positive and np.any(value <= 0)):
        raise ValueError(name + ' requires six finite' + (' positive' if positive else '') + ' values')
    return value


def valid_sha(value):
    return isinstance(value, str) and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def evaluate_point(r, v, config):
    """One unchanged native COM, one shared observer, all six residuals together.

    Sky residuals use a spherical logarithmic map at the observed direction;
    their norm equals the great-circle separation. They are not a display TAN.
    """
    from canonical_observer import production_observables
    r, v = np.asarray(r, float), np.asarray(v, float)
    if r.shape != (3,) or v.shape != (3,) or not np.isfinite(np.r_[r, v]).all():
        raise ValueError('A finite native position and velocity are required')
    cat = production_observables(r[None], v[None], np.array([0], dtype=np.uint64))
    if cat['canonical_observer_id'] != OBSERVER_ID or cat['observer_config_sha256'] != OBSERVER_SHA:
        raise ValueError('Shared observer identity changed')
    require_close(cat['_native']['physical_to_icrs'], config['physical_to_icrs'], 'feedback physical mapping')
    a, d = np.radians([config['target']['ra_deg'], config['target']['dec_deg']])
    target_ray = np.array([np.cos(d)*np.cos(a), np.cos(d)*np.sin(a), np.sin(d)])
    east = np.array([-np.sin(a), np.cos(a), 0.])
    north = np.cross(target_ray, east)
    am, dm = np.radians([cat['ra'][0], cat['dec'][0]])
    model_ray = np.array([np.cos(dm)*np.cos(am), np.cos(dm)*np.sin(am), np.sin(dm)])
    xy = np.array([model_ray @ east, model_ray @ north])
    sine, cosine = float(np.linalg.norm(xy)), float(model_ray @ target_ray)
    if cosine < 0 and sine < 1e-12:
        raise ValueError('Antipodal sky position has no unique logarithmic residual')
    rho = np.arctan2(sine, cosine)
    sky = np.degrees(rho) * xy / sine if sine > 1e-15 else np.zeros(2)
    observed = dict(zip(TARGET_KEYS, [float(cat[k][0]) for k in ['ra', 'dec', 'rh', 'vlos', 'pmra', 'pmdec']]))
    residual = np.r_[sky, [observed[k]-config['target'][k] for k in TARGET_KEYS[2:]]]
    h, u = r-SUN, v-SUN_VELOCITY
    direct_los = float(u @ h / np.linalg.norm(h))
    require_close(direct_los, observed['vlos_kms'], 'direct feedback LOS', 1e-8)
    target_r = np.asarray(config['target']['r_grf_native_kpc'])
    target_v = np.asarray(config['target']['v_grf_native_kms'])
    angle = lambda x, y: float(np.degrees(np.arctan2(np.linalg.norm(np.cross(x, y)), x @ y)))
    return dict(observables=observed, residual=dict(zip(RESIDUAL_KEYS, residual)),
                residual_vector=residual, sky_separation_deg=float(np.degrees(rho)),
                r_grf_kpc=r, v_grf_kms=v, native_position_residual_kpc=r-target_r,
                native_velocity_residual_kms=v-target_v,
                GC_speed_kms=float(np.linalg.norm(v)), target_GC_speed_kms=float(np.linalg.norm(target_v)),
                GC_velocity_direction_error_deg=angle(v, target_v),
                point_definition='same frozen old_full mass-mean position and velocity; not a Walker fit',
                final_state_registration_applied=False, solar_or_camera_fit_applied=False)


def load_history(index_path, model):
    """Verify an extraction index and old-star caches, without re-reading raw HDF5.

    Raw snapshot hashes are inherited from the frozen extractor. This does not
    certify a live raw file as unchanged; that is a separate extraction check.
    """
    index_path = Path(index_path).resolve()
    index = load(index_path)
    c = index['config']
    if c['canonical_observer_id'] != OBSERVER_ID or c['observer_config_sha256'] != OBSERVER_SHA:
        raise ValueError('Legacy or mixed observer history')
    require_close(c['sun_position_native_kpc'], SUN, 'history Sun')
    require_close(c['sun_velocity_native_kms'], SUN_VELOCITY, 'history Sun velocity')
    cases = [case for case in c['cases'] if case['model'] == model]
    if len(cases) != 1 or cases[0]['initial'] != 0:
        raise ValueError('Exactly one requested model with initial snapshot 0 is required')
    inputs = index['server_inputs']['models'][model]
    for item in inputs['files'].values():
        if hashlib.sha256(item['text'].encode()).hexdigest() != item['sha256']:
            raise ValueError('Frozen input text/hash mismatch')
    if 'case.json' not in inputs['files']:
        raise ValueError('Run identity absent: re-extract with captured case.json; do not relabel an old index')
    case = json.loads(inputs['files']['case.json']['text'])
    if not isinstance(case.get('run_id'), str) or not case['run_id'].strip():
        raise ValueError('An explicit captured run_id is required')
    galaxy = np.asarray(inputs['mergegalaxy.par'][1], float)
    cube = np.asarray(inputs['mergecube.par'][1], float)
    if 'parameters_12' in case.get('orbit_provenance', {}):
        require_close(galaxy, case['orbit_provenance']['parameters_12'], 'run/merge identity', 2e-6)
    parameters = {}
    for line in inputs['files']['GZWJL.PARAM']['text'].splitlines():
        fields = line.split('%', 1)[0].split()
        if len(fields) >= 2 and fields[0] not in {'InitCondFile', 'OutputDir'}:
            parameters[fields[0]] = fields[1:]
    context = dict(model=model, observer=OBSERVER_SHA, parameters=parameters,
                   body_angles=galaxy[6:].tolist(), cube_sha256=inputs['files']['mergecube.par']['sha256'],
                   fixed_input_sha256=inputs.get('fixed_input_sha256', {}))
    frames, sources, fixed_ids, fixed_cohort, extraction_source = [], [], None, None, None
    for number, reference in sorted(index['states'][model].items(), key=lambda item: int(item[0])):
        sp, cp = Path(reference['state_json']), Path(reference['cache_npz'])
        if not sp.is_absolute(): sp = index_path.parent/sp
        if not cp.is_absolute(): cp = index_path.parent/cp
        if not valid_sha(reference.get('state_sha256')) or sha(sp) != reference['state_sha256']:
            raise ValueError('State SHA absent or changed; create a newly verified extraction index')
        state = load(sp)
        if (state.get('schema_id') != 'canonical_orbit_state_v1' or state.get('state') != 'PASS'
                or state['model'] != model or state['numsp'] != int(number)):
            raise ValueError('Not a verified evolved extraction state')
        if (state['canonical_observer_id'] != OBSERVER_ID or state['observer_config_sha256'] != OBSERVER_SHA
                or state['host']['canonical_observer_id'] != OBSERVER_ID
                or state['host']['observer_config_sha256'] != OBSERVER_SHA):
            raise ValueError('State/host observer mismatch')
        for name in ['mergegalaxy.par', 'mergecube.par', 'GZWJL.PARAM']:
            if state['source_configuration_sha256'][name] != inputs['files'][name]['sha256']:
                raise ValueError('Run inputs do not match extracted snapshot: '+name)
        if not valid_sha(state['raw_original_snapshot_sha256']) or state['cache_sha256'] != sha(cp):
            raise ValueError('Missing raw identity or changed native cache')
        with np.load(cp, allow_pickle=False) as cache:
            all_ids, old = cache['id'], cache['is_old']
            if all_ids.dtype != np.uint64 or old.dtype != np.bool_ or old.shape != all_ids.shape:
                raise ValueError('Preserve uint64 IDs and the explicit old mask')
            ids = np.sort(all_ids[old])
            if len(ids) == 0 or len(np.unique(all_ids)) != len(all_ids):
                raise ValueError('Empty or duplicated member IDs')
            point = state['centers']['old_full']
            r, v = np.asarray(point['r_grf_kpc'], float), np.asarray(point['v_grf_kms'], float)
            require_close(cache['r_center_old_kpc'], r, 'cache/state COM position', 1e-8)
            require_close(cache['v_center_old_kms'], v, 'cache/state COM velocity', 1e-8)
            m = cache['mass_msun'][old].astype(float)
            dr, dv = cache['dr_native_kpc'][old], cache['dv_native_kms'][old]
            if not np.isfinite(m).all() or np.any(m <= 0):
                raise ValueError('Invalid old tracer mass')
            require_close(np.average(dr.astype(float)+r, weights=m, axis=0), r, 'cached old COM position', 2e-6)
            require_close(np.average(dv.astype(float)+v, weights=m, axis=0), v, 'cached old COM velocity', 2e-6)
        if point['N'] != len(ids):
            raise ValueError('COM/member count mismatch')
        id_sha = hashlib.sha256(ids.astype('<u8').tobytes()).hexdigest()
        cohort = state['cohort_file_sha256']
        if fixed_ids is None: fixed_ids, fixed_cohort = id_sha, cohort
        if (id_sha, cohort) != (fixed_ids, fixed_cohort):
            raise ValueError('Old tracer membership changed between epochs')
        method = dict(extractor_sha256=state['extractor_sha256'], source_modules=state['source_module_sha256'])
        if extraction_source is None: extraction_source = method
        if method != extraction_source:
            raise ValueError('Extraction/observer source changed between history frames')
        time = float(state['physical_time_Gyr'])
        if not np.isfinite(time) or not np.isfinite(np.r_[r, v]).all():
            raise ValueError('Invalid evolved native six-vector/time')
        frames.append(dict(snapshot=int(number), time_Gyr=time, r=r, v=v, state=state))
        sources.append(dict(snapshot=int(number), state_json=str(sp), state_sha256=sha(sp),
                            cache_sha256=state['cache_sha256'], raw_snapshot_sha256=state['raw_original_snapshot_sha256']))
    if len(frames) < 2 or frames[0]['snapshot'] != 0 or abs(frames[0]['time_Gyr']) > 1e-12:
        raise ValueError('History must include snapshot 0 and at least one evolved state')
    if any(b['time_Gyr'] <= a['time_Gyr'] for a, b in zip(frames, frames[1:])):
        raise ValueError('History times are not strictly increasing')
    context['fixed_old_ID_sha256'] = fixed_ids
    context['extraction_source'] = extraction_source
    if case.get('identity_kind') == 'content_addressed_legacy_run':
        if case.get('initial_raw_snapshot_sha256') != frames[0]['state']['raw_original_snapshot_sha256']:
            raise ValueError('Declared legacy run identity does not match actual snapshot 0')
    return dict(model=model, run_id=case['run_id'], index_path=str(index_path), index_sha256=sha(index_path),
                context=context, context_sha256=digest(context), frames=frames, inputs=inputs,
                initial=np.r_[frames[0]['r'], frames[0]['v']], sources=sources,
                physical_to_icrs=c['physical_to_icrs'],
                verification_scope='Frozen extracted state/cache SHA and same old IDs verified; raw HDF5 hashes inherited, not freshly re-read.')


def state_at(history, time):
    frames = history['frames']
    for frame in frames:
        if abs(time-frame['time_Gyr']) < 1e-12:
            return frame['r'], frame['v'], dict(snapshots=[frame['snapshot']], fraction=0.)
    for a, b in zip(frames, frames[1:]):
        if a['time_Gyr'] < time < b['time_Gyr']:
            if b['snapshot'] != a['snapshot']+1:
                raise ValueError('No adjacent snapshot support at the requested common time; extract both bracket endpoints')
            f = (time-a['time_Gyr'])/(b['time_Gyr']-a['time_Gyr'])
            return (1-f)*a['r']+f*b['r'], (1-f)*a['v']+f*b['v'], dict(snapshots=[a['snapshot'], b['snapshot']], fraction=f)
    raise ValueError('Requested time is outside actual evolution; extrapolation is forbidden')


def controls(config):
    f = config.get('feedback', {})
    if 'comparison_time_Gyr' not in f or 'search_time_window_Gyr' in f:
        raise ValueError('Freeze one comparison_time_Gyr for response plans/updates; time search is an evaluation stage')
    time = float(f['comparison_time_Gyr'])
    if not np.isfinite(time) or time <= 0:
        raise ValueError('Use one explicit positive comparison_time_Gyr for all six components and all trials')
    return f, time


def evaluate_history(history, config):
    f = config.get('feedback', {})
    if 'search_time_window_Gyr' in f:
        return search_time(history, config)
    f, time = controls(config)
    require_close(history['physical_to_icrs'], config['physical_to_icrs'], 'history/measurement mapping')
    result = dict(model=history['model'], run_id=history['run_id'], comparison_time_Gyr=time,
                  available_time_Gyr=[history['frames'][0]['time_Gyr'], history['frames'][-1]['time_Gyr']],
                  context_sha256=history['context_sha256'], index_sha256=history['index_sha256'],
                  evolved_6D_match_validated=False, verification_scope=history['verification_scope'])
    if not result['available_time_Gyr'][0] <= time <= result['available_time_Gyr'][1]:
        result['status'] = 'INSUFFICIENT_TIME_COVERAGE'
        return result
    r, v, support = state_at(history, time)
    result.update(evaluate_point(r, v, config), interpolation=support,
                  status='MEASURED_NO_ACCEPTANCE_TOLERANCES',
                  interpolation_scope='Linear native old_full COM r/v only between adjacent snapshots; interpolation error not estimated')
    if 'residual_scales' in f:
        scale = six(f['residual_scales'], 'residual_scales', True)
        result['scaled_residual_norm_squared'] = float(np.sum((result['residual_vector']/scale)**2))
        result['scaling_scope'] = 'Explicit engineering scales, not an observational covariance or statistical likelihood'
    if 'acceptance_tolerances' in f:
        tol = six(f['acceptance_tolerances'], 'acceptance_tolerances', True)
        passed = np.abs(result['residual_vector']) <= tol
        result.update(component_pass=dict(zip(RESIDUAL_KEYS, passed)),
                      status='ENGINEERING_6D_TOLERANCES_MET' if passed.all() else 'EVOLVED_6D_MISMATCH',
                      evolved_6D_match_validated=bool(passed.all()),
                      acceptance_scope='All six old_full COM components at one actual/interpolated time; PA and Walker fit not certified')
    return result


def search_time(history, config):
    """Select a single epoch by the complete six-residual engineering metric."""
    from copy import deepcopy
    from scipy.optimize import minimize_scalar
    f = config['feedback']
    if 'comparison_time_Gyr' in f:
        raise ValueError('Choose either a fixed comparison time or a six-dimensional time search')
    window = np.asarray(f['search_time_window_Gyr'], float)
    if window.shape != (2,) or not np.isfinite(window).all() or not 0 <= window[0] < window[1]:
        raise ValueError('Use an explicit increasing nonnegative search_time_window_Gyr')
    scale = six(f['residual_scales'], 'time-search residual_scales', True)
    frames = history['frames']
    if window[0] < frames[0]['time_Gyr'] or window[1] > frames[-1]['time_Gyr']:
        return dict(model=history['model'],run_id=history['run_id'],status='INSUFFICIENT_TIME_COVERAGE',
                    requested_search_time_window_Gyr=window,evolved_6D_match_validated=False,
                    available_time_Gyr=[frames[0]['time_Gyr'],frames[-1]['time_Gyr']])
    def score(t):
        r,v,_ = state_at(history, t)
        return float(np.sum((evaluate_point(r,v,config)['residual_vector']/scale)**2))
    candidates = [(score(float(t)),float(t)) for t in window]
    for a,b in zip(frames,frames[1:]):
        left,right = max(a['time_Gyr'],window[0]),min(b['time_Gyr'],window[1])
        if left >= right: continue
        if b['snapshot'] != a['snapshot']+1:
            raise ValueError('Time search cannot bridge missing adjacent snapshots')
        candidates.extend([(score(float(left)),float(left)),(score(float(right)),float(right))])
        fit=minimize_scalar(score,bounds=(left,right),method='bounded',options={'xatol':1e-10})
        if not fit.success: raise ValueError('Joint six-dimensional time search did not converge')
        candidates.append((float(fit.fun),float(fit.x)))
    _,time=min(candidates)
    chosen=deepcopy(config); del chosen['feedback']['search_time_window_Gyr']
    chosen['feedback']['comparison_time_Gyr']=time
    result=evaluate_history(history,chosen)
    result.update(time_selection='Joint six-residual engineering score over every supported adjacent interval',
                  requested_search_time_window_Gyr=window,
                  selected_time_on_window_boundary=bool(min(abs(time-window))<1e-8),
                  next_step='Freeze the selected common time explicitly before constructing response trials')
    return result


def require_feedback_support(history):
    identities = history['context']['fixed_input_sha256']
    if set(identities) != {'mw_ic', 'dwarf_ic', 'cube'} or not all(valid_sha(v) for v in identities.values()):
        raise ValueError('Feedback requires frozen hashes of the three fixed IC inputs; extract with feedback.freeze_physics_inputs=true')


def initial_merger(history, new_initial):
    """Only change initial placement/boost; body and cube transforms stay fixed."""
    from .solve_orbit_candidates import initial_biases
    inputs = history['inputs']
    galaxy, cube = np.asarray(inputs['mergegalaxy.par'][1], float), np.asarray(inputs['mergecube.par'][1], float)
    initial = history['frames'][0]
    biases = initial_biases(initial, galaxy, cube)
    values = galaxy.copy()
    values[:6] = cm.placement_for_target(new_initial[:3]/biases['length_to_kpc'], new_initial[3:]/biases['velocity_to_kms'],
        galaxy[9:12], cube[9:12], biases['cM'], biases['bM'], biases['cD'], biases['bD'])
    parsed = np.array([float(format(x, '.9g')) for x in values], dtype=np.float32)
    r, v = cm.initial_state_from_merge(parsed, cube, biases['cM'], biases['bM'], biases['cD'], biases['bD'])
    require_close(r*biases['length_to_kpc'], new_initial[:3], 'candidate linear initial position', 5e-4)
    require_close(v*biases['velocity_to_kms'], new_initial[3:], 'candidate linear initial velocity', 5e-4)
    return dict(parameters_12=values, input_names=inputs['mergegalaxy.par'][0],
                initial_r_grf_kpc=new_initial[:3], initial_v_grf_kms=new_initial[3:],
                note='Linear merger closure only; validate actual merged IC and new evolution before acceptance')


def signature(config):
    f, time = controls(config)
    return digest(dict(observer=OBSERVER_SHA, mapping=config['physical_to_icrs'],
                       target={k:config['target'][k] for k in TARGET_KEYS}, comparison_time_Gyr=time,
                       scales=f.get('residual_scales'), tolerances=f.get('acceptance_tolerances')))


def plan_feedback(history, config):
    require_feedback_support(history)
    f, _ = controls(config)
    steps = six(f['initial_steps'], 'initial_steps', True)
    six(f['residual_scales'], 'residual_scales', True)
    trials = []
    for axis, name in enumerate(CONTROL_KEYS):
        for sign in [-1, 1]:
            p = history['initial'].copy()
            p[axis] += sign*steps[axis]
            trials.append(dict(axis=axis, control=name, sign=sign, requested_initial=p,
                               merger=initial_merger(history, p), required_index=None))
    return dict(status='AWAITING_PERTURBED_EVOLUTIONS', model=history['model'],
                context_sha256=history['context_sha256'], feedback_signature=signature(config),
                baseline_index=history['index_path'], baseline_run_id=history['run_id'],
                baseline_initial=history['initial'], initial_steps=steps, trials=trials,
                simulations_submitted=False, evolved_6D_match_validated=False)


def propose_update(history, pairs, config):
    """Central responses from the *same model*, never from two different dwarfs."""
    require_feedback_support(history)
    f, _ = controls(config)
    scale = six(f['residual_scales'], 'residual_scales', True)
    steps = six(f['initial_steps'], 'initial_steps', True)
    damping, radius = float(f['damping']), float(f['trust_radius'])
    if not np.isfinite([damping, radius]).all() or damping < 0 or not 0 < radius <= 1:
        raise ValueError('Require damping>=0 and 0<trust_radius<=1 in units of the measured initial steps')
    baseline = evaluate_history(history, config)
    if 'residual_vector' not in baseline:
        raise ValueError('Baseline has no evolution at the requested common time')
    if len(pairs) != 6 or {p['axis'] for p in pairs} != set(range(6)):
        raise ValueError('Six distinct central-difference pairs are required; do not invent missing responses')
    dp, dy, provenance, used_runs = [], [], [], {history['run_id']}
    for pair in sorted(pairs, key=lambda p:p['axis']):
        axis, minus, plus = pair['axis'], pair['minus'], pair['plus']
        residuals = []
        for sign, trial in [(-1, minus), (1, plus)]:
            require_feedback_support(trial)
            if trial['context_sha256'] != history['context_sha256'] or trial['model'] != history['model']:
                raise ValueError('Response trials differ in model, members, body orientation, fixed inputs or physics')
            if trial['run_id'] in used_runs:
                raise ValueError('Distinct physical run IDs are required for baseline and each perturbation')
            used_runs.add(trial['run_id'])
            expected = history['initial'].copy(); expected[axis] += sign*steps[axis]
            require_close(trial['initial'], expected, 'actual perturbed initial six-vector', 5e-4)
            result = evaluate_history(trial, config)
            if 'residual_vector' not in result:
                raise ValueError('A perturbation has not evolved to the common comparison time')
            residuals.append(result['residual_vector'])
            provenance.append(dict(axis=axis, sign=sign, run_id=trial['run_id'], index_sha256=trial['index_sha256']))
        dp.append((plus['initial']-minus['initial'])/(2*steps))
        dy.append((residuals[1]-residuals[0])/(2*scale))
    # Account for actual parsed/measured IC perturbations, including small
    # cross-axis numerical offsets, rather than pretending the parser is exact.
    dp, dy = np.array(dp).T, np.array(dy).T
    if np.linalg.matrix_rank(dp) < 6:
        raise ValueError('Measured initial-condition perturbations are rank deficient')
    jacobian = np.linalg.solve(dp.T, dy.T).T
    singular = np.linalg.svd(jacobian, compute_uv=False)
    maximum_condition = float(f.get('maximum_condition_number', 1e8))
    if not np.isfinite(maximum_condition) or maximum_condition < 1:
        raise ValueError('Invalid maximum_condition_number')
    if singular[-1] <= singular[0]/maximum_condition or singular[0] == 0:
        raise ValueError('Evolved response cannot constrain all six IC controls; obtain independent evolutions')
    r0 = baseline['residual_vector']/scale
    system = np.vstack([jacobian, np.sqrt(damping)*np.eye(6)])
    rhs = np.r_[-r0, np.zeros(6)]
    delta = np.linalg.lstsq(system, rhs, rcond=None)[0]
    if np.linalg.norm(delta) > radius: delta *= radius/np.linalg.norm(delta)
    predicted = r0+jacobian@delta
    if np.linalg.norm(predicted) >= np.linalg.norm(r0):
        raise ValueError('No predicted local improvement; require new response measurements')
    initial = history['initial']+delta*steps
    return dict(status='AWAITING_FORWARD_EVOLUTION', model=history['model'],
                context_sha256=history['context_sha256'], feedback_signature=signature(config),
                baseline_index_sha256=history['index_sha256'], baseline_run_id=history['run_id'],
                comparison_time_Gyr=baseline['comparison_time_Gyr'], control_order=CONTROL_KEYS,
                residual_order=RESIDUAL_KEYS, residual_scales=scale, initial_steps=steps,
                baseline_initial=history['initial'], proposed_initial=initial, delta_initial=delta*steps,
                baseline_residual=baseline['residual_vector'], predicted_residual=predicted*scale,
                scaled_jacobian=jacobian, response_singular_values=singular,
                measured_initial_design=dp, perturbation_sources=provenance,
                trust_radius=radius, damping=damping, merger=initial_merger(history, initial),
                evolved_6D_match_validated=False, simulations_submitted=False,
                prediction_scope='Local response prediction, not an evolved measurement; Sun and final-state data unchanged')


def verify_proposal(history, measurement, proposal, config):
    if proposal.get('status') != 'AWAITING_FORWARD_EVOLUTION' or proposal['model'] != history['model']:
        raise ValueError('Not a forward-evolution proposal for this model')
    if (proposal['context_sha256'] != history['context_sha256']
            or proposal['feedback_signature'] != signature(config)):
        raise ValueError('Proposal/verification physics, observer, target, comparison time or scales changed')
    if history['run_id'] == proposal['baseline_run_id']:
        raise ValueError('Verification must use a new physical evolution')
    require_close(history['initial'], proposal['proposed_initial'], 're-evolved proposed initial six-vector', 5e-4)
    if 'residual_vector' not in measurement:
        return dict(status='INSUFFICIENT_TIME_COVERAGE', next_update_allowed=False)
    scale = six(proposal['residual_scales'], 'proposal scales', True)
    old = float(np.sum((six(proposal['baseline_residual'], 'baseline residual')/scale)**2))
    predicted = float(np.sum((six(proposal['predicted_residual'], 'predicted residual')/scale)**2))
    actual = float(np.sum((measurement['residual_vector']/scale)**2))
    if old <= predicted: raise ValueError('Proposal has no predicted improvement')
    ratio = (old-actual)/(old-predicted)
    return dict(status='ACTUAL_IMPROVEMENT' if actual < old else 'REJECTED_NO_ACTUAL_IMPROVEMENT',
                actual_to_predicted_reduction_ratio=ratio, previous_score=old, actual_score=actual,
                accept_for_next_iteration=bool(actual < old),
                recommended_next_trust_radius=proposal['trust_radius']*(.5 if ratio < .25 else 1.),
                evolved_6D_match_validated=measurement['evolved_6D_match_validated'])


def run_feedback(args, config):
    if args.index is None or not args.model or args.run_directory is None:
        raise ValueError('Feedback stages require --index, --model and a fresh --run-directory')
    if args.stage != 'evaluate' and not args.physical_candidates:
        raise ValueError('Initial-condition proposals require --physical-candidates')
    history = load_history(args.index, args.model)
    measurement = evaluate_history(history, config)
    if args.stage == 'evaluate':
        result = measurement
        if args.proposal:
            result['proposal_verification'] = verify_proposal(history, measurement, load(args.proposal), config)
    elif args.stage == 'plan-feedback':
        result = plan_feedback(history, config)
    else:
        if args.feedback_spec is None:
            raise ValueError('update-feedback requires --feedback-spec with six actual central-difference pairs')
        spec = load(args.feedback_spec)
        if spec.get('schema_id') != 'canonical_orbit_feedback_trials_v1':
            raise ValueError('Unknown feedback trial specification')
        def trial(ref):
            p = Path(ref['index'])
            if not p.is_absolute(): p = args.feedback_spec.resolve().parent/p
            h = load_history(p, args.model)
            if h['index_sha256'] != ref['sha256']:
                raise ValueError('Trial index SHA differs from feedback specification')
            return h
        pairs = [dict(axis=p['axis'], minus=trial(p['minus']), plus=trial(p['plus'])) for p in spec['pairs']]
        result = propose_update(history, pairs, config)
        result['trial_spec_sha256'] = sha(args.feedback_spec)
    result.update(schema_id='canonical_orbit_feedback_result_v1', canonical_observer_id=OBSERVER_ID,
                  observer_config_sha256=OBSERVER_SHA, input_index_sha256=history['index_sha256'],
                  frozen_state_sources=history['sources'], configuration_sha256=sha(args.config),
                  program_sha256=sha(__file__))
    out = args.run_directory.resolve()
    out.mkdir(parents=True, exist_ok=False)
    def save(name, value):
        (out/name).write_text(json.dumps(plain(value), ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    save('FEEDBACK_RESULT.json', result)
    rows = []
    for frame in history['frames']:
        point = evaluate_point(frame['r'], frame['v'], config)
        rows.append(dict(snapshot=frame['snapshot'], time_Gyr=frame['time_Gyr'],
                         **point['observables'], **{'delta_'+k:v for k,v in point['residual'].items()},
                         sky_separation_deg=point['sky_separation_deg'],
                         GC_velocity_direction_error_deg=point['GC_velocity_direction_error_deg']))
    with (out/'EVOLVED_6D_RESIDUALS.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    from .solve_orbit_candidates import write_par
    if args.stage == 'plan-feedback':
        for item in result['trials']:
            write_par(out/('axis%d_%s.CANDIDATE.par' % (item['axis'], 'plus' if item['sign']>0 else 'minus')),
                      item['merger']['input_names'], item['merger'])
    elif args.stage == 'update-feedback':
        write_par(out/'next_initial.CANDIDATE.par', result['merger']['input_names'], result['merger'])
    print(json.dumps(dict(status=result['status'], model=args.model, output=str(out),
                         evolved_6D_match_validated=result['evolved_6D_match_validated'])))
    return 0
