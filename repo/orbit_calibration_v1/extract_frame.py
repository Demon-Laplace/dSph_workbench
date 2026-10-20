"""Read one native snapshot through the deployed observing implementation."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import re
import sys

import numpy as np

SCHEMA_ID = 'canonical_orbit_state_v1'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def plain(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {k: plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [plain(v) for v in value]
    return value


def save(path, value):
    with Path(path).open('x') as f:
        json.dump(plain(value), f, indent=2)
        f.write('\n')


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--workspace-root', required=True, type=Path)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--tag', required=True)
    parser.add_argument('--model', required=True)
    parser.add_argument('--snapshot', required=True, type=int)
    args = parser.parse_args(argv)
    if not re.fullmatch(r'Fornax\d+(?:_[A-Za-z0-9_.-]+)?', args.model):
        raise ValueError('Use a complete Fornax model ID')
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', args.tag):
        raise ValueError('Use a simple, unique analysis tag')
    if args.snapshot < 0:
        raise ValueError('Snapshot number must be non-negative')
    ROOT, SOURCE, TAG = args.workspace_root.resolve(), args.source.resolve(), args.tag
    assert (SOURCE / 'canonical_observer.py').is_file()
    model = ROOT / args.model
    output = model / 'output/analysis' / ('snapshot_%03d' % args.snapshot) / TAG
    source_snapshot = model / 'output' / ('snapshot_%03d.hdf5' % args.snapshot)
    cohort = model / 'output/analysis/astropy_native_direct_v1/cohort_ids.npz'
    assert source_snapshot.is_file()
    expected = dict(schema_id=SCHEMA_ID, extractor_sha256=sha(__file__),
                    source_module_sha256={p.name: sha(p) for p in sorted(SOURCE.glob('*.py'))},
                    source_configuration_sha256={n: sha(model / n) for n in ['mergegalaxy.par', 'mergecube.par', 'GZWJL.PARAM']},
                    cohort_file_sha256=sha(cohort))
    if (output / 'state.json').exists():
        state = json.loads((output / 'state.json').read_text())
        assert state['state'] == 'PASS'
        for key, value in expected.items():
            assert state.get(key) == value, 'Stale extraction identity: ' + key
        assert state['raw_original_snapshot_sha256'] == sha(source_snapshot)
        assert state['cache_sha256'] == sha(output / 'stars.npz')
        print('REUSED VERIFIED EXTRACTION', args.model, args.snapshot, flush=True)
        raise SystemExit(0)
    output.mkdir(parents=True, exist_ok=False)
    for env, subdir in [('TMPDIR', '.tmp'), ('MPLCONFIGDIR', '.mpl'), ('XDG_CACHE_HOME', '.cache')]:
        p = output / subdir
        p.mkdir()
        os.environ[env] = str(p)
    os.environ['DSPH_CANONICAL_OUTPUT_ROOT'] = str(output / 'canonical')
    os.environ['DSPH_PREEXISTING_ID_PATH'] = str(cohort)
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    os.chdir(model)
    sys.path.insert(0, str(SOURCE))
    import variable
    from snapshot_context import prepare_snapshot_context
    from analysis_core import Analysis
    from canonical_observer import (CANONICAL_OBSERVER_ID, observer_contract_sha256,
                                    production_observables, production_sky)

    before = sha(source_snapshot)
    snapshot = prepare_snapshot_context(str(model / 'output'), args.snapshot,
        core_radius=variable.fornax_core_radius, include_mw_gas=False,
        include_dark_matter=False, include_star_birth=True)
    coords = snapshot['star_coords']
    native = coords['_native']
    old = np.asarray(snapshot['fixed_old_local_mask'], bool)
    mass = snapshot['df'].loc[snapshot['total_dw_star_mask'], 'm'].to_numpy(dtype=float)
    position = native['host_relative_position_kpc']
    velocity = native['grf_velocity_native_kms']
    assert len(position) == len(mass) == len(old)


    def center(mask):
        return dict(r_grf_kpc=np.average(position[mask], axis=0, weights=mass[mask]),
                    v_grf_kms=np.average(velocity[mask], axis=0, weights=mass[mask]),
                    N=int(mask.sum()), mass_msun=float(mass[mask].sum(dtype=np.float64)))


    old_center = center(old)
    r0, v0 = old_center['r_grf_kpc'], old_center['v_grf_kms']
    inside = np.linalg.norm(position - r0, axis=1) < 5
    centers = dict(old_full=old_center, current_full=center(np.ones(len(old), bool)),
                   old_core5=center(old & inside), current_core5=center(inside))
    x = snapshot['x_kpc'][old]
    y = snapshot['y_kpc'][old]
    shape = Analysis.calculate_projected_shape(x, y, mass=mass[old], n_neighbors=30)
    rhalf = Analysis.half_light_radius(x, y, mass[old], shape['eps'], pa=shape['pa'],
                                       center_x=shape['center_x'], center_y=shape['center_y'])
    shape_state = dict(PA_deg=float((90 - np.degrees(shape['pa'])) % 180),
                       q=float(shape['q']), pa_math_rad=float(shape['pa']),
                       rhalf_major_kpc=float(rhalf), center_x_kpc=shape['center_x'],
                       center_y_kpc=shape['center_y'], definition='Frozen old tracers; mass times 30-neighbour local-density weighted covariance, native canonical TAN.')
    projected = coords['_sky']
    baseline = dict(canonical_observer_id=CANONICAL_OBSERVER_ID,
                    observer_config_sha256=observer_contract_sha256(),
                    RA_deg=projected['center_RA_deg'], Dec_deg=projected['center_Dec_deg'],
                    D_hel_kpc=projected['center_distance_kpc'],
                    R_gc_kpc=float(np.linalg.norm(r0)),
                    Walker_R5={key: value for key, value in snapshot['walker_R5'].items() if key != 'fit_particle_id'},
                    motion_R5=snapshot['motion_R5'])
    # Store small offsets in float32, centres in float64; retain the original ID binding.
    dr = (position - r0).astype(np.float32)
    dv = (velocity - v0).astype(np.float32)
    cached_mass = mass.astype(np.float32)


    def check_cached_observation(dr_cache, dv_cache, mass_cache):
        p = dr_cache.astype(float) + r0
        v = dv_cache.astype(float) + v0
        c = production_observables(p, v, np.asarray(coords['id'], np.uint64))
        s = production_sky(c, r0 - native['observer']['sun_position_kpc'])
        xs, ys = s['east_kpc'][old], s['north_kpc'][old]
        sh = Analysis.calculate_projected_shape(xs, ys, mass=mass_cache[old], n_neighbors=30)
        half = Analysis.half_light_radius(xs, ys, mass_cache[old], sh['eps'], pa=sh['pa'],
                                         center_x=sh['center_x'], center_y=sh['center_y'])
        angle = (90 - np.degrees(sh['pa'])) % 180
        result = dict(position_max_error_kpc=float(np.max(abs(p - position))),
                      velocity_max_error_kms=float(np.max(abs(v - velocity))),
                      LOS_max_error_kms=float(np.max(abs(c['vlos'] - coords['vlos']))),
                      R5_ID_mask_differences=int(np.count_nonzero(s['fit_R5_mask'] != projected['fit_R5_mask'])),
                      PA_error_deg=float(abs((angle - shape_state['PA_deg'] + 90) % 180 - 90)),
                      q_error=float(abs(sh['q'] - shape['q'])),
                      rhalf_relative_error=float(abs(half / rhalf - 1)),
                      cached_old_center_error_kpc=float(np.max(abs(np.average(p[old], axis=0, weights=mass_cache[old]) - r0))))
        result['pass'] = (result['position_max_error_kpc'] < 1e-3 and result['velocity_max_error_kms'] < 1e-3
                          and result['LOS_max_error_kms'] < 1e-3 and result['R5_ID_mask_differences'] == 0
                          and result['PA_error_deg'] < 1e-3 and result['q_error'] < 1e-5
                          and result['rhalf_relative_error'] < 1e-5 and result['cached_old_center_error_kpc'] < 1e-6)
        return result


    cache_regression = check_cached_observation(dr, dv, cached_mass)
    float32_attempt = dict(cache_regression)
    if not cache_regression['pass']:
        # Preserve existing native precision only when the measured scientific gates require it.
        dr, dv, cached_mass = position - r0, velocity - v0, mass
        cache_regression = check_cached_observation(dr, dv, cached_mass)
    assert cache_regression['pass'], cache_regression
    cache = output / 'stars.npz'
    np.savez_compressed(cache, id=np.asarray(coords['id'], np.uint64), mass_msun=cached_mass,
                        is_old=old, dr_native_kpc=dr, dv_native_kms=dv,
                        r_center_old_kpc=np.asarray(r0, np.float64),
                        v_center_old_kms=np.asarray(v0, np.float64))
    params = {}
    for line in (model / 'GZWJL.PARAM').read_text().splitlines():
        fields = line.split('%', 1)[0].split()
        if len(fields) >= 2: params[fields[0]] = fields[1]
    seconds = float(params['UnitLength_in_cm']) / float(params['UnitVelocity_in_cm_per_s'])
    gyr_per_unit = seconds / (365.25 * 86400 * 1e9)
    state = dict(expected, state='PASS', model=args.model, numsp=args.snapshot,
                  header_time=float(snapshot['tsnap']), unit_time_to_Gyr=gyr_per_unit,
                  physical_time_Gyr=float(snapshot['tsnap']) * gyr_per_unit,
                  host=snapshot['observer_host_state'], centers=centers,
                  shape_old_full=shape_state, baseline_observation=baseline,
                  cache_storage=dict(offset_dtype=str(dr.dtype), mass_dtype=str(cached_mass.dtype),
                                      centre_dtype='float64', id_dtype='uint64'),
                  cache_observing_regression=cache_regression, float32_cache_attempt=float32_attempt,
                  cache_position_reconstruction_max_error_kpc=float(np.max(abs(dr.astype(float) + r0 - position))),
                  cache_velocity_reconstruction_max_error_kms=float(np.max(abs(dv.astype(float) + v0 - velocity))),
                  cache_sha256=sha(cache), raw_original_snapshot_sha256=before,
                  original_snapshot_path=str(source_snapshot),
                  selected_id_sha256=hashlib.sha256(np.asarray(coords['id'], dtype='<u8').tobytes()).hexdigest(),
                  canonical_observer_id=CANONICAL_OBSERVER_ID,
                  observer_config_sha256=observer_contract_sha256())
    assert sha(source_snapshot) == before
    save(output / 'state.json', state)
    print('EXTRACTION PASS', args.model, args.snapshot, 'PA', shape_state['PA_deg'], flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
