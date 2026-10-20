"""One explicit observing contract and one set of native six-vector inputs."""
from pathlib import Path
import argparse
import copy
import hashlib
import json
import re
import shlex
import subprocess
import sys

import numpy as np

from . import SCHEMA_ID
from .target_frame import INPUT_FRAME, resolve_target

OBSERVER_ID = 'astropy_native_direct_v1'
OBSERVER_SHA = 'e7fd3ff559ba72ba8fed4bc4b149a403cbdd137a4ed4a7b986e5ee357ae9fde2'
SUN = np.array([-np.sqrt(8.122 ** 2 - .0208 ** 2), 0., .0208])
SUN_VELOCITY = np.array([12.9, 245.6, 7.78])
MODEL = re.compile(r'Fornax\d+(?:_[A-Za-z0-9_.-]+)?\Z')
LABEL = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]*\Z')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def require_close(actual, expected, name, tolerance=1e-10):
    a, b = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    if a.shape != b.shape or not np.all(np.isfinite(a)) or np.max(np.abs(a - b)) > tolerance:
        raise ValueError('Incompatible canonical ' + name)


def validate_config(config, require_cases=True):
    """Reject legacy inputs rather than silently translating an old observer."""
    c = copy.deepcopy(config)
    if c.get('schema_id') != SCHEMA_ID:
        raise ValueError('Use canonical_orbit_recalibration_v1; old frame/camera/CSV configs require a new explicit configuration.')
    if c.get('canonical_observer_id') != OBSERVER_ID or c.get('observer_config_sha256') != OBSERVER_SHA:
        raise ValueError('Orbit proposals require the fixed astropy_native_direct_v1 observer identity.')
    require_close(c.get('sun_position_native_kpc'), SUN, 'Sun position')
    require_close(c.get('sun_velocity_native_kms'), SUN_VELOCITY, 'Sun velocity')
    matrix = np.asarray(c.get('physical_to_icrs'), dtype=float)
    if (matrix.shape != (3, 3) or not np.all(np.isfinite(matrix))
            or not np.allclose(matrix.T @ matrix, np.eye(3), atol=1e-11, rtol=0)
            or abs(np.linalg.det(matrix) - 1) > 1e-11):
        raise ValueError('physical_to_icrs must be the verified proper canonical mapping.')
    target = c.get('target', {})
    if target.get('input_frame') != INPUT_FRAME:
        raise ValueError('Declare target.input_frame=icrs_heliocentric; GC Table 2 values cannot be silently reinterpreted.')
    required = ['ra_deg', 'dec_deg', 'D_kpc', 'vlos_kms', 'mu_alpha_cosdec_masyr',
                'mu_delta_masyr', 'morphology_PA_E_of_N_deg']
    if any(key not in target or not np.isfinite(target[key]) for key in required):
        raise ValueError('An explicit finite ICRS six-vector and morphological PA target are required.')
    if not (0 <= target['ra_deg'] < 360 and -90 <= target['dec_deg'] <= 90 and target['D_kpc'] > 0
            and 0 <= target['morphology_PA_E_of_N_deg'] < 180):
        raise ValueError('Invalid target distance, declination or morphological PA.')
    if require_cases and ('cases' not in c or not c['cases']):
        raise ValueError('At least one explicit model case is required.')
    c.setdefault('cases',[])
    seen = set()
    for case in c['cases']:
        if not MODEL.fullmatch(case.get('model', '')) or case['model'] in seen:
            raise ValueError('Use unique full Fornax model IDs.')
        seen.add(case['model'])
        for key in ['initial', 'latest']:
            if not isinstance(case.get(key), int) or isinstance(case[key], bool) or case[key] < 0:
                raise ValueError('Explicit non-negative initial/latest snapshot numbers are required.')
        selected = {case['initial'], case['latest']}
        for key in ['GC_radius_bracket', 'crossing_bracket']:
            pair = case.get(key)
            if (not isinstance(pair, list) or len(pair) != 2
                    or any(not isinstance(n, int) or isinstance(n, bool) or n < 0 for n in pair)
                    or pair[1] != pair[0] + 1):
                raise ValueError('Use explicit adjacent snapshot brackets: ' + key)
            selected.update(pair)
        if case['initial'] >= min(case['GC_radius_bracket']) or case['latest'] < max(selected):
            raise ValueError('Initial/anchor/latest snapshots are out of order.')
        requested = case.get('extract_snapshots', sorted(selected))
        if any(not isinstance(n, int) or isinstance(n, bool) or n < 0 for n in requested):
            raise ValueError('Invalid extraction snapshot list.')
        case['extract_snapshots'] = sorted(selected | set(requested))
    return c


def verify_observation_source(source, config):
    """Cross-check the complete target through the actual shared Astropy entry."""
    sys.path.insert(0, str(source))
    from canonical_observer import production_observables, observer_contract_sha256
    if observer_contract_sha256() != OBSERVER_SHA:
        raise ValueError('Configured observing implementation has a different contract hash.')
    target = config['target']
    conversion = resolve_target(target)
    matrix = np.asarray(config['physical_to_icrs'])
    a, d = np.radians([target['ra_deg'], target['dec_deg']])
    ray = np.array([np.cos(a) * np.cos(d), np.sin(a) * np.cos(d), np.sin(d)])
    east = np.array([-np.sin(a), np.cos(a), 0.])
    north = np.cross(ray, east)
    r = np.array(conversion['native_state']['r_grf_native_kpc'])
    v = np.array(conversion['native_state']['v_grf_native_kms'])
    # Independent algebraic inversion checks the supplied mapping as well as
    # the full Astropy transformation. No model/observer rotation is fitted.
    require_close(r,SUN+target['D_kpc']*(matrix.T@ray),'ICRS inverse position',1e-8)
    require_close(v,SUN_VELOCITY+matrix.T@(target['vlos_kms']*ray+4.740470463533349*target['D_kpc']
        *(target['mu_alpha_cosdec_masyr']*east+target['mu_delta_masyr']*north)),'ICRS inverse velocity',1e-8)
    result = production_observables(r[None], v[None], np.array([1], np.uint64))
    require_close(result['_native']['physical_to_icrs'], matrix, 'physical mapping', 1e-10)
    for key, expected in [('ra', target['ra_deg']), ('dec', target['dec_deg']), ('rh', target['D_kpc']),
                          ('vlos', target['vlos_kms']), ('pmra', target['mu_alpha_cosdec_masyr']),
                          ('pmdec', target['mu_delta_masyr'])]:
        require_close(float(result[key][0]), expected, 'target ' + key, 1e-8)
    for key, value in [('r_grf_native_kpc', r), ('v_grf_native_kms', v)]:
        if key in target:
            require_close(target[key], value, key, 1e-8)
        target[key] = value.tolist()
    target.update(R_gc_kpc=float(np.linalg.norm(r)), V3D_kms=float(np.linalg.norm(v)),
                  Vr_gc_kms=float(r @ v / np.linalg.norm(r)),
                  Vt_gc_kms=float(np.linalg.norm(np.cross(r, v)) / np.linalg.norm(r)))
    config['target_conversion_contract']=conversion


def capture_inputs(root, cases):
    out = {'models': {}}
    for case in cases:
        model = case['model']
        directory = root / model
        record = {'files': {}}
        for name in ['mergegalaxy.par', 'mergecube.par', 'GZWJL.PARAM']:
            path = directory / name
            record['files'][name] = {'text': path.read_text(encoding='utf-8'), 'sha256': sha(path)}
        for name in ['mergegalaxy.par', 'mergecube.par']:
            rows = [shlex.split(line, comments=True) for line in record['files'][name]['text'].splitlines()]
            rows = [row for row in rows if row]
            if len(rows) < 2 or len(rows[0]) != 2 or len(rows[1]) != 12:
                raise ValueError('Invalid one-pass merge input: ' + str(directory / name))
            if not np.all(np.isfinite([float(x) for x in rows[1]])):
                raise ValueError('Non-finite merge transform')
            record[name] = rows
        out['models'][model] = record
    return out


def latest_complete_snapshot(directory):
    import h5py
    available = []
    for path in directory.glob('snapshot_*.hdf5'):
        match = re.fullmatch(r'snapshot_(\d+)\.hdf5', path.name)
        if match: available.append((int(match[1]), path))
    for number, path in sorted(available, reverse=True):
        try:
            before = path.stat()
            with h5py.File(path, 'r') as stream:
                counts = stream['Header'].attrs['NumPart_ThisFile']
                for part, count in enumerate(counts):
                    if count:
                        for key in ['Coordinates', 'Velocities', 'ParticleIDs']:
                            if len(stream['PartType%d' % part][key]) != count:
                                raise ValueError('Incomplete particle dataset')
            after = path.stat()
            if (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns):
                return number
        except (OSError, KeyError, ValueError):
            continue
    raise ValueError('No complete, closed single-file snapshot found: ' + str(directory))


def pericentre_status(index):
    """Report sampled states without extrapolating an unobserved turning point."""
    records = []
    for case in index['config']['cases']:
        rows = []
        for number, ref in sorted(index['states'][case['model']].items(), key=lambda item: int(item[0])):
            state = load(ref['state_json'])
            point = state['centers']['old_full']
            r, v = np.asarray(point['r_grf_kpc']), np.asarray(point['v_grf_kms'])
            rows.append(dict(snapshot=int(number), R_gc_kpc=float(np.linalg.norm(r)),
                             Vr_gc_kms=float(r @ v / np.linalg.norm(r)),
                             physical_time_Gyr=state['physical_time_Gyr']))
        turns = [[a['snapshot'], b['snapshot']] for a, b in zip(rows, rows[1:])
                 if a['Vr_gc_kms'] < 0 <= b['Vr_gc_kms']]
        last = rows[-1]
        result = dict(model=case['model'], sampled_states=rows, latest_complete=last,
                      latest_inward=last['Vr_gc_kms'] < 0,
                      sampled_radial_turn_brackets=turns,
                      pericentre_radius_measured=False,
                      status='RADIAL_TURN_NEEDS_REFINEMENT' if turns else
                             ('LATEST_COMPLETE_STILL_INFALLING' if last['Vr_gc_kms'] < 0 else 'NO_TURN_BRACKET'),
                      note='Sparse sampled history; the last sampled minimum is not a measured pericentre. No potential or future orbit was extrapolated.')
        records.append(result)
        directory = Path(index['states'][case['model']][str(last['snapshot'])]['state_json']).parent
        path = directory / 'PERICENTRE_STATUS.json'
        if path.exists():
            if load(path) != result: raise ValueError('Existing pericentre status differs')
        else:
            path.write_text(json.dumps(result, indent=2), encoding='utf-8')
    return records


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--workspace-root', type=Path)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--run-directory', type=Path)
    parser.add_argument('--tag')
    parser.add_argument('--stage', choices=['validate', 'frame-check', 'extract', 'solve', 'all'], default='all')
    parser.add_argument('--physical-candidates', action='store_true',
                        help='Explicitly request unevolved physical orbit/body candidates; never a reference-frame correction.')
    args = parser.parse_args(argv)
    config = validate_config(load(args.config),require_cases=args.stage not in ['validate','frame-check'])
    source = args.source.resolve() if args.source else Path(__file__).resolve().parent.parent
    verify_observation_source(source, config)
    if args.stage in ['validate','frame-check']:
        print(json.dumps(dict(state='PASS',canonical_observer_id=OBSERVER_ID,
                              target=config['target'],conversion=config['target_conversion_contract']),indent=2))
        return 0
    if args.stage in ['solve','all'] and not args.physical_candidates:
        parser.error('Candidate solving changes physical initial conditions. Use --physical-candidates explicitly; use --stage frame-check for a frame correction check.')
    if args.workspace_root is None or args.run_directory is None or not args.tag:
        parser.error('--workspace-root, --run-directory and --tag are required for execution')
    if not LABEL.fullmatch(args.tag):
        raise ValueError('Use a simple unique analysis tag')
    root, work = args.workspace_root.resolve(), args.run_directory.resolve()
    if not work.is_dir():
        work.mkdir(parents=True)
    config['analysis_source'] = str(source)
    config['workspace_root'] = str(root)
    config['run_tag'] = args.tag
    inputs = capture_inputs(root, config['cases'])
    package = Path(__file__).resolve().parent
    index_path = work / 'frame_state_index.json'
    if args.stage in ['extract', 'all']:
        index = dict(config=config, server_inputs=inputs, states={},
                     program_source_sha256={str(p.relative_to(package)): sha(p) for p in package.rglob('*.py')
                                            if '__pycache__' not in p.parts})
        for case in config['cases']:
            model = case['model']
            if config.get('inspect_latest_complete', True):
                number = latest_complete_snapshot(root / model / 'output')
                case['current_complete_snapshot'] = number
                case['extract_snapshots'] = sorted(set(case['extract_snapshots']) | {number, max(0, number - 1)})
            case['output_dir'] = str(root / model / 'output/analysis' / ('snapshot_%03d' % case['latest']) / args.tag)
            index['states'][model] = {}
            for number in case['extract_snapshots']:
                subprocess.run([sys.executable, '-B', '-m', 'orbit_calibration_v1.extract_frame', '--workspace-root', str(root),
                                '--source', str(source), '--tag', args.tag, '--model', model,
                                '--snapshot', str(number)], cwd=str(package.parent), check=True)
                directory = root / model / 'output/analysis' / ('snapshot_%03d' % number) / args.tag
                index['states'][model][str(number)] = dict(state_json=str(directory / 'state.json'), cache_npz=str(directory / 'stars.npz'))
                (work / 'progress.json').write_text(json.dumps(dict(model=model, snapshot=number, state='extracted'), indent=2))
        if index_path.exists():
            if load(index_path) != index:
                raise ValueError('Existing index differs; choose a new analysis tag/directory.')
        else:
            with index_path.open('x') as stream:
                json.dump(index, stream, indent=2)
    else:
        index = load(index_path)
        if index['config']['canonical_observer_id'] != OBSERVER_ID:
            raise ValueError('Existing index has a different observer')
        dynamic = {'analysis_source', 'workspace_root', 'run_tag'}
        original = copy.deepcopy({k: v for k, v in index['config'].items() if k not in dynamic})
        for case in original['cases']:
            for key in ['output_dir', 'current_complete_snapshot']:
                case.pop(key, None)
        for a, b in zip(original['cases'], config['cases']):
            a['extract_snapshots'] = b['extract_snapshots']
        requested = {k: v for k, v in config.items() if k not in dynamic}
        if original != requested:
            raise ValueError('Configuration differs from the frozen extraction index')
        program = {str(p.relative_to(package)): sha(p) for p in package.rglob('*.py') if '__pycache__' not in p.parts}
        if program != index['program_source_sha256']:
            raise ValueError('Orbit program changed since extraction; use a new tag')
        for case in config['cases']:
            for name, record in inputs['models'][case['model']]['files'].items():
                if record['sha256'] != index['server_inputs']['models'][case['model']]['files'][name]['sha256']:
                    raise ValueError('Source inputs changed since extraction')
        for refs in index['states'].values():
            for ref in refs.values():
                state = load(ref['state_json'])
                if state['cache_sha256'] != sha(ref['cache_npz']) or state['raw_original_snapshot_sha256'] != sha(state['original_snapshot_path']):
                    raise ValueError('Snapshot or frozen cache changed since extraction')
    statuses = pericentre_status(index)
    status_path = work / 'PERICENTRE_STATUS.json'
    if not status_path.exists():
        status_path.write_text(json.dumps(statuses, indent=2), encoding='utf-8')
    if args.stage in ['solve', 'all']:
        subprocess.run([sys.executable, '-B', '-m', 'orbit_calibration_v1.solve_orbit_candidates', '--index', str(index_path)],
                       cwd=str(package.parent), check=True)
    print('PASS: requested canonical orbit-correction stages complete; proposals are not re-evolved simulations.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
