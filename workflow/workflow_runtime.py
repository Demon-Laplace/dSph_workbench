"""Configured analysis runs, Slurm submission, and result verification."""
import ast
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone


def now():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')
    temporary.replace(path)


def confined(root, path):
    path = Path(path)
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f'Write would escape configured workspace: {path}')
    return path


def label(value):
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', value):
        raise ValueError(f'Invalid name: {value!r}')
    return value


def config(root):
    path = root / 'workflow.json'
    return json.loads(path.read_text()) if path.exists() else {}


def model_directory(root, model):
    if not re.fullmatch(r'(?:Fornax|MW|Sculptor)\d+[A-Za-z0-9_.-]*', model):
        raise ValueError('Use a full model ID, e.g. Fornax7003')
    base = (root / config(root).get('models_root', '.')).resolve()
    return confined(base, base / model)


def code_dir(args):
    setting = getattr(args, 'code', None) or config(args.root).get('analysis_code', 'tools/server-analysis')
    return (args.root / setting).resolve()


def run_dir(args):
    return confined(args.root, args.root / 'workspace' / 'analysis' / label(args.name))


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


ANALYSIS_SOURCE_PACKAGES = ('orbit_calibration_v1',)
SOURCE_PACKAGE_EXCLUDED_DIRS = {
    '__pycache__', 'cache', 'data', 'figures', 'logs', 'output', 'outputs', 'results',
}


def analysis_package_files(source):
    """Select small canonical package sources for both export and frozen runs.

    This is deliberately narrower than copying a package directory wholesale:
    scientific output, local caches, links and unversioned data never become
    executable frozen inputs. An absent optional package preserves older code.
    """
    source = Path(source)
    paths = []

    def linked(path):
        return path.is_symlink() or getattr(path, 'is_junction', lambda: False)()

    def visit(directory):
        for path in sorted(directory.iterdir()):
            if path.name.startswith('.') or path.name.lower() in SOURCE_PACKAGE_EXCLUDED_DIRS:
                continue
            if linked(path):
                raise ValueError('Linked analysis package source is not allowed: ' + str(path))
            if path.is_dir():
                visit(path)
                continue
            allowed = (path.suffix == '.py' or path.name.lower() in {'readme', 'readme.md'}
                       or (path.suffix == '.json' and path.name.startswith('config')))
            if not allowed or not path.is_file():
                continue
            if path.stat().st_size > 2 * 1024**2:
                raise ValueError('Unexpected large source file: ' + str(path))
            if path.suffix == '.py':
                ast.parse(path.read_bytes(), filename=str(path))
            elif path.suffix == '.json':
                json.loads(path.read_bytes())
            paths.append(path)

    for name in ANALYSIS_SOURCE_PACKAGES:
        directory = source / name
        if linked(directory):
            raise ValueError('Linked analysis package source is not allowed: ' + str(directory))
        if not directory.exists():
            continue
        if not directory.is_dir():
            raise ValueError('Analysis source package must be a directory: ' + str(directory))
        visit(directory)
    return sorted(paths)


CANONICAL_OBSERVER_ID = 'astropy_native_direct_v1'
CANONICAL_SCHEMA_ID = 'elinfo_astropy_native_direct_v1'
TAN_DISPLAY_CONTRACT_ID = 'heliocentric_tan_xi_eta_2p1_v1'
TAN_DISPLAY_RELATIVE_DIR = 'display/' + TAN_DISPLAY_CONTRACT_ID
CANONICAL_SUMMARY_COLUMNS = {
    'numsp', 'star_mass', 'distance', 'sigma', 'pmra', 'pmdec',
    'walker_R5_k_kms_per_deg', 'walker_R5_v0_kms',
    'walker_R5_sigma_raw_mle_kms', 'walker_R5_sigma_after_mle_kms',
    'walker_R5_nstar',
}


def code_observation_metadata(code):
    """Read the code identity without importing scientific dependencies."""
    module = code / 'canonical_observer.py'
    if not module.is_file() or module.is_symlink():
        raise ValueError('New analyses require the canonical observer release; select validated --code or activate it first')
    value = None
    for node in ast.parse(module.read_bytes(), filename=str(module)).body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'CANONICAL_OBSERVER_ID' for t in node.targets):
            value = ast.literal_eval(node.value)
    if value != CANONICAL_OBSERVER_ID:
        raise ValueError('Unsupported observer identity in analysis code')
    return {'canonical_observer_id': value, 'schema_id': CANONICAL_SCHEMA_ID,
            'analysis_relative_dir': 'output/analysis/' + value,
            'canonical_module_sha256': sha(module), 'cohort_definition': 'selected_then_frozen'}


def observation_output_directory(model_directory, observation):
    if observation is None:
        return model_directory
    if (observation.get('canonical_observer_id') != CANONICAL_OBSERVER_ID
            or observation.get('schema_id') != CANONICAL_SCHEMA_ID
            or observation.get('analysis_relative_dir') != 'output/analysis/' + CANONICAL_OBSERVER_ID):
        raise ValueError('Unknown or inconsistent observation manifest; prepare a new run')
    return confined(model_directory, model_directory / observation['analysis_relative_dir'])


def code_display_metadata(code):
    """Read the optional renderer contract; absent constants mean historical raw display."""
    module = code / 'canonical_display.py'
    if not module.is_file():
        return None
    if module.is_symlink():
        raise ValueError('Display module must be a regular pinned source file')
    values = {}
    def constant(node):
        if isinstance(node, ast.Name):
            return values[node.id]
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
            return constant(node.left) + constant(node.right)
        return ast.literal_eval(node)
    for node in ast.parse(module.read_bytes(), filename=str(module)).body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in ('TAN_DISPLAY_CONTRACT_ID', 'DISPLAY_RELATIVE_DIR'):
                    try:
                        values[target.id] = constant(node.value)
                    except (ValueError, KeyError, TypeError) as error:
                        raise ValueError('Display identity must use explicit string constants') from error
    if 'TAN_DISPLAY_CONTRACT_ID' not in values:
        return None
    display = {'display_contract_id': values['TAN_DISPLAY_CONTRACT_ID'],
               'relative_directory': values.get('DISPLAY_RELATIVE_DIR'),
               'display_module_sha256': sha(module)}
    validate_display_metadata(display)
    return display


def validate_display_metadata(display):
    if (display.get('display_contract_id') != TAN_DISPLAY_CONTRACT_ID
            or display.get('relative_directory') != TAN_DISPLAY_RELATIVE_DIR
            or not re.fullmatch(r'[0-9a-f]{64}', display.get('display_module_sha256', ''))):
        raise ValueError('Unknown or inconsistent display manifest; prepare a new run')


def display_output_directory(canonical_root, display):
    if display is None:
        return canonical_root
    validate_display_metadata(display)
    return confined(canonical_root, canonical_root / display['relative_directory'])


def verify_png(path):
    if path.is_symlink() or not path.is_file():
        raise ValueError(f'PNG must be a regular file: {path}')
    with path.open('rb') as image:
        if image.read(8) != b'\x89PNG\r\n\x1a\n':
            raise ValueError(f'Invalid PNG: {path}')
        image.seek(-12, 2)
        if image.read() != b'\x00\x00\x00\x00IEND\xaeB`\x82':
            raise ValueError(f'Incomplete PNG: {path}')


def verify_display_outputs(canonical_root, display, selected=None, *, summary_path, code):
    """Verify a renderer's explicit artifact list, never glob unproven videos."""
    root = display_output_directory(canonical_root, display)
    if display is None:
        raise ValueError('A display contract is required for versioned output verification')
    if sha(code / 'canonical_display.py') != display['display_module_sha256']:
        raise ValueError('Pinned display module changed')
    path = root / 'display_manifest.json'
    if path.is_symlink() or not path.is_file():
        raise ValueError('TAN display manifest is missing or not a regular file')
    manifest = json.loads(path.read_text(encoding='utf-8'))
    render_sources = {}
    for name in ('PlotFig.py', 'canonical_display.py', 'display_outputs.py'):
        if not (code / name).is_file():
            raise ValueError('Renderer source is missing: ' + name)
    for module in sorted(code.glob('*.py')):
        name = module.name
        if not re.fullmatch(r'[A-Za-z0-9_.-]+\.py', name):
            raise ValueError('Unsafe renderer source basename: ' + name)
        if module.is_symlink() or not module.is_file():
            raise ValueError('Renderer source must be a regular pinned file: ' + name)
        render_sources[name] = sha(module)
    expected_identity = {**display, **canonical_csv_identity(summary_path),
                         'render_source_sha256': render_sources}
    for key, expected in expected_identity.items():
        if manifest.get(key) != expected:
            raise ValueError('Display manifest identity mismatch: ' + key)
    frames, videos = manifest.get('frames'), manifest.get('videos')
    if not isinstance(frames, list) or not frames or not isinstance(videos, list):
        raise ValueError('Display manifest needs nonempty frames and an explicit videos list')
    numbers = []
    def artifact(record, relative):
        item = confined(root, root / relative)
        if item.is_symlink() or not item.is_file():
            raise ValueError('Display artifact is missing or not a regular file: ' + str(item))
        if item.stat().st_size != record.get('bytes') or sha(item) != record.get('sha256'):
            raise ValueError('Display artifact checksum mismatch: ' + str(item))
        return item
    for frame in frames:
        number = frame.get('numsp')
        if type(number) is not int or number < 0:
            raise ValueError('Invalid display frame snapshot number')
        numbers.append(number)
        image = artifact(frame, f'output_img/check_nature_{number:03d}.png')
        verify_png(image)
        sidecar = image.with_name(image.name + '.display.json')
        if sidecar.is_symlink() or not sidecar.is_file():
            raise ValueError('Frame display sidecar is missing or not a regular file')
        identity = json.loads(sidecar.read_text(encoding='utf-8'))
        for key, expected in {**expected_identity, 'png_sha256': frame['sha256'],
                              'numsp': number, 'bytes': frame['bytes']}.items():
            if identity.get(key) != expected:
                raise ValueError('Frame display sidecar mismatch: ' + key)
    if len(numbers) != len(set(numbers)) or (selected is not None and set(numbers) != set(selected)):
        raise ValueError('Display manifest has missing, extra or duplicate snapshots')
    names = set()
    for video in videos:
        name = video.get('name', '')
        if (not re.fullmatch(r'[A-Za-z0-9_.-]+\.mp4', name)
                or 'tan_xi_eta_2p1_v1' not in name or name in names):
            raise ValueError('Invalid or duplicate TAN video name')
        names.add(name)
        artifact(video, name)
    return manifest


def canonical_csv_identity(path):
    """Read an identity marker, never certify legacy data from a filename."""
    identities = []
    with Path(path).open(encoding='utf-8') as stream:
        for line in stream:
            if line.startswith('# canonical_observer_id='):
                fields = dict(item.split('=', 1) for item in line[2:].split() if '=' in item)
                identities.append(fields)
            elif line.strip() and not line.startswith('#'):
                break
    if len(identities) != 1:
        raise ValueError('Canonical summary requires exactly one observer identity marker')
    identity = identities[0]
    if identity.get('canonical_observer_id') != CANONICAL_OBSERVER_ID:
        raise ValueError('Legacy or different-observer summary is not canonical')
    for field in ('observer_config_sha256', 'cohort_file_sha256'):
        if not re.fullmatch(r'[0-9a-f]{64}', identity.get(field, '')):
            raise ValueError('Missing or malformed canonical identity field: ' + field)
    return identity


def verify_canonical_csv(path, observation, code, cohort_path):
    identity = canonical_csv_identity(path)
    if identity['cohort_file_sha256'] != observation['cohort_file_sha256'] or sha(cohort_path) != observation['cohort_file_sha256']:
        raise ValueError('Summary or frozen cohort identity differs from the prepared run')
    module = code / 'canonical_observer.py'
    if sha(module) != observation['canonical_module_sha256']:
        raise ValueError('Pinned observer module changed')
    environment = os.environ.copy()
    environment['PYTHONDONTWRITEBYTECODE'] = '1'
    environment['DSPH_PREEXISTING_ID_PATH'] = str(cohort_path)
    command = ('import sys; sys.path.insert(0,sys.argv[1]); '
               'from canonical_observer import require_canonical_csv; '
               'require_canonical_csv(sys.argv[2])')
    checked = subprocess.run([sys.executable, '-c', command, str(code), str(path)], env=environment,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if checked.returncode:
        raise ValueError('Pinned observer/cohort CSV guard rejected the summary: ' + checked.stderr.strip())
    return identity


def prepare(args):
    label(args.model)
    if not re.fullmatch(r'Fornax\d+(?:_.+)?', args.model):
        raise ValueError('Analysis currently supports Fornax production models')
    source = model_directory(args.root, args.model).resolve()
    if not source.is_dir():
        raise ValueError(f'Model missing or outside workspace: {source}')
    candidates = [p for p in (source, source / 'output') if list(p.glob('snapshot_*.hdf5'))]
    if args.snapshots:
        snapshot_dir = source if args.snapshots == 'model' else source / 'output'
    elif len(candidates) == 1:
        snapshot_dir = candidates[0]
    else:
        raise ValueError('Missing or ambiguous snapshots; select --snapshots model or output')
    snapshots = sorted(p for p in snapshot_dir.glob('snapshot_*.hdf5') if p.is_file())
    if not snapshots:
        raise ValueError(f'No snapshots: {snapshot_dir}')
    required = [source / n for n in ('mergegalaxy.par', 'GZWJL.PARAM', 'mergecube.par')]
    for path in required:
        if not path.is_file():
            raise ValueError(f'Missing input: {path}')
    ini = source / 'Fornax.ini'
    if not ini.is_file():
        matches = list(source.glob('IC_snapshot*_Fornax*.ini'))
        if len(matches) == 1:
            ini = matches[0]
        elif len(matches) > 1:
            raise ValueError('Multiple particle IC files; resolve the model inputs first')
        else:
            tokens = next((shlex.split(line, comments=True) for line in required[0].read_text().splitlines() if shlex.split(line, comments=True)), [])
            if len(tokens) < 2:
                raise ValueError('No dwarf IC reference in mergegalaxy.par')
            ini = (source / tokens[1]).resolve()
    if not ini.is_file():
        raise ValueError(f'Missing particle IC: {ini}')
    code = code_dir(args)
    for name in ('GetInfo.py', 'PlotFig.py', 'data_processing.py'):
        if not (code / name).is_file():
            raise ValueError(f'Missing production code: {code / name}')
    observation = code_observation_metadata(code)
    display = code_display_metadata(code)
    package_files = analysis_package_files(code)
    cohort_source = Path(os.environ.get('DSPH_PREEXISTING_ID_PATH') or
                         str(snapshot_dir / 'analysis' / CANONICAL_OBSERVER_ID / 'cohort_ids.npz')).resolve()
    if not cohort_source.is_file() or cohort_source.is_symlink() or cohort_source.suffix != '.npz':
        raise ValueError('A verified frozen uint64 particle_id cohort NPZ is required: ' + str(cohort_source))
    observation.update(cohort_source=str(cohort_source), cohort_file_sha256=sha(cohort_source),
                       cohort_relative_path='inputs/cohort_ids.npz')
    destination = run_dir(args)
    if destination.exists() or destination.is_symlink():
        raise ValueError(f'Run already exists; choose a different --name: {destination}')
    print(f'Source: {source}\nCode: {code}\nSnapshots: {snapshot_dir}\nRun: {destination}')
    if args.dry_run:
        return 0
    destination.mkdir(parents=True)
    runtime = destination / args.model
    runtime.mkdir()
    (runtime / 'output').mkdir()
    manifest = {'created_utc': now(), 'model': args.model, 'source': str(source),
                'analysis_source': str(code), 'code_sha256': {}, 'inputs': {}, 'snapshots': {},
                'adaptations': [], 'status': 'prepared'}
    manifest['observation'] = observation
    if display is not None:
        manifest['display'] = display
    # Copy the complete module set; no import from a model's obsolete local scripts.
    for path in sorted(code.glob('*.py')):
        target = destination / 'code' / path.name
        target.parent.mkdir(exist_ok=True)
        shutil.copy2(path, target)
        if path.name == 'snapshot_io.py':
            text = target.read_text()
            old = 'marker_path = f"/tmp/dsph_boxsize_{run_id}.flag"'
            if old in text:
                text = text.replace(old, 'marker_path = os.path.join(os.environ["TMPDIR"], f"dsph_boxsize_{run_id}.flag")')
                target.write_text(text)
                manifest['adaptations'].append('snapshot_io: temporary marker follows run TMPDIR')
        manifest['code_sha256'][str(target.relative_to(destination))] = sha(target)
    for path in package_files:
        target = destination / 'code' / path.relative_to(code)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        manifest['code_sha256'][str(target.relative_to(destination))] = sha(target)
    cohort_target = destination / observation['cohort_relative_path']
    cohort_target.parent.mkdir(exist_ok=True)
    shutil.copy2(cohort_source, cohort_target)
    if sha(cohort_target) != observation['cohort_file_sha256'] or sha(cohort_source) != observation['cohort_file_sha256']:
        raise ValueError('Frozen cohort changed while preparing the run')
    manifest['inputs'][observation['cohort_relative_path']] = {
        'source': str(cohort_source), 'sha256': observation['cohort_file_sha256']}
    for path in required:
        target = runtime / path.name
        shutil.copy2(path, target)
        manifest['inputs'][str(target.relative_to(destination))] = {'source': str(path), 'sha256': sha(target)}
    # A canonical Fornax.ini link avoids breaking mergegalaxy's old relative path.
    linked = [(ini, runtime / 'Fornax.ini')]
    linked += [(p, runtime / 'output' / p.name) for p in snapshots]
    for path, target in linked:
        target.symlink_to(os.path.relpath(path, target.parent))
        st = path.stat()
        record = {'source': str(path), 'size': st.st_size, 'mtime_ns': st.st_mtime_ns}
        field = 'snapshots' if target.parent.name == 'output' else 'inputs'
        manifest[field][str(target.relative_to(destination))] = record
    for name in ('SFH_Fornax.csv', 'Sample1_candidates.csv', 'Sample1_density.csv', 'Fornax_surface_brightness_profile.csv'):
        path = source / name if (source / name).is_file() else args.root / name
        if path.is_file():
            target = destination / name if name in ('SFH_Fornax.csv', 'Fornax_surface_brightness_profile.csv') else runtime / name
            shutil.copy2(path, target)
            manifest['inputs'][str(target.relative_to(destination))] = {'source': str(path), 'sha256': sha(target)}
    # Cube paths remain readable; analysis uses these as inputs only.
    for path in args.root.glob('Vube*.ini'):
        (destination / path.name).symlink_to(os.path.relpath(path, destination))
    write_json(destination / 'manifest.json', manifest)
    print('Prepared. Run statistics before plots; previous elinfo is preserved at the source.')
    return 0


def selected_snapshots(manifest, selection):
    match = re.fullmatch(r'(\d+),(\d+)', selection)
    if not match or int(match[1]) > int(match[2]):
        raise ValueError('--range must be START,END (inclusive)')
    available = {int(re.search(r'snapshot_(\d+)', p)[1]) for p in manifest['snapshots']}
    selected = sorted(n for n in available if int(match[1]) <= n <= int(match[2]))
    if not selected:
        raise ValueError('No available snapshots in requested range')
    return selected


def verify_inputs(directory, manifest):
    for name, expected in manifest['code_sha256'].items():
        if sha(directory / name) != expected:
            raise ValueError(f'Pinned code changed: {name}; prepare a new run')
    for name, entry in {**manifest['inputs'], **manifest['snapshots']}.items():
        path = directory / name
        if 'sha256' in entry:
            if sha(path) != entry['sha256']:
                raise ValueError(f'Input changed: {name}')
        else:
            st = path.stat()
            if st.st_size != entry['size'] or st.st_mtime_ns != entry['mtime_ns']:
                raise ValueError(f'Source changed since preparation: {name}; prepare a new run')


def verify_summary(path, selected, observation=None, code=None, cohort_path=None):
    if observation is not None:
        if code is None or cohort_path is None:
            raise ValueError('Canonical summary verification requires pinned code and cohort')
        verify_canonical_csv(path, observation, code, cohort_path)
    with path.open(encoding='utf-8') as stream:
        reader = csv.DictReader(line for line in stream if not line.startswith('#'))
        required = {'numsp', 'star_mass', 'distance', 'sigma'}
        if observation is not None:
            required = CANONICAL_SUMMARY_COLUMNS
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'Incomplete summary schema: {path}')
        numbers = [int(float(row['numsp'])) for row in reader]
    if observation is not None and not numbers:
        raise ValueError('Canonical summary contains no completed snapshot rows')
    if len(numbers) != len(set(numbers)) or not set(selected).issubset(numbers):
        raise ValueError('Summary has missing or duplicate snapshots')


def run(args):
    directory = run_dir(args)
    lock = (directory / '.run.lock').open('a')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close()
        raise ValueError('This analysis run is already active')
    manifest = json.loads((directory / 'manifest.json').read_text())
    selected = selected_snapshots(manifest, args.range)
    verify_inputs(directory, manifest)
    if args.processes < 1:
        raise ValueError('--processes must be positive')
    runtime = directory / manifest['model']
    observation = manifest.get('observation')
    display = manifest.get('display')
    output_root = observation_output_directory(runtime, observation)
    display_root = display_output_directory(output_root, display)
    summary_path = output_root / f"elinfo_{manifest['model']}.csv"
    cohort_path = confined(directory, directory / observation['cohort_relative_path']) if observation is not None else None
    summary_options = {'observation': observation, 'code': directory / 'code', 'cohort_path': cohort_path}
    environment = os.environ.copy()
    environment.update({'PYTHONDONTWRITEBYTECODE': '1', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
                        'MKL_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1', 'MPLBACKEND': 'Agg'})
    if observation is not None:
        environment['DSPH_CANONICAL_OUTPUT_ROOT'] = str(output_root)
        environment['DSPH_PREEXISTING_ID_PATH'] = str(cohort_path)
    for key, subdir in [('TMPDIR', '.tmp'), ('MPLCONFIGDIR', '.mpl'), ('XDG_CACHE_HOME', '.cache')]:
        path = directory / subdir
        path.mkdir(exist_ok=True)
        environment[key] = str(path)
    logdir = directory / 'logs'
    logdir.mkdir(exist_ok=True)
    status = {'started_utc': now(), 'status': 'running', 'stages': [], 'range': args.range, 'python': sys.executable}
    status['canonical_observer_id'] = observation['canonical_observer_id'] if observation is not None else 'legacy_unversioned'
    status['analysis_output_directory'] = str(output_root)
    if display is not None:
        status['display_contract_id'] = display['display_contract_id']
        status['display_output_directory'] = str(display_root)
    write_json(directory / 'status.json', status)
    stages = ['info', 'plot'] if args.stage == 'all' else [args.stage]
    try:
        for stage in stages:
            if stage == 'plot':
                verify_summary(summary_path, selected, **summary_options)
            script = 'GetInfo.py' if stage == 'info' else 'PlotFig.py'
            command = [sys.executable, str(directory / 'code' / script), '--range', args.range,
                       '--processes', str(args.processes), '--worker-threads', '1']
            record = {'stage': stage, 'command': command, 'started_utc': now()}
            status['stages'].append(record)
            write_json(directory / 'status.json', status)
            print(shlex.join(command), flush=True)
            with (logdir / f'{stage}.log').open('a') as log:
                log.write(f'\n=== {now()} ===\n'); log.flush()
                subprocess.run(command, cwd=runtime, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
            if stage == 'info':
                verify_summary(summary_path, selected, **summary_options)
            else:
                if display is not None:
                    verify_display_outputs(output_root, display, selected,
                                           summary_path=summary_path, code=directory / 'code')
                else:
                    for number in selected:
                        verify_png(display_root / 'output_img' / f'check_nature_{number:03d}.png')
            record.update(status='complete', finished_utc=now())
        status['status'] = 'complete'
    except BaseException as error:
        status.update(status='failed', error=str(error))
        raise
    finally:
        status['finished_utc'] = now()
        try:
            write_json(directory / 'status.json', status)
        finally:
            lock.close()
    return 0


def submit(args):
    directory = run_dir(args)
    manifest = json.loads((directory / 'manifest.json').read_text())
    selected_snapshots(manifest, args.range)
    verify_inputs(directory, manifest)
    settings = config(args.root)['slurm']
    if args.processes < 1:
        raise ValueError('--processes must be positive')
    receipt = directory / 'submission.json'
    if receipt.exists():
        if not getattr(args, 'resume', False):
            raise ValueError('This run was already submitted; use --resume after it finishes')
        previous = json.loads(receipt.read_text())
        job, _, cluster = previous['job'].partition(';')
        queued = subprocess.run(['squeue', '-M', cluster or settings['cluster'], '--noheader', '--jobs', job],
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if queued.returncode or queued.stdout.strip():
            raise ValueError('Previous job is still active or its state cannot be confirmed')
    script = directory / 'analysis.sbatch'
    command = [settings.get('python', 'python3'), str(args.root / 'workflow.py'), '--root', str(args.root),
               'run', '--name', args.name, '--range', args.range, '--processes', str(args.processes), '--stage', args.stage]
    text = '#!/bin/bash\nset -euo pipefail\n'
    text += f'export TMPDIR={shlex.quote(str(directory / ".tmp"))}\nmkdir -p "$TMPDIR"\n'
    text += 'export PYTHONDONTWRITEBYTECODE=1\n'
    text += '. ' + shlex.quote(settings['module_init']) + '\nmodule purge\nmodule load ' + shlex.join(settings['modules']) + '\n'
    text += 'exec ' + shlex.join(command) + '\n'
    script.write_text(text)
    command = ['sbatch', '--parsable', '--clusters', settings['cluster'], '--partition', settings['partition'],
               '--qos', settings['qos'], '--account', settings['account'], '--nodes', '1', '--ntasks', '1',
               '--cpus-per-task', str(args.processes), '--time', args.time or settings['time'],
               '--chdir', str(directory), '--output', str(directory / 'slurm-%j.log'), '--job-name', 'dsph-' + args.name]
    if args.test_only:
        command.append('--test-only')
    command.append(str(script))
    print(shlex.join(command), flush=True)
    submitted_utc = now()
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(result.stdout, end=''); print(result.stderr, end='', file=sys.stderr)
    if result.returncode:
        return result.returncode
    if not args.test_only:
        if receipt.exists():
            with (directory / 'submission_history.jsonl').open('a') as history:
                history.write(json.dumps(json.loads(receipt.read_text())) + '\n')
        write_json(receipt, {'submitted_utc': submitted_utc, 'job': result.stdout.strip(), 'command': command})
    return 0
