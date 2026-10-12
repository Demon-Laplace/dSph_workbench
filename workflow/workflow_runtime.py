"""Configured analysis runs, Slurm submission, and result verification."""
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


def code_dir(args):
    setting = getattr(args, 'code', None) or config(args.root).get('analysis_code', 'tools/server-analysis')
    return (args.root / setting).resolve()


def run_dir(args):
    return confined(args.root, args.root / 'workspace' / 'analysis' / label(args.name))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(args):
    label(args.model)
    if not re.fullmatch(r'Fornax\d+(?:_.+)?', args.model):
        raise ValueError('Analysis currently supports Fornax production models')
    source = (args.root / args.model).resolve()
    if not source.is_dir() or not source.is_relative_to(args.root):
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


def verify_summary(path, selected):
    with path.open() as stream:
        reader = csv.DictReader(line for line in stream if not line.startswith('#'))
        required = {'numsp', 'star_mass', 'distance', 'sigma'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f'Incomplete summary schema: {path}')
        numbers = [int(float(row['numsp'])) for row in reader]
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
    environment = os.environ.copy()
    environment.update({'PYTHONDONTWRITEBYTECODE': '1', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
                        'MKL_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1', 'MPLBACKEND': 'Agg'})
    for key, subdir in [('TMPDIR', '.tmp'), ('MPLCONFIGDIR', '.mpl'), ('XDG_CACHE_HOME', '.cache')]:
        path = directory / subdir
        path.mkdir(exist_ok=True)
        environment[key] = str(path)
    logdir = directory / 'logs'
    logdir.mkdir(exist_ok=True)
    status = {'started_utc': now(), 'status': 'running', 'stages': [], 'range': args.range, 'python': sys.executable}
    write_json(directory / 'status.json', status)
    stages = ['info', 'plot'] if args.stage == 'all' else [args.stage]
    try:
        for stage in stages:
            if stage == 'plot':
                verify_summary(runtime / f"elinfo_{manifest['model']}.csv", selected)
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
                verify_summary(runtime / f"elinfo_{manifest['model']}.csv", selected)
            else:
                for number in selected:
                    path = runtime / 'output_img' / f'check_nature_{number:03d}.png'
                    with path.open('rb') as image:
                        if image.read(8) != b'\x89PNG\r\n\x1a\n':
                            raise ValueError(f'Invalid PNG: {path}')
                        image.seek(-12, 2)
                        if image.read() != b'\x00\x00\x00\x00IEND\xaeB`\x82':
                            raise ValueError(f'Incomplete PNG: {path}')
            record.update(status='complete', finished_utc=now())
        status['status'] = 'complete'
    except BaseException as error:
        status.update(status='failed', error=str(error))
        raise
    finally:
        status['finished_utc'] = now()
        write_json(directory / 'status.json', status)
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
