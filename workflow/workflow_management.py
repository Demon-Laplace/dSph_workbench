"""Workspace catalogues and transport. All writes stay below configured roots."""
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import shutil
import time

import workflow_runtime as runtime


def link(root, target, source, replace=False):
    runtime.confined(root, target.parent)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.is_symlink():
        if target.resolve() != source.resolve():
            if not replace:
                raise ValueError(f'Existing link points elsewhere: {target}')
            temporary = target.with_name(target.name + '.next')
            temporary.symlink_to(os.path.relpath(source, target.parent), target_is_directory=source.is_dir())
            temporary.replace(target)
    elif target.exists():
        if target.resolve() != source.resolve():
            raise ValueError(f'Refusing to replace existing path: {target}')
    else:
        target.symlink_to(os.path.relpath(source, target.parent), target_is_directory=source.is_dir())


def organize(args):
    root = args.root
    base = root / 'workspace'
    for name in ('models/production', 'models/isolation', 'models/milky-way', 'code', 'inputs/observations',
                 'inputs/shared', 'inputs/parameters', 'templates', 'notebooks', 'studies', 'results/legacy',
                 'results/remote', 'analysis', 'catalog', 'archive'):
        runtime.confined(root, base / name).mkdir(parents=True, exist_ok=True)
    records = []
    for path in sorted(root.iterdir()):
        name = path.name
        if name.startswith('.') or name.startswith('workflow') or name in {'workspace', 'WORKFLOW.md', '__pycache__', 'tests'}:
            continue
        category = None
        if path.is_dir():
            if re.fullmatch(r'Fornax\d+(?:_.*)?', name) and not any(x in name for x in ('paper_figures', 'orientation_check')):
                category = 'models/production'
            elif re.match(r'MW\d+R', name) or re.fullmatch(r'MW\d+', name):
                category = 'models/milky-way'
            elif name == 'isoFornax':
                for model in sorted(path.iterdir()):
                    if model.is_dir() and re.fullmatch(r'Fornax\d+(?:_.*)?', model.name):
                        link(root, base / 'models/isolation' / model.name, model)
                category = 'models/isolation'
            elif name == 'observer_review_delete_20261006':
                category = 'archive'
            elif name in ('Fornax', 'MWbase'):
                category = 'templates'
            elif name.startswith('dSph_workbench') or name in ('tools', 'test_pyscript'):
                category = 'code'
            elif name == 'modelfile':
                category = 'inputs/parameters'
            else:
                category = 'studies'
        elif path.suffix == '.ipynb':
            category = 'notebooks'
        elif path.suffix in ('.png', '.mp4', '.pdf'):
            category = 'results/legacy'
        elif path.suffix in ('.py', '.sh'):
            category = 'code/legacy'
        elif name.startswith(('Sample1_', 'SFH_', 'Fornax_surface_brightness')):
            category = 'inputs/observations'
        elif path.suffix in ('.ini', '.npy', '.par'):
            category = 'inputs/shared'
        if category:
            target = base / category / name
            link(root, target, path)
            records.append({'name': name, 'category': category, 'source': str(path), 'entry': str(target.relative_to(root))})
    code = runtime.code_dir(args)
    if code.is_dir():
        link(root, base / 'code/production', code, replace=True)
    if args.repo.is_dir():
        link(root, base / 'code/github', args.repo, replace=True)
        link(root, base / 'inputs/parameters/github', args.repo / 'IniFile', replace=True)
    runtime.write_json(base / 'catalog/layout.json', {'created_utc': runtime.now(), 'entries': records})
    print(f'Organized {len(records)} entries under {base}; existing model data stay in place.')
    return 0


def ssh_config(args):
    remote = runtime.config(args.root)['remote']
    command = ['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=15',
               '-o', 'ServerAliveInterval=10', '-o', 'ServerAliveCountMax=2']
    if remote.get('proxy_jump'):
        command += ['-J', remote['proxy_jump']]
    socket = remote.get('control_socket')
    if socket and Path(socket).is_socket():
        command += ['-S', socket]
    return remote, command


def remote(args):
    return subprocess.run(remote_command(args, args.arguments)).returncode


def remote_command(args, arguments):
    settings, command = ssh_config(args)
    if arguments and arguments[0] == '--':
        arguments = arguments[1:]
    if not arguments or arguments[0] in ('remote', 'pull'):
        raise ValueError('Specify a remote workflow command, e.g. remote status')
    root = settings['root']
    cmd = ['python3', root + '/workflow.py', '--root', root, *arguments]
    setup = '. ' + shlex.quote(settings['module_init']) + ' && module load ' + shlex.join(settings['modules']) + ' && '
    return command + [settings['host'], setup + 'PYTHONDONTWRITEBYTECODE=1 ' + shlex.join(cmd)]


def production_results(root, model=None, name=None):
    """Small verified manifest of top-level production videos/tables, never frames."""
    import workflow_isolation as isolation
    expected_frames = None
    if name:
        runtime.label(name)
        run = runtime.confined(root, root / 'workspace/analysis' / name)
        manifest = json.loads((run / 'manifest.json').read_text())
        model = manifest['model']
        source = runtime.confined(run, run / model)
        state = json.loads((run / 'status.json').read_text()) if (run / 'status.json').exists() else {}
        if state.get('status') == 'complete' and state.get('range'):
            expected_frames = len(runtime.selected_snapshots(manifest, state['range']))
    else:
        source = runtime.model_directory(root, model)
    runtime.label(model)
    files = sorted(source.glob('*.mp4')) + list(source.glob('elinfo_' + model + '.csv'))
    if not files: raise ValueError('No production videos or summary table available')
    records = []
    for path in files:
        runtime.confined(source, path)
        if path.is_symlink() or not path.is_file() or not re.fullmatch(r'[A-Za-z0-9_.-]+', path.name):
            raise ValueError('Invalid production output: ' + str(path))
        records.append({'name': path.name, 'bytes': path.stat().st_size, 'sha256': isolation.sha(path)})
    return {'model': model, 'directory': str(source), 'files': records, 'expected_frames': expected_frames}


def fetch_production(args, model=None, name=None):
    import workflow_isolation as isolation
    settings, ssh = ssh_config(args)
    command = ('from pathlib import Path; import json,workflow_management as m; '
               'print(json.dumps(m.production_results(Path(' + repr(settings['root']) + '), '
               'model=' + repr(model) + ',name=' + repr(name) + ')))')
    response = isolation.remote_exec(args, ['python3', '-c', command], capture=True)
    result = json.loads(response.stdout)
    destination = runtime.model_directory(args.root, result['model'])
    filters = []
    for record in result['files']:
        filename = record['name']
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', filename) or not filename.endswith(('.mp4', '.csv')):
            raise ValueError('Invalid production result filename')
        target = runtime.confined(destination, destination / filename)
        if target.is_symlink() or (target.exists() and not target.is_file()):
            raise ValueError('Production result target is not a regular file')
        filters.append('--include=/' + filename)
    command = ['rsync', '-rt', '--checksum', '--itemize-changes', '-e', shlex.join(ssh), *filters,
               '--exclude=*', settings['host'] + ':' + result['directory'] + '/', str(destination) + '/']
    if args.dry_run: command.insert(1, '--dry-run')
    else: destination.mkdir(parents=True, exist_ok=True)
    print(shlex.join(command), flush=True)
    subprocess.run(command, check=True)
    if args.dry_run: return 0
    for record in result['files']:
        path = destination / record['name']
        if path.stat().st_size != record['bytes'] or isolation.sha(path) != record['sha256']:
            raise ValueError('Production output checksum mismatch; retry fetching')
        if path.suffix == '.mp4':
            info = json.loads(subprocess.check_output(['ffprobe', '-v', 'error', '-select_streams', 'v:0',
                '-show_entries', 'stream=nb_frames,duration', '-of', 'json', str(path)], text=True))
            streams = info.get('streams', [])
            frames = int(streams[0].get('nb_frames', 0)) if len(streams) == 1 else 0
            if frames < 1 or (result['expected_frames'] is not None and frames != result['expected_frames']):
                raise ValueError('Production video frame count mismatch')
            record['frames'] = frames
    receipt = runtime.confined(args.root, args.root / 'workspace/receipts/production' / result['model'] /
                              ((name or 'model') + '.json'))
    receipt.parent.mkdir(parents=True, exist_ok=True)
    runtime.write_json(receipt, dict(result, local_directory=str(destination), verified_utc=runtime.now()))
    print('Saved production video/summary: ' + str(destination))
    return 0


def pull(args):
    runtime.label(args.name)
    return fetch_production(args, name=args.name)


def pipeline(args):
    result = runtime.prepare(args)
    if result or args.dry_run:
        return result
    args.test_only = not args.submit
    return runtime.submit(args)


def status(args):
    base = args.root / 'workspace/analysis'
    paths = [runtime.run_dir(args)] if args.name else sorted(base.glob('*'))
    rows = []
    for path in paths:
        if not (path / 'manifest.json').is_file():
            continue
        manifest = json.loads((path / 'manifest.json').read_text())
        state = json.loads((path / 'status.json').read_text()) if (path / 'status.json').exists() else {'status': 'prepared'}
        submission = json.loads((path / 'submission.json').read_text()) if (path / 'submission.json').exists() else {}
        if submission and state.get('started_utc', '') < submission['submitted_utc']:
            state = {'status': 'queued'}
        row = {'name': path.name, 'model': manifest['model'], 'status': state['status'], 'job': submission.get('job', '-')}
        if submission and state['status'] != 'complete':
            job, _, cluster = submission['job'].partition(';')
            accounting = subprocess.run(['sacct', '-M', cluster, '-j', job, '--noheader', '--parsable2', '--format=JobIDRaw,State'],
                                        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if accounting.returncode == 0:
                states = [line.split('|')[1] for line in accounting.stdout.splitlines() if line.split('|')[0] == job]
                if states:
                    row['scheduler_state'] = states[0]
                    if states[0].split()[0].rstrip('+') in {'FAILED', 'CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY', 'NODE_FAIL', 'PREEMPTED'}:
                        row['status'] = 'failed'
        rows.append(row)
    if getattr(args, 'json', False):
        print(json.dumps(rows))
        return 0
    for row in rows:
        print(row['name'], row['model'], row['status'], row['job'], row.get('scheduler_state', ''))
    if 'slurm' in runtime.config(args.root):
        return subprocess.run(['squeue', '-M', 'all', '--me', '-o', '%.18i %.12P %.32j %.2t %.10M %R']).returncode
    return 0


def collect(args):
    deadline = time.monotonic() + args.timeout
    while True:
        result = subprocess.run(remote_command(args, ['status', '--name', args.name, '--json']),
                                text=True, stdout=subprocess.PIPE, check=True)
        rows = json.loads(result.stdout)
        if not rows:
            raise ValueError('Remote run does not exist')
        state = rows[0]['status']
        print(f'{args.name}: {state}', flush=True)
        if state == 'complete':
            args.dry_run = False
            return pull(args)
        if state == 'failed':
            raise ValueError('Remote analysis failed; inspect its status and stage logs before resuming')
        if not args.wait or time.monotonic() >= deadline:
            return 1
        time.sleep(min(15, max(0, deadline - time.monotonic())))


def parameter_path(args, number):
    if not number.isdigit():
        raise ValueError('Use numeric Fornax model IDs')
    name = f'IC_Fornax{number}.ini'
    if runtime.config(args.root).get('ic', {}).get('root'):
        from workflow_models import ic_paths
        root, _ = ic_paths(args)
        draft = root / 'drafts' / name
        return draft if draft.is_file() else root / 'IniFile' / name
    draft = args.root / 'workspace/inputs/parameters/drafts' / name
    return draft if draft.is_file() else args.repo / 'IniFile' / name


def new_ini(args):
    if not args.model.isdigit():
        raise ValueError('Use a numeric model ID')
    source = parameter_path(args, args.base)
    if runtime.config(args.root).get('ic', {}).get('root'):
        from workflow_models import ic_paths
        root, _ = ic_paths(args)
        target = runtime.confined(root, root / 'drafts' / f'IC_Fornax{args.model}.ini')
        if (root / 'IniFile' / target.name).exists():
            raise ValueError('Model ID already exists; choose a new ID')
    else:
        target = runtime.confined(args.root, args.root / 'workspace/inputs/parameters/drafts' / f'IC_Fornax{args.model}.ini')
    if target.exists() or (args.repo / 'IniFile' / target.name).exists():
        raise ValueError('Model ID already exists; choose a new ID')
    target.parent.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, str(args.repo / 'model_tool.py'), 'create', str(target), '--base', str(source)]
    for setting in args.settings:
        command += ['--set', setting]
    return subprocess.run(command).returncode


def verify_params(args):
    batch = runtime.confined(args.root, args.root / 'workspace/inputs/parameter_batches' / runtime.label(args.batch))
    manifest = json.loads((batch / 'manifest.json').read_text())
    for name, digest in manifest['files'].items():
        if not re.fullmatch(r'IC_Fornax\d+\.ini', name) or runtime.sha(batch / name) != digest:
            raise ValueError(f'Transfer checksum mismatch: {name}')
    runtime.write_json(batch / 'verified.json', {'verified_utc': runtime.now(), 'files': manifest['files']})
    print(f'Verified {len(manifest["files"])} parameter files in {batch}')
    return 0


def push_params(args):
    runtime.label(args.batch)
    paths = [parameter_path(args, n) for n in args.models]
    subprocess.run([sys.executable, str(args.repo / 'model_tool.py'), 'check', *map(str, paths)], check=True)
    staging = runtime.confined(args.root, args.root / 'workspace/transfers/parameters' / args.batch)
    if staging.exists():
        raise ValueError('Transfer batch already exists; use a new --batch name')
    staging.mkdir(parents=True)
    for path in paths:
        shutil.copy2(path, staging / path.name)
    runtime.write_json(staging / 'manifest.json', {'created_utc': runtime.now(), 'files': {p.name: runtime.sha(p) for p in paths}})
    settings, ssh = ssh_config(args)
    target = settings['root'] + '/workspace/inputs/parameter_batches/' + args.batch
    quoted = shlex.quote(target)
    subprocess.run(ssh + [settings['host'], f'test ! -e {quoted} && mkdir -p {quoted}'], check=True)
    subprocess.run(['rsync', '-rt', '--checksum', '-e', shlex.join(ssh), str(staging) + '/', settings['host'] + ':' + target + '/'], check=True)
    args.arguments = ['verify-params', '--batch', args.batch]
    return remote(args)


def fetch_model(args):
    runtime.label(args.model)
    return fetch_production(args, model=args.model)
