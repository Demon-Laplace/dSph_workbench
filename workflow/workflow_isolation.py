"""Resumable local IC transfer and isolated Slurm simulation without routine rereads."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile
from types import SimpleNamespace
import uuid

import workflow_runtime as runtime
import workflow_versions as versions


def sha(path):
    """Explicit diagnostic helper retained for legacy callers, not normal paths here."""
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def profile(root, name):
    profiles = runtime.config(root).get('isolation_profiles', {})
    if name not in profiles:
        raise ValueError(f'Configure isolation_profiles.{name} in the server workflow.json first')
    result = dict(profiles[name])
    template = runtime.confined(root, root / result['template']).resolve()
    if not (template / 'GZWJL.PARAM').is_file():
        raise ValueError('Isolation template is missing GZWJL.PARAM')
    for key in ('converter', 'solver'):
        if not isinstance(result.get(key), list) or not result[key] or not all(isinstance(s, str) for s in result[key]):
            raise ValueError(f'Profile {key} must be an explicit argument list')
    result['template'] = str(template)
    return result


def rewrite_parameters(source, target, values):
    return rewrite_parameter_text(source.read_text(), target, values)


def rewrite_parameter_text(text, target, values):
    remaining = set(values)
    seen = set()
    lines = []
    for line in text.splitlines():
        fields = line.split('%', 1)[0].split()
        if fields and fields[0] in values:
            key = fields[0]
            if key in seen:
                raise ValueError(f'Duplicate simulation parameter: {key}')
            seen.add(key)
            remaining.discard(key)
            line = key + ' ' + str(values[key])
        lines.append(line)
    if remaining:
        raise ValueError('Template lacks parameters: ' + ', '.join(sorted(remaining)))
    target.write_text('\n'.join(lines) + '\n')


def request_inputs(request):
    if set(request['files']) != {'Fornax.ini', 'parameters.ini'}:
        raise ValueError('Unexpected IC transfer file set')
    return request['files']


def verify_bundle(stage, request):
    """Legacy compatibility: retain the filename allowlist, not content checks."""
    request_inputs(request)
    return {'data_verification_performed': False}


def request_matches(saved, requested):
    """Keep run-name bindings without reading or comparing particle contents."""
    if any(saved.get(k) != requested.get(k) for k in ('model', 'name', 'profile')):
        return False
    if set(saved.get('files', {})) != set(requested.get('files', {})):
        return False
    if not versions.same_source(saved, requested):
        return False
    return all(not (saved.get(key) and requested.get(key)) or saved[key] == requested[key]
               for key in ('ic_version', 'parameter_version'))


def version_directory(root, model):
    base = runtime.confined(root, root / runtime.config(root).get('isolation_root', 'isoFornax'))
    runtime.model_directory(root, model)  # Validate the model label, not its contents.
    return runtime.confined(root, base / model)


def resolve_request(root, request):
    request_inputs(request)
    runtime.label(request['name'])
    runtime.label(request.get('requested_name', request['name']))
    settings = profile(root, request['profile'])
    entry = versions.resolve(
        version_directory(root, request['model']), request, settings,
        (Path(settings['template']) / 'GZWJL.PARAM').read_text(),
        lambda name: case_dir(root, name, request['model']))
    entry['existing'] = Path(entry['directory']).is_dir()
    return entry


def result_record(root, run, request, status, reused):
    directory = version_directory(root, request['model'])
    versions.update(directory, request['model'], run.name, status)
    entry = versions.lookup(directory, request['model'], run.name)
    actual = entry['request'] if entry else request
    return dict(status, directory=str(run), reused=reused,
                **{key: actual.get(key) for key in
                   ('model', 'name', 'requested_name', 'build_id', 'source_version',
                    'source_directory', 'profile', 'ic_version', 'parameter_version',
                    'version_contract')})


def case_dir(root, name, model=None):
    runtime.label(name)
    base = runtime.confined(root, root / runtime.config(root).get('isolation_root', 'isoFornax'))
    if model:
        if not re.fullmatch(r'(?:Fornax|MW|Sculptor)\d+[A-Za-z0-9_.-]*', model):
            raise ValueError('Use a full model ID')
        current = runtime.confined(root, base / model / 'runs' / name)
        previous = runtime.confined(root, root / model / 'isolation' / name)
        if previous.is_dir():
            if current.exists() and current.resolve() != previous.resolve():
                raise ValueError('Ambiguous existing isolation runs; reconcile before retrying')
            # Frozen completed/active cases retain their original absolute paths.
            return previous
        return current
    matches = list(base.glob('*/runs/' + name))
    matches += list(root.glob('*/isolation/' + name))
    legacy = root / 'workspace/models/isolation' / name
    if legacy.is_dir():
        matches.append(legacy)
    unique = {p.resolve() for p in matches if p.is_dir()}
    if len(unique) != 1:
        raise ValueError('Specify --model for an absent or ambiguous isolation run')
    return runtime.confined(root, unique.pop())


def verify_case(run):
    """Legacy compatibility: no routine file-content verification is performed."""
    return {'data_verification_performed': False}


def verify_completion(run):
    """Legacy callers may record process completion without reopening snapshots."""
    return {'scope': 'solver process exit only', 'data_verification_performed': False}


def discard_received_transfer(root, stage, run, request):
    """Release only the accepted upload staging copy, never the simulation case."""
    transfers = runtime.confined(root, root / 'workspace/transfers/isolation')
    if stage != run and stage.parent == transfers:
        shutil.rmtree(stage)


def receive(root, stage, submit=False):
    stage = runtime.confined(root, stage).resolve()
    request = json.loads((stage / 'request.json').read_text())
    request_inputs(request)
    settings = template_text = None
    run = case_dir(root, request['name'], request['model'])
    if request.get('version_contract') == versions.CONTRACT:
        entry, settings, template_text = versions.binding(version_directory(root, request['model']), request)
        request = entry['request']
        run = runtime.confined(root, Path(entry['directory']))
    elif stage != run or not (run / 'status.json').is_file():
        entry = resolve_request(root, request)
        request = entry['request']
        run = runtime.confined(root, Path(entry['directory']))
        if request.get('version_contract') == versions.CONTRACT:
            _, settings, template_text = versions.binding(version_directory(root, request['model']), request)
    removed = root / 'workspace/receipts/isolation-cleaned' / request['model'] / (request['name'] + '.json')
    if removed.exists():
        raise ValueError('This failed attempt was cleaned; use a new run name')
    run.parent.mkdir(parents=True, exist_ok=True)
    lock = runtime.confined(root, run.parent / (run.name + '.lock'))
    with lock.open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        reused = run.exists()
        if run.exists():
            saved = json.loads((run / 'request.json').read_text())
            if not request_matches(saved, request):
                # Untagged historical requests can still resume themselves.
                if not (stage == run and saved == request):
                    raise ValueError('Run name is already bound to different ICs or settings')
        else:
            if settings is None:
                raise ValueError('A new run requires a reserved version binding')
            # Assemble privately, then atomically publish a complete case.
            with tempfile.TemporaryDirectory(prefix='isolation-', dir=run.parent) as temporary:
                prepared = Path(temporary) / 'case'
                prepared.mkdir()
                for name in request['files']:
                    shutil.copy2(stage / name, prepared / name)
                for name in ('workflow_isolation.py', 'workflow_runtime.py', 'workflow_versions.py'):
                    shutil.copy2(Path(__file__).with_name(name), prepared / name)
                runtime.write_json(prepared / 'request.json', request)
                runtime.write_json(prepared / 'profile.json', settings)
                template = Path(settings['template'])
                values = dict(settings.get('parameter_overrides', {}))
                values.update(InitCondFile=str(run / 'Fornax.IC_000'), OutputDir=str(run / 'output'), ResubmitOn='0')
                rewrite_parameter_text(template_text, prepared / 'GZWJL.PARAM', values)
                cooling = template / 'TREECOOL'
                if settings.get('cooling_required') and not cooling.is_file():
                    raise ValueError('Missing required TREECOOL in isolation template')
                if cooling.is_file():
                    shutil.copy2(cooling, prepared / 'TREECOOL')
                (prepared / 'output').mkdir()
                command = [settings.get('python', 'python3'), str(run / 'workflow_isolation.py'),
                           'run-job', '--root', str(root), '--name', run.name, '--model', request['model']]
                script = '#!/bin/bash\nset -euo pipefail\n'
                script += 'cd ' + shlex.quote(str(run)) + '\n'
                script += 'export TMPDIR="$PWD/.tmp"\nmkdir -p "$TMPDIR"\n'
                script += settings.get('setup', '') + '\n'
                script += 'exec ' + shlex.join(command) + '\n'
                (prepared / 'run.sh').write_text(script)
                runtime.write_json(prepared / 'status.json', {'state': 'prepared', 'created_utc': runtime.now()})
                prepared.rename(run)
        status = json.loads((run / 'status.json').read_text())
        if status['state'] in ('submitted', 'running', 'complete', 'failed'):
            discard_received_transfer(root, stage, run, request)
            return result_record(root, run, request, status, reused)
        if status['state'] == 'submitting':
            raise ValueError('Previous submission outcome is uncertain; reconcile Slurm before retrying, to prevent duplicate jobs')
        # Do not overwrite an IC or submit into an existing simulation output.
        if any((run / 'output').iterdir()) or (run / 'Fornax.IC_000').exists():
            raise ValueError('Run already contains simulation data; new submission refused')
        settings = json.loads((run / 'profile.json').read_text())
        if not submit:
            discard_received_transfer(root, stage, run, request)
            return result_record(root, run, request, dict(status, submitted=False), reused)
        options = settings.get('sbatch', {})
        allowed = {'clusters', 'qos', 'partition', 'account', 'nodes', 'ntasks', 'cpus-per-task', 'time', 'mem', 'constraint'}
        if set(options) - allowed:
            raise ValueError('Unsupported Slurm profile option')
        command = ['sbatch', '--parsable', '--chdir', str(run), '--job-name', 'iso-' + run.name,
                   '--output', str(run / 'slurm-%j.out'), '--error', str(run / 'slurm-%j.out')]
        for key, value in options.items():
            command += ['--' + key, str(value)]
        command.append(str(run / 'run.sh'))
        if submit:
            status.update(state='submitting', submission_started_utc=runtime.now())
            runtime.write_json(run / 'status.json', status)
            versions.update(version_directory(root, request['model']), request['model'], run.name, status)
        result = subprocess.run(command, capture_output=True, text=True)
        runtime.write_json(run / 'submission.json', {'command': command, 'returncode': result.returncode,
                           'stdout': result.stdout, 'stderr': result.stderr})
        if result.returncode:
            # A nonzero sbatch exit is a definitive rejection, not a lost acknowledgement.
            status.update(state='prepared', last_error=result.stderr)
            runtime.write_json(run / 'status.json', status)
            versions.update(version_directory(root, request['model']), request['model'], run.name, status)
            raise ValueError('Slurm rejected the request: ' + result.stderr.strip())
        if submit:
            job = result.stdout.strip().split(';')[0]
            if not job.isdigit():
                raise ValueError('Unrecognized Slurm job ID; reconcile submission before retrying')
            status.update(state='submitted', job_id=job)
        runtime.write_json(run / 'status.json', status)
        discard_received_transfer(root, stage, run, request)
        return result_record(root, run, request, status, reused)


def resume(root, name, model, time_max, submit=False):
    """Extend TimeMax in place; never copy checkpoints or reconvert the IC."""
    run = case_dir(root, name, model)
    with (run.parent / (run.name + '.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        status = json.loads((run / 'status.json').read_text())
        settings = json.loads((run / 'profile.json').read_text())
        request = json.loads((run / 'request.json').read_text())
        original = (run / 'GZWJL.PARAM').read_text()
        values = versions.parameter_values(original)
        old_end = float(values['TimeMax'])
        if time_max == old_end and status.get('continuation') and status['state'] in ('submitted', 'running', 'complete'):
            return {'directory': str(run), 'TimeMax': time_max, 'rerun': '1',
                    'submitted': True, 'reused': True, 'job_id': status['job_id']}
        retry = status['state'] == 'failed' and time_max == old_end
        if time_max is None or not (retry or old_end < time_max < float('inf')):
            raise ValueError('The continuation TimeMax must extend the existing endpoint')
        plan = {'directory': str(run), 'TimeMax': time_max, 'rerun': '1',
                'output_directory_unchanged': True, 'submitted': False}
        if not submit:
            return plan
        if status['state'] not in ('complete', 'failed'):
            raise ValueError('Stop/wait for the current solver before resuming it')
        queue = subprocess.run(['squeue', '-M', 'all', '--me', '-h', '-o', '%Z'],
                               capture_output=True, text=True, check=True)
        if any(Path(line.strip()) == run or run in Path(line.strip()).parents
               for line in queue.stdout.splitlines() if line.strip().startswith('/')):
            raise ValueError('An active job still uses this run')
        output_match = re.search(r'^\s*OutputDir\s+(\S+)', original, re.M)
        output = run / output_match[1] if output_match else run / 'output'
        output = runtime.confined(root, output)
        checkpoints = output / 'restartfiles'
        for rank in range(int(settings['sbatch']['ntasks'])):
            if not (checkpoints / (values.get('RestartFile', 'restart') + '.' + str(rank))).is_file():
                raise ValueError('Missing binary rank checkpoint: ' + str(rank))
        number = len(list(run.glob('.resume-*'))) + 1
        attempt = run / ('.resume-%03d' % number)
        attempt.mkdir()
        for source in ('workflow_isolation.py', 'workflow_runtime.py', 'workflow_versions.py'):
            shutil.copy2(Path(__file__).with_name(source), attempt / source)
        shutil.copy2(run / 'GZWJL.PARAM', attempt / 'previous.PARAM')
        rewrite_parameter_text(original, attempt / 'GZWJL.PARAM', {'TimeMax': str(time_max)})
        for suffix in ('json', 'csv'):
            report = run / ('rotation_evolution.' + suffix)
            if report.is_file():
                shutil.copy2(report, attempt / report.name)
        settings = dict(settings, parameter_overrides=dict(
            settings.get('parameter_overrides', {}), TimeMax=str(time_max)))
        signature_settings = dict(settings, rerun='1')
        restart_time = float(values['TimeMax'])
        info = output / 'info.txt'
        if info.is_file():
            with info.open('rb') as stream:
                stream.seek(max(0, info.stat().st_size - 16384))
                times = re.findall(r'Sync-Point\s+\d+,\s*Time:\s*([0-9.eE+-]+)',
                                   stream.read().decode(errors='replace'))
            if times:
                restart_time = float(times[-1])
        if '--run-solver' in settings['solver']:
            checkout = root / runtime.config(root)['github']['checkout']
            shutil.copy2(checkout / 'local_tools/ic_generation/rotation_diagnostics.py',
                         attempt / 'rotation_diagnostics.py')
            settings['solver'] = list(settings['solver'])
            settings['solver'][1] = str(attempt / 'rotation_diagnostics.py')
            settings['solver'] += ['--rerun', '1', '--restart-time', str(restart_time)]
        else:
            settings['solver'] = list(settings['solver'])
            if settings['solver'][-1] not in ('1', '{rerun}'):
                settings['solver'].append('1')
        parameter = versions.continuation(version_directory(root, model), request, signature_settings,
                                         (attempt / 'GZWJL.PARAM').read_text(), dict(status))
        runtime.write_json(attempt / 'profile.json', settings)
        command = [settings.get('python', 'python3'), str(attempt / 'workflow_isolation.py'),
                   'run-job', '--root', str(root), '--model', model, '--name', name,
                   '--resume', str(attempt)]
        script = '#!/bin/bash\nset -euo pipefail\ncd ' + shlex.quote(str(run)) + '\n'
        script += 'export TMPDIR="$PWD/.tmp"\nmkdir -p "$TMPDIR"\n'
        script += settings.get('setup', '') + '\nrerun="1"\nexec ' + shlex.join(command) + '\n'
        (attempt / 'run.sh').write_text(script)
        batch = ['sbatch', '--parsable', '--chdir', str(run), '--job-name', 'resume-' + name,
                 '--output', str(run / 'slurm-%j.out')]
        for key, value in settings['sbatch'].items():
            batch += ['--' + key, str(value)]
        batch.append(str(attempt / 'run.sh'))
        for key in ('started_utc', 'finished_utc', 'error', 'completion'):
            status.pop(key, None)
        status.update(state='submitting', continuation=str(attempt), parameter_version=parameter)
        runtime.write_json(run / 'status.json', status)
        result = subprocess.run(batch, capture_output=True, text=True)
        runtime.write_json(attempt / 'submission.json', {
            'command': batch, 'returncode': result.returncode, 'stdout': result.stdout, 'stderr': result.stderr})
        if result.returncode:
            status.update(state='failed', error=result.stderr)
        else:
            job = result.stdout.strip().split(';')[0]
            if not job.isdigit():
                raise ValueError('Uncertain resume acknowledgement; reconcile before retrying')
            shutil.copy2(attempt / 'GZWJL.PARAM', run / 'GZWJL.PARAM')
            status.update(state='submitted', job_id=job)
        runtime.write_json(run / 'status.json', status)
        versions.update(version_directory(root, model), model, name, status)
        if result.returncode:
            raise ValueError('Slurm rejected continuation: ' + result.stderr)
        return dict(plan, submitted=True, job_id=job, parameter_version=parameter)


def resume_client(args):
    remote = runtime.config(args.root)['remote']
    command = ['python3', remote['root'] + '/workflow_isolation.py', 'resume',
               '--root', remote['root'], '--model', args.model, '--name', args.name,
               '--time-max', str(args.time_max)]
    if args.submit:
        command.append('--submit')
    return remote_exec(args, command).returncode


def run_job(root, name, model=None, resume=None):
    run = case_dir(root, name, model)
    model = json.loads((run / 'request.json').read_text())['model']
    # Same lock as submitter: a fast-starting job must not lose its status update.
    with (run.parent / (run.name + '.lock')).open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        status = json.loads((run / 'status.json').read_text())
        if status['state'] != 'submitted' or status.get('job_id') != os.environ.get('SLURM_JOB_ID'):
            raise ValueError('Run is not awaiting this Slurm job; duplicate execution refused')
        status.update(state='running', started_utc=runtime.now())
        runtime.write_json(run / 'status.json', status)
        versions.update(version_directory(root, model), model, run.name, status)
    try:
        control = runtime.confined(run, Path(resume)) if resume else run
        settings = json.loads((control / 'profile.json').read_text())
        substitutions = {'ascii': str(run / 'Fornax.ini'), 'binary': str(run / 'Fornax.IC_000'),
                         'prefix': str(run / 'Fornax'), 'parameters': 'GZWJL.PARAM' if resume else str(run / 'GZWJL.PARAM'),
                         'rerun': '1' if resume else '0'}
        for stage in (('solver',) if resume else ('converter', 'solver')):
            command = [arg.format(**substitutions) for arg in settings[stage]]
            with (run / (stage + '.log')).open('a' if resume else 'w') as log:
                subprocess.run(command, cwd=run, stdout=log, stderr=subprocess.STDOUT, check=True)
            if stage == 'converter' and not (run / 'Fornax.IC_000').is_file():
                raise ValueError('Converter did not produce the expected binary IC')
        status.update(state='complete', completion='solver process returned successfully')
    except Exception as error:
        status.update(state='failed', error=str(error))
        runtime.write_json(run / 'status.json', status)
        versions.update(version_directory(root, model), model, run.name, status)
        raise
    status['finished_utc'] = runtime.now()
    runtime.write_json(run / 'status.json', status)
    versions.update(version_directory(root, model), model, run.name, status)


def prepare_legacy(root, directory, profile_name='default'):
    """Prepare a new, manually created flat isolation model; never overwrite data."""
    root = root.resolve()
    run = runtime.confined(root, directory.resolve())
    base = runtime.confined(root, root / runtime.config(root).get('isolation_root', 'isoFornax'))
    if run.parent != base or not re.fullmatch(r'Fornax\d+[A-Za-z0-9_.-]*', run.name):
        raise ValueError('Manual preparation requires isoFornax/FornaxNUMBER as working directory')
    settings = profile(root, profile_name)
    source = run / (run.name + '.ini')
    parameters = {}
    for line in (run / 'GZWJL.PARAM').read_text().splitlines():
        fields = line.split('%', 1)[0].split()
        if len(fields) >= 2:
            parameters[fields[0]] = fields[1]
    if (run / parameters.get('InitCondFile', '')).resolve() != run / 'Fornax.IC_000':
        raise ValueError('GZWJL.PARAM must use Fornax.IC_000 for this workflow')
    output = runtime.confined(run, run / parameters['OutputDir'])
    if (output.exists() and any(output.iterdir())) or any(run.glob('snapshot*.hdf5')):
        raise ValueError('Model already has simulation output; preparation refused')
    lock = runtime.confined(run, run / '.prepare-ic.lock')
    with lock.open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        record = run / 'ic-prepared.json'
        binary = run / 'Fornax.IC_000'
        alias = run / 'Fornax.ini'
        tree = run / 'TREECOOL'
        for path in (record, binary):
            if path.is_symlink():
                raise ValueError(f'Refusing symlink: {path.name}')
        if binary.exists():
            if not binary.is_file():
                raise ValueError('Existing IC is not a regular file')
            output.mkdir(exist_ok=True)
            return {'state': 'prepared', 'model': run.name, 'directory': str(run),
                    'input': source.name, 'reused_existing_binary': True}
        if not source.is_file() or source.resolve().parent != run:
            raise ValueError(f'Upload the numbered particle file {source.name} into {run.name}')
        cooling = Path(settings['template']) / 'TREECOOL'
        if not tree.exists() and not cooling.is_file():
            raise ValueError('Missing TREECOOL in configured isolation template')
        if alias.exists() or alias.is_symlink():
            if not alias.is_file() or alias.resolve().parent != run:
                raise ValueError('Fornax.ini must name a regular input within this model')
        if tree.exists() and not tree.is_file():
            raise ValueError('Existing TREECOOL is not a regular file')
        if not alias.exists():
            alias.symlink_to(source.name)
        if not tree.exists():
            temporary_tree = run / '.TREECOOL.preparing'
            with temporary_tree.open('xb') as target, cooling.open('rb') as src:
                shutil.copyfileobj(src, target)
            os.replace(temporary_tree, tree)
        with tempfile.TemporaryDirectory(prefix='.ic-convert-', dir=run) as temporary:
            converted = Path(temporary) / 'Fornax.IC_000'
            substitutions = {'ascii': str(source), 'binary': str(converted),
                             'prefix': str(Path(temporary) / 'Fornax'), 'parameters': str(run / 'GZWJL.PARAM')}
            command = [arg.format(**substitutions) for arg in settings['converter']]
            with (run / 'converter.log').open('w') as log:
                subprocess.run(command, cwd=run, check=True, stdout=log, stderr=subprocess.STDOUT)
            if not converted.is_file():
                raise ValueError('Converter did not produce the expected IC')
            if binary.exists() or binary.is_symlink():
                raise ValueError('Existing binary IC is protected')
            converted.rename(binary)
        runtime.write_json(record, {'model': run.name, 'input': source.name,
                                   'converter': settings['converter'], 'created_utc': runtime.now()})
        output.mkdir(exist_ok=True)
        return {'state': 'prepared', 'model': run.name, 'directory': str(run), 'input': source.name}


def remote_exec(args, command, capture=False):
    import workflow_management as management
    settings, ssh = management.ssh_config(args)
    setup = '. ' + shlex.quote(settings['module_init']) + ' && module load ' + shlex.join(settings['modules']) + ' && '
    return subprocess.run(ssh + [settings['host'], setup + shlex.join(command)], check=True,
                          text=True, stdout=subprocess.PIPE if capture else None)


def pipeline(args):
    import workflow_models as models
    selected = [args.model] if isinstance(args.model, str) else args.model
    if len(set(selected)) != len(selected):
        raise ValueError('Duplicate models in isolation batch')
    for model in selected:
        models.model_name(model)
    for model in selected:
        result = pipeline_one(SimpleNamespace(**{**vars(args), 'model': model}))
        if result:
            return result
    return 0


def client_directory(root, model, name):
    # Both local entry workspaces share the configured model storage and receipts.
    base = runtime.model_directory(root, model).parent
    return runtime.confined(base, base / 'workspace/receipts/isolation' / model / runtime.label(name))


def pipeline_one(args):
    import workflow_models as models
    import workflow_management as management
    root, _ = models.ic_paths(args)
    models.model_name(args.model)
    runtime.label(args.name)
    remote = runtime.config(args.root)['remote']
    if getattr(args, 'existing', False) and args.dry_run:
        print(f'Existing IC: {runtime.model_directory(args.root, args.model) / "isolation/Fornax.ini"} (reuse; do not generate)\n'
              f'Server: {remote["root"]}/{remote.get("isolation_root", "isoFornax")}/{args.model}/runs/{args.name}; submit={args.submit}')
        return 0
    if args.dry_run:
        print(f'IC: {models.find_initial(args, args.ic_run, args.model) if args.ic_run else models.initial_directory(args, args.model, args.name)}\n'
              f'Server: {remote["root"]}/{remote.get("isolation_root", "isoFornax")}/{args.model}/runs/{args.name}\n'
              f'Profile: {args.profile}; submit={args.submit}')
        return 0
    model_root = runtime.model_directory(args.root, args.model)
    # Transport bookkeeping belongs outside the browsable model results.
    local = client_directory(args.root, args.model, args.name)
    runtime.confined(model_root, model_root / 'isolation').mkdir(parents=True, exist_ok=True)
    local.mkdir(parents=True, exist_ok=True)
    with (local / 'pipeline.lock').open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        if getattr(args, 'existing', False):
            run = model_root / 'isolation'
            input_paths = {'Fornax.ini': run / 'Fornax.ini', 'parameters.ini': run / ('IC_' + args.model + '.ini')}
            receipt = {'state': 'existing', 'build_id': 'existing:' + str(run.resolve()), 'model': args.model,
                       'source_directory': str(run),
                       'inputs': {name: str(path) for name, path in input_paths.items()}}
            print(f"{args.model}: reuse existing IC from {run}", flush=True)
        else:
            run = models.find_initial(args, args.ic_run, args.model) if args.ic_run else models.initial_directory(args, args.model, args.name)
            if not run.exists():
                if args.ic_run:
                    raise ValueError('Requested existing IC run does not exist')
                result = models.build_ic(SimpleNamespace(**{**vars(args), 'ini': None}))
                if result:
                    return result
            receipt = json.loads((run / 'build.json').read_text())
            if receipt['state'] != 'complete' or Path(receipt['source_ini']).stem != 'IC_' + args.model:
                raise ValueError('Existing run is incomplete or belongs to a different model')
            if not args.ic_run:
                if receipt.get('seed') != args.seed:
                    raise ValueError('Requested seed differs from the named IC run; choose its matching run or seed')
            input_paths = {name: run / name for name in ('Fornax.ini', 'parameters.ini')}
        request = {'name': args.name, 'model': args.model, 'profile': args.profile,
                   'build_id': receipt.get('build_id', run.name),
                   'source_directory': str(run),
                   'files': {name: path.name for name, path in input_paths.items()}}
        if getattr(args, 'ic_version', None):
            request['source_version'] = runtime.label(args.ic_version)
        saved = local / 'request.json'
        resolution = remote_exec(args, ['python3', remote['root'] + '/workflow_isolation.py',
                                       'resolve', '--root', remote['root'],
                                       '--request-json', json.dumps(request)], capture=True)
        selected = json.loads(resolution.stdout)
        request = selected['request']
        remote_case = selected['directory']
        print(f"{args.model}: {request['name']} / "
              f"{request.get('ic_version') or 'legacy-unversioned'} / "
              f"{request.get('parameter_version') or 'legacy-unversioned'}", flush=True)
        if selected['existing']:
            command = ['python3', remote['root'] + '/workflow_isolation.py', 'receive',
                       '--root', remote['root'], '--stage', remote_case]
            if args.submit:
                command.append('--submit')
            result = remote_exec(args, command, capture=True)
            status = json.loads(result.stdout)
            runtime.write_json(saved, request)
            runtime.write_json(local / 'server.json', status)
            print(json.dumps(status, indent=2))
            return 0
        # Unique staging avoids modifying a previously submitted run on retries.
        remote_stage = remote['root'] + '/workspace/transfers/isolation/' + args.model + '-' + args.name + '-' + uuid.uuid4().hex[:12]
        create = ('from pathlib import Path; r=Path(' + repr(remote['root']) + ').resolve(); '
                  'p=Path(' + repr(remote_stage) + '); assert p.resolve().is_relative_to(r); '
                  'p.mkdir(parents=True,exist_ok=False)')
        remote_exec(args, ['python3', '-c', create])
        with tempfile.TemporaryDirectory(prefix='upload-', dir=local) as temporary:
            stage = Path(temporary)
            runtime.write_json(stage / 'request.json', request)
            for name in request['files']:
                shutil.copy2(input_paths[name], stage / name)
            for name in ('workflow_isolation.py', 'workflow_runtime.py', 'workflow_versions.py'):
                shutil.copy2(Path(__file__).with_name(name), stage / name)
            _, ssh = management.ssh_config(args)
            subprocess.run(['rsync', '-rt', '-e', shlex.join(ssh), str(stage) + '/',
                            remote['host'] + ':' + remote_stage + '/'], check=True)
        command = ['python3', remote['root'] + '/workflow_isolation.py', 'receive',
                   '--root', remote['root'], '--stage', remote_stage]
        if args.submit:
            command.append('--submit')
        result = remote_exec(args, command, capture=True)
        status = json.loads(result.stdout)
        runtime.write_json(saved, request)
        runtime.write_json(local / 'server.json', status)
        print(json.dumps(status, indent=2))
        return 0


def status(args):
    remote = runtime.config(args.root)['remote']
    runtime.label(args.name)
    model = getattr(args, 'model', None)
    if model:
        runtime.model_directory(args.root, model)
        saved = client_directory(args.root, model, args.name) / 'server.json'
        if saved.is_file():
            selected = json.loads(saved.read_text())
            args = SimpleNamespace(**{**vars(args), 'name': selected.get('name') or args.name})
    command = ['python3', remote['root'] + '/workflow_isolation.py', 'status',
               '--root', remote['root'], '--name', args.name]
    if model:
        command += ['--model', model]
    return remote_exec(args, command).returncode


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=('receive', 'run-job', 'resume', 'locate', 'resolve', 'status', 'prepare-legacy'))
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--stage', type=Path)
    parser.add_argument('--directory', type=Path, default=Path.cwd())
    parser.add_argument('--profile', default='default')
    parser.add_argument('--name')
    parser.add_argument('--model')
    parser.add_argument('--submit', action='store_true')
    parser.add_argument('--resume', type=Path)
    parser.add_argument('--time-max', type=float)
    parser.add_argument('--request-json', help='declared IC build and requested run; no particle metadata')
    args = parser.parse_args()
    if args.operation == 'prepare-legacy':
        print(json.dumps(prepare_legacy(args.root.resolve(), args.directory, args.profile)))
    elif args.operation == 'receive':
        print(json.dumps(receive(args.root.resolve(), args.stage, args.submit)))
    elif args.operation == 'run-job':
        run_job(args.root.resolve(), args.name, args.model, args.resume)
    elif args.operation == 'resume':
        print(json.dumps(resume(args.root.resolve(), args.name, args.model, args.time_max, args.submit)))
    elif args.operation == 'resolve':
        print(json.dumps(resolve_request(args.root.resolve(), json.loads(args.request_json))))
    else:
        run = case_dir(args.root.resolve(), args.name, args.model)
        if args.operation == 'locate':
            print(json.dumps({'directory': str(run)}))
        else:
            print((run / 'status.json').read_text())


if __name__ == '__main__':
    main()
