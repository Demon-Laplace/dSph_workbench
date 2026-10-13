"""Resumable local IC -> verified transfer -> isolated Slurm simulation."""
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
import sys
import tempfile
from types import SimpleNamespace
import uuid

import workflow_runtime as runtime


def sha(path):
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
    remaining = set(values)
    seen = set()
    lines = []
    for line in source.read_text().splitlines():
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


def verify_bundle(stage, request):
    if set(request['files']) != {'Fornax.ini', 'parameters.ini'}:
        raise ValueError('Unexpected IC transfer file set')
    for name, digest in request['files'].items():
        path = stage / name
        if path.is_symlink() or sha(path) != digest:
            raise ValueError(f'IC transfer checksum mismatch: {name}')


def case_dir(root, name, model=None):
    runtime.label(name)
    if model:
        # Receiver is confined to data even when local IC tools use models_root.
        return runtime.confined(root, runtime.model_directory(root, model) / 'isolation' / name)
    legacy = runtime.confined(root, root / 'workspace/models/isolation' / name)
    matches = [p for p in root.glob('*/isolation/' + name) if p.is_dir()]
    if legacy.is_dir():
        matches.append(legacy)
    unique = {p.resolve() for p in matches}
    if len(unique) != 1:
        raise ValueError('Specify --model for an absent or ambiguous isolation run')
    return runtime.confined(root, unique.pop())


def verify_case(run):
    for name, digest in json.loads((run / 'prepared.json').read_text()).items():
        path = runtime.confined(run, run / name)
        if path.is_symlink() or sha(path) != digest:
            raise ValueError(f'Prepared case changed: {name}')


def verify_completion(run):
    import h5py
    values = {}
    for line in (run / 'GZWJL.PARAM').read_text().splitlines():
        fields = line.split('%', 1)[0].split()
        if len(fields) == 2:
            values[fields[0]] = fields[1]
    times = []
    for path in (run / 'output').glob('snapshot*.hdf5'):
        with h5py.File(path, 'r') as snapshot:
            times.append(float(snapshot['Header'].attrs['Time']))
    target = float(values['TimeMax'])
    if not times or max(times) < target - max(abs(target) * 1e-6, 1e-12):
        raise ValueError('No HDF5 snapshot reached TimeMax; inspect solver.log (time limit or incomplete run)')
    return {'final_snapshot_time': max(times), 'target_time': target,
            'scope': 'time reached; equilibrium and convergence not evaluated'}


def receive(root, stage, submit=False):
    stage = runtime.confined(root, stage).resolve()
    request = json.loads((stage / 'request.json').read_text())
    verify_bundle(stage, request)
    run = case_dir(root, request['name'], request['model'])
    run.parent.mkdir(parents=True, exist_ok=True)
    lock = runtime.confined(root, run.parent / (run.name + '.lock'))
    with lock.open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        if run.exists():
            saved = json.loads((run / 'request.json').read_text())
            if saved != request:
                raise ValueError('Run name is already bound to different ICs or settings; use a new name')
        else:
            settings = profile(root, request['profile'])
            # Assemble privately, then atomically publish a complete case.
            with tempfile.TemporaryDirectory(prefix='isolation-', dir=run.parent) as temporary:
                prepared = Path(temporary) / 'case'
                prepared.mkdir()
                for name in request['files']:
                    shutil.copy2(stage / name, prepared / name)
                for name in ('workflow_isolation.py', 'workflow_runtime.py'):
                    shutil.copy2(stage / name, prepared / name)
                runtime.write_json(prepared / 'request.json', request)
                runtime.write_json(prepared / 'profile.json', settings)
                template = Path(settings['template'])
                values = dict(settings.get('parameter_overrides', {}))
                values.update(InitCondFile=str(run / 'Fornax.IC_000'), OutputDir=str(run / 'output'), ResubmitOn='0')
                rewrite_parameters(template / 'GZWJL.PARAM', prepared / 'GZWJL.PARAM', values)
                cooling = template / 'TREECOOL'
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
                pinned = {p.name: sha(p) for p in prepared.iterdir() if p.is_file()}
                runtime.write_json(prepared / 'prepared.json', pinned)
                runtime.write_json(prepared / 'status.json', {'state': 'prepared', 'created_utc': runtime.now()})
                prepared.rename(run)
        status = json.loads((run / 'status.json').read_text())
        if status['state'] in ('submitted', 'running', 'complete', 'failed'):
            return dict(status, directory=str(run))
        if status['state'] == 'submitting':
            raise ValueError('Previous submission outcome is uncertain; reconcile Slurm before retrying, to prevent duplicate jobs')
        # Do not overwrite an IC or submit into an existing simulation output.
        verify_bundle(run, request)
        verify_case(run)
        if any((run / 'output').iterdir()) or (run / 'Fornax.IC_000').exists():
            raise ValueError('Run already contains simulation data; new submission refused')
        settings = json.loads((run / 'profile.json').read_text())
        options = settings.get('sbatch', {})
        allowed = {'clusters', 'qos', 'partition', 'account', 'nodes', 'ntasks', 'cpus-per-task', 'time', 'mem', 'constraint'}
        if set(options) - allowed:
            raise ValueError('Unsupported Slurm profile option')
        command = ['sbatch', '--parsable', '--chdir', str(run), '--job-name', 'iso-' + run.name,
                   '--output', str(run / 'slurm-%j.out'), '--error', str(run / 'slurm-%j.out')]
        for key, value in options.items():
            command += ['--' + key, str(value)]
        if not submit:
            command.append('--test-only')
        command.append(str(run / 'run.sh'))
        if submit:
            status.update(state='submitting', submission_started_utc=runtime.now())
            runtime.write_json(run / 'status.json', status)
        result = subprocess.run(command, capture_output=True, text=True)
        runtime.write_json(run / 'submission.json', {'command': command, 'returncode': result.returncode,
                           'stdout': result.stdout, 'stderr': result.stderr})
        if result.returncode:
            # A nonzero sbatch exit is a definitive rejection, not a lost acknowledgement.
            status.update(state='prepared', last_error=result.stderr)
            runtime.write_json(run / 'status.json', status)
            raise ValueError('Slurm rejected the request: ' + result.stderr.strip())
        if submit:
            job = result.stdout.strip().split(';')[0]
            if not job.isdigit():
                raise ValueError('Unrecognized Slurm job ID; reconcile submission before retrying')
            status.update(state='submitted', job_id=job)
        else:
            status.update(state='prepared', scheduler_preflight='passed')
        runtime.write_json(run / 'status.json', status)
        return dict(status, directory=str(run))


def run_job(root, name, model=None):
    run = case_dir(root, name, model)
    # Same lock as submitter: a fast-starting job must not lose its status update.
    with (run.parent / (run.name + '.lock')).open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        status = json.loads((run / 'status.json').read_text())
        if status['state'] != 'submitted' or status.get('job_id') != os.environ.get('SLURM_JOB_ID'):
            raise ValueError('Run is not awaiting this Slurm job; duplicate execution refused')
        status.update(state='running', started_utc=runtime.now())
        runtime.write_json(run / 'status.json', status)
    try:
        verify_bundle(run, json.loads((run / 'request.json').read_text()))
        verify_case(run)
        settings = json.loads((run / 'profile.json').read_text())
        substitutions = {'ascii': str(run / 'Fornax.ini'), 'binary': str(run / 'Fornax.IC_000'),
                         'prefix': str(run / 'Fornax'), 'parameters': str(run / 'GZWJL.PARAM')}
        for stage in ('converter', 'solver'):
            command = [arg.format(**substitutions) for arg in settings[stage]]
            with (run / (stage + '.log')).open('w') as log:
                subprocess.run(command, cwd=run, stdout=log, stderr=subprocess.STDOUT, check=True)
            if stage == 'converter' and not (run / 'Fornax.IC_000').is_file():
                raise ValueError('Converter did not produce the expected binary IC')
        status.update(state='complete', verification=verify_completion(run))
    except Exception as error:
        status.update(state='failed', error=str(error))
        runtime.write_json(run / 'status.json', status)
        raise
    status['finished_utc'] = runtime.now()
    runtime.write_json(run / 'status.json', status)


def remote_exec(args, command, capture=False):
    import workflow_management as management
    settings, ssh = management.ssh_config(args)
    setup = '. ' + shlex.quote(settings['module_init']) + ' && module load ' + shlex.join(settings['modules']) + ' && '
    return subprocess.run(ssh + [settings['host'], setup + shlex.join(command)], check=True,
                          text=True, stdout=subprocess.PIPE if capture else None)


def pipeline(args):
    import workflow_models as models
    import workflow_management as management
    root, _ = models.ic_paths(args)
    models.model_name(args.model)
    runtime.label(args.name)
    remote = runtime.config(args.root)['remote']
    if args.dry_run:
        print(f'IC: {models.find_initial(args, args.ic_run, args.model) if args.ic_run else models.initial_directory(args, args.model, args.name)}\n'
              f'Server: {remote["root"]}/{args.model}/isolation/{args.name}\n'
              f'Profile: {args.profile}; submit={args.submit}')
        return 0
    model_root = runtime.model_directory(args.root, args.model)
    local = runtime.confined(model_root, model_root / 'isolation' / args.name / 'client')
    local.mkdir(parents=True, exist_ok=True)
    with (local / 'pipeline.lock').open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        # Verify server configuration before starting an expensive build.
        probe = ('import json,pathlib; p=pathlib.Path(' + repr(remote['root']) + '); '
                 's=json.loads((p/"workflow.json").read_text()); '
                 'assert ' + repr(args.profile) + ' in s.get("isolation_profiles",{}), "Configure the server isolation profile first"')
        remote_exec(args, ['python3', '-c', probe])
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
        if sha(run / 'parameters.ini') != receipt['ini_sha256'] or sha(run / 'Fornax.ini') != receipt['particle_sha256']:
            raise ValueError('IC data changed after generation')
        if not args.ic_run:
            current = root / 'IniFile' / ('IC_' + args.model + '.ini')
            if sha(current) != receipt['ini_sha256'] or receipt['seed'] != args.seed:
                raise ValueError('Parameters or seed changed; choose a new run name')
        if receipt.get('catalog_sync', {}).get('state') != 'pushed':
            result = models.ic_catalog(SimpleNamespace(**{**vars(args), 'run': run.name}))
            if result:
                return result
        request = {'name': args.name, 'model': args.model, 'profile': args.profile,
                   'build_id': receipt.get('build_id', run.name),
                   'files': {'Fornax.ini': receipt['particle_sha256'], 'parameters.ini': receipt['ini_sha256']}}
        saved = local / 'request.json'
        if saved.exists() and json.loads(saved.read_text()) != request:
            raise ValueError('Local run name is bound to different ICs or settings')
        runtime.write_json(saved, request)
        remote_case = remote['root'] + '/' + args.model + '/isolation/' + args.name
        existing = ('import json; from pathlib import Path; r=Path(' + repr(remote['root']) + ').resolve(); '
                    'p=Path(' + repr(remote_case) + '); assert p.resolve().is_relative_to(r); '
                    'expected=json.loads(' + repr(json.dumps(request)) + '); '
                    'assert not p.exists() or json.loads((p/"request.json").read_text())==expected, "Run name already bound to different ICs"; '
                    'print("existing" if p.exists() else "missing")')
        if remote_exec(args, ['python3', '-c', existing], capture=True).stdout.strip() == 'existing':
            command = ['python3', remote_case + '/workflow_isolation.py', 'receive',
                       '--root', remote['root'], '--stage', remote_case]
            if args.submit:
                command.append('--submit')
            result = remote_exec(args, command, capture=True)
            status = json.loads(result.stdout)
            runtime.write_json(local / 'server.json', status)
            print(json.dumps(status, indent=2))
            return 0
        # Unique staging avoids modifying a previously submitted run on retries.
        remote_stage = remote['root'] + '/workspace/transfers/isolation/' + args.name + '-' + uuid.uuid4().hex[:12]
        create = ('from pathlib import Path; r=Path(' + repr(remote['root']) + ').resolve(); '
                  'p=Path(' + repr(remote_stage) + '); assert p.resolve().is_relative_to(r); '
                  'p.mkdir(parents=True,exist_ok=False)')
        remote_exec(args, ['python3', '-c', create])
        with tempfile.TemporaryDirectory(prefix='upload-', dir=local) as temporary:
            stage = Path(temporary)
            shutil.copy2(saved, stage / 'request.json')
            for name in request['files']:
                shutil.copy2(run / name, stage / name)
            for name in ('workflow_isolation.py', 'workflow_runtime.py'):
                shutil.copy2(Path(__file__).with_name(name), stage / name)
            _, ssh = management.ssh_config(args)
            subprocess.run(['rsync', '-rt', '--checksum', '-e', shlex.join(ssh), str(stage) + '/',
                            remote['host'] + ':' + remote_stage + '/'], check=True)
        command = ['python3', remote_stage + '/workflow_isolation.py', 'receive',
                   '--root', remote['root'], '--stage', remote_stage]
        if args.submit:
            command.append('--submit')
        result = remote_exec(args, command, capture=True)
        status = json.loads(result.stdout)
        runtime.write_json(local / 'server.json', status)
        print(json.dumps(status, indent=2))
        return 0


def status(args):
    remote = runtime.config(args.root)['remote']
    runtime.label(args.name)
    model = getattr(args, 'model', None)
    if model:
        runtime.model_directory(args.root, model)
    relative = (model + '/isolation/' + args.name) if model else None
    script = ('from pathlib import Path; r=Path(' + repr(remote['root']) + ').resolve(); '
              'items=[r/' + repr(relative) + '] if ' + repr(relative) + ' else list(r.glob("*/isolation/' + args.name + '"))+list(r.glob("workspace/models/isolation/' + args.name + '")); '
              'items=list({p.resolve() for p in items if p.is_dir()}); '
              'assert len(items)==1, "Specify --model for an absent or ambiguous run"; '
              'p=items[0]/"status.json"; assert p.resolve().is_relative_to(r); print(p.read_text())')
    return remote_exec(args, ['python3', '-c', script]).returncode


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=('receive', 'run-job'))
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--stage', type=Path)
    parser.add_argument('--name')
    parser.add_argument('--model')
    parser.add_argument('--submit', action='store_true')
    args = parser.parse_args()
    if args.operation == 'receive':
        print(json.dumps(receive(args.root.resolve(), args.stage, args.submit)))
    else:
        run_job(args.root.resolve(), args.name, args.model)


if __name__ == '__main__':
    main()
