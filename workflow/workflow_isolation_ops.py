"""Post-processing and conservative removal of failed isolation startup attempts."""
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

import workflow_runtime as runtime
import workflow_isolation as isolation


def accounting(job, cluster):
    if not str(job).isdigit():
        raise ValueError('Invalid Slurm job ID')
    result = subprocess.run(['sacct', '-M', cluster, '-j', str(job), '-n', '-P',
                             '--format=JobIDRaw,State,ElapsedRaw,ExitCode'],
                            capture_output=True, text=True, check=True)
    rows = [line.split('|') for line in result.stdout.splitlines()]
    rows = [r for r in rows if r[0] == str(job)]
    if len(rows) != 1:
        raise ValueError('Slurm accounting is unavailable or ambiguous; retry later')
    row = rows[0]
    return {'state': row[1].split()[0], 'elapsed_seconds': int(row[2]), 'exit_code': row[3]}


def active_directories():
    result = subprocess.run(['squeue', '-M', 'all', '--me', '--noheader', '--format=%Z'],
                            capture_output=True, text=True, check=True)
    return [Path(line.strip()).resolve() for line in result.stdout.splitlines()
            if line.strip().startswith('/')]


def tree_bytes(path):
    return sum(p.stat().st_size for base, _, files in os.walk(path)
               for p in (Path(base) / name for name in files) if not p.is_symlink())


def cleanup(root, model, keep, dry_run=False):
    """Keep a verified replacement; remove only short, terminal memory/segfault failures."""
    survivor = isolation.case_dir(root, keep, model)
    request = json.loads((survivor / 'request.json').read_text())
    status = json.loads((survivor / 'status.json').read_text())
    if status['state'] not in ('submitted', 'running', 'complete'):
        raise ValueError('Replacement must be submitted, running or complete')
    isolation.verify_bundle(survivor, request)
    isolation.verify_case(survivor)
    active = active_directories()  # A failed query aborts deletion.
    report = {'kept': str(survivor), 'removed': [], 'skipped': [], 'bytes': 0, 'dry_run': dry_run}
    if status['state'] != 'complete':
        settings = json.loads((survivor / 'profile.json').read_text())
        record = accounting(status['job_id'], settings['sbatch']['clusters'])
        progressing = False
        if (survivor / 'solver.log').exists():
            with (survivor / 'solver.log').open(errors='replace') as log:
                progressing = any(re.search(r'Sync-Point [1-9]\d*, Time:', line) for line in log)
        if record['state'] != 'RUNNING' or record['elapsed_seconds'] < 140 or not progressing:
            return dict(report, deferred='Waiting for 140 seconds of running and integration progress')
    for run in sorted(survivor.parent.iterdir()):
        if run == survivor or not run.is_dir() or run.is_symlink():
            continue
        runtime.confined(root, run)
        lock = run.parent / (run.name + '.lock')
        with lock.open('a') as stream:
            fcntl.flock(stream, fcntl.LOCK_EX)
            try:
                old = json.loads((run / 'request.json').read_text())
                state = json.loads((run / 'status.json').read_text())
                if old['model'] != model or old['files'] != request['files']:
                    continue
                if state['state'] != 'failed':
                    continue
                if any(p == run.resolve() or run.resolve() in p.parents for p in active):
                    raise ValueError('An active Slurm job uses this directory')
                settings = json.loads((run / 'profile.json').read_text())
                record = accounting(state['job_id'], settings['sbatch']['clusters'])
                if record['state'] not in ('FAILED', 'OUT_OF_MEMORY') or record['elapsed_seconds'] > 140:
                    raise ValueError('Not a terminal startup failure within 140 seconds')
                log = (run / 'solver.log').read_text(errors='replace') if (run / 'solver.log').exists() else ''
                memory_failure = record['state'] == 'OUT_OF_MEMORY' or re.search(
                    r'EXIT CODE: 11|out of memory|cannot allocate memory|failed to allocate|not enough memory', log, re.I)
                if not memory_failure:
                    raise ValueError('No memory/segfault failure evidence')
                isolation.verify_bundle(run, old)
                size = tree_bytes(run)
                tombstone = runtime.confined(root, root / 'workspace/receipts/isolation-cleaned' / model / (run.name + '.json'))
                if not dry_run:
                    tombstone.parent.mkdir(parents=True, exist_ok=True)
                    runtime.write_json(tombstone, {'removed_utc': runtime.now(), 'request': old,
                        'status': state, 'accounting': record, 'replacement': str(survivor), 'bytes': size})
                    shutil.rmtree(run)
                report['removed'].append(str(run)); report['bytes'] += size
            except (ValueError, KeyError, OSError) as error:
                report['skipped'].append({'directory': str(run), 'reason': str(error)})
    # Completed uploads duplicate the exact IC already verified in the retained case.
    transfers = runtime.confined(root, root / 'workspace/transfers/isolation')
    if transfers.exists():
        for stage in transfers.iterdir():
            if not stage.is_dir() or stage.is_symlink() or not (stage / 'request.json').is_file():
                continue
            old = json.loads((stage / 'request.json').read_text())
            if old.get('model') != model or old.get('files') != request['files']:
                continue
            # Never remove a stage for a different pending/prepared attempt.
            target = isolation.case_dir(root, old['name'], model)
            tombstone = root / 'workspace/receipts/isolation-cleaned' / model / (old['name'] + '.json')
            if target != survivor and str(target) not in report['removed'] and not tombstone.is_file():
                continue
            isolation.verify_bundle(stage, old)
            size = tree_bytes(stage)
            if not dry_run:
                shutil.rmtree(stage)
            report['removed'].append(str(stage)); report['bytes'] += size
    return report


def density_plan(root, model, name, selection):
    run = isolation.case_dir(root, name, model)
    if json.loads((run / 'status.json').read_text())['state'] != 'complete':
        raise ValueError('Density checks require a completed isolation run')
    isolation.verify_completion(run)
    settings = runtime.config(root)['isolation_check']
    script = runtime.confined(root, root / settings['script']).resolve()
    snapshots = sorted(int(p.stem.split('_')[1]) for p in (run / 'output').glob('snapshot_*.hdf5'))
    if selection == 'quick':
        selection = ','.join(str(n) for n in sorted({snapshots[0], snapshots[len(snapshots)//2], snapshots[-1]}))
    if selection != 'all' and not re.fullmatch(r'\d+(?:-\d+|(?:,\d+)*)', selection):
        raise ValueError('Use quick, all, 0-600 or 0,300,600')
    directory = runtime.confined(root, root / 'workspace/receipts/isolation-check' / model / name)
    directory.mkdir(parents=True, exist_ok=True)
    identity = {'selection': selection, 'script_sha256': isolation.sha(script)}
    token = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:12]
    output = runtime.confined(run, run / 'density_output')
    command = ['python3', str(script), '--model-dir', str(run / 'output'),
               '--output-dir', str(output), '--processes', str(settings.get('processes', 16))]
    if selection != 'all': command += ['--range', selection, '--no-video']
    shell = 'set -eu\ncd ' + shlex.quote(str(run)) + '\n' + settings['setup'] + '\n'
    shell += 'export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1\n'
    shell += shlex.join(command)
    if selection == 'all':
        verify_video = ('from pathlib import Path; import json,subprocess; p=Path(' + repr(str(output)) +
            '); videos=list(p.glob("*.mp4")); assert videos, "Density video was not generated"; '
            'v=max(videos,key=lambda v:v.stat().st_mtime); '
            'info=json.loads(subprocess.check_output(["ffprobe","-v","error","-select_streams","v:0",'
            '"-show_entries","stream=nb_frames","-of","json",str(v)])); '
            'assert any(int(s.get("nb_frames",0)) == ' + str(len(snapshots)) +
            ' for s in info["streams"]), "Density video frame count disagrees with snapshots"')
        shell += '\n' + shlex.join(['python3', '-c', verify_video])
    return {'model': model, 'run': run, 'directory': directory, 'receipt': directory / (token + '.json'),
            'identity': identity, 'shell': shell}


def batch_time(value, count):
    match = re.fullmatch(r'(?:(\d+)-)?(\d+):(\d+):(\d+)', str(value))
    if not match: raise ValueError('Density time must use HH:MM:SS or D-HH:MM:SS')
    days, hours, minutes, seconds = [int(x or 0) for x in match.groups()]
    seconds = (((days * 24 + hours) * 60 + minutes) * 60 + seconds) * count
    return f'{seconds // 3600:02}:{seconds // 60 % 60:02}:{seconds % 60:02}'


def check_many(root, models, name, selection='all'):
    from contextlib import ExitStack
    if not models or len(models) != len(set(models)):
        raise ValueError('Provide distinct model IDs')
    plans = [density_plan(root, model, name, selection) for model in models]
    settings = runtime.config(root)['isolation_check']
    cluster = settings['sbatch']['clusters']
    existing, pending = [], []
    with ExitStack() as stack:
        # Lock all models in a stable order, including against single-model requests.
        for plan in sorted(plans, key=lambda p: p['model']):
            stream = stack.enter_context((plan['directory'] / 'check.lock').open('a'))
            fcntl.flock(stream, fcntl.LOCK_EX)
        active = active_directories()
        for plan in plans:
            receipt, run = plan['receipt'], plan['run']
            if receipt.exists():
                previous = json.loads(receipt.read_text())
                if not previous.get('job_id'):
                    raise ValueError('Uncertain density submission; reconcile Slurm before retrying')
                previous['scheduler'] = accounting(previous['job_id'], cluster)
                previous.setdefault('model', plan['model'])
                existing.append(previous)
                continue
            if any(p == run.resolve() or run.resolve() in p.parents for p in active):
                raise ValueError('A job is already using this run; wait before density submission')
            # Batched jobs run in a shared directory, so also inspect each model's receipts.
            for old_receipt in plan['directory'].glob('*.json'):
                old = json.loads(old_receipt.read_text())
                if not old.get('job_id'):
                    raise ValueError('Uncertain density submission; reconcile Slurm before retrying')
                state = accounting(old['job_id'], cluster)['state']
                if state not in {'COMPLETED', 'FAILED', 'CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY', 'NODE_FAIL', 'PREEMPTED', 'BOOT_FAIL', 'DEADLINE'}:
                    raise ValueError('An active batch already includes this model; wait before submitting another selection')
            pending.append(plan)
        if not pending:
            return {'job_id': None, 'submitted_models': [], 'existing': existing}
        options = dict(settings['sbatch'])
        allowed = {'clusters', 'qos', 'partition', 'account', 'nodes', 'ntasks', 'cpus-per-task', 'time', 'mem'}
        if set(options) - allowed: raise ValueError('Unsupported density Slurm option')
        if 'time' in options: options['time'] = batch_time(options['time'], len(pending))
        if len(pending) == 1:
            work = pending[0]['run']
            shell = pending[0]['shell']
        else:
            token = hashlib.sha256(json.dumps([(p['model'], name, p['identity']) for p in pending], sort_keys=True).encode()).hexdigest()[:16]
            work = runtime.confined(root, root / 'workspace/receipts/isolation-check/batches' / token)
            work.mkdir(parents=True, exist_ok=True)
            shell = 'set -eu\n'
            for plan in pending:
                logfile = shlex.quote(str(plan['run'] / 'density-')) + '"${SLURM_JOB_ID}".log'
                shell += '(\n' + plan['shell'] + '\n) > ' + logfile + ' 2>&1\n'
        job_name = 'density-' + pending[0]['model'] if len(pending) == 1 else 'density-batch-' + str(len(pending))
        batch = ['sbatch', '--parsable', '--job-name', job_name, '--chdir', str(work), '--output', str(work / 'density-%j.log')]
        for key, value in options.items(): batch += ['--' + key, str(value)]
        batch += ['--wrap', shell]
        records = []
        for plan in pending:
            result = dict(plan['identity'], model=plan['model'], state='submitting', created_utc=runtime.now(),
                          directory=str(plan['run']), batch_directory=str(work), command=batch)
            runtime.write_json(plan['receipt'], result)
            records.append(result)
        submitted = subprocess.run(batch, capture_output=True, text=True)
        if submitted.returncode:
            for plan in pending: plan['receipt'].unlink()
            raise ValueError(submitted.stderr.strip())
        job = submitted.stdout.strip().split(';')[0]
        if not job.isdigit(): raise ValueError('Unrecognized density job ID; reconcile Slurm')
        for plan, result in zip(pending, records):
            result.update(state='submitted', job_id=job)
            runtime.write_json(plan['receipt'], result)
        return {'job_id': job, 'submitted_models': [p['model'] for p in pending], 'runs': records, 'existing': existing}


def check(root, model, name, selection='all'):
    result = check_many(root, [model], name, selection)
    return (result.get('runs') or result['existing'])[0]


def results(root, model, name):
    run = isolation.case_dir(root, name, model)
    if any(p == run.resolve() or run.resolve() in p.parents for p in active_directories()):
        raise ValueError('Run still has an active job; wait before fetching density results')
    output = runtime.confined(run, run / 'density_output')
    if not list(output.glob('snapshot_*_density_profile.png')):
        raise ValueError('No density figures available; run isolation-check first')
    directory = root / 'workspace/receipts/isolation-check' / model / name
    receipts = sorted(directory.glob('*.json'), key=lambda p: p.stat().st_mtime)
    if not receipts: raise ValueError('No density job receipt')
    latest = json.loads(receipts[-1].read_text())
    settings = runtime.config(root)['isolation_check']
    latest['scheduler'] = accounting(latest['job_id'], settings['sbatch']['clusters'])
    if latest['scheduler']['state'] != 'COMPLETED':
        raise ValueError('Density analysis has not completed successfully: ' + str(latest['scheduler']))
    return {'directory': str(output), 'analysis': latest}


def client(args):
    import workflow_management as management
    from workflow_models import model_name
    runtime.label(args.name)
    remote, ssh = management.ssh_config(args)
    if args.command == 'isolation-check':
        for model in args.model:
            model_name(model)
        command = ['python3', remote['root'] + '/workflow_isolation_ops.py', 'check-batch',
                   '--root', remote['root'], '--models', *args.model, '--name', args.name, '--range', args.range]
        response = isolation.remote_exec(args, command, capture=True)
        result = json.loads(response.stdout)
        summary = {key: result[key] for key in ('job_id', 'submitted_models')}
        summary['existing'] = [{key: row.get(key) for key in ('model', 'job_id', 'scheduler')} for row in result['existing']]
        print(json.dumps(summary, indent=2), flush=True)
        return 0
    for model in args.model:
        model_name(model)
        operation = {'isolation-clean-failed': 'cleanup', 'isolation-check': 'check', 'isolation-fetch': 'results'}[args.command]
        command = ['python3', remote['root'] + '/workflow_isolation_ops.py', operation,
                   '--root', remote['root'], '--model', model, '--name', args.name]
        if operation == 'cleanup' and args.dry_run: command += ['--dry-run']
        if operation == 'check': command += ['--range', args.range]
        response = isolation.remote_exec(args, command, capture=True)
        result = json.loads(response.stdout)
        print(model + ': ' + json.dumps(result, indent=2), flush=True)
        if operation == 'results':
            model_root = runtime.model_directory(args.root, model)
            destination = runtime.confined(model_root, model_root / 'isolation/density_output')
            source_receipt = isolation.client_directory(args.root, model, args.name).parent / 'density-source.json'
            identity = {'model': model, 'run': args.name, 'directory': result['directory']}
            if source_receipt.exists() and json.loads(source_receipt.read_text()) != identity:
                raise ValueError('Local density results belong to another run; move them before fetching a different run')
            if not source_receipt.exists() and destination.exists() and any(destination.iterdir()):
                raise ValueError('Local density results have unknown provenance; move them before fetching')
            receipt = isolation.client_directory(args.root, model, args.name) / 'density-result.json'
            receipt.parent.mkdir(parents=True, exist_ok=True)
            destination.mkdir(parents=True, exist_ok=True)
            runtime.write_json(source_receipt, identity)
            subprocess.run(['rsync', '-rt', '--checksum', '-e', shlex.join(ssh), '--include=*/',
                            '--include=*.png', '--include=*.mp4', '--include=*.csv', '--include=*.json',
                            '--include=*.txt', '--exclude=*', remote['host'] + ':' + result['directory'] + '/',
                            str(destination) + '/'], check=True)
            runtime.write_json(receipt, result)
            print('Saved density results: ' + str(destination))
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation', choices=['cleanup', 'check', 'check-batch', 'results'])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--model')
    parser.add_argument('--models', nargs='+')
    parser.add_argument('--name', required=True)
    parser.add_argument('--range', default='all')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    root = args.root.resolve()
    if args.operation == 'check-batch' and not args.models: parser.error('--models is required')
    if args.operation != 'check-batch' and not args.model: parser.error('--model is required')
    if args.operation == 'check-batch': value = check_many(root, args.models, args.name, args.range)
    elif args.operation == 'cleanup': value = cleanup(root, args.model, args.name, args.dry_run)
    elif args.operation == 'check': value = check(root, args.model, args.name, args.range)
    else: value = results(root, args.model, args.name)
    print(json.dumps(value))


if __name__ == '__main__':
    main()
