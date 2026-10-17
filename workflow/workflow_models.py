"""Prepare production cases and submit existing simulation launchers explicitly."""
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time

import workflow_runtime as runtime


def ic_paths(args):
    settings = runtime.config(args.root).get('ic', {})
    if not settings.get('root'):
        raise ValueError('Configure ic.root in workflow.json to the separate IC workspace')
    root = (args.root / settings['root']).resolve()
    if not root.is_dir():
        raise ValueError(f'Missing configured IC workspace: {root}')
    repo = (args.root / settings.get('repository', settings['root'])).resolve()
    return root, repo


def build_ic(args):
    root, repo = ic_paths(args)
    source = resolve_input(args.root, args.ini) if args.ini else resolve_input(root, f'IniFile/IC_{model_name(args.model)}.ini')
    model = model_name(args.model or source.stem.removeprefix('IC_'))
    destination = initial_directory(args, model, args.name)
    if destination.exists() or destination.is_symlink():
        raise ValueError('Initial-condition run already exists; refusing to overwrite')
    builder = args.repo / 'local_tools/ic_generation/build_ic.py'
    if not builder.is_file():
        raise ValueError(f'Missing versioned IC builder: {builder}')
    python = args.python or runtime.config(args.root).get('ic_python') or sys.executable
    return subprocess.run([python, str(builder), '--ini', str(source), '--output', str(destination),
                           '--model', model, '--threads', str(args.threads), '--seed', str(args.seed),
                           '--catalog-repo', str(repo)]).returncode


def initial_directory(args, model, name):
    directory = runtime.model_directory(args.root, model)
    return runtime.confined(directory, directory / 'ic' / runtime.label(name))


def find_initial(args, name, model=None):
    runtime.label(name)
    if model:
        current = initial_directory(args, model, name)
        if current.exists():
            return current
    else:
        base = (args.root / runtime.config(args.root).get('models_root', '.')).resolve()
        matches = [p for p in base.glob('*/ic/' + name) if p.is_dir()]
        if len(matches) > 1:
            raise ValueError('Multiple models have this run name; specify --model')
        if matches:
            return runtime.confined(base, matches[0])
    root, _ = ic_paths(args)
    return runtime.confined(root, root / 'generated' / name)


def ic_catalog(args):
    root, repo = ic_paths(args)
    command = [sys.executable, str(args.repo / 'local_tools/ic_generation/catalog.py'), '--repo', str(repo)]
    if args.run:
        command += ['--run', str(find_initial(args, args.run, getattr(args, 'model', None)))]
    else:
        command += ['--history-root', str(root)]
    return subprocess.run(command).returncode


def model_name(value):
    runtime.label(value)
    if not re.fullmatch(r'(?:Fornax|MW|Sculptor)\d+[A-Za-z0-9_.-]*', value):
        raise ValueError('Use a model name beginning with Fornax, MW or Sculptor and a number')
    return value


def resolve_input(root, value):
    path = (root / value).resolve()
    if not path.is_file():
        raise ValueError(f'Missing input: {path}')
    return path


def rewrite_first_line(path, transform):
    lines = path.read_text().splitlines(keepends=True)
    for index, line in enumerate(lines):
        tokens = shlex.split(line, comments=True)
        if tokens:
            lines[index] = ' '.join(transform(tokens)) + '\n'
            path.write_text(''.join(lines))
            return
    raise ValueError(f'No input references in {path}')


def new_case(args):
    name = model_name(args.model)
    if not name.startswith('Fornax'):
        raise ValueError('New-case preparation currently supports Fornax production models')
    template = runtime.model_directory(args.root, model_name(args.template)).resolve()
    dwarf = resolve_input(args.root, args.dwarf_ic)
    if any(c.isspace() for c in str(dwarf)):
        raise ValueError('Legacy IC merger requires paths without whitespace')
    destination = runtime.model_directory(args.root, name)
    if destination.exists() and (not destination.is_dir() or any(p.name not in {'ic', 'isolation'} for p in destination.iterdir())):
        raise ValueError('Model already exists; refusing to overwrite')
    required = [template / n for n in ('GZWJL.PARAM', 'mergegalaxy.par', 'mergecube.par')]
    if not all(p.is_file() for p in required):
        raise ValueError('Template lacks production parameters')
    # Resolve all references before creating anything.
    first = next(shlex.split(line, comments=True) for line in required[1].read_text().splitlines() if shlex.split(line, comments=True))
    if len(first) < 2:
        raise ValueError('Invalid mergegalaxy input references')
    mw = resolve_input(args.root, args.mw_ic) if args.mw_ic else resolve_input(template, first[0])
    cube_tokens = next(shlex.split(line, comments=True) for line in required[2].read_text().splitlines() if shlex.split(line, comments=True))
    cube = resolve_input(template, cube_tokens[0])
    destination.mkdir(parents=True, exist_ok=True)
    for path in required:
        shutil.copy2(path, destination / path.name)
    rewrite_first_line(destination / 'mergegalaxy.par', lambda tokens: [str(mw), str(dwarf), *tokens[2:]])
    rewrite_first_line(destination / 'mergecube.par', lambda tokens: [str(cube), *tokens[1:]])
    (destination / 'Fornax.ini').symlink_to(os.path.relpath(dwarf, destination))
    for name in ('subastrothin.sh', 'subastrofat.sh', 'submeso1.sh'):
        path = template / name
        if path.is_file():
            shutil.copy2(path, destination / path.name)
    if (template / 'TREECOOL').exists():
        (destination / 'TREECOOL').symlink_to((template / 'TREECOOL').resolve())
    runtime.write_json(destination / 'case.json', {'created_utc': runtime.now(), 'template': str(template),
                      'dwarf_ic': str(dwarf), 'mw_ic': str(mw), 'cube': str(cube),
                      'state': 'inputs_prepared', 'note': 'Orbit, orientation and physics copied from template. IC merging and simulation have not run.'})
    print(f'Prepared production inputs in {destination}; IC and isolation stages are preserved.')
    return 0


def simulate(args):
    if getattr(args, 'startup_timeout', 7200) <= 0:
        raise ValueError('--startup-timeout must be positive')
    model = runtime.model_directory(args.root, model_name(args.model)).resolve()
    launcher = model / args.launcher
    if args.launcher not in ('subastrothin.sh', 'subastrofat.sh', 'submeso1.sh') or not launcher.is_file():
        raise ValueError('Select an existing supported simulation launcher')
    parameters = {}
    for line in (model / 'GZWJL.PARAM').read_text().splitlines():
        fields = line.split('%', 1)[0].split()
        if len(fields) >= 2:
            parameters[fields[0]] = fields[1]
    initial = (model / parameters['InitCondFile']).resolve()
    output = runtime.confined(args.root, model / parameters['OutputDir'])
    if not initial.is_file():
        raise ValueError(f'IC missing: {initial}; run merge-ic {args.model} --name NAME --submit first')
    if list(output.glob('snapshot*')) or list(output.glob('restart*')):
        raise ValueError('Simulation already contains data; use its established restart procedure')
    if args.submit:
        from workflow_cleanup import active_jobs
        if any(directory == model for _, directory in active_jobs()):
            raise ValueError('An active Slurm job already uses this model')
    directory = runtime.confined(args.root, args.root / 'workspace/jobs' / runtime.label(args.name))
    if directory.exists():
        raise ValueError('Job name already used; choose a new --name')
    directory.mkdir(parents=True)
    script = directory / launcher.name
    shutil.copy2(launcher, script)
    command = ['sbatch', '--parsable', '--chdir', str(model), '--output', str(directory / 'slurm-%j.out')]
    if not args.submit:
        command.append('--test-only')
    command.append(str(script))
    print(shlex.join(command), flush=True)
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    runtime.write_json(directory / 'submission.json', {'time': runtime.now(), 'model': str(model),
                       'command': command, 'submitted': args.submit and result.returncode == 0,
                       'stdout': result.stdout, 'stderr': result.stderr, 'returncode': result.returncode,
                       'launcher_sha256': runtime.sha(script), 'parameter_sha256': runtime.sha(model / 'GZWJL.PARAM')})
    print(result.stdout, end='', flush=True); print(result.stderr, end='', flush=True)
    if args.submit and result.returncode == 0 and not getattr(args, 'no_wait_startup', False):
        return simulation_startup(args)
    return result.returncode


def simulation_startup(args):
    """Stop monitoring after 140 seconds of running and nonzero integration progress."""
    from workflow_isolation_ops import accounting
    timeout = getattr(args, 'startup_timeout', getattr(args, 'timeout', 7200))
    if timeout <= 0:
        raise ValueError('Startup timeout must be positive')
    model = runtime.model_directory(args.root, model_name(args.model)).resolve()
    directory = runtime.confined(args.root, args.root / 'workspace/jobs' / runtime.label(args.name))
    submitted = json.loads((directory / 'submission.json').read_text())
    if not submitted.get('submitted') or Path(submitted['model']).resolve() != model:
        raise ValueError('Startup monitoring requires a submitted job for this model')
    match = re.fullmatch(r'(\d+)(?:;([A-Za-z0-9_.-]+))?', submitted['stdout'].strip())
    if not match:
        raise ValueError('Unrecognized submitted job ID; reconcile Slurm before monitoring')
    job, cluster = match.group(1), match.group(2) or 'all'
    parameter = model / 'GZWJL.PARAM'
    if runtime.sha(parameter) != submitted['parameter_sha256']:
        raise ValueError('Simulation parameters changed since submission')
    parameters = {}
    for line in parameter.read_text().splitlines():
        fields = line.split('%', 1)[0].split()
        if len(fields) >= 2:
            parameters[fields[0]] = fields[1]
    output = runtime.confined(args.root, model / parameters['OutputDir'])
    deadline = time.monotonic() + timeout
    print(f'Automatically checking job {job}; exit monitoring after 140 running seconds and integration progress.', flush=True)
    terminal = {'COMPLETED', 'FAILED', 'CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY',
                'NODE_FAIL', 'PREEMPTED', 'BOOT_FAIL', 'DEADLINE', 'REVOKED'}
    while True:
        try:
            scheduler = accounting(job, cluster)
        except ValueError as error:
            if str(error) != 'Slurm accounting is unavailable or ambiguous; retry later':
                raise
            scheduler = {'state': 'ACCOUNTING_PENDING', 'elapsed_seconds': 0}
        progress = None
        info = output / 'info.txt'
        if info.is_file():
            for point in re.finditer(r'Sync-Point\s+(\d+),\s*Time:\s*([0-9.eE+-]+)', info.read_text(errors='replace')):
                if int(point.group(1)) > 0 and float(point.group(2)) > 0:
                    progress = {'sync_point': int(point.group(1)), 'simulation_time': float(point.group(2))}
        state = scheduler['state']
        record = {'job_id': job, 'cluster': cluster, 'model': str(model),
                  'checked_utc': runtime.now(), 'scheduler': scheduler, 'progress': progress,
                  'required_running_seconds': 140, 'state': 'waiting'}
        if state in terminal and (state != 'COMPLETED' or scheduler.get('exit_code') != '0:0'):
            record['state'] = 'failed'
            runtime.write_json(directory / 'startup.json', record)
            raise ValueError(f'Simulation job {job} failed during startup: {scheduler}')
        if state in {'RUNNING', 'COMPLETED'} and scheduler['elapsed_seconds'] >= 140 and progress:
            record['state'] = 'startup_verified'
            runtime.write_json(directory / 'startup.json', record)
            print(json.dumps(record, indent=2), flush=True)
            return 0
        if state == 'COMPLETED':
            record['state'] = 'failed'
            runtime.write_json(directory / 'startup.json', record)
            raise ValueError(f'Simulation job {job} exited before startup could be verified')
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            record['state'] = 'timed_out'
            runtime.write_json(directory / 'startup.json', record)
            raise ValueError('Startup check timed out; the simulation was not cancelled. Resume with simulation-startup.')
        runtime.write_json(directory / 'startup.json', record)
        time.sleep(min(15, remaining))
