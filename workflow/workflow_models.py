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
    monitored = prepare_startup_monitor(args, model, directory, script)
    command = ['sbatch', '--parsable', '--chdir', str(model), '--output', str(directory / 'slurm-%j.out')]
    if not args.submit:
        command.append('--test-only')
    command.append(str(monitored))
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


def prepare_startup_monitor(args, model, directory, launcher):
    """Freeze a batch wrapper; leave the original model launcher unchanged."""
    sources = ('workflow_models.py', 'workflow_runtime.py', 'workflow_merge.py')
    for name in sources:
        shutil.copy2(Path(__file__).with_name(name), directory / name)
    request = {'model': str(model), 'name': args.name,
               'parameter_sha256': runtime.sha(model / 'GZWJL.PARAM'),
               'launcher_sha256': runtime.sha(launcher),
               'launcher': str(launcher), 'timeout': getattr(args, 'startup_timeout', 7200),
               'monitor_sources': {name: runtime.sha(directory / name) for name in sources}}
    cluster = re.search(r'^#SBATCH\s+--clusters(?:=|\s+)(\S+)', launcher.read_text(), re.M)
    request['cluster'] = cluster.group(1) if cluster else 'all'
    command = [sys.executable, str(directory / 'workflow_models.py'), '--in-job',
               '--root', str(args.root), '--model', args.model, '--name', args.name]
    wrapper = write_monitored_launcher(directory, launcher, command)
    request['wrapper_sha256'] = runtime.sha(wrapper)
    runtime.write_json(directory / 'startup-request.json', request)
    return wrapper


def write_monitored_launcher(directory, launcher, monitor):
    text = launcher.read_text()
    interpreter = shlex.split(text.splitlines()[0][2:]) if text.startswith('#!') else ['/bin/sh']
    headers = [line for line in text.splitlines() if line.startswith('#SBATCH')]
    wrapper = directory / 'run-monitored.sh'
    wrapper.write_text('#!/bin/bash\n' + '\n'.join(headers) + '\n' +
        'env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 ' +
        shlex.join(monitor) + ' > ' + shlex.quote(str(directory / 'startup-monitor.log')) + ' 2>&1 &\n' +
        'monitor_pid=$!\n' +
        'finish() {\n  result=$?\n  trap - EXIT\n' +
        '  printf "%s\\n" "$result" > ' + shlex.quote(str(directory / 'solver-exit-status')) + '\n' +
        '  kill "$monitor_pid" 2>/dev/null || true\n' +
        '  wait "$monitor_pid" 2>/dev/null || true\n  exit "$result"\n}\n' +
        'trap finish EXIT\ntrap "exit 143" TERM\ntrap "exit 130" INT\n' +
        shlex.join([*interpreter, str(launcher)]) + '\nexit $?\n')
    return wrapper


def integration_progress(output):
    progress = None
    info = output / 'info.txt'
    if info.is_file():
        for point in re.finditer(r'Sync-Point\s+(\d+),\s*Time:\s*([0-9.eE+-]+)', info.read_text(errors='replace')):
            if int(point.group(1)) > 0 and float(point.group(2)) > 0:
                progress = {'sync_point': int(point.group(1)), 'simulation_time': float(point.group(2))}
    return progress


def finish_startup_report(args, model, directory, record, parameters):
    details = simulation_details(args.root, model, directory, record, parameters)
    report = directory / 'simulation-details.json'
    runtime.write_json(report, details)
    record['details_sha256'] = runtime.sha(report)
    runtime.write_json(directory / 'startup.json', record)
    print_startup_report(args, record)
    return 0


def print_startup_report(args, record):
    print(json.dumps(record, indent=2), flush=True)
    print('DSPH_SIMULATION_DETAILS ' + json.dumps({'model': args.model, 'name': args.name,
        'sha256': record['details_sha256']}), flush=True)


def in_job_startup(args):
    """A sleeping batch child: no extra Slurm allocation, MPI task or reserved CPU."""
    model = runtime.model_directory(args.root, model_name(args.model)).resolve()
    directory = runtime.confined(args.root, args.root / 'workspace/jobs' / runtime.label(args.name))
    request = json.loads((directory / 'startup-request.json').read_text())
    assert Path(request['model']).resolve() == model
    for name, expected in request['monitor_sources'].items():
        if runtime.sha(directory / name) != expected:
            raise ValueError('Frozen startup monitor changed: ' + name)
    if runtime.sha(Path(request['launcher'])) != request['launcher_sha256'] or \
            runtime.sha(directory / 'run-monitored.sh') != request['wrapper_sha256']:
        raise ValueError('Frozen startup launcher changed')
    parameter = model / 'GZWJL.PARAM'
    if runtime.sha(parameter) != request['parameter_sha256']:
        raise ValueError('Parameters changed after monitor preparation')
    parameters = {}
    for line in parameter.read_text().splitlines():
        fields = line.split('%', 1)[0].split()
        if len(fields) >= 2: parameters[fields[0]] = fields[1]
    output = runtime.confined(args.root, model / parameters['OutputDir'])
    started = time.monotonic()
    while True:
        elapsed = time.monotonic() - started
        record = {'job_id': os.environ['SLURM_JOB_ID'], 'cluster': os.environ.get('SLURM_CLUSTER_NAME', request['cluster']),
                  'model': str(model), 'checked_utc': runtime.now(), 'state': 'waiting',
                  'scheduler': {'state': 'RUNNING', 'elapsed_seconds': int(elapsed), 'exit_code': '0:0'},
                  'progress': integration_progress(output), 'required_running_seconds': 140,
                  'monitor_location': 'inside_batch_allocation'}
        if elapsed >= 140 and record['progress']:
            record['state'] = 'startup_verified'
            return finish_startup_report(args, model, directory, record, parameters)
        if elapsed >= request['timeout']:
            record['state'] = 'timed_out'
            runtime.write_json(directory / 'startup.json', record)
            return 2
        runtime.write_json(directory / 'startup.json', record)
        time.sleep(min(15, 140 - elapsed) if elapsed < 140 else 15)


def await_startup_report(args, directory, job, cluster, timeout):
    """Wait on a small batch-generated receipt; query Slurm at most once a minute."""
    from workflow_isolation_ops import accounting
    deadline = time.monotonic() + timeout
    next_accounting = 0
    while True:
        receipt = directory / 'startup.json'
        if receipt.is_file():
            record = json.loads(receipt.read_text())
            if record.get('job_id') != job:
                raise ValueError('Startup receipt belongs to another job')
            if record['state'] == 'startup_verified':
                if runtime.sha(directory / 'simulation-details.json') != record['details_sha256']:
                    raise ValueError('Startup report checksum mismatch')
                print_startup_report(args, record)
                return 0
            if record['state'] in {'failed', 'timed_out'}:
                raise ValueError('Batch startup monitor failed: ' + record['state'])
        if (directory / 'solver-exit-status').is_file():
            raise ValueError('Solver exited before startup report was available')
        now = time.monotonic()
        if now >= next_accounting:
            try:
                state = accounting(job, cluster)['state']
            except ValueError as error:
                if str(error) != 'Slurm accounting is unavailable or ambiguous; retry later': raise
                state = 'ACCOUNTING_PENDING'
            if state not in {'PENDING', 'RUNNING', 'CONFIGURING', 'COMPLETING', 'SUSPENDED',
                             'REQUEUED', 'REQUEUE_FED', 'REQUEUE_HOLD', 'RESIZING', 'ACCOUNTING_PENDING'}:
                raise ValueError('Simulation ended before startup report was available: ' + state)
            next_accounting = now + 60
        if now >= deadline:
            raise ValueError('Startup check timed out; simulation and batch monitor were not cancelled. Resume with simulation-startup.')
        time.sleep(min(15, deadline - now))


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
    if (directory / 'startup-request.json').is_file():
        return await_startup_report(args, directory, job, cluster, timeout)
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
        progress = integration_progress(output)
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
            return finish_startup_report(args, model, directory, record, parameters)
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


def simulation_details(root, model, directory, startup, parameters):
    """Summarize actual per-component particle masses from the merged binary."""
    import configparser
    import struct
    import numpy as np
    from workflow_merge import read_par
    galaxy_files = [p for p in (model / 'mergegalaxy.par', model / 'mergegalaxye.par') if p.is_file()]
    if len(galaxy_files) != 1:
        raise ValueError('Require one galaxy merge parameter file for the simulation report')
    galaxy = read_par(galaxy_files[0])
    cube = read_par(model / 'mergecube.par')
    components = []
    for role, filename in [('background_gas', cube[0][0]), ('milky_way', galaxy[0][0]), ('dwarf', galaxy[0][1])]:
        path = (model / filename).resolve()
        with path.open() as stream:
            counts = [int(n) for n in stream.readline().split()]
        if len(counts) != 4 or min(counts) < 0:
            raise ValueError('Invalid component particle header: ' + str(path))
        components.append({'role': role, 'model': path.stem, 'input_file': str(path),
                           'particle_counts': counts, 'total_particles': sum(counts), 'resolution': {}})
    initial = (model / parameters['InitCondFile']).resolve()
    total = sum(c['total_particles'] for c in components)
    with initial.open('rb') as stream:
        if stream.read(4) != struct.pack('<I', 256):
            raise ValueError('Invalid binary IC header for simulation report')
        header = stream.read(256)
        expected = tuple(sum(c['particle_counts'][tp] for c in components) for tp in range(4)) + (0, 0)
        if struct.unpack_from('<6I', header) != expected or any(struct.unpack_from('<6d', header, 24)):
            raise ValueError('Binary IC counts or mass table disagree with component inputs')
        if stream.read(4) != struct.pack('<I', 256):
            raise ValueError('Invalid binary IC header trailer')
        for size in [12 * total, 12 * total, 4 * total]:
            if stream.read(4) != struct.pack('<I', size):
                raise ValueError('Invalid binary IC particle block')
            stream.seek(size, 1)
            if stream.read(4) != struct.pack('<I', size):
                raise ValueError('Invalid binary IC particle block trailer')
        if stream.read(4) != struct.pack('<I', 4 * total):
            raise ValueError('Invalid binary IC mass block')
        offset = stream.tell()
    masses = np.memmap(initial, dtype='<f4', mode='r', offset=offset, shape=(total,))
    mass_unit = float(parameters['UnitMass_in_g']) / 1.989e33
    cursor = 0
    for tp in range(4):
        # The validated merger orders cube, MW, dwarf within each particle type.
        for component in components:
            count = component['particle_counts'][tp]
            if count:
                values = masses[cursor:cursor + count]
                if not np.isfinite(values).all() or np.any(values <= 0):
                    raise ValueError('Invalid IC masses in simulation report')
                component['resolution'][str(tp)] = {'count': count,
                    'particle_mass_min_msun': float(values.min()) * mass_unit,
                    'particle_mass_max_msun': float(values.max()) * mass_unit,
                    'particle_mass_mean_msun': float(values.mean(dtype=np.float64)) * mass_unit,
                    'total_mass_msun': float(values.sum(dtype=np.float64)) * mass_unit}
            cursor += count
    case = json.loads((model / 'case.json').read_text()) if (model / 'case.json').is_file() else {}
    source = {}
    if case.get('export_receipt'):
        source = json.loads(Path(case['export_receipt']).read_text())
    initial_parameters = {}
    if source.get('snapshot'):
        path = Path(source['snapshot']).parent.parent / 'parameters.ini'
        if path.is_file():
            config = configparser.ConfigParser(interpolation=None)
            config.optionxform = str
            config.read(path)
            initial_parameters = {section: dict(config[section]) for section in config.sections()}
    resources = {}
    if (directory / 'startup-request.json').is_file():
        launcher = Path(json.loads((directory / 'startup-request.json').read_text())['launcher'])
    else:
        submitted = json.loads((directory / 'submission.json').read_text())
        launcher = Path(submitted['command'][-1]) if submitted.get('command') else None
    if launcher and launcher.is_file():
        for option in re.finditer(r'^#SBATCH\s+--(nodes|ntasks|cpus-per-task|time)(?:=|\s+)(\S+)', launcher.read_text(), re.M):
            resources[option.group(1)] = option.group(2)
    return {'created_utc': runtime.now(), 'model': model.name, 'startup': dict(startup),
        'requested_slurm_resources': resources, 'total_particles': total,
        'particle_type_order': ['gas', 'halo', 'disk_stars', 'bulge_stars'],
        'mass_unit_msun': mass_unit, 'components': components,
        'isolation_source': source, 'dwarf_initial_parameters': initial_parameters,
        'orbit_reference': case.get('orbit_reference'), 'orientation_reference': case.get('orientation_reference'),
        'orbit': {'position_kpc': [float(x) for x in galaxy[1][:3]],
                  'velocity_km_s': [float(x) for x in galaxy[1][3:6]],
                  'dwarf_rotation_deg': [float(x) for x in galaxy[1][9:12]],
                  'box_shift_kpc': [float(x) for x in cube[1][:3]]},
        'parameters': parameters, 'binary_ic': str(initial), 'binary_ic_sha256': runtime.sha(initial),
        'parameter_sha256': runtime.sha(model / 'GZWJL.PARAM')}


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Internal production startup monitor')
    parser.add_argument('--in-job', action='store_true', required=True)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--model', required=True)
    parser.add_argument('--name', required=True)
    arguments = parser.parse_args()
    arguments.root = arguments.root.resolve()
    try:
        result = in_job_startup(arguments)
    except Exception as error:
        directory = runtime.confined(arguments.root, arguments.root / 'workspace/jobs' / runtime.label(arguments.name))
        runtime.write_json(directory / 'startup.json', {'job_id': os.environ.get('SLURM_JOB_ID'),
            'state': 'failed', 'checked_utc': runtime.now(), 'error': str(error)})
        print('Startup monitor failed: ' + str(error), file=sys.stderr)
        result = 2
    raise SystemExit(result)
