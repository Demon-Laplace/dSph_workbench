"""Prepare production cases and submit existing simulation launchers explicitly."""
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

import workflow_runtime as runtime


def build_ic(args):
    source = resolve_input(args.root, args.ini)
    destination = runtime.confined(args.root, args.root / 'workspace/models/initial' / runtime.label(args.name))
    if destination.exists() or destination.is_symlink():
        raise ValueError('Initial-condition run already exists; refusing to overwrite')
    builder = args.repo / 'local_tools/ic_generation/build_ic.py'
    if not builder.is_file():
        raise ValueError(f'Missing versioned IC builder: {builder}')
    python = args.python or runtime.config(args.root).get('ic_python') or sys.executable
    return subprocess.run([python, str(builder), '--ini', str(source), '--output', str(destination),
                           '--threads', str(args.threads), '--seed', str(args.seed)]).returncode


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
    template = (args.root / model_name(args.template)).resolve()
    dwarf = resolve_input(args.root, args.dwarf_ic)
    if any(c.isspace() for c in str(dwarf)):
        raise ValueError('Legacy IC merger requires paths without whitespace')
    destination = runtime.confined(args.root, args.root / 'workspace/models/production' / name)
    alias = args.root / name
    if destination.exists() or alias.exists() or alias.is_symlink():
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
    destination.mkdir(parents=True)
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
    alias.symlink_to(os.path.relpath(destination, args.root), target_is_directory=True)
    print(f'Created {destination}; legacy path {alias} is compatible.')
    return 0


def simulate(args):
    model = runtime.confined(args.root, args.root / model_name(args.model)).resolve()
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
        raise ValueError(f'IC merging must finish first: {initial}')
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
    print(result.stdout, end=''); print(result.stderr, end='')
    return result.returncode
