"""One-pass production IC merge, legacy binary conversion and verified cleanup."""
import argparse
import fcntl
import hashlib
import itertools
import json
import os
from pathlib import Path
import shlex
import shutil
import struct
import subprocess
import sys

import workflow_runtime as runtime


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_par(path):
    lines = [shlex.split(line, comments=True) for line in path.read_text().splitlines()]
    lines = [line for line in lines if line]
    if len(lines) < 2 or len(lines[0]) != 2 or len(lines[1]) != 12:
        raise ValueError(f'Expected two filenames and twelve transform values: {path}')
    import math
    if not all(math.isfinite(float(x)) for x in lines[1]):
        raise ValueError(f'Invalid transform: {path}')
    if any(any(c.isspace() for c in x) for x in lines[0]):
        raise ValueError('Legacy particle filenames must not contain whitespace')
    return lines[:2]


def plan(root, model):
    directory = runtime.model_directory(root, model).resolve()
    galaxy = [directory / n for n in ('mergegalaxy.par', 'mergegalaxye.par') if (directory / n).is_file()]
    if len(galaxy) != 1:
        raise ValueError('Require exactly one mergegalaxy.par or mergegalaxye.par')
    cube = directory / 'mergecube.par'
    g, c = read_par(galaxy[0]), read_par(cube)
    # The second input must mean the first merge output, not unrelated particles.
    if Path(c[0][1]) not in (Path('mergegalaxy.dat'), Path('mergegalaxye.dat')):
        raise ValueError('mergecube.par must consume mergegalaxy.dat / mergegalaxye.dat')
    inputs = [(directory / x).resolve() for x in (g[0][0], g[0][1], c[0][0])]
    if not all(p.is_file() for p in inputs):
        raise ValueError('Missing particle inputs: ' + ', '.join(str(p) for p in inputs if not p.is_file()))
    parameters = {}
    param = directory / 'GZWJL.PARAM'
    for line in param.read_text().splitlines():
        fields = line.split('%', 1)[0].split()
        if len(fields) >= 2:
            parameters[fields[0]] = fields[1]
    target = runtime.confined(directory, directory / parameters['InitCondFile'])
    if target.parent != directory or not target.name.endswith('IC_000'):
        raise ValueError('InitCondFile must be an IC_000 file directly inside the model directory')
    if target.exists() or target.is_symlink():
        raise ValueError(f'Existing IC is protected: {target}')
    output = runtime.confined(directory, directory / parameters['OutputDir'])
    if list(output.glob('snapshot*')) or list(output.glob('restart*')):
        raise ValueError('Model already has simulation data; refusing to rebuild IC')
    files = list(dict.fromkeys([galaxy[0], cube, param, *inputs]))
    identity = {str(p): {'size': p.stat().st_size, 'mtime_ns': p.stat().st_mtime_ns,
                         'inode': p.stat().st_ino} for p in files}
    return directory, galaxy[0], cube, target, files, identity


def validate_binary(ascii_path, binary):
    """Validate every block and every particle against the final ASCII, in chunks."""
    import numpy as np
    with ascii_path.open() as f:
        counts = [int(x) for x in f.readline().split()]
    if len(counts) != 4 or min(counts) < 0:
        raise ValueError('Invalid merged header')
    n, gas = sum(counts), counts[0]
    offsets = []
    sizes = [256, 12*n, 12*n, 4*n, 4*n] + ([4*gas, 4*gas] if gas else [])
    with binary.open('rb') as f:
        for size in sizes:
            if f.read(4) != struct.pack('<I', size):
                raise ValueError('Missing/invalid binary block marker')
            offsets.append(f.tell())
            if size == 256 and len(offsets) == 1:
                header = f.read(size)
                if (struct.unpack_from('<6I', header) != tuple(counts+[0, 0]) or
                    struct.unpack_from('<6I', header, 96) != tuple(counts+[0, 0]) or
                    any(struct.unpack_from('<6d', header, 24))):
                    raise ValueError('Binary header counts/mass table mismatch')
            else:
                f.seek(size, 1)
            if f.read(4) != struct.pack('<I', size):
                raise ValueError('Truncated/mismatched binary block')
        if f.read(1):
            raise ValueError('Unexpected extra binary blocks')
    arrays = [np.memmap(binary, dtype=dtype, mode='r', offset=offset, shape=(size//4,))
              for offset, size, dtype in zip(offsets[1:], sizes[1:], ['<f4', '<f4', '<i4', '<f4', '<f4', '<f4'])]
    with ascii_path.open() as f:
        next(f)
        start = 0
        while True:
            lines = list(itertools.islice(f, 65536))
            if not lines:
                break
            data = np.loadtxt(lines, ndmin=2)
            end = start + len(data)
            if data.shape[1] != 10 or end > n or not np.isfinite(data).all():
                raise ValueError('Invalid merged particle rows')
            expected = [data[:, 1:4].astype('<f4').ravel(), data[:, 4:7].astype('<f4').ravel(),
                        data[:, 0].astype('<i4'), data[:, 7].astype('<f4')]
            for index, value in enumerate(expected):
                scale = 3 if index < 2 else 1
                if not np.array_equal(arrays[index][start*scale:end*scale], value):
                    raise ValueError(f'Binary particle block {index} differs at rows {start}:{end}')
            stop = min(end, gas)
            if start < stop:
                if (not np.array_equal(arrays[4][start:stop], data[:stop-start, 8].astype('<f4')) or
                        np.any(arrays[5][start:stop] != 0)):
                    raise ValueError('Gas energy/density blocks disagree')
            start = end
        if start != n:
            raise ValueError('Merged particle count mismatch')
    return {'particles': n, 'counts': counts, 'blocks': sizes, 'all_particle_fields_verified': True}


def settings(root):
    cfg = runtime.config(root)
    result = dict(cfg.get('ic_merge', {}))
    if 'converter' not in result:
        result['converter'] = cfg.get('isolation_profiles', {}).get('default', {}).get('converter')
    if not result.get('converter') or not isinstance(result['converter'], list):
        raise ValueError('Configure ic_merge.converter as [executable, "{ascii}", "{binary}"]')
    return result


def prepare(root, model, work, expected=None, frozen_settings=None):
    directory, galaxy, cube, target, files, identity = plan(root, model)
    if expected is not None and identity != expected:
        raise ValueError('Inputs changed after submission; submit a new preparation')
    cfg = frozen_settings if frozen_settings is not None else settings(root)
    hashes = {str(p): sha(p) for p in files}
    source = Path(__file__).with_name('ic_merge.c')
    executable = work / 'ic_merge'
    subprocess.run([cfg.get('compiler', 'gcc'), '-O2', '-std=c99', '-ffp-contract=off', '-fno-fast-math',
                    '-Wall', '-Wextra', str(source), '-lm', '-o', str(executable)], check=True)
    ascii_path, binary = work / 'merged.dat', work / target.name
    for path in (ascii_path, binary):
        if path.exists() or path.is_symlink():
            raise ValueError(f'Previous preparation output exists: {path}')
    with (work / 'merge.log').open('w') as log:
        subprocess.run([str(executable), str(galaxy), str(cube), str(ascii_path)], cwd=directory,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    command = [token.format(ascii=str(ascii_path), binary=str(binary)) for token in cfg['converter']]
    with (work / 'convert.log').open('w') as log:
        subprocess.run(command, cwd=directory, stdout=log, stderr=subprocess.STDOUT, check=True)
    validation = validate_binary(ascii_path, binary)
    if plan(root, model)[-1] != identity:
        raise ValueError('Inputs changed during merge; refusing to publish IC')
    receipt = {'created_utc': runtime.now(), 'model': model, 'inputs_sha256': hashes,
               'source_sha256': sha(source), 'converter_command': command,
               'converter_sha256': sha(Path(shutil.which(command[0]) or command[0])), 'binary': str(target), 'binary_sha256': sha(binary),
               'validation': validation, 'intermediate_bytes': ascii_path.stat().st_size}
    runtime.write_json(work / 'validated.json', receipt)
    # prepare runs under the model preparation lock. This filesystem forbids
    # hard links; rename publishes the complete file atomically on the same FS.
    # Recheck after validation; never replace an existing model IC.
    if target.exists() or target.is_symlink():
        raise ValueError(f'Existing IC is protected: {target}')
    binary.rename(target)
    ascii_path.unlink()  # Only this invocation's verified intermediate is removed.
    receipt['intermediate_deleted'] = True
    runtime.write_json(work / 'completed.json', receipt)
    print(json.dumps(receipt, indent=2))
    return receipt


def command(args):
    cfg = settings(args.root)
    work = runtime.confined(args.root, args.root / 'workspace/merge-ic' / args.model / runtime.label(args.name))
    work.parent.mkdir(parents=True, exist_ok=True)
    with (work.parent / '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if (work / 'submission.json').is_file():
            result = work / ('completed.json' if (work / 'completed.json').is_file() else 'submission.json')
            print(result.read_text())
            print('Existing submission: no duplicate job submitted.')
            return 0
        if work.exists():
            raise ValueError('Preparation name already exists; inspect it and choose a new --name')
        info = plan(args.root, args.model)
        from workflow_cleanup import active_jobs
        if any(d == info[0] or d.is_relative_to(work.parent) for _, d in active_jobs()):
            raise ValueError('An active Slurm job already uses this model/preparation')
        options = cfg.get('sbatch')
        if not options:
            raise ValueError('Configure ic_merge.sbatch for a compute-node preparation')
        work.mkdir()
        runtime.write_json(work / 'request.json', {'model': args.model, 'identity': info[-1], 'settings': cfg})
        for name in ('workflow_merge.py', 'ic_merge.c'):
            shutil.copy2(Path(__file__).with_name(name), work / name)
        script = work / 'prepare.sh'
        script.write_text('#!/bin/bash\nset -eu\n' + cfg.get('setup', '') + '\n' +
                          'export PYTHONPATH=' + shlex.quote(str(args.root)) + '\n' +
                          shlex.join([cfg.get('python', 'python3'), str(work / 'workflow_merge.py'),
                                      '--root', str(args.root), '--work', str(work), args.model]) + '\n')
        batch = ['sbatch', '--parsable', '--chdir', str(work), '--job-name', 'ic-'+args.model,
                 '--output', str(work / 'slurm-%j.out')]
        batch += [f'--{key}={value}' for key, value in options.items()]
        if not args.submit:
            batch.append('--test-only')
        batch.append(str(script))
        result = subprocess.run(batch, text=True, capture_output=True)
        record = {'command': batch, 'submitted': args.submit and result.returncode == 0,
                  'stdout': result.stdout, 'stderr': result.stderr, 'returncode': result.returncode}
        runtime.write_json(work / ('submission.json' if record['submitted'] else 'preflight.json'), record)
        print(json.dumps(record, indent=2))
        # Preflight is disposable; permit submitting the same requested name later.
        if not args.submit and result.returncode == 0:
            shutil.rmtree(work)
        return result.returncode


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('model')
    a = p.parse_args()
    a.root = a.root.resolve()
    a.work = runtime.confined(a.root, a.work).resolve()
    with (a.work.parent / '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        request = json.loads((a.work / 'request.json').read_text())
        prepare(a.root, a.model, a.work, request['identity'], request['settings'])
