#!/usr/bin/env python3
"""Build particle ICs in a new directory using the existing Schwarzschild solver."""
import argparse
import configparser
import hashlib
import importlib.util
import json
import os
import re
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime, timezone
import uuid

import catalog


def component_type(name):
    if 'gas' in name:
        return 0
    if 'stellar' in name or 'star' in name:
        return 2
    if 'dark' in name or 'dm' in name:
        return 1
    raise ValueError(f'Unknown particle type for component {name!r}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ini', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--model', help='model ID used for the exported particle filename')
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--catalog-repo', type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    source, output = args.ini.resolve(), args.output.absolute()
    model = args.model or source.stem.removeprefix('IC_')
    if not re.fullmatch(r'(?:Fornax|MW|Sculptor)\d+[A-Za-z0-9_.-]*', model):
        parser.error('Use --model with a full model ID')
    if not source.is_file() or output.exists() or output.is_symlink():
        parser.error('INI must exist and output must be a new directory; overwriting is forbidden')
    if args.threads < 1 or not 0 <= args.seed < 2**31:
        parser.error('threads must be positive and seed must be in [0, 2**31)')
    if importlib.util.find_spec('agama') is None:
        parser.error('This Python has no AGAMA; select an AGAMA environment with --python in workflow.py')
    import numpy as np
    config = configparser.RawConfigParser()
    config.read(source)
    components = []
    for section in config.sections():
        if section.lower().startswith('component_') and config.has_option(section, 'nbody'):
            kind = component_type(section.lower()[10:])
            count = config.getint(section, 'nbody')
            if count <= 0:
                parser.error('Component particle counts must be positive')
            components.append((kind, section.lower(), count))
    if not components:
        parser.error('No particle components configured')
    solver = Path(__file__).with_name('schwarzschild.py')
    output.mkdir(parents=True)
    params = output / 'parameters.ini'
    shutil.copy2(source, params)
    manifest = {'state': 'building', 'build_id': str(uuid.uuid4()),
                'seed': args.seed, 'threads': args.threads,
                'python': sys.executable, 'source_ini': str(source),
                'ini_sha256': hashlib.sha256(params.read_bytes()).hexdigest(),
                'solver_sha256': hashlib.sha256(solver.read_bytes()).hexdigest()}
    receipt = output / 'build.json'
    receipt.write_text(json.dumps(manifest, indent=2) + '\n')
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS='1')
    worker = ('import sys,runpy,agama,numpy; '
              'numpy.random.seed(int(sys.argv[3])); agama.setRandomSeed(int(sys.argv[3])); '
              'sys.argv=sys.argv[1:3]; runpy.run_path(sys.argv[0],run_name="__main__")')
    try:
        with (output / 'build.log').open('w') as log:
            subprocess.run([sys.executable, '-c', worker, str(solver), str(params), str(args.seed)],
                           cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        counts = [0, 0, 0, 0]
        masses = {}
        arrays = []
        for kind, section, expected in sorted(components):
            data = np.loadtxt(output / f'model_{section}.nbody', ndmin=2)
            if data.shape != (expected, 7) or not np.isfinite(data).all() or not (data[:, 6] > 0).all():
                raise ValueError(f'Invalid particles exported for {section}')
            counts[kind] += len(data)
            masses[section] = float(data[:, 6].sum())
            # Match the existing funcini exporter: DM temp=1, gas/stars temp=0.
            arrays.append(np.column_stack((data, np.full(len(data), int(kind == 1)), np.full(len(data), kind))))
        data = np.vstack(arrays)
        destination = output / (model + '.ini')
        with destination.open('w') as stream:
            stream.write(' '.join(map(str, counts)) + '\n')
            np.savetxt(stream, np.column_stack((np.arange(len(data)), data)),
                       fmt=['%d'] + ['%.17g'] * 7 + ['%d', '%d'])
        # Keep the historical converter input name without duplicating particle data.
        (output / 'Fornax.ini').symlink_to(destination.name)
        manifest.update(state='complete', generated_utc=datetime.now(timezone.utc).isoformat(),
                        counts=counts, component_masses=masses,
                        particle_ic=destination.name,
                        particle_sha256=hashlib.sha256(destination.read_bytes()).hexdigest())
    except Exception as error:
        manifest.update(state='failed', error=str(error))
        receipt.write_text(json.dumps(manifest, indent=2) + '\n')
        raise
    receipt.write_text(json.dumps(manifest, indent=2) + '\n')
    print(f'Built {sum(counts)} particles: {destination}')
    try:
        catalog.register(output, args.catalog_repo)
    except Exception as error:
        print(f'IC generation succeeded; catalogue sync pending: {error}\n'
              f'Retry with workflow.py ic-catalog --run {output.name}', file=sys.stderr)
        return 3
    print('IC catalogue updated and pushed to GitHub main.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
