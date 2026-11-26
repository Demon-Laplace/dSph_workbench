#!/usr/bin/env python3
"""Build particle ICs in a new directory using the existing Schwarzschild solver."""
import argparse
import configparser
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


def generate_default_plots(data, model, build, directory=None):
    path = Path(__file__).with_name('plot_ic.py')
    spec = importlib.util.spec_from_file_location('ic_default_plots', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.generate(data, model, directory or module.default_directory(build), build=build)


def publish_flat(output, destination):
    path = Path(__file__).with_name('flat_layout.py')
    spec = importlib.util.spec_from_file_location('ic_flat_layout', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.publish(output, destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ini', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--model', help='model ID used for the exported particle filename')
    parser.add_argument('--threads', type=int,
                        help='optional OpenMP thread limit; default: inherit the AGAMA/OpenMP environment')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--plots-dir', type=Path, help='default: model/isolation for a standard ic/ build')
    parser.add_argument('--flat-output', type=Path,
                        help='publish directly in the model isolation folder; workspace retains only tool links')
    parser.add_argument('--no-plots', action='store_true', help='explicitly skip default IC figures')
    parser.add_argument('--no-catalog', action='store_true',
                        help='local build only: skip catalogue registration and GitHub synchronization')
    parser.add_argument('--catalog-repo', type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    source, output = args.ini.resolve(), args.output.absolute()
    model = args.model or source.stem.removeprefix('IC_')
    standard_model = re.fullmatch(r'(?:Fornax|MW|Sculptor)\d+[A-Za-z0-9_.-]*', model)
    local_model = args.no_catalog and re.fullmatch(r'X\d+[A-Za-z0-9_.-]*', model)
    if not (standard_model or local_model):
        parser.error('Use a full model ID; experimental X IDs require --no-catalog')
    if not source.is_file() or output.exists() or output.is_symlink():
        parser.error('INI must exist and output must be a new directory; overwriting is forbidden')
    if (args.threads is not None and args.threads < 1) or not 0 <= args.seed < 2**31:
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
    two_beta = any(config.get(section, 'beta', fallback='').strip().lower() == 'two_zone'
                   for section in config.sections() if section.lower().startswith('component_'))
    solver = Path(__file__).with_name('schwarzschild_two_beta.py' if two_beta else 'schwarzschild.py')
    output.mkdir(parents=True)
    params = output / 'parameters.ini'
    shutil.copy2(source, params)
    manifest = {'state': 'building', 'build_id': str(uuid.uuid4()),
                'seed': args.seed, 'threads': args.threads,
                'thread_policy': 'inherited' if args.threads is None else 'explicit',
                'python': sys.executable, 'source_ini': str(source),
                'solver_file': str(solver), 'parameter_snapshot': 'parameters.ini',
                'ini_sha256': None, 'solver_sha256': None,
                'version_tracking': 'declared_build_id_v1',
                'data_verification_performed': False,
                'sampling_recipe': {
                    section: {key: config.get(section, key, fallback=None)
                              for key in ('icbeta', 'ickappa', 'beta', 'numorbits', 'nbody')}
                     for section in config.sections() if section.lower().startswith('component_')}}
    if args.flat_output:
        flat = args.flat_output.absolute()
        model_directory = flat.parent if flat.name == 'isolation' else flat
        manifest.update(local_model_directory=str(model_directory),
                        local_output_directory=str(args.flat_output.absolute()))
    try:
        manifest['generator_revision'] = subprocess.check_output(
            ['git', '-C', str(Path(__file__).parent), 'rev-parse', 'HEAD'],
            text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        manifest['generator_revision'] = None
    if two_beta:
        base_solver = solver.with_name('schwarzschild.py')
        manifest['base_solver_file'] = str(base_solver)
        manifest['base_solver_sha256'] = None
    receipt = output / 'build.json'
    receipt.write_text(json.dumps(manifest, indent=2) + '\n')
    env = dict(os.environ)
    if args.threads is not None:
        env.update(OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS='1')
    worker = ('import sys,runpy,agama,numpy; '
              'numpy.random.seed(int(sys.argv[3])); agama.setRandomSeed(int(sys.argv[3])); '
              'sys.argv=sys.argv[1:3]; ')
    if two_beta:
        worker += ('reader=runpy.run_path(sys.argv[0]); '
                   'reader["build_components"](reader["createModel"](sys.argv[1]))')
    else:
        worker += 'runpy.run_path(sys.argv[0],run_name="__main__")'
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
                        particle_sha256=None)
    except Exception as error:
        manifest.update(state='failed', error=str(error))
        receipt.write_text(json.dumps(manifest, indent=2) + '\n')
        if args.flat_output:
            publish_flat(output, args.flat_output)
        raise
    receipt.write_text(json.dumps(manifest, indent=2) + '\n')
    if args.no_plots:
        manifest['plots'] = {'state': 'disabled'}
    else:
        try:
            manifest['plots'] = generate_default_plots(data, model, output, args.plots_dir)
            print('Initial model figures: '+manifest['plots']['directory'])
        except Exception as error:
            manifest['plots'] = {'state': 'failed', 'error': str(error)}
            print(f'IC generation succeeded, but default figures failed: {error}', file=sys.stderr)
    receipt.write_text(json.dumps(manifest, indent=2) + '\n')
    if args.flat_output:
        publication = publish_flat(output, args.flat_output)
        directory, prefix = Path(publication['directory']), publication['prefix']
        print(f'Built {sum(counts)} particles: {directory / (prefix + destination.name)}')
        print(f'IC parameter snapshot: {directory / (prefix + params.name)}')
    else:
        print(f'Built {sum(counts)} particles: {destination}')
    if args.no_catalog:
        print('Local IC build; catalogue registration and GitHub synchronization disabled.')
        return 0
    try:
        catalog.register(output, args.catalog_repo)
    except Exception as error:
        print(f'IC generation succeeded; catalogue sync pending: {error}\n'
              f'Retry with workflow.py ic-catalog --run {output.name}', file=sys.stderr)
        return 3
    finally:
        if args.flat_output:
            publish_flat(output, args.flat_output)
    print('IC catalogue updated and pushed to GitHub main.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
