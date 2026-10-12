"""Small, versioned IC parameter catalogue. Particle data never enter Git."""
import configparser
import csv
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

TABLE = 'IC_MODELS.csv'
FIELDS = ['model', 'status', 'run', 'gas_mass_Msun', 'stellar_mass_Msun',
          'gas_particles', 'stellar_particles', 'gas_scaleRadius_kpc', 'gas_scaleHeight_kpc',
          'stellar_scaleRadius_kpc', 'stellar_sersicIndex', 'stellar_axisRatioZ',
          'stellar_beta', 'gas_kappa', 'stellar_kappa', 'other_components',
          'seed', 'generated_utc', 'parameter_file', 'parameter_sha256',
          'particle_sha256', 'record_id']


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def git(repo, *args):
    return subprocess.check_output(['git', '-C', str(repo), *args], text=True).strip()


def parameters(path, model):
    parser = configparser.RawConfigParser()
    parser.read(path)
    sections = {name.lower(): dict(parser.items(name)) for name in parser.sections()}
    row = dict.fromkeys(FIELDS, '')
    row.update(model=model, parameter_file=Path(path).name, parameter_sha256=sha(path))
    known = set()
    for group in ('gas', 'stellar'):
        densities = [(name, values) for name, values in sections.items()
                     if name == 'density_' + group or name.startswith('density_' + group + '_')]
        components = [(name, values) for name, values in sections.items()
                      if name == 'component_' + group or name.startswith('component_' + group + '_')]
        known.update(name for name, _ in densities + components)
        masses = [float(values['mass']) for _, values in densities]
        if not all(math.isfinite(mass) and mass > 0 for mass in masses):
            raise ValueError(f'Invalid density mass in {path}')
        row[group + '_mass_Msun'] = format(sum(masses) * 1e10, '.12g')
        row[group + '_particles'] = str(sum(int(values.get('nbody', 0)) for _, values in components))
        def values_for(items, key):
            values = [(name.split('_', 1)[1], values[key]) for name, values in items if key in values]
            return values[0][1] if len(values) == 1 else '; '.join(name + '=' + value for name, value in values)
        row[group + '_scaleRadius_kpc'] = values_for(densities, 'scaleradius')
        if group == 'gas':
            row['gas_scaleHeight_kpc'] = values_for(densities, 'scaleheight')
        else:
            row['stellar_sersicIndex'] = values_for(densities, 'sersicindex')
            row['stellar_axisRatioZ'] = values_for(densities, 'axisratioz')
            row['stellar_beta'] = values_for(components, 'beta')
        row[group + '_kappa'] = values_for(components, 'ickappa')
    # Do not silently classify an unfamiliar historical component as stars or DM.
    others = {name: values for name, values in sections.items()
              if name.startswith(('density_', 'component_')) and name not in known}
    row['other_components'] = json.dumps(others, sort_keys=True) if others else ''
    return row


def historical_rows(root):
    root = Path(root)
    rows = []
    for ini in sorted((root / 'IniFile').glob('IC_*.ini')):
        model = ini.stem[3:]
        directory = root / model
        frozen = directory / ini.name
        particle = directory / 'Fornax.ini'
        current = parameters(ini, model)
        if frozen.is_file() and particle.is_file():
            row = parameters(frozen, model)
            with particle.open() as stream:
                header = stream.readline().split()
            try:
                counts = [int(n) for n in header]
                valid = len(counts) == 4 and min(counts) >= 0 and sum(counts) > 0
                valid = valid and counts[0] == int(row['gas_particles']) and counts[2] == int(row['stellar_particles'])
            except ValueError:
                valid = False
            row.update(status='existing_ic' if valid else 'existing_ic_check_needed', run=model,
                       record_id='legacy:' + model + ':' + row['parameter_sha256'])
            rows.append(row)
            if current['parameter_sha256'] == row['parameter_sha256']:
                continue
        current.update(status='parameters_only', record_id='parameters:' + model)
        rows.append(current)
    return rows


def generated_row(directory):
    directory = Path(directory)
    build = json.loads((directory / 'build.json').read_text())
    if build['state'] != 'complete':
        raise ValueError('Only successful IC builds can enter the catalogue')
    params = directory / 'parameters.ini'
    particle = directory / build['particle_ic']
    if params.parent != directory or particle.resolve().parent != directory.resolve():
        raise ValueError('Build receipt points outside the IC run')
    if sha(params) != build['ini_sha256'] or sha(particle) != build['particle_sha256']:
        raise ValueError('IC files changed after generation; catalogue registration refused')
    model = Path(build['source_ini']).stem.removeprefix('IC_')
    row = parameters(params, model)
    row.update(status='generated', run=directory.name, parameter_file=Path(build['source_ini']).name,
               seed=str(build['seed']), generated_utc=build.get('generated_utc', ''),
               particle_sha256=build['particle_sha256'],
               record_id='build:' + build.get('build_id', directory.name + ':' + build['particle_sha256']))
    return row


def update(repo, additions, push=True):
    repo = Path(repo).resolve()
    lock = Path(git(repo, 'rev-parse', '--git-path', 'ic-catalog.lock'))
    if not lock.is_absolute():
        lock = repo / lock
    with lock.open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        if push:
            if git(repo, 'branch', '--show-current') != 'main':
                raise ValueError('IC catalogue synchronization requires the main branch')
            if git(repo, 'status', '--porcelain'):
                raise ValueError('Catalogue checkout has uncommitted edits; commit them before retrying synchronization')
            git(repo, 'fetch', 'origin', 'main')
            ancestor = subprocess.run(['git', '-C', str(repo), 'merge-base', '--is-ancestor', 'origin/main', 'HEAD']).returncode
            if ancestor:
                git(repo, 'merge', '--ff-only', 'origin/main')
        table = repo / TABLE
        rows = {}
        if table.exists():
            with table.open(newline='') as source:
                reader = csv.DictReader(source)
                if reader.fieldnames != FIELDS:
                    raise ValueError('Unexpected IC catalogue schema')
                for row in reader:
                    if row['record_id'] in rows:
                        raise ValueError('Duplicate catalogue record ID')
                    rows[row['record_id']] = row
        for row in additions:
            rows[row['record_id']] = row
        with tempfile.NamedTemporaryFile(mode='w', newline='', dir=repo, delete=False) as target:
            writer = csv.DictWriter(target, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(sorted(rows.values(), key=lambda row: (row['model'], row['status'], row['run'], row['record_id'])))
            temporary = Path(target.name)
        os.replace(temporary, table)
        if push:
            git(repo, 'add', '--', TABLE)
            changed = subprocess.run(['git', '-C', str(repo), 'diff', '--cached', '--quiet', '--', TABLE])
            if changed.returncode not in (0, 1):
                changed.check_returncode()
            if changed.returncode == 1:
                git(repo, 'commit', '--only', '-m', 'Update IC model catalogue', '--', TABLE)
            git(repo, 'push', 'origin', 'HEAD:refs/heads/main')
        return len(rows)


def register(directory, repo):
    directory = Path(directory)
    receipt = directory / 'build.json'
    build = json.loads(receipt.read_text())
    try:
        row = generated_row(directory)
        (directory / 'catalog-entry.json').write_text(json.dumps(row, indent=2) + '\n')
        update(repo, [row])
        build['catalog_sync'] = {'state': 'pushed', 'commit': git(repo, 'rev-parse', 'HEAD')}
    except Exception as error:
        build['catalog_sync'] = {'state': 'pending', 'error': str(error)}
        receipt.write_text(json.dumps(build, indent=2) + '\n')
        raise
    receipt.write_text(json.dumps(build, indent=2) + '\n')


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--history-root', type=Path)
    source.add_argument('--run', type=Path)
    parser.add_argument('--no-push', action='store_true', help='stage a historical import locally')
    args = parser.parse_args()
    if args.run:
        if args.no_push:
            parser.error('--no-push only applies to historical imports')
        register(args.run.resolve(), args.repo)
        print('IC catalogue synchronized.')
    else:
        print(f'{update(args.repo, historical_rows(args.history_root), push=not args.no_push)} IC catalogue records')


if __name__ == '__main__':
    main()
