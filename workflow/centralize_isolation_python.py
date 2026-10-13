"""Archive duplicate legacy Python files, retaining entries required by user shells.

The inventory is JSONL with model and files[name].sha256/mtime_ns fields.
Only audited top-level .py files in numeric Fornax isolation models are replaced.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_shells(model, files):
    names = {p.name for p in model.iterdir() if p.suffix in ('.sh', '.slurm')}
    expected = {n for n in files if n.endswith(('.sh', '.slurm'))}
    if names != expected:
        raise ValueError(f'Shell file set changed: {model}')
    for name in names:
        p = model / name
        if digest(p) != files[name]['sha256'] or p.stat().st_mtime_ns != files[name]['mtime_ns']:
            raise ValueError(f'Shell changed: {p}')


def centralize(root, rows, source, archive, verify=False):
    root = root.resolve()
    iso = root / 'isoFornax'
    shared = iso / 'python'
    if iso.is_symlink() or shared.is_symlink():
        raise ValueError('Isolation/shared directory must not be a symlink')
    archive = archive.resolve()
    archive.relative_to(root / 'observer_review_delete_20261006')
    for row in rows:
        if not re.fullmatch(r'Fornax\d+', row['model']):
            raise ValueError('Only numeric legacy models are supported')
        model = iso / row['model']
        if model.is_symlink():
            raise ValueError('Model directory must not be a symlink')
        check_shells(model, row['files'])
        for name, info in row['files'].items():
            if not name.endswith('.py'):
                continue
            if Path(name).name != name:
                raise ValueError('Invalid filename')
            original = model / name
            target = shared / name
            backup = archive / row['model'] / name
            expected = info['sha256']
            if digest(source / name) != expected:
                raise ValueError(f'Different Python version: {original}')
            if target.exists() and (target.is_symlink() or digest(target) != expected):
                raise ValueError(f'Shared source differs: {target}')
            if original.is_symlink():
                if original.resolve() != target or digest(original) != expected or digest(backup) != expected:
                    raise ValueError(f'Invalid migrated entry: {original}')
                continue
            if verify:
                raise ValueError(f'Unmigrated entry: {original}')
            if original.exists() and digest(original) != expected:
                raise ValueError(f'Original source changed: {original}')
            if backup.exists() and digest(backup) != expected:
                raise ValueError(f'Archive collision: {backup}')
            if not original.exists() and not backup.exists():
                raise ValueError(f'Missing original and backup: {original}')
            shared.mkdir(exist_ok=True)
            if not target.exists():
                shutil.copy2(source / name, target)
            backup.parent.mkdir(parents=True, exist_ok=True)
            # Rename preserves the original metadata; interrupted runs resume here.
            if original.exists():
                if backup.exists():
                    raise ValueError(f'Original and archive both exist: {original}')
                original.rename(backup)
            temporary = model / ('.' + name + '.centralize')
            if temporary.is_symlink():
                temporary.unlink()
            os.symlink(os.path.relpath(target, model), temporary)
            os.replace(temporary, original)
        check_shells(model, row['files'])
        print(json.dumps({'model': row['model'], 'verified' if verify else 'centralized': True}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--inventory', required=True, type=Path)
    parser.add_argument('--source', type=Path, default=Path(__file__).parent / 'isolation_python')
    parser.add_argument('--archive', required=True, type=Path)
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--stop', type=int)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    rows = [json.loads(line) for line in args.inventory.read_text().splitlines()]
    if not args.verify:
        queue = subprocess.run(['squeue', '-M', 'all', '--me', '--noheader', '--format=%i|%T|%Z'],
                               check=True, capture_output=True, text=True)
        # Some Slurm versions print CLUSTER: headings even with --noheader.
        jobs = [line for line in queue.stdout.splitlines() if line.strip() and not line.startswith('CLUSTER:')]
        if jobs:
            raise SystemExit('Refusing migration while user jobs are queued or running')
    centralize(args.root, rows[args.start:args.stop], args.source, args.archive, args.verify)


if __name__ == '__main__':
    main()
