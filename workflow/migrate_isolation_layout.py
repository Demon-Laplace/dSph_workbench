"""Move an explicitly listed legacy cohort, preserving files and original IC workspace paths.

Run with a JSON plan containing root, source_root, models and receipt. Server
migration is refused: server isolation must stay separate. This does not infer
that a same-number production uses the IC.
"""
import argparse
import fcntl
import json
import os
from pathlib import Path
import re


def confined(root, path):
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError('Path escapes configured root: ' + str(path))
    return path


def links_in(directory):
    with os.scandir(directory) as entries:
        for entry in entries:
            path = Path(entry.path)
            if entry.is_symlink():
                yield path
            elif entry.is_dir(follow_symlinks=False):
                yield from links_in(path)


def migrate(plan, server=False):
    root, source_root = Path(plan['root']).resolve(), Path(plan['source_root']).resolve()
    if server:
        raise ValueError('Server isolation directories must remain separate; migration is local-only')
    receipt = confined(root, Path(plan['receipt']))
    receipt.parent.mkdir(parents=True, exist_ok=True)
    with receipt.with_suffix('.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        for model in plan['models']:
            if not re.fullmatch(r'Fornax\d+', model):
                raise ValueError('Invalid model: ' + model)
            source = source_root / model
            target = confined(root, root / model / 'isolation')
            if source.is_symlink():
                if source.resolve() != target.resolve() or not target.is_dir():
                    raise ValueError('Unexpected existing source link: ' + str(source))
                print(model + ': already migrated', flush=True)
                continue
            if not source.is_dir() or target.exists() or target.is_symlink():
                raise ValueError('Source absent or destination already occupied: ' + model)
            if source.stat().st_dev != root.stat().st_dev:
                raise ValueError('Migration requires same filesystem (atomic rename)')
            identity = (source.stat().st_dev, source.stat().st_ino)
            # stat the top-level inputs, directories and products without reading snapshots.
            before = {p.name: (p.lstat().st_dev, p.lstat().st_ino, p.lstat().st_size)
                      for p in source.iterdir() if not p.is_symlink()}
            links = []
            for path in links_in(source):
                raw = os.readlink(path)
                if not os.path.isabs(raw):
                    absolute = Path(os.path.abspath(path.parent / raw))
                    if absolute.is_relative_to(source):
                        absolute = target / absolute.relative_to(source)
                    new_path = target / path.relative_to(source)
                    links.append((str(path.relative_to(source)), raw, os.path.relpath(absolute, new_path.parent)))
            target.parent.mkdir(parents=True, exist_ok=True)
            record = {'model': model, 'source': str(source), 'target': str(target),
                      'directory_identity': identity, 'relative_links': links, 'state': 'planned'}
            with receipt.open('a') as log:
                log.write(json.dumps(record) + '\n'); log.flush(); os.fsync(log.fileno())
            source.rename(target)
            try:
                source.symlink_to(os.path.relpath(target, source.parent), target_is_directory=True)
                for relative, old, new in links:
                    path = target / relative
                    if old != new:
                        path.unlink(); path.symlink_to(new)
                after = {p.name: (p.lstat().st_dev, p.lstat().st_ino, p.lstat().st_size)
                         for p in target.iterdir() if not p.is_symlink()}
                if before != after or identity != (target.stat().st_dev, target.stat().st_ino):
                    raise ValueError('Directory identity changed unexpectedly')
                if source.resolve() != target:
                    raise ValueError('Compatibility alias failed')
            except BaseException:
                for relative, old, new in links:
                    path = target / relative
                    if path.is_symlink():
                        path.unlink(); path.symlink_to(old)
                if source.is_symlink(): source.unlink()
                target.rename(source)
                raise
            record.update(state='verified', preserved_top_level_entries=len(before))
            with receipt.open('a') as log:
                log.write(json.dumps(record) + '\n'); log.flush(); os.fsync(log.fileno())
            print(model + ': moved and verified', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('plan', type=Path)
    parser.add_argument('--server', action='store_true')
    args = parser.parse_args()
    migrate(json.loads(args.plan.read_text()), args.server)
