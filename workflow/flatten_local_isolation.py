"""Apply an audited local isolation flattening plan without copying scientific data."""
import argparse
from contextlib import ExitStack
import fcntl
import json
import os
from pathlib import Path


def identity(path):
    stat = path.stat()
    return [stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns]


def apply(plan, journal):
    root = Path(plan['root']).resolve()
    source_root = Path(plan['source_root']).resolve()
    config = json.loads((root / 'workflow.json').read_text())
    if not config.get('ic', {}).get('root') or (root / config['ic']['root']).resolve() != source_root:
        raise ValueError('This plan requires the configured local IC workspace; server migration is forbidden')

    def inside(value):
        path = Path(value)
        if not path.resolve().is_relative_to(root):
            raise ValueError('Path escapes local workspace: ' + value)
        return path

    with ExitStack() as stack:
        for value in plan['locks']:
            path = inside(value)
            if path.exists():
                lock = stack.enter_context(path.open('a'))
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        # Check every collision and identity before moving the first file.
        for entry in plan['moves']:
            src, dst = inside(entry['source']), inside(entry['target'])
            if src.is_symlink() or identity(src) != entry['identity']:
                raise ValueError('Source changed: ' + str(src))
            if dst.exists() or dst.is_symlink():
                raise ValueError('Target occupied: ' + str(dst))
            if src.suffix == '.sh' or (src.is_dir() and any(src.rglob('*.sh'))):
                raise ValueError('Shell files are protected')
        for entry in plan['links']:
            path = Path(entry['path'])
            if not any(path.parent.resolve().is_relative_to(base) for base in (root, source_root)):
                raise ValueError('Alias outside allowed workspaces')
            if not path.is_symlink() or os.readlink(path) != entry['old']:
                raise ValueError('Alias changed: ' + str(path))
        for entry in plan['remove_links']:
            path = Path(entry['path'])
            if not path.parent.resolve().is_relative_to(root) or not path.is_symlink() or os.readlink(path) != entry['old']:
                raise ValueError('Source link changed: ' + str(path))
        log = stack.enter_context(Path(journal).open('x'))

        def record(value):
            log.write(json.dumps(value) + '\n')
            log.flush()
            os.fsync(log.fileno())

        for entry in plan['moves']:
            src, dst = Path(entry['source']), Path(entry['target'])
            dst.parent.mkdir(parents=True, exist_ok=True)
            src.rename(dst)
            if identity(dst) != entry['identity']:
                raise RuntimeError('Rename did not preserve file identity')
            record({'move': entry})
        for entry in plan['links']:
            path = Path(entry['path'])
            pending = path.with_name('.' + path.name + '.flatten')
            pending.symlink_to(entry['new'])
            pending.replace(path)
            assert path.exists(), path
            record({'alias': entry})
        for entry in plan['metadata_updates']:
            path = inside(entry['path'])
            assert path.read_text() == entry['original']
            path.write_text(entry['content'])
            record({'metadata': entry['path']})
        for entry in plan['remove_links']:
            Path(entry['path']).unlink()
            record({'removed_link': entry})
        removed = []
        for value in plan['remove_dirs']:
            path = inside(value)
            if path.exists() and not any(path.iterdir()):
                path.rmdir()
                removed.append(value)
                record({'removed_empty_directory': value})
        metadata = {entry['path'] for entry in plan['metadata_updates']}
        for entry in plan['moves']:
            if entry['target'] not in metadata and identity(Path(entry['target'])) != entry['identity']:
                raise RuntimeError('Final data identity verification failed')
        return {'models': plan['models'], 'moved_entries': len(plan['moves']),
                'aliases_updated': len(plan['links']), 'source_records': len(metadata),
                'empty_directories_removed': len(removed), 'shell_files_changed': 0,
                'verification': 'device, inode, size and mtime preserved for data; source metadata paths updated'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('plan', type=Path)
    parser.add_argument('--journal', required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(apply(json.loads(args.plan.read_text()), args.journal), indent=2))
