#!/usr/bin/env python3
"""Apply or undo a reviewed list of reversible moves inside one workspace."""
import argparse
import json
import os
from pathlib import Path

from workflow_runtime import confined, write_json, now


def execute(root, plan_path, restore=False):
    root = root.resolve()
    plan = json.loads(plan_path.read_text())
    records = plan['moves']
    pending = []
    for entry in records:
        source = confined(root, root / entry['destination' if restore else 'source'])
        target = confined(root, root / entry['source' if restore else 'destination'])
        if source.is_symlink():
            raise ValueError(f'Refusing to move a symlink: {source}')
        if not source.exists():
            if target.exists():
                continue  # Resume an interrupted, already-applied move.
            raise ValueError(f'Missing source: {source}')
        if target.exists() or target.is_symlink():
            raise ValueError(f'Destination exists: {target}')
        if source.stat().st_dev != root.stat().st_dev:
            raise ValueError('Archive moves must be atomic on the same filesystem')
        pending.append((source, target))
    journal = plan_path.with_suffix('.journal.jsonl')
    with journal.open('a') as log:
        for source, target in pending:
            target.parent.mkdir(parents=True, exist_ok=True)
            source.rename(target)
            log.write(json.dumps({'time': now(), 'source': str(source.relative_to(root)),
                                  'destination': str(target.relative_to(root)), 'restore': restore}) + '\n')
            log.flush()
    print(f'{len(pending)} moves completed; journal: {journal}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--restore', action='store_true')
    args = parser.parse_args()
    execute(args.root, confined(args.root, args.plan), args.restore)


if __name__ == '__main__':
    main()
