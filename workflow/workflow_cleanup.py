#!/usr/bin/env python3
"""Inventory or delete inactive Slurm stdout/stderr under one explicit data root."""
import argparse
import fnmatch
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shlex
import stat
import subprocess
import time

LOG = re.compile(r'slurm[-_](\d+)(?:_\d+)?\.(?:out|err|log)$')


def declared_logs(base, files):
    patterns = []
    for name in files:
        if not name.endswith(('.sh', '.slurm', '.sbatch')):
            continue
        path = Path(base) / name
        if path.is_symlink() or path.stat().st_size > 1000000:
            continue
        for line in path.read_text(errors='replace').splitlines():
            if not line.startswith('#SBATCH '):
                continue
            tokens = shlex.split(line[len('#SBATCH '):])
            for index, token in enumerate(tokens):
                value = None
                if token in ('-o', '-e', '--output', '--error') and index + 1 < len(tokens):
                    value = tokens[index + 1]
                elif token.startswith(('--output=', '--error=')):
                    value = token.split('=', 1)[1]
                if value and Path(value).parent == Path('.'):
                    patterns.append(re.sub(r'%[A-Za-z]', '*', value))
    return patterns


def active_jobs():
    result = subprocess.run(['squeue', '-M', 'all', '--me', '--noheader', '--format=%i|%Z'],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True)
    jobs = []
    for line in result.stdout.splitlines():
        if not line.strip() or line.startswith('CLUSTER:'):
            continue
        job, directory = line.split('|', 1)
        jobs.append((job.strip().split('_')[0], Path(directory.strip()).resolve()))
    return jobs


def candidates(root):
    for base, dirs, files in os.walk(str(root)):
        dirs[:] = [d for d in dirs if d not in ('.git', '__pycache__', '.cache', '.mpl', '.tmp')
                   and not (Path(base) / d).is_symlink()]
        patterns = declared_logs(base, files)
        for name in files:
            match = LOG.fullmatch(name)
            if match or (name.endswith(('.out', '.err', '.log')) and any(fnmatch.fnmatchcase(name, p) for p in patterns)):
                path = Path(base) / name
                info = path.lstat()
                if stat.S_ISREG(info.st_mode):
                    yield path, match[1] if match else '', info


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--apply', action='store_true')
    parser.add_argument('--queue-snapshot', type=Path, help='scan only: queue snapshot from login node')
    parser.add_argument('--plan', type=Path, help='apply an existing inventory without rescanning')
    args = parser.parse_args()
    root = args.root.resolve()
    receipt = args.receipt.resolve()
    if root == Path('/') or not receipt.is_relative_to(root):
        parser.error('Receipt must be inside a non-root data directory')
    if args.apply and args.queue_snapshot:
        parser.error('Deletion requires live queue queries, not a saved snapshot')
    jobs = [(job, Path(path)) for job, path in json.loads(args.queue_snapshot.read_text())] if args.queue_snapshot else active_jobs()
    print('Active jobs:', len(jobs), flush=True)
    receipt.mkdir(parents=True, exist_ok=True)
    records = []
    total = 0
    for path, job, info in ([] if args.plan else candidates(root)):
        if any(job == active or path.is_relative_to(directory) for active, directory in jobs):
            continue
        if time.time() - info.st_mtime < 300:
            continue
        records.append({'path': str(path.relative_to(root)), 'job_id': job, 'bytes': info.st_size,
                        'mtime_ns': info.st_mtime_ns, 'inode': info.st_ino, 'device': info.st_dev})
        total += info.st_size
    if args.plan:
        stored = json.loads(args.plan.read_text())
        if Path(stored['root']).resolve() != root:
            parser.error('Inventory belongs to a different root')
        records = stored['files']
        total = sum(entry['bytes'] for entry in records)
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    plan = receipt / ('slurm-plan-' + stamp + '.json')
    plan.write_text(json.dumps({'root': str(root), 'files': records}, indent=2) + '\n')
    print(f'Candidates: {len(records)} files, {total / 2**30:.3f} GiB; {plan}', flush=True)
    if not args.apply:
        return
    jobs = active_jobs()
    deleted = skipped = freed = 0
    with (receipt / ('slurm-deleted-' + stamp + '.jsonl')).open('x') as report:
        for index, entry in enumerate(records):
            if index and index % 200 == 0:
                jobs = active_jobs()
            path = root / entry['path']
            if not path.resolve().is_relative_to(root) or path.is_symlink():
                skipped += 1; continue
            if any(path.is_relative_to(directory) or entry['job_id'] == job for job, directory in jobs):
                skipped += 1; continue
            try:
                current = path.lstat()
                # Lustre mount device numbers differ between login/compute nodes.
                # The inode, size and mtime identify the file across those mounts.
                if (current.st_size, current.st_mtime_ns, current.st_ino) != (entry['bytes'], entry['mtime_ns'], entry['inode']):
                    skipped += 1; continue
                path.unlink()
                report.write(json.dumps(entry) + '\n'); report.flush()
                deleted += 1; freed += entry['bytes']
            except FileNotFoundError:
                skipped += 1
    summary = {'deleted': deleted, 'skipped': skipped, 'logical_bytes_removed': freed, 'finished_utc': stamp}
    (receipt / ('slurm-summary-' + stamp + '.json')).write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
