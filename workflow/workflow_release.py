"""Release analysis code from Git; machine settings and simulation data stay outside Git."""
import json
from pathlib import Path, PurePosixPath
import shlex
import subprocess

import workflow_runtime as runtime
import workflow_management as management


def git(repo, *arguments):
    return subprocess.run(['git', '-C', str(repo), *arguments], check=True,
                          stdout=subprocess.PIPE, text=True).stdout.strip()


def build(args):
    directory = runtime.confined(args.root, args.root / 'workspace/releases' / runtime.label(args.name))
    if directory.exists():
        raise ValueError('Release name already exists')
    if git(args.repo, 'status', '--porcelain'):
        raise ValueError('Commit or resolve checkout changes before creating a release')
    commit = git(args.repo, 'rev-parse', 'HEAD')
    entries = git(args.repo, 'ls-tree', '-r', commit, '--', 'repo', 'workflow', 'IniFile').splitlines()
    blobs = []
    for entry in entries:
        metadata, name = entry.split('\t', 1)
        mode, kind, blob = metadata.split()
        if mode not in ('100644', '100755') or kind != 'blob':
            raise ValueError(f'Only regular tracked release files are supported: {name}')
        if '..' in PurePosixPath(name).parts or PurePosixPath(name).is_absolute():
            raise ValueError('Invalid release path')
        size = int(git(args.repo, 'cat-file', '-s', blob))
        if size > 2 * 1024**2:
            raise ValueError(f'Unexpected large source file: {name}')
        blobs.append((name, blob))
    if not any(name == 'repo/GetInfo.py' for name, _ in blobs):
        raise ValueError('Commit the production analysis baseline before releasing')
    directory.mkdir(parents=True)
    hashes = {}
    for name, blob in blobs:
        target = directory / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(subprocess.run(['git', '-C', str(args.repo), 'cat-file', 'blob', blob],
                                         check=True, stdout=subprocess.PIPE).stdout)
        hashes[name] = runtime.sha(target)
    runtime.write_json(directory / 'release.json', {'created_utc': runtime.now(), 'git_commit': commit,
                                                   'name': args.name, 'files': hashes})
    print(f'Release {args.name}: {commit}, {len(hashes)} files')
    return 0


def verify(directory):
    manifest = json.loads((directory / 'release.json').read_text())
    for name, expected in manifest['files'].items():
        path = runtime.confined(directory, directory / name)
        if path.is_symlink() or not path.is_file() or runtime.sha(path) != expected:
            raise ValueError(f'Release checksum mismatch: {name}')
    actual = {str(p.relative_to(directory)) for p in directory.rglob('*') if p.is_file()}
    if actual != set(manifest['files']) | {'release.json'}:
        raise ValueError('Unexpected files in immutable release')
    if 'repo/GetInfo.py' not in manifest['files']:
        raise ValueError('Release lacks production analysis')
    return manifest


def activate(args):
    directory = runtime.confined(args.root, args.root / 'workspace/releases' / runtime.label(args.name))
    manifest = verify(directory)
    config_path = runtime.confined(args.root, args.root / 'workflow.json')
    settings = runtime.config(args.root)
    old_code = settings.get('analysis_code')
    settings['analysis_code'] = str((directory / 'repo').relative_to(args.root))
    settings['analysis_release'] = {'name': args.name, 'git_commit': manifest['git_commit']}
    journal = runtime.confined(args.root, args.root / 'workspace/catalog')
    journal.mkdir(parents=True, exist_ok=True)
    with (journal / 'activations.jsonl').open('a') as stream:
        stream.write(json.dumps({'time': runtime.now(), 'previous_code': old_code,
                                 'release': args.name, 'git_commit': manifest['git_commit']}) + '\n')
    runtime.write_json(config_path, settings)
    # Update the navigation entry atomically; existing prepared runs have pinned copies.
    link = runtime.confined(args.root, args.root / 'workspace/code/production')
    link.parent.mkdir(parents=True, exist_ok=True)
    if link.is_symlink() or not link.exists():
        temporary = link.with_name('production.next')
        if temporary.is_symlink():
            temporary.unlink()
        temporary.symlink_to(directory / 'repo', target_is_directory=True)
        temporary.replace(link)
    print(f'Activated analysis release {args.name} ({manifest["git_commit"]})')
    return 0


def deploy(args):
    runtime.label(args.name)
    source = runtime.confined(args.root, args.root / 'workspace/releases' / args.name)
    verify(source)
    settings, ssh = management.ssh_config(args)
    target = settings['root'] + '/workspace/releases/' + args.name
    quoted = shlex.quote(target)
    subprocess.run(ssh + [settings['host'], f'test ! -e {quoted} && mkdir -p {quoted}'], check=True)
    subprocess.run(['rsync', '-rt', '--checksum', '-e', shlex.join(ssh), str(source) + '/',
                    settings['host'] + ':' + target + '/'], check=True)
    return subprocess.run(management.remote_command(args, ['activate-release', '--name', args.name])).returncode
