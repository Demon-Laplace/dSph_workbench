"""Publish a scoped production-source snapshot to a GitHub main branch."""
import ast
import json
from pathlib import Path
import shutil
import subprocess
import shlex
import hashlib

import workflow_runtime as runtime
from workflow_release import git
import workflow_management as management


def source_bytes(path):
    return path.read_bytes().replace(b'\r\n', b'\n')


def source_sha(path):
    return hashlib.sha256(source_bytes(path)).hexdigest()


def source_baseline(root, repo=None):
    path = root / 'workspace/catalog/github-source-baseline.json'
    if path.exists():
        return json.loads(path.read_text())['files']
    name = runtime.config(root).get('analysis_release', {}).get('name')
    if name:
        path = root / 'workspace/releases' / runtime.label(name) / 'workflow/source_manifest.json'
    elif repo is not None:
        path = repo / 'workflow/source_manifest.json'
    return json.loads(path.read_text()).get('files', {}) if path.exists() else {}


def source_files(source):
    paths = [p for pattern in ('*.py', '*.sh', '*.slurm') for p in source.glob(pattern) if not p.is_symlink()]
    paths += [source / n for n in ('AGENTS.md', 'SERVER_WORKFLOW.md')
              if (source / n).is_file() and not (source / n).is_symlink()]
    for folder in ('tests', 'build_production_ic', 'check_density_profile', 'evolution_pipeline', 'paper_figures'):
        directory = source / folder
        if directory.is_dir() and not directory.is_symlink():
            paths += [p for p in directory.rglob('*') if p.is_file() and not p.is_symlink()
                      and p.suffix in {'.py', '.md', '.json', '.sh', '.slurm'}
                      and (p.suffix != '.json' or p.name.startswith('config'))
                      and not any(part.startswith('.') or part == '__pycache__' for part in p.relative_to(source).parts)
                      and not any(parent.is_symlink() for parent in p.parents if parent != source)]
    # Canonical orbit tools use the same package policy as analysis prepare;
    # exporting a release must not include data or omit frozen dependencies.
    paths += runtime.analysis_package_files(source)
    if not (source / 'GetInfo.py').is_file():
        raise ValueError('Expected the production analysis source directory')
    for path in paths:
        if path.stat().st_size > 2 * 1024**2:
            raise ValueError(f'Unexpected large source file: {path}')
        if path.suffix == '.py':
            ast.parse(path.read_bytes(), filename=str(path))
        elif path.suffix in {'.sh', '.slurm'}:
            subprocess.run(['bash', '-n', str(path)], check=True)
    return sorted(paths)


def sync(args):
    settings = runtime.config(args.root).get('github', {})
    repo = runtime.confined(args.root, args.root / settings.get('checkout', 'tools/dSph_workbench'))
    source = Path(getattr(args, 'source', None) or (args.root / settings.get('source', 'tools/server-analysis'))).resolve()
    branch = settings.get('branch', 'main')
    url = settings.get('url', 'https://github.com/Demon-Laplace/dSph_workbench.git')
    files = source_files(source)
    if not repo.exists():
        repo.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(['git', 'clone', url, str(repo)], check=True)
    if git(repo, 'remote', 'get-url', 'origin') != url:
        raise ValueError('Checkout origin differs from configured GitHub URL')
    if git(repo, 'status', '--porcelain'):
        raise ValueError('Sync checkout has uncommitted changes; resolve them without discarding edits')
    baseline = getattr(args, 'source_baseline', None)
    if baseline is None:
        baseline = source_baseline(args.root, repo)
    subprocess.run(['git', '-C', str(repo), 'fetch', 'origin'], check=True)
    exists = subprocess.run(['git', '-C', str(repo), 'show-ref', '--verify', '--quiet', 'refs/heads/' + branch]).returncode == 0
    remote_exists = subprocess.run(['git', '-C', str(repo), 'show-ref', '--verify', '--quiet', 'refs/remotes/origin/' + branch]).returncode == 0
    if exists:
        git(repo, 'switch', branch)
    elif remote_exists:
        git(repo, 'switch', '--track', '-c', branch, 'origin/' + branch)
    else:
        git(repo, 'switch', '-c', branch)
    if remote_exists:
        # A failed prior push can leave the local branch ahead. Preserve those commits.
        remote_is_ancestor = subprocess.run(['git', '-C', str(repo), 'merge-base', '--is-ancestor', 'origin/' + branch, 'HEAD']).returncode == 0
        if not remote_is_ancestor:
            git(repo, 'merge', '--ff-only', 'origin/' + branch)
    destination = repo / 'repo'
    destination.mkdir(exist_ok=True)
    names = {str(p.relative_to(source)) for p in files}
    for path in files:
        name = str(path.relative_to(source))
        current = destination / name
        if name in baseline and current.is_file():
            latest = source_sha(current)
            if latest != baseline[name] and source_sha(path) != latest:
                raise ValueError(f'GitHub and source both changed: {name}; reconcile before synchronizing')
    manifest_path = repo / 'workflow/source_manifest.json'
    for relative, previous_hash in baseline.items():
        if relative not in names:
            target = runtime.confined(destination, destination / relative)
            if target.is_file() and source_sha(target) != previous_hash:
                raise ValueError(f'GitHub changed a source-deleted file: {relative}; reconcile before synchronizing')
            target.unlink(missing_ok=True)
    for path in files:
        target = destination / path.relative_to(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(source_bytes(path))
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    runtime.write_json(manifest_path, {'files': {str(p.relative_to(source)): source_sha(p) for p in files}})
    git(repo, 'add', '-u', '--', 'repo')
    git(repo, 'add', '-f', '--', 'workflow/source_manifest.json', *['repo/' + name for name in sorted(names)])
    changed = subprocess.run(['git', '-C', str(repo), 'diff', '--cached', '--quiet']).returncode
    if changed == 1:
        git(repo, 'commit', '-m', args.message)
    elif changed != 0:
        raise ValueError('Unable to inspect staged source changes')
    # Never force-push. Authentication or divergence failures leave the commit intact.
    subprocess.run(['git', '-C', str(repo), 'push', '--set-upstream', 'origin', 'HEAD:refs/heads/' + branch], check=True)
    commit = git(repo, 'rev-parse', 'HEAD')
    receipt = runtime.confined(args.root, args.root / 'workspace/catalog')
    receipt.mkdir(parents=True, exist_ok=True)
    runtime.write_json(receipt / 'github-sync.json', {'time': runtime.now(), 'branch': branch,
                       'commit': commit, 'source_files': len(files), 'pushed': True})
    if not getattr(args, 'source', None):
        runtime.write_json(receipt / 'github-source-baseline.json', {'commit': commit,
                           'files': {str(p.relative_to(source)): source_sha(p) for p in files}})
    print(f'GitHub synchronized: {branch} {commit}')
    return 0


def export_source(args):
    settings = runtime.config(args.root).get('github', {})
    source = (args.root / settings.get('source', 'tools/server-analysis')).resolve()
    directory = runtime.confined(args.root, args.root / 'workspace/transfers/source' / runtime.label(args.name))
    if directory.exists():
        raise ValueError('Source export name already exists')
    files = source_files(source)
    (directory / 'code').mkdir(parents=True)
    hashes = {}
    for path in files:
        name = str(path.relative_to(source))
        target = directory / 'code' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        hashes[name] = runtime.sha(target)
        if runtime.sha(path) != hashes[name]:
            raise ValueError(f'Source changed during export: {name}')
    runtime.write_json(directory / 'manifest.json', {'created_utc': runtime.now(), 'files': hashes,
                       'baseline_files': source_baseline(args.root, args.repo)})
    print(f'Exported {len(files)} source files: {args.name}')
    return 0


def record_source_sync(args):
    directory = runtime.confined(args.root, args.root / 'workspace/transfers/source' / runtime.label(args.name))
    manifest = json.loads((directory / 'manifest.json').read_text())
    files = {}
    for name, expected in manifest['files'].items():
        path = runtime.confined(directory / 'code', directory / 'code' / name)
        if path.is_symlink() or runtime.sha(path) != expected:
            raise ValueError('Source export changed before recording synchronization')
        files[name] = source_sha(path)
    catalog = runtime.confined(args.root, args.root / 'workspace/catalog')
    catalog.mkdir(parents=True, exist_ok=True)
    runtime.write_json(catalog / 'github-source-baseline.json', {'commit': args.commit, 'files': files})
    return 0


def sync_server(args):
    name = 'github-' + runtime.now().replace(':', '').replace('+', '-').replace('.', '-')
    subprocess.run(management.remote_command(args, ['export-source', '--name', name]), check=True)
    settings, ssh = management.ssh_config(args)
    directory = runtime.confined(args.root, args.root / 'workspace/transfers/source' / name)
    directory.mkdir(parents=True)
    remote = settings['root'] + '/workspace/transfers/source/' + name + '/'
    subprocess.run(['rsync', '-rt', '--checksum', '-e', shlex.join(ssh),
                    settings['host'] + ':' + remote, str(directory) + '/'], check=True)
    manifest = json.loads((directory / 'manifest.json').read_text())
    for filename, digest in manifest['files'].items():
        path = runtime.confined(directory / 'code', directory / 'code' / filename)
        if path.is_symlink() or runtime.sha(path) != digest:
            raise ValueError('Source transfer checksum mismatch: ' + filename)
    actual = {str(p.relative_to(directory / 'code')) for p in (directory / 'code').rglob('*') if p.is_file()}
    if actual != set(manifest['files']):
        raise ValueError('Unexpected files in source export')
    args.source = directory / 'code'
    args.source_baseline = manifest['baseline_files']
    result = sync(args)
    receipt = json.loads((args.root / 'workspace/catalog/github-sync.json').read_text())
    subprocess.run(management.remote_command(args, ['record-source-sync', '--name', name,
                   '--commit', receipt['commit']]), check=True)
    return result
