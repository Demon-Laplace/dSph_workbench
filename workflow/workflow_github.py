"""Publish a scoped production-source snapshot to a GitHub main branch."""
import ast
from pathlib import Path
import shutil
import subprocess

import workflow_runtime as runtime
from workflow_release import git


def source_files(source):
    paths = [p for p in source.glob('*.py') if not p.is_symlink()]
    paths += [source / n for n in ('AGENTS.md', 'SERVER_WORKFLOW.md') if (source / n).is_file()]
    if (source / 'tests').is_dir() and not (source / 'tests').is_symlink():
        paths += [p for p in (source / 'tests').rglob('*.py') if not p.is_symlink()
                  and not any(parent.is_symlink() for parent in p.parents if parent != source)]
    if not (source / 'GetInfo.py').is_file():
        raise ValueError('Expected the production analysis source directory')
    for path in paths:
        if path.stat().st_size > 2 * 1024**2:
            raise ValueError(f'Unexpected large source file: {path}')
        if path.suffix == '.py':
            ast.parse(path.read_bytes(), filename=str(path))
    return sorted(paths)


def sync(args):
    settings = runtime.config(args.root).get('github', {})
    repo = runtime.confined(args.root, args.root / settings.get('checkout', 'tools/dSph_workbench'))
    source = (args.root / settings.get('source', 'tools/server-analysis')).resolve()
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
    destination = repo / 'production_analysis'
    destination.mkdir(exist_ok=True)
    names = {str(p.relative_to(source)) for p in files}
    # Only remove previously tracked source snapshots that no longer exist at source.
    for name in git(repo, 'ls-files', 'production_analysis').splitlines():
        relative = str(Path(name).relative_to('production_analysis'))
        if relative not in names:
            (repo / name).unlink(missing_ok=True)
    for path in files:
        target = destination / path.relative_to(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
    git(repo, 'add', '-u', '--', 'production_analysis')
    git(repo, 'add', '-f', '--', *['production_analysis/' + name for name in sorted(names)])
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
    print(f'GitHub synchronized: {branch} {commit}')
    return 0
