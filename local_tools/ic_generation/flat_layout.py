"""Merge IC files into the flat isolation folder; tool views live in workspace."""
import json
import os
from pathlib import Path


def link_target(path):
    return Path(os.path.abspath(path.parent / os.readlink(path)))


def publish(directory, destination, prefix=None):
    """Rename files, not copy particles; retain a fixed-name compatibility view."""
    directory, destination = Path(directory).absolute(), Path(destination).absolute()
    if directory == destination or directory in destination.parents:
        raise ValueError('The model destination must be outside the working view')
    destination.mkdir(parents=True, exist_ok=True)
    receipt = directory / 'build.json'
    build = json.loads(receipt.read_text()) if receipt.is_file() else {}
    entries = list(directory.iterdir())
    if any(path.is_dir() and not path.is_symlink() for path in entries):
        raise ValueError('An IC working view must contain files, not subdirectories')
    if prefix is None:
        prefix = build.get('local_file_prefix')
    if prefix is None:
        prefix = ''
        if any((destination / path.name).exists() or (destination / path.name).is_symlink()
               for path in entries):
            prefix = directory.name + '__'
    if '/' in prefix or prefix in ('.', '..'):
        raise ValueError('Invalid IC filename prefix')
    pending = []
    for source in entries:
        target = destination / (prefix + source.name)
        if source.is_symlink() and link_target(source) == target:
            continue
        if target.exists() or target.is_symlink():
            raise FileExistsError('Existing model file is preserved: ' + str(target))
        pending.append((source, target, link_target(source) if source.is_symlink() else None))
    if receipt.is_file():
        model = destination.parent if destination.name == 'isolation' else destination
        build.update(local_model_directory=str(model),
                     local_output_directory=str(destination), local_file_prefix=prefix)
        receipt.write_text(json.dumps(build, indent=2) + '\n')
    targets = {source: target for source, target, _ in pending}
    for source, target, original_link in pending:
        if original_link is None:
            source.rename(target)
            source.symlink_to(os.path.relpath(target, source.parent))
    for source, target, original_link in pending:
        if original_link is not None:
            actual = targets.get(original_link, original_link)
            target.symlink_to(os.path.relpath(actual, target.parent))
            source.unlink()
            source.symlink_to(os.path.relpath(target, source.parent))
    return {'directory': str(destination), 'prefix': prefix,
            'moved_files': len(pending)}
