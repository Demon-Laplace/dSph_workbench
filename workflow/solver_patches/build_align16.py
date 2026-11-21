#!/usr/bin/env python3
"""Build an independent legacy GIZMO with a 16-byte allocator, never in place."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess


def aligned_source(text):
    old = '''  if((n % 8) > 0)
    n = (n / 8 + 1) * 8;

  if(n < 8)
    n = 8;'''
    new = '''  /* Preserve SIMD alignment through allocations, resizing and relocation. */
  if((n % 16) > 0)
    n = (n / 16 + 1) * 16;

  if(n < 16)
    n = 16;'''
    if text.count(old) != 4:
        raise ValueError('Expected exactly four original allocator rounding blocks')
    base = 'if(!(Base = malloc(n)))'
    if text.count(base) != 1:
        raise ValueError('Expected the original memory-pool allocation')
    return text.replace(old, new).replace(base, 'if(posix_memalign(&Base, 16, n) != 0)')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--destination', required=True, type=Path)
    parser.add_argument('--jobs', type=int, default=4)
    args = parser.parse_args()
    source, destination = args.source.resolve(), args.destination.absolute()
    if destination.exists() or source == destination or source in destination.parents:
        raise ValueError('Choose a new independent destination, outside the source tree')
    original = (source / 'system/mymalloc.c').read_text()
    patched = aligned_source(original)
    def ignore(directory, names):
        return [name for name in names if name.endswith('.o') or name == '.git'
                or (name.startswith('GIZMO') and name != 'GIZMO_config.h')]
    shutil.copytree(source, destination, ignore=ignore)
    (destination / 'system/mymalloc.c').write_text(patched)
    # Preserve the exact pre-existing compile-time settings and optimization flags.
    configuration = (source / 'GIZMO_config.h').read_text()
    command = ['make', '-j' + str(args.jobs), 'EXEC=GIZMO_align16']
    record = {'state': 'building', 'source': str(source), 'directory': str(destination),
              'created_utc': datetime.now(timezone.utc).isoformat(),
              'command': command, 'alignment_bytes': 16,
              'changes': ['pool base via posix_memalign(16)',
                          'all four allocation/reallocation lengths rounded to 16'],
              'physics_configuration': configuration}
    manifest = destination / 'aligned-build.json'
    manifest.write_text(json.dumps(record, indent=2) + '\n')
    try:
        with (destination / 'build.log').open('w') as log:
            subprocess.run(command, cwd=destination, stdout=log, stderr=subprocess.STDOUT, check=True)
        if (destination / 'GIZMO_config.h').read_text() != configuration:
            raise ValueError('Compile-time configuration changed; do not activate this build')
        record.update(state='built', executable=str(destination / 'GIZMO_align16'))
    except Exception as error:
        record.update(state='failed', error=str(error))
        raise
    finally:
        manifest.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2), flush=True)


if __name__ == '__main__':
    main()
