"""Keep historical notebook imports and data paths explicit after code relocation."""
import os
from pathlib import Path
import sys


def configure(workspace, context='.'):
    workspace = Path(workspace).resolve()
    directory = (workspace / context).resolve()
    if not directory.is_relative_to(workspace) or not directory.is_dir():
        raise ValueError('Notebook context must be an existing directory inside the workspace')
    source = Path(__file__).resolve().parent / 'analysis'
    for name in ('basefunc', 'variable', 'readsnap', 'zim', 'mysubs'):
        loaded = sys.modules.get(name)
        if loaded is not None and Path(loaded.__file__).resolve().parent != source:
            raise RuntimeError('Restart this notebook kernel to load the reorganized legacy helpers')
    os.chdir(directory)
    os.environ['DSPH_LEGACY_WORKSPACE'] = str(directory)
    if str(source) not in sys.path:
        sys.path.insert(0, str(source))
    return directory
