"""Recognize a model from its directory or a local isolation/output subdirectory."""
from pathlib import Path
import re


def model_id(directory=None):
    directory = Path.cwd() if directory is None else Path(directory)
    directory = directory.resolve()
    for candidate in (directory, *directory.parents):
        if re.fullmatch(r'(?:Fornax|MW|Sculptor)\d+[A-Za-z0-9_.-]*', candidate.name):
            return candidate.name
    return directory.name
