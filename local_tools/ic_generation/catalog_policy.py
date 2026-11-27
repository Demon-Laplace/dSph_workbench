"""Catalogue publication defaults shared by IC build entry points."""
from pathlib import Path


def skip_catalog_for_paths(paths, explicit=None):
    """Skip only the X-series experiment tree or an explicit local-only build."""
    if explicit is not None:
        return explicit
    for path in paths:
        for candidate in (str(path), str(Path(path).resolve())):
            parts = candidate.replace('\\', '/').casefold().split('/')
            if any(parts[index:index + 2] == ['analysis', 'x_series_ic_experiments']
                   for index in range(len(parts) - 1)):
                return True
    return False
