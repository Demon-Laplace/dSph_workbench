#!/usr/bin/env python3
"""Verify complete small GIZMO runs, including actual time advancement.

Example: python verify_simulation.py CASE/output --counts 512,0,1024,0,0,0
         --time-max 0.00001 --report validation.json
This checks numerical/file integrity, not physical equilibrium or convergence.
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np


def verify(directory, counts, time_max):
    paths = sorted(Path(directory).glob('snapshot_*.hdf5'))
    if len(paths) < 2:
        raise ValueError('At least two snapshots are required to verify time advancement')
    rows = []
    first_ids = None
    first_mass = None
    for path in paths:
        with h5py.File(path, 'r') as snapshot:
            header = snapshot['Header'].attrs
            actual = header['NumPart_Total'].tolist()
            if actual != counts or int(header['NumFilesPerSnapshot']) != 1:
                raise ValueError(f'{path}: particle counts or snapshot parts differ from expected')
            ids = []
            masses = []
            for kind, count in enumerate(counts):
                if not count:
                    continue
                group = snapshot[f'PartType{kind}']
                for name in ('Coordinates', 'Velocities', 'ParticleIDs'):
                    if len(group[name]) != count:
                        raise ValueError(f'{path}: invalid {name} length')
                for name, dataset in group.items():
                    values = dataset[:]
                    if values.dtype.kind == 'f' and not np.isfinite(values).all():
                        raise ValueError(f'{path}: non-finite {name}')
                if 'Masses' in group:
                    mass = group['Masses'][:]
                    if len(mass) != count or not (mass > 0).all():
                        raise ValueError(f'{path}: invalid masses')
                    # Only scalar sums use double precision; bulk arrays retain stored dtype.
                    masses.append(float(mass.sum(dtype=np.float64)))
                else:
                    constant_mass = float(header['MassTable'][kind])
                    if not np.isfinite(constant_mass) or constant_mass <= 0:
                        raise ValueError(f'{path}: invalid constant mass')
                    masses.append(constant_mass * count)
                ids.extend(group['ParticleIDs'][:].tolist())
                if kind == 0:
                    if not (group['Density'][:] > 0).all() or not (group['InternalEnergy'][:] >= 0).all():
                        raise ValueError(f'{path}: invalid gas density or internal energy')
            if len(set(ids)) != sum(counts):
                raise ValueError(f'{path}: duplicate particle IDs')
            if first_ids is None:
                first_ids, first_mass = set(ids), masses
            elif set(ids) != first_ids or not np.allclose(masses, first_mass, rtol=2e-6, atol=0):
                raise ValueError(f'{path}: particle IDs or component masses changed')
            rows.append({'file': path.name, 'time': float(header['Time']),
                         'counts': actual, 'component_masses': masses})
    times = np.array([row['time'] for row in rows])
    if not np.isfinite(times).all() or times[0] != 0 or not (np.diff(times) > 0).all():
        raise ValueError('Snapshots must advance monotonically from time zero')
    if not np.isclose(times[-1], time_max, rtol=1e-10, atol=0):
        raise ValueError(f'Final snapshot time {times[-1]} differs from requested {time_max}')
    return {'passed': True, 'scope': 'fixed-particle-count smoke test; not a convergence test',
            'finite_fields': True, 'unique_ids': True, 'mass_conserved': True, 'snapshots': rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--counts', required=True, help='six comma-separated particle counts')
    parser.add_argument('--time-max', required=True, type=float)
    parser.add_argument('--report', required=True, type=Path)
    args = parser.parse_args()
    counts = [int(value) for value in args.counts.split(',')]
    if len(counts) != 6 or min(counts) < 0 or sum(counts) == 0 or not args.time_max > 0:
        parser.error('Expected six nonnegative counts and a positive final time')
    report = verify(args.directory, counts, args.time_max)
    args.report.write_text(json.dumps(report, indent=2) + '\n')
    print(f'Verified {len(report["snapshots"])} snapshots: {args.directory}')


if __name__ == '__main__':
    main()
