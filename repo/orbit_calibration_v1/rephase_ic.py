"""Translate/boost a complete dwarf in a NEW format-1 IC; never rotate it.

The target is the absolute physical COM of an explicitly frozen old-star
cohort. All dwarf particles receive the same displacement and velocity boost.
The old run's initial snapshot must first agree by ID with the source binary.
This utility builds inputs only; submission and startup use the existing workflow.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import struct
import time

import numpy as np


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024**2), b''):
            h.update(b)
    return h.hexdigest()


def binary_arrays(path):
    """Strict supported merger layout: HEADER, POS, VEL, ID, MASS, U, RHO."""
    path = Path(path)
    with path.open('rb') as f:
        if f.read(4) != struct.pack('<I', 256):
            raise ValueError('Require little-endian single-file format-1 IC')
        header = f.read(256)
    counts = struct.unpack_from('<6I', header)
    n, ng = sum(counts), counts[0]
    if (any(struct.unpack_from('<6d', header, 24)) or
            struct.unpack_from('<6I', header, 96) != counts or
            struct.unpack_from('<d', header, 72)[0] != 0):
        raise ValueError('Require explicit per-particle masses and noncosmological IC')
    sizes = [256, 12*n, 12*n, 4*n, 4*n] + ([4*ng, 4*ng] if ng else [])
    offsets = []
    with path.open('rb') as f:
        for size in sizes:
            if f.read(4) != struct.pack('<I', size):
                raise ValueError('Invalid binary block marker')
            offsets.append(f.tell())
            f.seek(size, 1)
            if f.read(4) != struct.pack('<I', size):
                raise ValueError('Invalid binary block trailer')
        if f.read(1):
            raise ValueError('Unsupported additional binary fields')
    arrays = {}
    names = ['position', 'velocity', 'id', 'mass', 'gas_u', 'gas_rho']
    for name, off, size in zip(names, offsets[1:], sizes[1:]):
        dtype = '<u4' if name == 'id' else '<f4'
        shape = (n, 3) if name in ('position', 'velocity') else (size // 4,)
        arrays[name] = np.memmap(path, mode='r', dtype=dtype, offset=off, shape=shape)
    return counts, sizes, offsets, arrays


def cohort_ids(path):
    with np.load(path, allow_pickle=False) as z:
        if 'particle_id' in z:
            ids = z['particle_id'].astype(np.uint64)
        else:
            first = int(z['first'].reshape(-1)[0])
            ids = np.r_[np.uint64(first), np.uint64(first) + np.cumsum(z['delta'], dtype=np.uint64)]
    if ids.ndim != 1 or not len(ids) or len(np.unique(ids)) != len(ids):
        raise ValueError('Cohort must contain unique uint64 IDs')
    return ids


def weighted_com(p, v, m, rows, L, V):
    weights = np.asarray(m[rows], dtype=np.float64)
    return np.r_[np.average(p[rows], axis=0, weights=weights) * L,
                 np.average(v[rows], axis=0, weights=weights) * V]


def join_ids(sorted_ids, sorted_rows, wanted):
    wanted = np.asarray(wanted, dtype=np.uint64)
    at = np.searchsorted(sorted_ids, wanted)
    if np.any(at >= len(sorted_ids)) or not np.array_equal(sorted_ids[np.minimum(at, len(sorted_ids)-1)], wanted):
        raise ValueError('Particle ID missing from initial binary')
    return sorted_rows[at]


def validate_snapshot(arrays, counts, snapshot, box_size, gas_velocity_atol=0., periodic_position_ulps=0.):
    """Full ID join; explicit audited initial-output allowances, exact mass/stars."""
    import h5py
    n = sum(counts)
    ids = np.asarray(arrays['id'])
    # Isolation star formation already creates generation bits. Full IDs remain
    # unique but need not be contiguous, nor bounded by the particle count.
    order = np.argsort(ids)
    sorted_ids = ids[order]
    if np.any(sorted_ids[1:] == sorted_ids[:-1]):
        raise ValueError('Duplicate full particle IDs in initial binary')
    result = {}
    with h5py.File(snapshot, 'r') as f:
        if float(f['Header'].attrs['Time']) != 0 or tuple(f['Header'].attrs['NumPart_Total']) != tuple(counts):
            raise ValueError('Initial snapshot time/counts mismatch')
        for tp, count in enumerate(counts):
            if not count:
                continue
            g = f['PartType%d' % tp]
            actual_ids = g['ParticleIDs'][:]
            if len(actual_ids) != count:
                raise ValueError('Invalid snapshot IDs')
            rows = join_ids(sorted_ids, order, actual_ids)
            lo, hi = sum(counts[:tp]), sum(counts[:tp+1])
            if len(np.unique(actual_ids)) != count or np.any(rows < lo) or np.any(rows >= hi):
                raise ValueError('Snapshot type/ID identity mismatch')
            result[str(tp)] = {}
            for name, key in [('position', 'Coordinates'), ('velocity', 'Velocities'), ('mass', 'Masses')]:
                worst = 0.0
                for start in range(0, count, 262144):
                    stop = min(count, start + 262144)
                    a = arrays[name][rows[start:stop]]
                    b = g[key][start:stop]
                    if not np.isfinite(a).all() or (name == 'mass' and np.any(a <= 0)):
                        raise ValueError('Invalid initial phase/mass values')
                    delta = b.astype(np.float64) - a
                    allowance = 0.
                    if name == 'position':
                        delta -= box_size * np.rint(delta / box_size)
                        allowance = periodic_position_ulps * np.abs(np.spacing(np.maximum(np.abs(a), np.abs(b)).astype(np.float32)))
                    elif name == 'velocity' and tp == 0:
                        allowance = gas_velocity_atol
                    # Initial gas output includes synchronized hydro kicks;
                    # periodic coordinates can acquire one float32 rounding.
                    # Such allowances must be audited separately and explicit.
                    if np.any(np.abs(delta) > allowance):
                        raise ValueError('IC disagrees with original snapshot: type%d/%s' % (tp, key))
                    worst = max(worst, float(np.max(np.abs(delta))))
                result[str(tp)][key + '_max_abs_error'] = worst
    return (sorted_ids, order), result


def build(config):
    start = time.monotonic()
    cfg = dict(config)
    source, dest = Path(cfg['source_binary']).resolve(), Path(cfg['destination_binary']).resolve()
    receipt = Path(cfg['receipt']).resolve()
    if source == dest or source.parent == dest.parent or dest.exists() or receipt.exists():
        raise ValueError('Require a fresh destination directory and fresh IC/report paths')
    if not dest.parent.is_dir() or not receipt.parent.is_dir():
        raise ValueError('Destination parents must already exist')
    if cfg['canonical_observer_id'] != 'astropy_native_direct_v1':
        raise ValueError('Unsupported observer')
    inputs = {name: Path(cfg[name]).resolve() for name in
              ['source_binary', 'source_snapshot', 'cohort_file', 'candidate_file', 'initial_state_file', 'dwarf_ascii']}
    for name, expected in cfg['expected_input_sha256'].items():
        if sha(inputs[name]) != expected:
            raise ValueError('Frozen input checksum mismatch: ' + name)
    candidate = json.loads(inputs['candidate_file'].read_text())
    if candidate['canonical_observer_id'] != cfg['canonical_observer_id']:
        raise ValueError('Candidate observer mismatch')
    target = np.asarray(candidate['initial_absolute_6D_seed'], dtype=np.float64)
    L, V = float(cfg['length_to_kpc']), float(cfg['velocity_to_kms'])
    if target.shape != (6,) or not np.isfinite(target).all() or not (L > 0 and V > 0):
        raise ValueError('Invalid target/units')
    counts, sizes, offsets, old = binary_arrays(source)
    if list(counts) != cfg['expected_counts']:
        raise ValueError('IC particle counts changed')
    n = sum(counts)
    ids = cohort_ids(inputs['cohort_file'])
    if len(ids) != cfg['expected_cohort_N']:
        raise ValueError('Cohort count changed')
    # Read the actual source IDs; never infer identity from row count or strip
    # generation bits. The historical merger adds one component ID offset.
    with inputs['dwarf_ascii'].open() as f:
        dwarf_counts = [int(x) for x in f.readline().split()]
    if dwarf_counts != cfg['expected_dwarf_counts']:
        raise ValueError('Dwarf source counts changed')
    dwarf_ids = np.loadtxt(inputs['dwarf_ascii'], skiprows=1, usecols=0, dtype=np.uint64, ndmin=1)
    dwarf_ids += np.uint64(cfg['dwarf_id_offset'])
    if len(dwarf_ids) != sum(dwarf_counts) or len(np.unique(dwarf_ids)) != len(dwarf_ids):
        raise ValueError('Invalid complete dwarf source IDs')
    if not np.isin(ids, dwarf_ids).all():
        raise ValueError('Frozen cohort is outside the source dwarf ID set')
    gas_atol = float(cfg.get('initial_snapshot_gas_velocity_atol_code', 0.))
    pos_ulps = float(cfg.get('initial_snapshot_periodic_position_ulps', 0.))
    if not (0 <= gas_atol <= 2e-5 and 0 <= pos_ulps <= 1):
        raise ValueError('Initial-output allowance exceeds independently audited supported bounds')
    mapping, original_check = validate_snapshot(old, counts, inputs['source_snapshot'],
                                                float(cfg['box_size_code']), gas_atol, pos_ulps)
    rows = join_ids(*mapping, ids)
    dwarf_rows = join_ids(*mapping, dwarf_ids)
    selected_dwarf = np.zeros(n, dtype=bool)
    selected_dwarf[dwarf_rows] = True
    actual_dwarf_counts = [int(selected_dwarf[sum(counts[:tp]):sum(counts[:tp+1])].sum()) for tp in range(4)]
    if actual_dwarf_counts != dwarf_counts:
        raise ValueError('Dwarf source IDs/type counts disagree with initial binary')
    original = weighted_com(old['position'], old['velocity'], old['mass'], rows, L, V)
    state = json.loads(inputs['initial_state_file'].read_text())
    expected = np.r_[state['old_P_absolute_kpc'], state['old_V_absolute_kms']]
    if state['old_N'] != len(ids) or not np.allclose(original, expected, rtol=0, atol=1e-8):
        raise ValueError('Original frozen-cohort COM differs from fitted initial state')
    delta = (target - original) / np.r_[[L]*3, [V]*3]
    provenance = {name: {'path': str(p), 'size': p.stat().st_size, 'sha256': sha(p)} for name, p in inputs.items()}
    # No hard links: the original file is never opened writable.
    with dest.open('xb') as f, source.open('rb') as original_stream:
        shutil.copyfileobj(original_stream, f, 8 * 1024**2)
    changed = {}
    for key, offset, shift in [('position', offsets[1], delta[:3]), ('velocity', offsets[2], delta[3:])]:
        writable = np.memmap(dest, mode='r+', dtype='<f4', offset=offset, shape=(n, 3))
        requested = old[key][dwarf_rows].astype(np.float64) + shift
        if key == 'position' and (np.any(requested < 0) or np.any(requested >= cfg['box_size_code'])):
            raise ValueError('Dwarf leaves initial periodic box; wrapping requires separate explicit support')
        writable[dwarf_rows] = requested
        writable.flush()
        changed[key] = {'max_float32_rounding_error_code': float(np.max(np.abs(writable[dwarf_rows] - requested)))}
        del writable
    _, _, _, new = binary_arrays(dest)
    for name in old:
        if name in ('position', 'velocity'):
            shift = delta[:3] if name == 'position' else delta[3:]
            for a in range(0, n, 262144):
                b = min(n, a + 262144)
                expected_block = np.array(old[name][a:b], copy=True)
                selected = selected_dwarf[a:b]
                expected_block[selected] = expected_block[selected].astype(np.float64) + shift
                if not np.array_equal(expected_block, new[name][a:b]):
                    raise ValueError('Unexpected position/velocity change')
        elif not np.array_equal(old[name], new[name]):
            raise ValueError('Protected field changed: ' + name)
    # Also verify exact block payload hashes for every unmodified block/header.
    with source.open('rb') as a, dest.open('rb') as b:
        for index, (offset, size) in enumerate(zip(offsets, sizes)):
            if index in (1, 2):
                continue
            a.seek(offset-4); b.seek(offset-4)
            if a.read(size+8) != b.read(size+8):
                raise ValueError('Protected binary block changed')
    final = weighted_com(new['position'], new['velocity'], new['mass'], rows, L, V)
    if not np.allclose(final, target, rtol=0, atol=1e-5):
        raise ValueError('New frozen-cohort COM misses target beyond float32 tolerance')
    if sha(source) != provenance['source_binary']['sha256']:
        raise ValueError('Original IC changed during preparation')
    result = {'state': 'VALIDATED_NEW_ORBIT_ONLY_IC', 'canonical_observer_id': cfg['canonical_observer_id'],
              'provenance': provenance, 'counts': list(counts), 'dwarf_N': len(dwarf_rows), 'cohort_N': len(ids),
              'dwarf_counts': actual_dwarf_counts, 'dwarf_selection': 'complete source ASCII full IDs plus historical merger offset',
              'full_ID_range': [int(old['id'].min()), int(old['id'].max())],
              'original_frozen_COM_absolute_physical': original.tolist(),
              'requested_frozen_COM_absolute_physical': target.tolist(),
              'actual_frozen_COM_absolute_physical': final.tolist(), 'COM_rounding_error': (final-target).tolist(),
              'coherent_delta_code_position_velocity': delta.tolist(),
              'original_snapshot_full_ID_phase_mass_check': original_check,
              'initial_output_comparison_allowances': {'gas_velocity_absolute_code': gas_atol,
                  'periodic_position_float32_ulps': pos_ulps, 'stellar_velocity_and_all_masses_exact': True},
              'protected_host_cube_and_nonphase_fields_verified': True, 'internal_rotation_applied': False,
              'quantization': changed, 'new_binary': str(dest), 'new_binary_sha256': sha(dest),
              'source_code_sha256': sha(__file__), 'wall_seconds': time.monotonic()-start,
              'warning': candidate.get('warning', 'Initial candidate still requires independent live evolution')}
    with receipt.open('x') as f:
        json.dump(result, f, indent=2); f.write('\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    args = parser.parse_args()
    result = build(json.loads(Path(args.config).read_text()))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
