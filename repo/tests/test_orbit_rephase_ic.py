"""Physical preservation and fail-closed tests for orbit-only IC translation."""
import json
from pathlib import Path
import struct
import tempfile
import unittest

import h5py
import numpy as np

from orbit_calibration_v1 import rephase_ic as ic


class RephaseTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root/'old').mkdir(); (self.root/'new').mkdir()
        self.source = self.root/'old/source.IC_000'
        counts = [3, 1, 3, 0, 0, 0]
        # IDs 4..6: one dwarf gas particle and two dwarf stars.
        ids = np.array([0, 1, 4, 2, 3, 5+(1<<25), 6+(3<<25)], dtype='<u4')
        p = (np.arange(21).reshape(7, 3) / 8 + 520).astype('<f4')
        v = (np.arange(21).reshape(7, 3) / 10 - 2).astype('<f4')
        m = np.array([10, 11, 2, 30, 40, 3, 4], dtype='<f4')
        header = bytearray(256)
        struct.pack_into('<6I', header, 0, *counts)
        struct.pack_into('<6I', header, 96, *counts)
        struct.pack_into('<i', header, 124, 1)
        blocks = [header, p.tobytes(), v.tobytes(), ids.tobytes(), m.tobytes(),
                  np.array([7, 8, 9], dtype='<f4').tobytes(), np.zeros(3, dtype='<f4').tobytes()]
        with self.source.open('wb') as f:
            for block in blocks:
                f.write(struct.pack('<I', len(block))); f.write(block); f.write(struct.pack('<I', len(block)))
        snap = self.root/'old/snapshot_000.hdf5'
        with h5py.File(snap, 'w') as f:
            h = f.create_group('Header'); h.attrs['Time'] = 0.; h.attrs['NumPart_Total'] = counts
            cursor = 0
            for tp, n in enumerate(counts):
                if n:
                    g = f.create_group('PartType%d' % tp)
                    # Intentionally permute rows so an order-based comparison is wrong.
                    rows = np.arange(cursor, cursor+n)[::-1]
                    for name, data in [('ParticleIDs', ids), ('Coordinates', p), ('Velocities', v), ('Masses', m)]:
                        g.create_dataset(name, data=data[rows])
                cursor += n
        dwarf = self.root/'dwarf.ini'
        dwarf.write_text('1 0 2 0\n0 0\n%d 0\n%d 0\n' % (1+(1<<25), 2+(3<<25)))
        cohort = self.root/'cohort.npz'; np.savez(cohort, particle_id=ids[-2:].astype(np.uint64))
        old = ic.weighted_com(p, v, m, np.array([5, 6]), 1.0001, 2.)
        target = old + [13, 17, -4, -23, 19, 5]
        candidate = self.root/'candidate.json'
        candidate.write_text(json.dumps({'canonical_observer_id': 'astropy_native_direct_v1',
                                         'initial_absolute_6D_seed': target.tolist()}))
        state = self.root/'state.json'
        state.write_text(json.dumps({'old_N': 2, 'old_P_absolute_kpc': old[:3].tolist(),
                                     'old_V_absolute_kms': old[3:].tolist()}))
        self.cfg = {'source_binary': str(self.source), 'source_snapshot': str(snap), 'cohort_file': str(cohort),
                    'candidate_file': str(candidate), 'initial_state_file': str(state),
                    'destination_binary': str(self.root/'new/new.IC_000'), 'receipt': str(self.root/'new/receipt.json'),
                    'canonical_observer_id': 'astropy_native_direct_v1', 'length_to_kpc': 1.0001,
                    'velocity_to_kms': 2., 'expected_counts': counts, 'expected_cohort_N': 2,
                    'dwarf_ascii':str(dwarf),'expected_dwarf_counts':[1,0,2,0], 'dwarf_id_offset':4,'box_size_code': 1000.,
                    'expected_input_sha256': {}}
        self.cfg['expected_input_sha256'] = {name: ic.sha(self.cfg[name]) for name in
                                            ['source_binary', 'cohort_file', 'candidate_file', 'initial_state_file']}

    def test_complete_dwarf_shift_protects_host_mass_ids_and_gas(self):
        compressed = self.root/'delta_cohort.npz'
        np.savez(compressed, first=np.array([5+(1<<25)], dtype=np.uint64), delta=np.array([1+(2<<25)], dtype=np.uint64))
        np.testing.assert_array_equal(ic.cohort_ids(compressed), ic.cohort_ids(self.cfg['cohort_file']))
        self.cfg['cohort_file'] = str(compressed)
        self.cfg['expected_input_sha256']['cohort_file'] = ic.sha(compressed)
        old_sha = ic.sha(self.source)
        result = ic.build(self.cfg)
        self.assertEqual(ic.sha(self.source), old_sha)
        self.assertEqual(result['dwarf_N'], 3)  # Gas is moved even though cohort contains only stars.
        self.assertEqual(result['state'], 'VALIDATED_NEW_ORBIT_ONLY_IC')
        self.assertFalse(result['internal_rotation_applied'])
        _, _, _, a = ic.binary_arrays(self.source)
        _, _, _, b = ic.binary_arrays(self.cfg['destination_binary'])
        selected = a['id'] >= 4
        for key in ['position', 'velocity']:
            np.testing.assert_array_equal(a[key][~selected], b[key][~selected])
            # Internal pair separations/relative velocities retain float32 precision.
            np.testing.assert_allclose(np.diff(a[key][selected], axis=0), np.diff(b[key][selected], axis=0), atol=4e-6)
        for key in ['id', 'mass', 'gas_u', 'gas_rho']:
            np.testing.assert_array_equal(a[key], b[key])

    def test_refuses_existing_destination(self):
        Path(self.cfg['destination_binary']).write_bytes(b'protected')
        with self.assertRaises(ValueError): ic.build(self.cfg)
        self.assertEqual(Path(self.cfg['destination_binary']).read_bytes(), b'protected')

    def test_refuses_wrong_snapshot_phase(self):
        with h5py.File(self.cfg['source_snapshot'], 'r+') as f:
            f['PartType0/Coordinates'][0, 0] += 1
        with self.assertRaisesRegex(ValueError, 'disagrees'): ic.build(self.cfg)
        self.assertFalse(Path(self.cfg['destination_binary']).exists())

    def test_refuses_changed_candidate_sha(self):
        with Path(self.cfg['candidate_file']).open('a') as f: f.write(' ')
        with self.assertRaisesRegex(ValueError, 'checksum'): ic.build(self.cfg)

    def test_refuses_wrong_initial_cohort_COM(self):
        self.cfg['expected_input_sha256'].pop('initial_state_file')
        state = Path(self.cfg['initial_state_file'])
        x = json.loads(state.read_text()); x['old_V_absolute_kms'][0] += 1; state.write_text(json.dumps(x))
        with self.assertRaisesRegex(ValueError, 'COM differs'): ic.build(self.cfg)

    def test_refuses_invalid_origin_offset(self):
        self.cfg['dwarf_id_offset'] = 6
        with self.assertRaisesRegex(ValueError, 'outside'): ic.build(self.cfg)

    def test_refuses_extra_format_blocks(self):
        with self.source.open('ab') as f: f.write(struct.pack('<I', 0)*2)
        with self.assertRaisesRegex(ValueError, 'additional'): ic.binary_arrays(self.source)

    def test_documented_gas_readout_and_periodic_coordinates(self):
        self.cfg.update(initial_snapshot_gas_velocity_atol_code=2e-5, initial_snapshot_periodic_position_ulps=1.)
        with h5py.File(self.cfg['source_snapshot'], 'r+') as f:
            f['PartType0/Velocities'][0, 0] += 1e-5
            f['PartType0/Coordinates'][0, 0] += self.cfg['box_size_code']
        result = ic.build(self.cfg)
        self.assertTrue(result['initial_output_comparison_allowances']['stellar_velocity_and_all_masses_exact'])

    def test_gas_readout_allowance_does_not_relax_stellar_velocity(self):
        self.cfg.update(initial_snapshot_gas_velocity_atol_code=2e-5, initial_snapshot_periodic_position_ulps=1.)
        with h5py.File(self.cfg['source_snapshot'], 'r+') as f:
            f['PartType2/Velocities'][0, 0] += 1e-5
        with self.assertRaisesRegex(ValueError, 'disagrees'): ic.build(self.cfg)


if __name__ == '__main__':
    unittest.main()
