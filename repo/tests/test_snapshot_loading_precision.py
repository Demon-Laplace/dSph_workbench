from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

from snapshot_io import SnapshotLoaderMixin


class SnapshotVectorLoadingTests(unittest.TestCase):
    def test_all_vector_columns_and_uint64_ids_across_storage_layouts(self):
        for dtype in (np.float32, np.float64):
            for compressed in (False, True):
                for parts in (1, 2):
                    with self.subTest(dtype=dtype, compressed=compressed, parts=parts), tempfile.TemporaryDirectory() as tmp:
                        base = str(Path(tmp) / 'snapshot_000')
                        expected = {tp: [] for tp in range(5)}
                        for part in range(parts):
                            filename = base + ('.'+str(part) if parts > 1 else '') + '.hdf5'
                            with h5py.File(filename, 'w') as stream:
                                header = stream.create_group('Header')
                                for key, value in dict(Time=1.25, NumPart_Total=[parts]*5+[0],
                                        NumPart_ThisFile=[1]*5+[0], MassTable=[0, .0002, 0, 0, 0, 0],
                                        NumFilesPerSnapshot=parts, Flag_Cooling=0).items():
                                    header.attrs[key] = value
                                for tp in range(5):
                                    group = stream.create_group('PartType'+str(tp))
                                    position = np.array([[tp+.123456789+part, -tp-.987654321, part+.333333333]], dtype=dtype)
                                    velocity = np.array([[-17.123456789-tp, 28.987654321+part, -39.333333333]], dtype=dtype)
                                    options = dict(compression='gzip') if compressed else {}
                                    group.create_dataset('Coordinates', data=position, **options)
                                    group.create_dataset('Velocities', data=velocity, **options)
                                    identifier = np.uint64(2**60+100*tp+part)
                                    group['ParticleIDs'] = np.array([identifier], dtype=np.uint64)
                                    group['Masses'] = np.array([.0001], dtype=dtype)
                                    if tp == 4:
                                        group['StellarFormationTime'] = np.array([.25], dtype=dtype)
                                    expected[tp].append((position[0], velocity[0], identifier))
                        first = base + ('.0' if parts > 1 else '') + '.hdf5'
                        class Loader(SnapshotLoaderMixin):
                            boxsize = 0
                            def _resolve_snapshot_hdf5(self):
                                return first, base, '.hdf5'
                        frame, age = Loader()._load_snapshot_direct_hdf5(dftype=3,
                            box=False, gas_fields=(), star_fields=('birth',), particle_types=tuple(range(5)))
                        rows = [row for tp in range(5) for row in expected[tp]]
                        positions = np.array([row[0] for row in rows], dtype=np.float32)
                        velocities = np.array([row[1] for row in rows], dtype=np.float32)
                        np.testing.assert_array_equal(frame[['x', 'y', 'z']].to_numpy(), positions)
                        np.testing.assert_array_equal(frame[['vx', 'vy', 'vz']].to_numpy(), velocities)
                        np.testing.assert_array_equal(frame['id'].to_numpy(), np.array([row[2] for row in rows], dtype=np.uint64))
                        self.assertEqual(frame['x'].dtype, np.dtype('float32'))
                        self.assertEqual(frame['vx'].dtype, np.dtype('float32'))
                        self.assertEqual(frame['id'].dtype, np.dtype('uint64'))
                        self.assertEqual(frame['m'].dtype, np.dtype('float64'))
                        self.assertEqual(age, 1.25)


if __name__ == '__main__':
    unittest.main()
