import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np

PATH = Path(__file__).resolve().parents[2]/'local_tools/ic_generation/rotation_diagnostics.py'
SPEC = importlib.util.spec_from_file_location('rotation_diagnostics_test', PATH)
DIAG = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DIAG)


class RotationMeasurements(unittest.TestCase):
    def ring(self):
        phi = np.linspace(0, 2*np.pi, 1000, endpoint=False)
        x = np.column_stack((np.cos(phi), np.sin(phi), .05*np.sin(3*phi)))
        v = np.column_stack((-6*np.sin(phi), 6*np.cos(phi), np.zeros(len(phi))))
        return x, v, np.ones(len(phi))

    def test_streaming_is_signed_not_absolute_speed(self):
        x, v, m = self.ring()
        v[::2] *= -1
        result = DIAG.moments(x, v, m, np.array([0., 0., 1.]), np.ones(len(m), bool))
        self.assertAlmostEqual(result['mean_vphi_kms'], 0., places=12)
        self.assertAlmostEqual(result['mean_abs_vphi_kms'], 6., places=12)
        self.assertAlmostEqual(result['sigma_R_phi_z_kms'][1], 6., places=12)
        self.assertAlmostEqual(result['prograde_mass_fraction'], .5, places=12)
        self.assertAlmostEqual(result['V_over_sigma_1d_cyl'], 0., places=10)

    def test_translation_and_bulk_do_not_create_or_remove_rotation(self):
        x, v, m = self.ring()
        result = DIAG.diagnose({'old': (x+[10, -3, 6], v+[20, -30, 4], m)}, 'test', 0.)
        old = result['populations']['old']['3D5']
        self.assertAlmostEqual(old['mean_vphi_kms'], 6., places=10)
        self.assertAlmostEqual(old['prograde_mass_fraction'], 1.)
        np.testing.assert_allclose(result['density_short_axis'], [0., 0., 1.], atol=1e-12)

    def test_spatial_membership_does_not_use_particle_ids(self):
        x, v, m = self.ring()
        x = np.vstack((x, [20., 0., 0.]))
        v = np.vstack((v, [0., -100., 0.]))
        m = np.append(m, .1)
        result = DIAG.diagnose({'old': (x, v, m)}, 'test', 0.)
        self.assertEqual(result['populations']['old']['3D5']['N'], 1000)
        self.assertAlmostEqual(result['populations']['old']['3D5']['mean_vphi_kms'], 6., places=10)

    def test_startup_needs_both_elapsed_threshold_and_nonzero_integration(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory/'output').mkdir()
            info = directory/'output/info.txt'
            info.write_text('Sync-Point 0, Time: 0, Systemstep: 0\n')
            self.assertIsNone(DIAG.startup_evidence(directory, 200))
            info.write_text('Sync-Point 4, Time: 0.0002, Systemstep: 0.00005\n')
            self.assertIsNone(DIAG.startup_evidence(directory, 139))
            evidence = DIAG.startup_evidence(directory, 155)
            self.assertEqual(evidence['state'], 'started')
            self.assertAlmostEqual(evidence['time_code'], .0002)


if __name__ == '__main__':
    unittest.main()
