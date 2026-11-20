import importlib.util
from pathlib import Path
import unittest

import numpy as np

PATH = Path(__file__).resolve().parents[2]/'local_tools/ic_generation/disk_diagnostics.py'
SPEC = importlib.util.spec_from_file_location('disk_diagnostics_test', PATH)
DIAG = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DIAG)


class Kepler:
    def potential(self, x):
        return -1/np.linalg.norm(x, axis=1)

    def force(self, x):
        return -x/np.linalg.norm(x, axis=1)[:, None]**3

    def Rcirc(self, E):
        return -.5/E


class DiskMeasurements(unittest.TestCase):
    def test_circularity_distinguishes_circular_retrograde_and_eccentric(self):
        x = np.array([[1., 0., 0.]]*3)
        v = np.array([[0., 1., 0.], [0., -1., 0.], [.7, .7, 0.]])
        epsilon, energy = DIAG.circularity(Kepler(), x, v)
        self.assertAlmostEqual(epsilon[0], 1.)
        self.assertAlmostEqual(epsilon[1], -1.)
        self.assertLess(epsilon[2], .8)
        self.assertTrue(np.all(energy < 0))

    def test_all_prograde_need_not_mean_radially_cold(self):
        phi = np.linspace(0., 2*np.pi, 1000, endpoint=False)
        x = np.column_stack((np.cos(phi), np.sin(phi), np.zeros(len(phi))))
        radial = np.where(np.arange(len(phi)) % 2, 10., -10.)
        v = np.column_stack((radial*np.cos(phi)-3*np.sin(phi),
                             radial*np.sin(phi)+3*np.cos(phi), np.zeros(len(phi))))
        result = DIAG.kinematics(x, v, np.ones(len(phi)), np.ones(len(phi), bool))
        self.assertAlmostEqual(result['prograde_mass_fraction'], 1.)
        self.assertAlmostEqual(result['V_over_sigma_R'], .3)

    def test_toomre_units_and_invalid_energy_do_not_hide_particles(self):
        self.assertAlmostEqual(DIAG.toomre_q(5., 30., 10., 2.), 150/67.2)
        epsilon, _ = DIAG.circularity(Kepler(), np.array([[1., 0., 0.]]),
                                     np.array([[0., 2., 0.]]))
        result = DIAG.epsilon_summary(epsilon, np.ones(1), np.ones(1, bool))
        self.assertEqual(result['mass_fraction_without_bound_circular_reference'], 1.)
        self.assertEqual(result['mass_fraction_epsilon_gt_0p8'], 0.)


if __name__ == '__main__':
    unittest.main()
