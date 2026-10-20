import unittest

import numpy as np

from basefunc import Analysis
from snapshot_metrics import (detrended_dispersion_in_aperture, stellar_mass_components,
    member_star_projected_kinematics, old_star_projected_kinematics)


class SigmaReApertureTests(unittest.TestCase):
    def test_new_member_star_changes_dispersion_without_changing_old_light_sample(self):
        snapshot = dict(x_kpc=np.array([-.2, .2, 0., 0., 0.]),
            y_kpc=np.array([0., 0., -.2, .2, 0.]),
            star_coords={'vlos': np.array([0., 0., 0., 0., 10.])},
            fixed_old_local_mask=np.array([True, True, True, True, False]))
        members = member_star_projected_kinematics(snapshot)
        old = old_star_projected_kinematics(snapshot)
        result = detrended_dispersion_in_aperture(
            members['x_kpc'], members['y_kpc'], members['vlos'], .5)
        self.assertEqual(result['nstar'], 5)
        self.assertAlmostEqual(result['sigma'], 4.)
        self.assertEqual(len(old['vlos']), 4)
        self.assertAlmostEqual(np.std(old['vlos_detrended']), 0.)

    def test_total_stellar_mass_includes_new_stars(self):
        result = stellar_mass_components(
            masses=np.array([1.0, 2.0, 3.0, 4.0]),
            aperture_mask=np.array([True, True, True, False]),
            old_star_mask=np.array([True, True, False, False]),
        )
        self.assertEqual(result["total"], 6.0)
        self.assertEqual(result["old_tracer"], 3.0)
        self.assertEqual(result["new"], 3.0)

    def test_circular_aperture_and_gradient_removal(self):
        x = np.array([0.0, 0.4, -0.4, 0.0, 0.0, 1.2])
        y = np.array([0.0, 0.0, 0.0, 0.4, -0.4, 0.0])
        residual = np.array([2.0, -0.5, -0.5, -0.5, -0.5, 50.0])
        vlos = 12.0 + 2.0 * x - 3.0 * y + residual

        result = detrended_dispersion_in_aperture(
            x,
            y,
            vlos,
            radius_kpc=0.5,
            circular=True,
        )

        self.assertEqual(result["nstar"], 5)
        self.assertTrue(np.isfinite(result["sigma"]))
        self.assertLess(result["sigma"], 2.0)
        self.assertAlmostEqual(result["gradient"]["a"], 2.0, places=10)
        self.assertAlmostEqual(result["gradient"]["b"], -3.0, places=10)

    def test_half_light_radius_ignores_invalid_weights(self):
        x = np.array([0.0, 1.0, 2.0, 100.0])
        y = np.zeros_like(x)
        mass = np.array([1.0, 1.0, 2.0, np.nan])

        radius = Analysis.half_light_radius(
            x,
            y,
            mass,
            ep=0.0,
        )

        self.assertAlmostEqual(radius, 1.0)


if __name__ == "__main__":
    unittest.main()
