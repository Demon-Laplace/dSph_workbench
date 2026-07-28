import unittest

import numpy as np

from basefunc import Analysis
from snapshot_metrics import detrended_dispersion_in_aperture


class SigmaReApertureTests(unittest.TestCase):
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
