import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "evolution_pipeline" / "fornax_evolution_pipeline.py"
SPEC = importlib.util.spec_from_file_location("fornax_evolution_pipeline", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
from coordinate_transform import CoordinateTransform


class EvolutionPipelineTests(unittest.TestCase):
    def test_empty_knn_cgm_shell_is_a_valid_missing_measurement(self):
        frame = pd.DataFrame(
            {
                "tp": [0],
                "x": [0.0],
                "y": [0.0],
                "z": [0.0],
                "vx": [0.0],
                "vy": [0.0],
                "vz": [0.0],
                "m": [1.0],
                "temp": [100.0],
                "nh": [1.0],
            }
        )
        result = MODULE.local_cgm_measurement(
            frame,
            np.zeros(3),
            np.zeros(3),
            {"method": "knn_shell", "minimum_particles": 1},
        )
        self.assertEqual(result["cgm_particle_count"], 0.0)
        self.assertTrue(np.isnan(result["ram_pressure_dyn_cm2"]))

    def test_heliocentric_velocity_transform_is_astropy_version_compatible(self):
        transformed = CoordinateTransform.to_heliocentric(
            np.array([10.0]),
            np.array([0.0]),
            np.array([0.0]),
            np.array([100.0]),
            np.array([200.0]),
            np.array([10.0]),
        )
        self.assertEqual(len(transformed), 6)
        self.assertTrue(all(np.all(np.isfinite(values)) for values in transformed))

    def test_actual_crossing_requires_tolerance(self):
        distance = np.array([150.0, 142.0, 139.3, 136.0])
        self.assertEqual(
            MODULE.actual_crossing_index(distance, 139.4, 0.25, "first_crossing"), 2
        )
        self.assertIsNone(
            MODULE.actual_crossing_index(distance, 139.4, 0.05, "first_crossing")
        )

    def test_pericentre_requires_post_minimum_rise(self):
        descending = np.array([200.0, 180.0, 160.0, 140.0, 130.0])
        self.assertIsNone(
            MODULE.reached_pericentre_index(
                descending,
                {"enabled": True, "minimum_post_points": 2, "minimum_rise_kpc": 0.5},
            )
        )
        crossed = np.array([200.0, 170.0, 140.0, 130.0, 130.3, 131.0, 132.0])
        self.assertEqual(
            MODULE.reached_pericentre_index(
                crossed,
                {"enabled": True, "minimum_post_points": 3, "minimum_rise_kpc": 0.5},
            ),
            3,
        )

    def test_smoothing_scale_is_configuration_driven(self):
        values = np.array([10.0, 9.0, 8.5, 8.0, 7.5, 7.0, 6.5])
        smoothed = MODULE.smooth_series(
            values, {"method": "savgol", "window_snapshots": 5, "polyorder": 2}
        )
        self.assertEqual(smoothed.shape, values.shape)
        self.assertTrue(np.all(np.isfinite(smoothed)))

    def test_derived_timescales_and_gas_fraction(self):
        frame = pd.DataFrame(
            {
                "analysis_config_sha256": ["x"] * 5,
                "snapshot": np.arange(5),
                "time_gyr": np.arange(5, dtype=float),
                "gas_mass_msun": np.array([10.0, 8.0, 6.0, 4.0, 2.0]) * 1.0e6,
                "re_major_kpc": np.full(5, 0.5),
                "sigma_los_kms": np.full(5, 10.0),
                "distance_heliocentric_kpc": np.array([160.0, 150.0, 140.0, 139.3, 138.0]),
                "distance_galactocentric_kpc": np.array([158.0, 148.0, 138.0, 137.0, 136.0]),
            }
        )
        config = {
            "smoothing": {
                "method": "savgol",
                "window_snapshots": 5,
                "polyorder": 2,
                "minimum_abs_rate_msun_per_gyr": 0.0,
            },
            "comparison_epoch": {
                "method": "heliocentric_distance",
                "target_heliocentric_distance_kpc": 139.4,
                "tolerance_kpc": 0.25,
                "branch": "first_crossing",
            },
            "pericentre_detection": {"enabled": True, "minimum_post_points": 2},
        }
        derived = MODULE.add_derived_columns(frame, config)
        self.assertTrue(
            np.allclose(
                derived["stellar_dynamical_time_gyr"], 0.9777922216807892 * 0.05
            )
        )
        self.assertTrue(np.all(np.isfinite(derived["tau_gas_smoothed_gyr"])))
        self.assertEqual(MODULE.gas_fraction(3.0, 1.0), 0.25)
        self.assertEqual(int(derived["is_comparison_epoch"].sum()), 1)
        self.assertEqual(int(derived["is_pericentre"].sum()), 0)


if __name__ == "__main__":
    unittest.main()
