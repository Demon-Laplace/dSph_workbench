"""Full Astropy six-dimensional frame and portable orbit-entry regressions."""
from contextlib import redirect_stdout
from copy import deepcopy
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

SOURCE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SOURCE))
from orbit_calibration_v1 import cli
from orbit_calibration_v1.target_frame import FRAME, resolve_target


class FrameCorrections(unittest.TestCase):
    def config(self):
        return json.loads((SOURCE / 'orbit_calibration_v1/config.frame_check.example.json')
                          .read_text(encoding='utf-8'))

    def test_adopted_target_full_roundtrip_through_shared_observer(self):
        config = self.config()
        cli.verify_observation_source(SOURCE, config)
        contract = config['target_conversion_contract']
        self.assertLess(max(abs(v) for v in contract['roundtrip_error'].values()), 1e-9)
        self.assertAlmostEqual(config['target']['vlos_kms'], 55.3)
        self.assertAlmostEqual(config['target']['Vr_gc_kms'], -40.168618654688, places=8)
        self.assertFalse(contract['physical_model_rotation_applied'])

    def test_GC_spherical_frame_is_not_silently_accepted(self):
        config = self.config()
        config['target']['input_frame'] = 'li2021_galactocentric_spherical'
        with self.assertRaises(ValueError):
            cli.validate_config(config, require_cases=False)

    def test_frame_label_is_required(self):
        target = self.config()['target']
        del target['input_frame']
        with self.assertRaises(ValueError):
            resolve_target(target)

    def test_GC_spherical_columns_rejected_even_with_ICRS_label(self):
        target = self.config()['target']
        target['vtheta_kms'] = -100
        with self.assertRaises(ValueError):
            resolve_target(target)

    def test_inconsistent_derived_GC_target_is_rejected(self):
        target = self.config()['target']
        target['Vr_gc_kms'] = -36.8
        with self.assertRaises(ValueError):
            resolve_target(target)

    def test_positive_sun_side_is_not_inferred(self):
        config = self.config()
        config['sun_position_native_kpc'][0] *= -1
        with self.assertRaises(ValueError):
            cli.validate_config(config, require_cases=False)

    def test_heliocentric_zero_motion_stays_zero(self):
        target = self.config()['target']
        target.update(vlos_kms=0, mu_alpha_cosdec_masyr=0, mu_delta_masyr=0)
        restored = resolve_target(target)['restored_observables']
        for key in ('vlos_kms', 'mu_alpha_cosdec_masyr', 'mu_delta_masyr'):
            self.assertAlmostEqual(restored[key], 0, places=10)

    def test_random_sky_nonzero_PM_and_RA_wraparound(self):
        rng = np.random.default_rng(6020)
        for _ in range(24):
            target = self.config()['target']
            target.update(ra_deg=float(rng.uniform(0, 360)), dec_deg=float(rng.uniform(-85, 85)),
                          D_kpc=float(rng.uniform(5, 300)), vlos_kms=float(rng.uniform(-250, 250)),
                          mu_alpha_cosdec_masyr=float(rng.uniform(-1, 1)),
                          mu_delta_masyr=float(rng.uniform(-1, 1)))
            self.assertLess(max(abs(v) for v in resolve_target(target)['roundtrip_error'].values()), 1e-8)

    def test_frame_check_needs_no_simulation_cases(self):
        config = cli.validate_config(self.config(), require_cases=False)
        self.assertEqual(config['cases'], [])
        with self.assertRaises(ValueError):
            cli.validate_config(self.config())

    def test_comoving_native_probe_has_zero_LOS(self):
        from canonical_observer import astropy_native_direct_observables
        sun = np.array([-np.sqrt(8.122**2 - .0208**2), 0, .0208])
        velocity = np.array([12.9, 245.6, 7.78])
        result = astropy_native_direct_observables(
            np.array([sun + [4, 2, 7]]), velocity[None], np.array([1], np.uint64),
            host_center_raw=[0, 0, 0], host_bulk_raw=[0, 0, 0],
            length_to_kpc=1, velocity_to_kms=1, **FRAME)
        self.assertAlmostEqual(result['vlos_native_kms'][0], 0, places=11)

    def test_fixed_physical_mapping_cannot_be_replaced_by_camera_roll(self):
        config = self.config()
        config['physical_to_icrs'] = np.eye(3).tolist()
        cli.validate_config(config, require_cases=False)
        with self.assertRaisesRegex(ValueError, 'ICRS inverse'):
            cli.verify_observation_source(SOURCE, config)

    def test_module_and_all_compatibility_entries_run_outside_checkout(self):
        config = SOURCE / 'orbit_calibration_v1/config.frame_check.example.json'
        commands = [[str(SOURCE / 'calibrate_orbit.py'), '--config', str(config)],
                    [str(SOURCE / 'analyze_orbit_orientation.py'), '--config', str(config)],
                    [str(SOURCE / 'run_calibration_from_config.py'), str(config)]]
        with tempfile.TemporaryDirectory() as directory:
            for command in commands:
                with self.subTest(command=command):
                    result = subprocess.run([sys.executable, '-B', *command, '--stage', 'frame-check'],
                                            cwd=directory, capture_output=True, text=True, timeout=30)
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    record = json.loads(result.stdout)
                    self.assertEqual(record['state'], 'PASS')
                    self.assertFalse(record['conversion']['physical_model_rotation_applied'])
            self.assertEqual(list(Path(directory).iterdir()), [])

    def test_candidate_solver_requires_explicit_opt_in_before_model_reads(self):
        config = SOURCE / 'orbit_calibration_v1/config.example.json'
        result = subprocess.run([sys.executable, '-B', str(SOURCE / 'calibrate_orbit.py'),
                                 '--config', str(config), '--stage', 'solve'],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 2)
        self.assertIn('Use --physical-candidates explicitly', result.stderr)

    def test_positional_config_entry_rejects_old_CSV_schema_before_launch(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'old.json'
            path.write_text('{"schema_id": "legacy_csv_response"}', encoding='utf-8')
            result = subprocess.run([sys.executable, '-B', str(SOURCE / 'run_calibration_from_config.py'),
                                     str(path), '--stage', 'frame-check'],
                                    capture_output=True, text=True, timeout=30)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('Old CSV response coefficients', result.stderr)


if __name__ == '__main__':
    unittest.main(verbosity=2)
