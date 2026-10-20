"""Small native-vector regressions; no Astropy, snapshots, or Slurm required."""
from __future__ import annotations

import argparse
import ast
import copy
import inspect
import json
from pathlib import Path
import subprocess
import sys
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
PARSER = argparse.ArgumentParser(add_help=False)
PARSER.add_argument('--source-dir', type=Path, default=HERE.parent)
ARGS, REMAINING = PARSER.parse_known_args()
SOURCE = ARGS.source_dir.resolve()
sys.path.insert(0, str(SOURCE))

from orbit_calibration_v1 import cli
from orbit_calibration_v1 import solve_orbit_candidates as solver
from orbit_calibration_v1.math import candidate_math as math


def example():
    return json.loads((SOURCE / 'orbit_calibration_v1/config.example.json').read_text(encoding='utf-8'))


def frame(position, velocity, number, mass=None):
    mass = np.array([1., 2., 3.] if mass is None else mass, dtype=float)
    dr = np.array([[-3., 0., .5], [0., 3., .5], [1., -2., -.5]])
    dv = np.array([[0., -3., .5], [3., 0., .5], [-2., 1., -.5]])
    dr -= np.average(dr, axis=0, weights=mass)
    dv -= np.average(dv, axis=0, weights=mass)
    return dict(id=np.array([2**53 + 17, 2**53 + 101, 2**63 + 7], dtype=np.uint64),
                mass=mass, dr=dr, dv=dv, r=np.asarray(position, float),
                v=np.asarray(velocity, float), state=dict(numsp=number,
                header_time=1 + .01 * (number - 10), unit_time_to_Gyr=.9777929880599285,
                cache_sha256=('%064x' % number)))


class CanonicalOrbitTests(unittest.TestCase):
    def test_cli_worker_modules_start_in_real_subprocesses_without_math_shadowing(self):
        expected = {'orbit_calibration_v1.extract_frame',
                    'orbit_calibration_v1.solve_orbit_candidates'}
        dispatches = []
        for node in ast.walk(ast.parse(inspect.getsource(cli.main))):
            if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                    and isinstance(node.func.value, ast.Name)
                    and node.func.value.id == 'subprocess' and node.func.attr == 'run'):
                command = node.args[0]
                self.assertIsInstance(command, ast.List)
                self.assertEqual(ast.unparse(command.elts[0]), 'sys.executable')
                prefix = [ast.literal_eval(item) for item in command.elts[1:4]]
                self.assertEqual(prefix[:2], ['-B', '-m'])
                dispatches.append(prefix[2])
        self.assertEqual(set(dispatches), expected)
        self.assertEqual(len(dispatches), len(expected))
        for module in sorted(expected):
            with self.subTest(module=module):
                result = subprocess.run([sys.executable, '-B', '-m', module, '--help'],
                                        cwd=str(SOURCE), capture_output=True, text=True,
                                        timeout=30, check=False)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn('usage:', result.stdout.lower())

    def test_example_is_explicit_and_validation_does_not_rewrite_input(self):
        config = example()
        before = copy.deepcopy(config)
        checked = cli.validate_config(config)
        self.assertEqual(config, before)
        self.assertEqual(checked['schema_id'], 'canonical_orbit_recalibration_v1')
        self.assertEqual(checked['target']['morphology_PA_E_of_N_deg'], 43.)
        self.assertEqual([case['model'] for case in checked['cases']], ['Fornax7004', 'Fornax7005'])

    def test_legacy_schema_observer_and_contract_sha_are_rejected(self):
        for key, value in [('schema_id', 'legacy_orbit_v2'),
                           ('canonical_observer_id', 'legacy_positive_x'),
                           ('observer_config_sha256', '0' * 64)]:
            with self.subTest(key=key):
                config = example()
                config[key] = value
                with self.assertRaises(ValueError):
                    cli.validate_config(config)

    def test_positive_x_sun_and_opposite_solar_rotation_are_rejected(self):
        for key, axis in [('sun_position_native_kpc', 0), ('sun_velocity_native_kms', 1)]:
            with self.subTest(key=key):
                config = example()
                config[key][axis] *= -1
                with self.assertRaises(ValueError):
                    cli.validate_config(config)

    def test_reflection_and_nonorthogonal_physical_maps_are_rejected(self):
        for matrix in [np.diag([-1., 1., 1.]), np.diag([1., 1., 1.01]), np.eye(2)]:
            config = example()
            config['physical_to_icrs'] = matrix.tolist()
            with self.subTest(matrix=matrix.tolist()), self.assertRaises(ValueError):
                cli.validate_config(config)

    def test_invalid_or_incomplete_targets_are_rejected(self):
        invalid = [('ra_deg', -1.), ('ra_deg', 360.), ('dec_deg', -91.),
                   ('dec_deg', 91.), ('D_kpc', 0.), ('D_kpc', -139.4),
                   ('morphology_PA_E_of_N_deg', -1.), ('morphology_PA_E_of_N_deg', 180.)]
        invalid += [(key, float('nan')) for key in example()['target']]
        invalid += [('vlos_kms', float('inf')), ('mu_delta_masyr', float('-inf'))]
        for key, value in invalid:
            config = example()
            config['target'][key] = value
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                cli.validate_config(config)
        config = example()
        del config['target']['mu_alpha_cosdec_masyr']
        with self.assertRaises(ValueError):
            cli.validate_config(config)

    def test_model_ids_duplicate_models_and_bad_snapshot_brackets_are_rejected(self):
        bad = []
        config = example(); config['cases'][0]['model'] = '../Fornax7004'; bad.append(config)
        config = example(); config['cases'][1]['model'] = 'Fornax7004'; bad.append(config)
        config = example(); config['cases'][0]['initial'] = True; bad.append(config)
        config = example(); config['cases'][0]['latest'] = 0; bad.append(config)
        config = example(); config['cases'][0]['GC_radius_bracket'] = [275, 277]; bad.append(config)
        for config in bad:
            with self.subTest(config=config), self.assertRaises(ValueError):
                cli.validate_config(config)

    def test_adopted_native_target_recovers_all_six_heliocentric_components(self):
        config = example()
        # Fixed native regression fixture, independently recorded for this target.
        r = np.array([-39.66138111144952, -48.25711878005904, -126.89994789458696])
        v = np.array([17.932795557157682, -101.37607711932961, 77.71745722001768])
        observed = solver.independent_point_observables(r, v, config)
        target = config['target']
        for actual, expected in [('RA_deg', 'ra_deg'), ('Dec_deg', 'dec_deg'),
                                 ('D_kpc', 'D_kpc'), ('vlos_kms', 'vlos_kms'),
                                 ('mu_alpha_cosdec_masyr', 'mu_alpha_cosdec_masyr'),
                                 ('mu_delta_masyr', 'mu_delta_masyr')]:
            self.assertAlmostEqual(observed[actual], target[expected], delta=1e-9)
        self.assertAlmostEqual(observed['R_gc_kpc'], 141.44034587156966, delta=1e-9)
        self.assertGreater(abs(observed['R_gc_kpc'] - target['D_kpc']), 2.)
        self.assertAlmostEqual(observed['Vr_gc_kms'], -40.168618654688025, delta=1e-9)

    def test_empirical_energy_and_angular_momentum_recover_known_seed(self):
        r0, v0 = np.array([100., 0., 0.]), np.array([-30., 40., 0.])
        rf, vf = np.array([50., 0., 0.]), np.array([-40., 60., 0.])
        coefficient = math.empirical_coefficients(r0, v0, rf, vf)
        self.assertAlmostEqual(coefficient['C_E'], 2700., delta=1e-10)
        self.assertAlmostEqual(coefficient['k_L'], .75, delta=1e-14)
        solved = math.solve_initial_velocity(r0, v0, rf, vf, 2700., .75, radial_sign=-1)
        np.testing.assert_allclose(solved['velocity'], v0, rtol=0, atol=1e-12)
        self.assertFalse(solved['pericentre_constrained'])

    def test_infeasible_energy_and_radial_orbit_plane_are_rejected(self):
        with self.assertRaises(ValueError):
            math.solve_initial_velocity([100, 0, 0], [-30, 40, 0], [50, 0, 0], [-40, 60, 0], 6000., .75)
        with self.assertRaises(ValueError):
            math.solve_initial_velocity([100, 0, 0], [-30, 40, 0], [50, 0, 0], [-40, 60, 0], 2700., 0.)
        with self.assertRaises(ValueError):
            math.decompose_state([100, 0, 0], [-50, 0, 0])

    def test_orbit_frame_registration_preserves_directed_basis_but_not_speed_claim(self):
        rotation = math.orbit_plane_rotation([50, 0, 0], [-40, 60, 0], [0, 80, 0], [-30, -10, 0])
        expected = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        np.testing.assert_allclose(rotation, expected, rtol=0, atol=1e-14)
        self.assertAlmostEqual(np.linalg.det(rotation), 1., delta=1e-14)
        self.assertNotAlmostEqual(np.linalg.norm(rotation @ [-40, 60, 0]), np.linalg.norm([-30, -10, 0]))

    def test_morphological_pa_roll_has_north_to_east_sign(self):
        east, north, los = np.eye(3)
        axis = np.sin(np.radians(30.)) * east + np.cos(np.radians(30.)) * north
        rotation = math.pa_roll_rotation(los, 13.)
        self.assertAlmostEqual(math.axis_pa_deg(rotation @ axis, east, north), 43., delta=1e-12)
        np.testing.assert_allclose(rotation @ los, los, rtol=0, atol=1e-14)
        self.assertAlmostEqual(math.pa_difference_deg(1., 179.), 2., delta=1e-14)

    def test_one_pass_placement_keeps_all_input_centre_and_bulk_offsets(self):
        r, v = np.array([100., 200., -50.]), np.array([-10., 20., 30.])
        cm, bm = [1., 2., 3.], [4., -5., 6.]
        cd, bd = [-.2, .4, .6], [-1., 2., -3.]
        placement = math.placement_for_target(r, v, [0, 0, 0], [0, 0, 0], cm, bm, cd, bd)
        np.testing.assert_allclose(placement, [101.2, 201.6, -47.6, -5., 13., 39.], atol=1e-12, rtol=0)
        galaxy = np.r_[placement, [90., 0., 0.], [0., 0., 0.]]
        cube = np.array([500., 500., 500., -20., 11., 13., 41., 72., -15., 0., 0., 0.])
        actual_r, actual_v = math.initial_state_from_merge(galaxy, cube, cm, bm, cd, bd)
        np.testing.assert_allclose(actual_r, r, atol=1e-12, rtol=0)
        np.testing.assert_allclose(actual_v, v, atol=1e-12, rtol=0)

    def test_nonzero_cube_rotation_conjugates_internal_correction(self):
        old, cube = [103.42, 35.84, -1.8], [21., -12., 37.]
        internal = math.pa_roll_rotation([.2, -.4, .8], 13.)
        plane = math.axis_angle_rotation([1., 2., 3.], 17.)
        result = math.corrected_internal_angles(old, cube, internal, plane)
        q = math.legacy_matrix(cube)
        expected = internal @ plane @ q @ math.legacy_matrix(old)
        np.testing.assert_allclose(q @ result['matrix'], expected, rtol=0, atol=4e-7)
        self.assertLess(result['native_matrix_error'], 4e-7)

    def test_crossing_solves_interpolated_native_vector_not_scalar_radii(self):
        left, right = frame([8, 6, 0], [-3, 5, 1], 10), frame([6, 0, 0], [-1, 4, 2], 11)
        anchor = solver.radial_anchor(left, right, 8., np.zeros(3))
        fraction = (104. - np.sqrt(5056.)) / 80.
        self.assertAlmostEqual(anchor['fraction'], fraction, delta=1e-12)
        self.assertAlmostEqual(np.linalg.norm(anchor['r']), 8., delta=2e-12)
        self.assertGreater(abs(anchor['fraction'] - .5), .05)
        sun = np.array([-8.1219733661223, 0., .0208])
        left['r'] += sun; right['r'] += sun
        heliocentric = solver.radial_anchor(left, right, 8., sun)
        self.assertAlmostEqual(heliocentric['fraction'], fraction, delta=1e-12)
        self.assertAlmostEqual(np.linalg.norm(heliocentric['r'] - sun), 8., delta=2e-12)
        self.assertGreater(abs(np.linalg.norm(heliocentric['r']) - 8.), 1.)

    def test_id_bound_interpolation_records_and_removes_cache_first_moments(self):
        left = frame([8, 6, 0], [-3, 5, 1], 10)
        right = frame([6, 0, 0], [-1, 4, 2], 11, [1.0000001, 1.9999998, 3.0000003])
        right['dr'] += [1e-7, -2e-7, 3e-7]
        right['dv'] += [-1e-7, 2e-7, -3e-7]
        anchor = solver.interpolate_frames(left, right, .3)
        np.testing.assert_allclose(np.average(anchor['dr'], axis=0, weights=anchor['mass']), 0, atol=1e-14)
        np.testing.assert_allclose(np.average(anchor['dv'], axis=0, weights=anchor['mass']), 0, atol=1e-14)
        np.testing.assert_allclose(anchor['r'], .7 * left['r'] + .3 * right['r'], rtol=0, atol=1e-14)
        self.assertGreater(np.max(abs(anchor['dr_recentering_kpc'])), 0.)
        np.testing.assert_array_equal(anchor['id'], left['id'])
        self.assertEqual(anchor['id'].dtype, np.dtype('uint64'))
        right['id'] = right['id'][::-1]
        with self.assertRaises(AssertionError):
            solver.interpolate_frames(left, right, .3)


if __name__ == '__main__':
    unittest.main(argv=[sys.argv[0]] + REMAINING)
