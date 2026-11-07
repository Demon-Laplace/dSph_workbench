"""Physical and search regressions for perspective-preserving static poses."""
import os
os.environ.setdefault('DSPH_SKIP_CUBE_PATH', '1')
import unittest
import numpy as np
from scipy.spatial.transform import Rotation
from orbit_calibration_v1 import orientation_search as search


class OrientationTests(unittest.TestCase):
    def test_reflections_and_invalid_rotations_are_rejected(self):
        for matrix in (np.diag([-1., 1., 1.]), np.diag([1., 1., 1.001]),
                       np.full((3, 3), np.nan)):
            with self.assertRaises(ValueError):
                search.project_tan(np.zeros((3, 3)), matrix, np.eye(3), 100.)

    def test_depth_changes_overlap_and_matches_independent_rays(self):
        offset = np.array([[1., 0., -30.], [1., 0., 30.], [-2., 0., 0.]])
        q = Rotation.from_rotvec([.015, -.025, .01]).as_matrix()
        xy, distance = search.project_tan(offset, q, np.eye(3), 100., with_distance=True)
        h = offset @ q.T + [0., 0., 100.]
        np.testing.assert_allclose(xy, np.rad2deg(h[:, :2] / h[:, 2, None]), atol=1e-12)
        np.testing.assert_allclose(distance, np.linalg.norm(h, axis=1), atol=1e-12)
        baseline = search.project_tan(offset, np.eye(3), np.eye(3), 100.)
        self.assertGreater(baseline[0, 0], baseline[1, 0])
        self.assertGreater(np.max(abs(xy - baseline)), .2)

    def test_float32_projection_matches_canonical_astropy(self):
        from canonical_observer import production_observables, production_sky
        from canonical_display import tangent_coordinates
        rng = np.random.default_rng(54)
        dr = (rng.normal(size=(600, 3)) * [2., 1., 8.]).astype(np.float32)
        r = np.array([-18., -70., -115.]); v = np.array([50., -20., 80.])
        point = production_observables(r[None], v[None], np.array([0], dtype=np.uint64))
        sun = point['_native']['observer']['sun_position_kpc']
        sky = production_sky(point, r - sun)
        basis = np.array([sky['east_native'], sky['north_native'], sky['LOS_native']])
        q = Rotation.from_euler('xyz', [43., -71., 129.], degrees=True).as_matrix()
        coords = production_observables(dr @ q.T + r, np.tile(v, (len(dr), 1)),
                                        np.arange(len(dr), dtype=np.uint64))
        actual = np.column_stack(tangent_coordinates(production_sky(coords, r - sun)))
        fast, distance = search.project_tan(dr, q, basis, sky['center_distance_kpc'], with_distance=True)
        np.testing.assert_allclose(fast, actual, atol=2e-6, rtol=0)
        np.testing.assert_allclose(distance, coords['rh'], atol=5e-6, rtol=0)

    def test_rotating_offsets_and_velocities_preserves_com_and_internal_norms(self):
        rng = np.random.default_rng(3); mass = rng.uniform(.5, 2., 200)
        r = rng.normal(size=(200, 3)); v = rng.normal(size=(200, 3))
        r -= np.average(r, axis=0, weights=mass); v -= np.average(v, axis=0, weights=mass)
        q = search.sobol_rotations(3, 4)[4]
        np.testing.assert_allclose(np.average(r @ q.T, axis=0, weights=mass), 0., atol=1e-15)
        np.testing.assert_allclose(np.average(v @ q.T, axis=0, weights=mass), 0., atol=1e-15)
        np.testing.assert_allclose(np.linalg.norm(v @ q.T, axis=1), np.linalg.norm(v, axis=1))
        np.testing.assert_allclose(np.linalg.norm((r[1:] - r[:-1]) @ q.T, axis=1),
                                   np.linalg.norm(r[1:] - r[:-1], axis=1))

    def test_so3_sampling_is_broad_and_reproducible(self):
        q = search.sobol_rotations(10, 12)
        np.testing.assert_array_equal(q, search.sobol_rotations(10, 12))
        np.testing.assert_allclose(np.linalg.det(q), 1., atol=2e-15)
        np.testing.assert_allclose(np.mean(q[:, :, 0], axis=0), 0., atol=.015)
        np.testing.assert_allclose(np.mean(q[:, :, 0] ** 2, axis=0), 1/3, atol=.015)
        self.assertGreater(np.rad2deg(Rotation.from_matrix(q).magnitude()).max(), 175.)

    def test_light_leaving_window_is_not_renormalized_away(self):
        xy = np.array([[0., 0.], [.05, .02], [-.05, -.02]])
        image = search.cic_image(xy, np.ones(3)/3, .2, 32)
        outside = search.cic_image(xy + [1., 0.], np.ones(3)/3, .2, 32)
        self.assertAlmostEqual(image.sum(), 1.)
        self.assertEqual(outside.sum(), 0.)

    def test_sources_behind_observer_are_not_silently_removed(self):
        with self.assertRaises(ValueError):
            search.project_tan(np.array([[0., 0., -150.], [1., 1., 0.]]),
                               np.eye(3), np.eye(3), 100.)

    def test_refinement_improves_an_off_grid_tilt(self):
        rng = np.random.default_rng(78)
        dr = rng.normal(size=(1800, 3)) * [1.5, .7, .25]
        dr[:300] += [1.8, -.3, .8]
        q = Rotation.from_euler('xyz', [32., 45., -71.], degrees=True).as_matrix()
        seed = Rotation.from_rotvec(np.deg2rad([3., -4., 2.])).as_matrix() @ q
        reference = search.project_tan(dr, q, np.eye(3), 100.)
        objective = search.DensityObjective(dr, np.ones(len(dr)), np.eye(3), 100.,
                    reference, np.ones(len(dr)), [.4, 1.4], npix=24)
        fit = search.refine(objective, [seed], radius_deg=8., max_evaluations=280)[0]
        self.assertLess(fit['score'], objective(seed) * .01)
        self.assertLess(search.rotation_separation_deg(fit['matrix'], q), 1.)

    def test_nontrivial_3d_pose_recovered_beyond_pa_only(self):
        rng = np.random.default_rng(8)
        dr = rng.normal(size=(1600, 3)) * [1.6, .5, .15]
        dr[:250] += [2., .4, 1.]
        q = search.sobol_rotations(7, 91)[37]
        reference = search.project_tan(dr, q, np.eye(3), 110.)
        objective = search.DensityObjective(dr, np.ones(len(dr)), np.eye(3), 110.,
                            reference, np.ones(len(dr)), [.5, 1.5], npix=24)
        matrices, scores = search.global_search(objective, power=7, seed=91)
        self.assertLess(scores.min(), 1e-20)
        self.assertGreater(objective(np.eye(3)), .05)
        selected = search.distinct_candidates(matrices, scores, 3, 15.)
        self.assertLess(search.rotation_separation_deg(matrices[selected[0]], q), 1e-10)
        refined = search.refine(objective, matrices[selected[:1]], radius_deg=5., max_evaluations=30)
        self.assertLessEqual(refined[0]['score'], scores[selected[0]] + 1e-20)


if __name__ == '__main__':
    unittest.main()
