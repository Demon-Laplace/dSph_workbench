"""Independent ray and solid-angle checks for the production TAN renderer.

Only display helpers are imported.  No snapshot, observer, cohort, or production
output is loaded or changed.  Run directly, or via unittest discovery.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
PARSER = argparse.ArgumentParser()
PARSER.add_argument('--source-dir', type=Path, default=HERE.parent)
PARSER.add_argument('--output-json', type=Path)
ARGS, _ = PARSER.parse_known_args()
SOURCE = ARGS.source_dir.resolve()
sys.path.insert(0, str(SOURCE))
import canonical_display as display


def basis(ra_deg, dec_deg):
    """Independent Cartesian ICRS basis, ordered centre/east/north."""
    a, d = np.deg2rad([ra_deg, dec_deg])
    return (np.array([np.cos(d)*np.cos(a), np.cos(d)*np.sin(a), np.sin(d)]),
            np.array([-np.sin(a), np.cos(a), 0.]),
            np.array([-np.sin(d)*np.cos(a), -np.sin(d)*np.sin(a), np.cos(d)]))


def absolute_rays(ra_deg, dec_deg):
    a, d = np.deg2rad(ra_deg), np.deg2rad(dec_deg)
    return np.stack([np.cos(d)*np.cos(a), np.cos(d)*np.sin(a), np.sin(d)], axis=-1)


def triangle_solid_angle(a, b, c):
    """Spherical triangle area from three unit vectors, no TAN antiderivative."""
    return 2*np.arctan2(abs(np.linalg.det(np.stack([a, b, c]))),
                        1 + a@b + b@c + c@a)


def rectangle_solid_angle(xlo, xhi, ylo, yhi):
    corners = np.array([[xlo, ylo, 1.], [xhi, ylo, 1.],
                        [xhi, yhi, 1.], [xlo, yhi, 1.]])
    corners /= np.linalg.norm(corners, axis=1)[:, None]
    return (triangle_solid_angle(*corners[[0, 1, 2]])
            + triangle_solid_angle(*corners[[0, 2, 3]]))


class TanDisplayTests(unittest.TestCase):
    def test_raw_points_equal_direct_ray_projection_at_wrap_and_high_declination(self):
        for ra0, dec0 in [(359.8, -34.5), (22., 80.), (280., -85.)]:
            dra = np.array([0., .9, -1.7, 180., np.nan])
            ddec = np.array([0., .7, -1.1, 0., 0.])
            rays = absolute_rays(ra0 + dra, dec0 + ddec)
            n, east, north = basis(ra0, dec0)
            denominator = rays @ n
            valid = np.isfinite(dra) & (denominator > 0)
            expected = np.full((len(dra), 2), np.nan)
            expected[valid] = np.rad2deg(np.column_stack([
                rays[valid] @ east, rays[valid] @ north]) / denominator[valid, None])
            xi, eta = display.raw_offsets_to_tangent(
                {'center_RA_deg': ra0, 'center_Dec_deg': dec0}, dra, ddec)
            np.testing.assert_allclose(np.column_stack([xi, eta]), expected,
                                       rtol=0, atol=3e-12, equal_nan=True)

    def test_physical_gnomonic_units_are_plane_degrees_not_axis_angles(self):
        distance = 147.2
        plane = np.array([[-2.1, .3], [0., 0.], [.9, 2.1]])
        n, east, north = basis(41., -34.)
        ratios = np.deg2rad(plane)
        rays = n + ratios[:, 0, None]*east + ratios[:, 1, None]*north
        rays /= np.linalg.norm(rays, axis=1)[:, None]
        denominator = rays @ n
        projected = dict(center_distance_kpc=distance,
                         east_kpc=distance*(rays@east)/denominator,
                         north_kpc=distance*(rays@north)/denominator)
        xi, eta = display.tangent_coordinates(projected)
        np.testing.assert_allclose(np.column_stack([xi, eta]), plane, rtol=0, atol=1e-12)
        self.assertGreater(abs(eta[-1] - np.rad2deg(np.arctan(ratios[-1, 1]))), 5e-4)
        gx, gy = np.meshgrid([-2., 1.], [-1., 0., 2.])
        x, y = display.tangent_kpc_to_display(projected, gx, gy)
        np.testing.assert_allclose(x, np.rad2deg(gx/distance), rtol=0, atol=0)
        np.testing.assert_allclose(y, np.rad2deg(gy/distance), rtol=0, atol=0)

    def test_grid_preserves_ids_los_and_number_weighted_cells(self):
        distance = 1.
        plane = np.array([[-2.09, -2.08], [.010, .011], [.020, .021],
                          [2.09, 2.08], [2.11, 0.], [0., 0.], [np.nan, 1.]])
        ids = np.arange(2**63 + 100, 2**63 + 107, dtype=np.uint64)
        velocity = np.array([-4., 7., 13., 30., 1e5, -1e5, 1.])
        projected = dict(particle_id=ids, center_distance_kpc=distance,
                         east_kpc=np.deg2rad(plane[:, 0])*distance,
                         north_kpc=np.deg2rad(plane[:, 1])*distance,
                         front_mask=np.array([True, True, True, True, True, False, True]))
        frozen = {name: value.copy() if isinstance(value, np.ndarray) else value
                  for name, value in projected.items()}
        frozen_velocity = velocity.copy()
        result = display.tangent_velocity_grid40(projected, velocity, display_v0_kms=3.)
        expected_count = np.zeros((40, 40), dtype=np.int64)
        expected_total = np.zeros((40, 40))
        edges = np.linspace(-2.1, 2.1, 41)
        for row in range(4):
            i, j = np.searchsorted(edges, plane[row], side='right') - 1
            expected_count[i, j] += 1
            expected_total[i, j] += velocity[row] - 3.
        expected_mean = np.divide(expected_total, expected_count,
                                  out=np.full((40, 40), np.nan), where=expected_count > 0)
        np.testing.assert_array_equal(result['count'], expected_count)
        np.testing.assert_allclose(result['mean_HRF_minus_v0'], expected_mean,
                                   rtol=0, atol=0, equal_nan=True)
        np.testing.assert_array_equal(result['display_particle_id'], ids[:4])
        np.testing.assert_array_equal(result['display_mask'], [True]*4 + [False]*3)
        self.assertEqual(result['display_contract_id'], display.TAN_DISPLAY_CONTRACT_ID)
        np.testing.assert_array_equal(result['xi_vertices_deg'][0], edges)
        np.testing.assert_array_equal(result['eta_vertices_deg'][:, 0], edges)
        self.assertEqual(result['display_color_limits_kms'], (-12., 12.))
        for name in projected:
            np.testing.assert_array_equal(projected[name], frozen[name])
        np.testing.assert_array_equal(velocity, frozen_velocity)

    def test_empty_grid_and_invalid_id_velocity_contracts(self):
        projected = dict(particle_id=np.array([], dtype=np.uint64),
                         center_distance_kpc=100., east_kpc=np.array([]),
                         north_kpc=np.array([]), front_mask=np.array([], dtype=bool))
        grid = display.tangent_velocity_grid40(projected, np.array([]), display_v0_kms=0.)
        self.assertEqual(grid['count'].sum(), 0)
        self.assertTrue(np.isnan(grid['mean_HRF_minus_v0']).all())
        for bad in [np.array([1., 2.]), np.array([1, 1], dtype=np.uint64)]:
            with self.assertRaises(ValueError):
                display.tangent_velocity_grid40(dict(projected, particle_id=bad),
                                                 np.zeros(2), display_v0_kms=0.)
        with self.assertRaises(ValueError):
            display.tangent_velocity_grid40(projected, np.array([]), display_v0_kms=np.nan)

    def test_arrow_is_direct_geodesic_tan_endpoint_and_no_cosdec_rescale(self):
        direction = np.array([3., -4.])
        arrow = display.tangent_motion_arrow(direction, length_deg=.65)
        np.testing.assert_array_equal(arrow['origin_xi_eta_deg'], [0., 0.])
        for ra0, dec0 in [(41., -34.), (359.9, 83.)]:
            n, east, north = basis(ra0, dec0)
            unit = direction/np.linalg.norm(direction)
            psi = np.deg2rad(.65)
            ray = np.cos(psi)*n + np.sin(psi)*(unit[0]*east + unit[1]*north)
            expected = np.rad2deg(np.array([ray@east, ray@north])/(ray@n))
            np.testing.assert_allclose(arrow['tip_xi_eta_deg'], expected, rtol=0, atol=5e-13)
            self.assertAlmostEqual(np.rad2deg(np.arccos(ray@n)), .65, places=11)
        for direction, length in [([0., 0.], .65), ([1., np.nan], .65), ([1., 0.], 0.),
                                  ([1., 0.], 90.), ([1., 0.], np.nan)]:
            with self.assertRaises(ValueError):
                display.tangent_motion_arrow(direction, length_deg=length)

    def test_pixel_solid_angle_matches_independent_spherical_polygons(self):
        conversion = (180./np.pi*3600.)**2
        for npix in [1, 7, 120]:
            area = display.tangent_pixel_area_arcsec2(npix)
            self.assertEqual(area.shape, (npix, npix))
            self.assertTrue(np.all(area > 0))
            edges = np.deg2rad(np.linspace(-2.1, 2.1, npix + 1))
            total = rectangle_solid_angle(edges[0], edges[-1], edges[0], edges[-1])
            self.assertAlmostEqual(area.sum()/conversion, total, places=15)
            for i, j in [(0, 0), (npix//2, npix//2), (npix-1, 0)]:
                expected = rectangle_solid_angle(edges[i], edges[i+1], edges[j], edges[j+1])
                self.assertAlmostEqual(area[i, j]/conversion, expected, places=17)
            np.testing.assert_allclose(area, area[::-1, :], rtol=2e-12, atol=1e-8)
        for bad in [0, -2, 1.2, True]:
            with self.assertRaises(ValueError):
                display.tangent_pixel_area_arcsec2(bad)

    def test_legacy_raw_helpers_keep_their_contract(self):
        projected = dict(particle_id=np.array([1], dtype=np.uint64),
                         delta_RA_deg=np.array([.1]), delta_Dec_deg=np.array([.2]),
                         front_mask=np.array([True]), center_RA_deg=41., center_Dec_deg=-34.)
        grid = display.angular_velocity_grid40(projected, np.array([8.]), display_v0_kms=3.)
        self.assertEqual(grid['display_contract_id'], 'heliocentric_delta_radec_2p1_v1')
        self.assertEqual(grid['mean_HRF_minus_v0'][np.isfinite(grid['mean_HRF_minus_v0'])][0], 5.)
        self.assertIn('tip_delta_RA_Dec_deg', display.center_motion_arrow(projected, [1., 0.]))
        self.assertEqual(display.DISPLAY_RELATIVE_DIR,
                         'display/heliocentric_tan_xi_eta_2p1_v1')


if __name__ == '__main__':
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(TanDisplayTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    document = dict(state='PASS' if result.wasSuccessful() else 'FAIL',
                    source_dir=str(SOURCE), source_sha256=sha(SOURCE/'canonical_display.py'),
                    test_source_sha256=sha(Path(__file__)), tests_run=result.testsRun,
                    errors=[str(t) for t, _ in result.errors],
                    failures=[str(t) for t, _ in result.failures],
                    reference='independent ICRS ray dot products and spherical triangle solid angles',
                    snapshot_reads=False, production_mutations=False)
    if ARGS.output_json:
        if ARGS.output_json.exists():
            raise FileExistsError(ARGS.output_json)
        ARGS.output_json.parent.mkdir(parents=True, exist_ok=True)
        ARGS.output_json.write_text(json.dumps(document, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(document))
    raise SystemExit(not result.wasSuccessful())
