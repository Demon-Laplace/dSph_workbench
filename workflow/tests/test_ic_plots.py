import importlib.util
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import re
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

PATH = Path(__file__).resolve().parents[2]/'local_tools/ic_generation/plot_ic.py'
SPEC = importlib.util.spec_from_file_location('ic_plots_test', PATH)
PLOT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PLOT)


class ICPlotInputs(unittest.TestCase):
    def test_velocity_panels_plot_every_particle_without_bin_averaging(self):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        populations, summaries = {}, {}
        for name, count in (('Stars', 53), ('Gas', 7)):
            x = np.column_stack((np.arange(count, dtype=float), np.zeros(count), np.ones(count)))
            v = np.tile([3., 4., 0.], (count, 1))
            mass = np.ones(count)
            populations[name] = [x, v, mass]
            summaries[name] = {'All': PLOT.kinematic_summary(x, v, mass)}
        fig = PLOT._kinematics_page(populations, 'Fornax9001', summaries, plt)
        try:
            speed_ax, rotation_ax = fig.axes[1:3]
            self.assertEqual(len(speed_ax.lines), 0)
            self.assertEqual([len(points.get_offsets()) for points in speed_ax.collections], [53, 7])
            self.assertEqual([len(points.get_offsets()) for points in rotation_ax.collections], [52, 6])
            for points, name in zip(speed_ax.collections, ('Stars', 'Gas')):
                np.testing.assert_allclose(points.get_offsets()[:, 0], np.linalg.norm(populations[name][0], axis=1))
                np.testing.assert_allclose(points.get_offsets()[:, 1], 5.)
                self.assertTrue(points.get_rasterized())
            self.assertGreater(speed_ax.get_xlim()[1], 52.)
        finally:
            plt.close(fig)

    def test_rotation_summary_separates_speed_from_net_rotation_and_uses_mass_weights(self):
        x = np.array([[1., 0., 0.], [1., 0., 0.]])
        v = np.array([[0., 10., 0.], [0., -10., 0.]])
        equal = PLOT.kinematic_summary(x, v, np.ones(2))
        self.assertEqual(equal['mean_speed_3d_kms'], 10.)
        self.assertEqual(equal['mean_vphi_kms'], 0.)
        self.assertEqual(equal['prograde_mass_fraction'], .5)
        weighted = PLOT.kinematic_summary(x, v, np.array([1., 3.]))
        self.assertEqual(weighted['mean_vphi_kms'], -5.)
        self.assertEqual(weighted['prograde_mass_fraction'], .25)
        self.assertEqual(weighted['retrograde_mass_fraction'], .75)

    def test_rotation_summary_keeps_native_axis_and_handles_axis_particles(self):
        x = np.array([[1., 0., 0.], [0., 0., 1.]])
        v = np.array([[0., -10., 0.], [0., 3., 0.]])
        summary = PLOT.kinematic_summary(x, v, np.ones(2))
        self.assertEqual(summary['mean_vphi_kms'], -10.)
        self.assertEqual(summary['prograde_mass_fraction'], 0.)
        self.assertEqual(summary['retrograde_mass_fraction'], .5)
        self.assertEqual(summary['zero_Lz_mass_fraction'], .5)
        empty = PLOT.kinematic_summary(np.empty((0, 3)), np.empty((0, 3)), np.empty(0))
        self.assertEqual(empty['N'], 0)
        self.assertIsNone(empty['mean_vphi_kms'])

    def test_radial_speed_profile_uses_3d_shells(self):
        x = np.tile([0., 0., 1.5], (20, 1))
        v = np.tile([3., 4., 0.], (20, 1))
        values = PLOT.radial_kinematics(x, v, np.ones(20), np.array([.5, 1., 2.]))
        self.assertEqual(values.shape, (1, 7))
        self.assertAlmostEqual(values[0, 0], np.sqrt(2.))
        self.assertEqual(values[0, 1], 5.)
        self.assertTrue(np.isnan(values[0, 2]))

    def test_report_writes_one_three_page_pdf_and_failed_refresh_preserves_it(self):
        munoz_spec = importlib.util.spec_from_file_location('pdf_test_munoz', PATH.with_name('munoz2018.py'))
        munoz = importlib.util.module_from_spec(munoz_spec)
        munoz_spec.loader.exec_module(munoz)
        rng = np.random.default_rng(17)
        n = 256
        data = np.column_stack((rng.normal(size=(n, 3))*.3, rng.normal(size=(n, 3))*5,
                                np.full(n, 1e-6), np.zeros(n), np.r_[np.zeros(n//2), np.full(n//2, 2)]))
        original = data.copy()
        with tempfile.TemporaryDirectory() as temporary, patch.dict(sys.modules, munoz2018=munoz):
            directory = Path(temporary)
            result = PLOT.generate(data, 'Fornax9001', directory)
            self.assertEqual(result['files'], ['Fornax9001_initial_diagnostics.pdf'])
            self.assertEqual(result['page_count'], 3)
            report = directory / result['files'][0]
            content = report.read_bytes()
            self.assertTrue(content.startswith(b'%PDF-'))
            self.assertEqual(len(re.findall(rb'/Type /Page\b', content)), 3)
            self.assertEqual(list(directory.glob('*.png')), [])
            np.testing.assert_array_equal(data, original)
            profiles = PLOT.generate(data, 'Fornax9001', directory, profile_only=True)
            self.assertEqual(profiles['files'], ['Fornax9001_initial_profiles.pdf'])
            self.assertEqual(profiles['page_count'], 1)
            self.assertEqual(len(re.findall(rb'/Type /Page\b', (directory/profiles['files'][0]).read_bytes())), 1)
            self.assertEqual(report.read_bytes(), content)
            with self.assertRaises(FileExistsError):
                PLOT.generate(data, 'Fornax9001', directory)
            with patch.object(PLOT, '_profiles_page', side_effect=RuntimeError('render failed')):
                with self.assertRaisesRegex(RuntimeError, 'render failed'):
                    PLOT.generate(data, 'Fornax9001', directory, overwrite=True)
            self.assertEqual(report.read_bytes(), content)
            self.assertEqual(list(directory.glob('.*.pdf')), [])

    def test_in_place_refresh_reuses_selected_build_plot_prefix_and_preserves_build_record(self):
        with tempfile.TemporaryDirectory() as temporary:
            build = Path(temporary)
            receipt = {'build_id': 'second-build', 'source_ini': 'IC_X001.ini',
                       'plots': {'prefix': 'X001_second-build'}}
            record = build / 'build.json'
            record.write_text(json.dumps(receipt))
            original = record.read_bytes()
            result = {'state': 'complete', 'directory': str(build),
                      'prefix': 'X001_second-build', 'files': ['X001_second-build_initial_diagnostics.pdf']}
            legacy = build / 'X001_second-build_kinematics.png'
            legacy.write_bytes(b'old derived figure')
            unrelated = build / 'X001_kinematics.png'
            unrelated.write_bytes(b'another build figure')
            with patch.object(sys, 'argv', ['plot_ic.py', '--build', str(build), '--overwrite']), \
                 patch.object(PLOT.np, 'loadtxt', return_value=np.ones((2, 7))), \
                 patch.object(PLOT, 'generate', return_value=result) as generate, \
                 redirect_stdout(io.StringIO()):
                PLOT.main()
            self.assertTrue(generate.call_args.kwargs['overwrite'])
            self.assertEqual(generate.call_args.kwargs['prefix'], 'X001_second-build')
            self.assertEqual(record.read_bytes(), original)
            self.assertFalse(legacy.exists())
            self.assertTrue(unrelated.exists())
            refresh = json.loads((build / 'initial_figures.json').read_text())
            self.assertEqual(refresh['source_build_id'], 'second-build')
            self.assertFalse(refresh['particles_regenerated'])

    def test_volume_density_uses_shell_volume_and_mass_units(self):
        x = np.array([[.5, 0., 0.], [.75, 0., 0.], [1.5, 0., 0.]])
        mass = np.array([1., 3., 2.])*1e-10
        edges = np.array([0., 1., 2., 3.])
        radius, rho, error, counts = PLOT.spherical_density_profile(x, mass, edges)
        volumes = 4*np.pi/3*np.diff(edges**3)*1e9
        np.testing.assert_allclose(rho[:2]*volumes[:2], [4., 2.])
        np.testing.assert_allclose(error[:2]*volumes[:2], [np.sqrt(10.), 2.])
        np.testing.assert_array_equal(counts, [2, 1, 0])
        self.assertTrue(np.isnan(rho[2]))
        self.assertTrue(np.isnan(error[2]))

    def test_volume_density_uses_3d_radius_not_cylindrical_radius(self):
        x = np.array([[0., 0., 1.5]])
        _, rho, _, counts = PLOT.spherical_density_profile(x, np.array([1e-10]), np.array([0., 1., 2.]))
        np.testing.assert_array_equal(counts, [0, 1])
        self.assertTrue(np.isnan(rho[0]))
        self.assertGreater(rho[1], 0.)

    def test_cold_rotating_ring_distinguishes_cylindrical_and_xyz_scatter(self):
        phi = np.linspace(0., 2*np.pi, 40, endpoint=False)
        x = np.column_stack((np.cos(phi), np.sin(phi), np.zeros_like(phi)))
        v = np.column_stack((-10*np.sin(phi), 10*np.cos(phi), np.zeros_like(phi)))
        values = PLOT.curve(x, v, np.ones(len(phi)), np.array([.5, 1.5]))
        self.assertEqual(values.shape, (1, 9))
        self.assertAlmostEqual(values[0, 1], 10.)
        np.testing.assert_allclose(values[0, 2:5], 0., atol=1e-13)
        np.testing.assert_allclose(values[0, 6:9], [10/np.sqrt(2), 10/np.sqrt(2), 0.])

    def test_xyz_sigmas_are_mass_weighted_and_bulk_velocity_invariant(self):
        x = np.tile([1., 0., 0.], (20, 1))
        v = np.tile([[2., 5., 1.], [6., 9., 3.]], (10, 1))
        mass = np.tile([1., 3.], 10)
        values = PLOT.curve(x, v, mass, np.array([.5, 1.5]))
        shifted = PLOT.curve(x, v+[100., -200., 50.], mass, np.array([.5, 1.5]))
        np.testing.assert_allclose(values[0, 6:9], [np.sqrt(3), np.sqrt(3), np.sqrt(.75)])
        np.testing.assert_allclose(shifted[0, 6:9], values[0, 6:9])
        self.assertAlmostEqual(values[0, 4], values[0, 8])

    def test_dispersion_bins_below_particle_threshold_are_omitted(self):
        x = np.tile([1., 0., 0.], (19, 1))
        values = PLOT.curve(x, np.zeros_like(x), np.ones(len(x)), np.array([.5, 1.5]))
        self.assertEqual(values.shape, (0, 9))

    def test_standard_build_uses_shallow_model_isolation_directory(self):
        self.assertEqual(PLOT.default_directory(Path('/data/Fornax8002/ic/cold-v1')),
                         Path('/data/Fornax8002/isolation'))
        self.assertEqual(PLOT.default_directory(Path('/data/custom-build')),
                         Path('/data/custom-build'))

    def test_flat_build_uses_model_root_and_keeps_project_observation_config(self):
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            model = project / 'Fornax8003'
            build = project / 'workspace/ic-builds/Fornax8003/v1'
            build.mkdir(parents=True)
            (build / 'build.json').write_text(json.dumps(
                {'local_output_directory': str(model), 'local_model_directory': str(model)}))
            config = {'ic_plotting': {'observation_profile': 'profile.csv'}}
            (project / 'workflow.json').write_text(json.dumps(config))
            self.assertEqual(PLOT.default_directory(build), model)
            self.assertEqual(PLOT.project_settings(build), (project, config))

    def test_density_map_is_mass_weighted_not_particle_count_weighted(self):
        image = PLOT.surface_image(np.array([-.5, .5]), np.array([0., 0.]),
                                   np.array([1., 3.]), extent=1., bins=2)
        self.assertAlmostEqual(float(image.sum()), 4e4)
        self.assertAlmostEqual(float(image.max()/image[image > 0].min()), 3.)

    def test_fornax_reference_preserves_prepared_magnitudes_and_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root/'Fornax_surface_brightness_profile.csv'
            original = ('r_mid_kpc,mu_v_mag_arcsec2,mu_v_err_total,calibration_name\n'
                        '0.1,23.5,0.2,existing\n0.2,24.1,0.3,existing\n')
            path.write_text(original)
            data, metadata = PLOT.load_surface_brightness_observation(observation_profile=path)
            np.testing.assert_allclose(data['radius_kpc'], [.1, .2])
            np.testing.assert_allclose(data['mu_v'], [23.5, 24.1])
            np.testing.assert_allclose(data['mu_v_error'], [.2, .3])
            self.assertEqual(path.read_text(), original)
            self.assertEqual(metadata['state'], 'loaded')
            self.assertFalse(metadata['include_WLM'])
            self.assertFalse(metadata['observation_recalibrated'])

    def test_fornax_reference_uses_configured_profile_not_working_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            project, observations = root/'analysis', root/'observations'
            build = project/'Fornax8002/ic/cold-v1'
            build.mkdir(parents=True)
            observations.mkdir()
            path = observations/'profile.csv'
            path.write_text('r_mid_kpc,mu_v_mag_arcsec2,mu_v_err_total\n0.2,24,0.1\n')
            (project/'workflow.json').write_text(json.dumps(
                {'ic_plotting': {'observation_profile': '../observations/profile.csv'}}))
            np.save(project/'radiusbin.npy', [999.])
            np.save(project/'flux.npy', [999.])
            data, metadata = PLOT.load_surface_brightness_observation(build=build)
            np.testing.assert_allclose(data['mu_v'], [24.])
            self.assertEqual(metadata['source_file'], str(path))
            self.assertFalse(metadata['include_WLM'])

    def test_stellar_light_conversion_matches_existing_apparent_magnitude_formula(self):
        surface = np.array([1., 10., 50.])
        for distance_kpc in (100., 139.4, 139.6):
            # One physical pc^2: apparent flux and angular-area distance factors cancel.
            magnitude = (-2.5*np.log10(surface/2.6)+5*np.log10(distance_kpc*1000)
                         +4.83-5-.02)
            area_arcsec2 = (206265/(distance_kpc*1000))**2
            expected = magnitude+2.5*np.log10(area_arcsec2)
            np.testing.assert_allclose(PLOT.surface_mass_to_mu_v(surface, 2.6), expected)

    def test_stellar_light_uses_annular_area_and_leaves_empty_bins_undefined(self):
        x = np.array([[.5, 0., 0.]])
        mass = np.array([np.pi*1e6/1e10])
        radius, mu = PLOT.stellar_brightness_profile(
            x, mass, np.array([0., 1.]), np.array([1., 2.]), 2.6)
        self.assertAlmostEqual(mu[0], float(PLOT.surface_mass_to_mu_v(np.array([1.]), 2.6)[0]))
        self.assertTrue(np.isnan(mu[1]))

    def test_photometry_reads_active_source_without_executing_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            project = Path(tmp)
            build = project/'Fornax8002/ic/cold-v1'
            build.mkdir(parents=True)
            active = project/'active/repo'
            active.mkdir(parents=True)
            (active/'variable.py').write_text('sim_plot_mass_to_light_ratio=2.6\nd_today=139.4\n'
                                            'raise RuntimeError("must not execute")\n')
            (project/'workflow.json').write_text(json.dumps({'analysis_code': 'active/repo'}))
            settings = PLOT.photometric_settings(build=build)
            self.assertEqual(settings['mass_to_light_V'], 2.6)
            self.assertEqual(settings['mass_to_light_source'], str(active/'variable.py'))
            self.assertEqual(settings['reference_display_distance_kpc'], 139.4)

    def test_angular_observation_radii_use_common_distance_without_changing_input(self):
        data = {'radius_kpc': np.array([1.]), 'inner_kpc': np.array([.9]),
                'outer_kpc': np.array([1.1]), 'radius_deg': np.array([.4]),
                'inner_deg': np.array([.3]), 'outer_deg': np.array([.5]),
                'mu_v': np.array([24.]), 'mu_v_error': np.array([.2])}
        displayed = PLOT.observation_at_distance(data, 139.4)
        np.testing.assert_allclose(displayed['radius_kpc'], np.deg2rad([.4])*139.4)
        np.testing.assert_allclose(data['radius_kpc'], [1.])
        np.testing.assert_allclose(displayed['mu_v'], data['mu_v'])
