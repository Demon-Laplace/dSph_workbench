import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

PATH = Path(__file__).resolve().parents[2]/'local_tools/ic_generation/plot_ic.py'
SPEC = importlib.util.spec_from_file_location('ic_plots_test', PATH)
PLOT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PLOT)


class ICPlotInputs(unittest.TestCase):
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
