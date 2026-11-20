import importlib.util
from pathlib import Path
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

    def test_density_map_is_mass_weighted_not_particle_count_weighted(self):
        image = PLOT.surface_image(np.array([-.5, .5]), np.array([0., 0.]),
                                   np.array([1., 3.]), extent=1., bins=2)
        self.assertAlmostEqual(float(image.sum()), 4e4)
        self.assertAlmostEqual(float(image.max()/image[image > 0].min()), 3.)
