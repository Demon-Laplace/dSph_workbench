import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np

PATH = Path(__file__).resolve().parents[2]/'local_tools/ic_generation/munoz2018.py'
SPEC = importlib.util.spec_from_file_location('munoz2018_test', PATH)
REF = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REF)


class PublishedMunozReference(unittest.TestCase):
    def test_central_and_effective_values_match_the_published_rounded_tables(self):
        distance = 147.
        re = np.deg2rad(18.4/60)*distance
        mu = REF.sersic_mu_v(np.array([0., re]), distance)
        self.assertAlmostEqual(mu[0], 23.59)
        self.assertLess(abs(mu[1]-24.77), .02)
        self.assertLess(abs(re*1000-787.), 1.)

    def test_angular_profile_is_distance_independent(self):
        theta = np.array([0., 1., 18.4, 35.])
        a = REF.sersic_mu_v(np.deg2rad(theta/60)*139.4, 139.4)
        b = REF.sersic_mu_v(np.deg2rad(theta/60)*147., 147.)
        np.testing.assert_allclose(a, b)

    def test_reference_is_never_claimed_to_be_new_observed_bins(self):
        metadata = REF.metadata(139.4)
        self.assertEqual(metadata['data_kind'], 'published_sersic_fit_reconstruction')
        self.assertEqual(metadata['new_independent_radial_measurements'], 0)
        self.assertIn('not a full', metadata['uncertainty_shown'])
        self.assertFalse(metadata['independent_absolute_calibration_of_existing_yang_profile'])

    def test_export_labels_sampled_fit_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'published_fit.csv'
            REF.export_curve(path, 139.4, points=8)
            with path.open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 8)
            self.assertTrue(all(r['data_kind']=='published_fit_evaluation_not_observed_bin' for r in rows))
