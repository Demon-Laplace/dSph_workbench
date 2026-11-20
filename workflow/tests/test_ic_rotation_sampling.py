"""No dynamics: exercise the actual wrapper's argument selection."""
import contextlib
import configparser
import importlib.util
import io
from pathlib import Path
from types import ModuleType
import unittest
from unittest.mock import Mock, patch
import sys


class SamplingTests(unittest.TestCase):
    def setUp(self):
        path = Path(__file__).resolve().parents[2] / 'local_tools/ic_generation/schwarzschild.py'
        agama = ModuleType('agama')
        agama.setUnits = Mock()
        spec = importlib.util.spec_from_file_location('sampling_wrapper_test', path)
        self.module = importlib.util.module_from_spec(spec)
        with patch.dict(sys.modules, agama=agama):
            spec.loader.exec_module(self.module)
        self.potential = object()

    def options(self, values):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            options = self.module.orbit_sampling_options(values, self.potential)
        return options, output.getvalue()

    def test_final_beta_never_silently_substitutes_for_icbeta(self):
        options, log = self.options({'numorbits': '30000', 'beta': '-0.1', 'ickappa': '.7'})
        self.assertNotIn('beta', options)
        self.assertIn('WARNING', log)
        self.assertIn('Eddington', log)

    def test_explicit_icbeta_activates_axisymmetric_rotation_sampling(self):
        options, log = self.options({'numorbits': '30000', 'icbeta': '0', 'beta': '-.1', 'ickappa': '.7'})
        self.assertEqual(options['beta'], 0.)
        self.assertEqual(options['kappa'], .7)
        self.assertIs(options['potential'], self.potential)
        self.assertIn('axisymmetric_Jeans', log)

    def test_intentional_nonrotating_recipes_remain_nonrotating(self):
        options, log = self.options({'numorbits': '100'})
        self.assertNotIn('kappa', options)
        self.assertNotIn('WARNING', log)
        self.assertIn('Eddington', log)

    def test_8001_only_changes_the_explicit_stellar_seed_beta(self):
        root = Path(__file__).resolve().parents[2]
        a, b = configparser.RawConfigParser(), configparser.RawConfigParser()
        a.read(root / 'IniFile/IC_Fornax5028.ini')
        b.read(root / 'IniFile/IC_Fornax8001.ini')
        original = {s: dict(a.items(s)) for s in a.sections()}
        candidate = {s: dict(b.items(s)) for s in b.sections()}
        self.assertEqual(candidate['Component_stellar'].pop('icbeta'), '0.0')
        self.assertEqual(candidate, original)

    def test_8002_is_an_explicit_thin_disk_not_a_silent_change_to_8001(self):
        root = Path(__file__).resolve().parents[2]
        parent, cold = configparser.RawConfigParser(), configparser.RawConfigParser()
        parent.read(root / 'IniFile/IC_Fornax8001.ini')
        cold.read(root / 'IniFile/IC_Fornax8002.ini')
        self.assertEqual(dict(parent.items('Density_gas')), dict(cold.items('Density_gas')))
        self.assertEqual(dict(parent.items('Component_gas')), dict(cold.items('Component_gas')))
        for field in ('mass', 'sersicindex', 'scaleradius'):
            self.assertEqual(parent.get('Density_stellar', field), cold.get('Density_stellar', field))
        self.assertEqual(cold.getfloat('Density_stellar', 'axisratioz'), .12)
        self.assertEqual(cold.getfloat('Component_stellar', 'axisratioz'), .12)
        self.assertEqual(cold.getfloat('Component_stellar', 'icbeta'), .65)
        self.assertEqual(cold.getfloat('Component_stellar', 'ickappa'), 1.)
        self.assertFalse(cold.has_option('Component_stellar', 'beta'))
        self.assertFalse(cold.has_option('Component_stellar', 'kinemgrid'))
        self.assertEqual(parent.getfloat('Component_stellar', 'beta'), -.1)


if __name__ == '__main__':
    unittest.main()
