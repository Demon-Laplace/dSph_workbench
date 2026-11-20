"""Regression checks for routing dual-beta ICs through the normal build entry."""
import importlib.util
import json
from pathlib import Path
import runpy
import sys
import tempfile
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import numpy as np

SOURCE = Path(__file__).resolve().parents[2] / 'local_tools/ic_generation'


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, SOURCE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TwoBetaBuildTests(unittest.TestCase):
    def test_both_beta_regions_have_unit_independent_constraint_scales(self):
        reader = load('two_beta_scaling', 'schwarzschild_two_beta.py')
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / 'parameters.ini'
            source.write_text('[Component_stellar]\nbeta=two_zone\n'
                'beta_inner=-0.2\nkinemgrid_inner=0,.1,.2,.35,.5\n'
                'beta_outer=0.8\nkinemgrid_outer=1,1.4,2,3\n')
            _, specs = reader.prepare_config(source, SimpleNamespace())
        matrix = np.array([[2., 3., 1., 5.], [4., 1., 2., 3.]])
        for region in specs['component_stellar']:
            self.assertTrue(region['normalize'])
            args = (.006, region['beta'], region['normalize'])
            reference = reader.moment_constraint(matrix, *args)
            rescaled = reader.moment_constraint(matrix * 1e-9, *args)
            np.testing.assert_allclose(reference, rescaled, rtol=1e-14)

    def test_normal_builder_runs_dual_reader_and_exports_numbered_ic(self):
        self.check_builder(fail=False)

    def test_failed_constraints_do_not_export_or_register_ic(self):
        self.check_builder(fail=True)

    def check_builder(self, fail):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / 'IC_Fornax5029.ini'
            source.write_text('[Component_stellar]\nnbody=3\nbeta=two_zone\n'
                              'beta_inner=-0.2\nbeta_outer=0.8\n')
            output = root / 'build'
            catalog = ModuleType('catalog')
            catalog.register = Mock()
            with patch.dict(sys.modules, catalog=catalog):
                builder = load('two_beta_builder', 'build_ic.py')
            agama = ModuleType('agama')
            agama.setRandomSeed = Mock()
            model = object()
            create = Mock(return_value=model)

            def build_components(value):
                self.assertIs(value, model)
                if fail:
                    raise RuntimeError('Two-beta constraints did not pass')
                particles = np.zeros((3, 7))
                particles[:, 6] = .002
                np.savetxt(output / 'model_component_stellar.nbody', particles)

            def worker(command, **kwargs):
                self.assertEqual(Path(command[3]).name, 'schwarzschild_two_beta.py')
                self.assertEqual(Path(kwargs['cwd']), output)
                reader = dict(createModel=create, build_components=build_components)
                with patch.dict(sys.modules, agama=agama), \
                     patch.object(sys, 'argv', ['-c'] + command[3:]), \
                     patch.object(runpy, 'run_path', return_value=reader) as entry:
                    exec(command[2], {})
                entry.assert_called_once_with(command[3])

            original_find = importlib.util.find_spec
            argv = ['build_ic.py', '--ini', str(source), '--output', str(output)]
            with patch.object(sys, 'argv', argv), \
                 patch.object(builder.subprocess, 'check_output', return_value='test-generator-revision\n'), \
                 patch.object(builder.subprocess, 'run', side_effect=worker), \
                 patch.object(builder.importlib.util, 'find_spec',
                     side_effect=lambda name: True if name == 'agama' else original_find(name)):
                if fail:
                    with self.assertRaisesRegex(RuntimeError, 'Two-beta constraints'):
                        builder.main()
                else:
                    self.assertEqual(builder.main(), 0)
            create.assert_called_once_with(str(output / 'parameters.ini'))
            receipt = json.loads((output / 'build.json').read_text())
            self.assertIn('base_solver_sha256', receipt)
            self.assertEqual(receipt['state'], 'failed' if fail else 'complete')
            self.assertEqual((output / 'Fornax5029.ini').exists(), not fail)
            if fail:
                catalog.register.assert_not_called()
            else:
                catalog.register.assert_called_once()
                self.assertEqual(receipt['counts'], [0, 0, 3, 0])
                self.assertEqual((output / 'Fornax.ini').resolve(), output / 'Fornax5029.ini')

    def test_failed_region_diagnostics_survive_component_failure(self):
        reader = load('two_beta_reader', 'schwarzschild_two_beta.py')
        comp = SimpleNamespace(nbody=3, two_beta_validation={'regions': ['failed inner']})
        model = SimpleNamespace(components={'component_stellar': comp},
                                potential=object(), original_solver=object())
        agama = ModuleType('agama')
        agama.writeSnapshot = Mock()
        with patch.dict(sys.modules, agama=agama), \
             patch.object(reader, 'runComponent', side_effect=RuntimeError('invalid beta')), \
             patch.object(reader.Path, 'write_text') as write:
            with self.assertRaisesRegex(RuntimeError, 'invalid beta'):
                reader.build_components(model)
        self.assertEqual(json.loads(write.call_args.args[0]), comp.two_beta_validation)
        agama.writeSnapshot.assert_not_called()


if __name__ == '__main__':
    unittest.main()
