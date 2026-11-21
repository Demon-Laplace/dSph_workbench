import importlib.util
from pathlib import Path
import unittest
import json
import tempfile

import workflow_models as models
import workflow_runtime as runtime

FILE = Path(__file__).resolve().parents[1] / 'solver_patches/build_align16.py'
SPEC = importlib.util.spec_from_file_location('solver_align16_build', FILE)
BUILD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BUILD)


class AllocatorPatchTests(unittest.TestCase):
    def test_all_four_rounding_blocks_and_pool_base_are_fixed(self):
        block = '  if((n % 8) > 0)\n    n = (n / 8 + 1) * 8;\n\n  if(n < 8)\n    n = 8;'
        source = 'if(!(Base = malloc(n)))\n' + '\n'.join([block] * 4)
        patched = BUILD.aligned_source(source)
        self.assertEqual(patched.count('if((n % 16) > 0)'), 4)
        self.assertIn('posix_memalign(&Base, 16, n)', patched)
        self.assertNotIn('(n % 8)', patched)

    def test_unrecognized_and_already_patched_sources_are_not_silently_rebuilt(self):
        with self.assertRaises(ValueError):
            BUILD.aligned_source('some other allocator')


class SolverSelectionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.binary = self.root / 'software/GIZMO_align16'
        self.binary.parent.mkdir()
        self.binary.write_text('#!/bin/sh\nexit 0\n')
        self.binary.chmod(0o755)
        runtime.write_json(self.root / 'build.json', {'state': 'built', 'alignment_bytes': 16,
                           'executable': str(self.binary)})
        runtime.write_json(self.root / 'validation.json', {'state': 'passed'})
        runtime.write_json(self.root / 'workflow.json', {'production_solver': {
            'executable': str(self.binary), 'build_manifest': 'build.json',
            'validation': 'validation.json',
            'replaces': ['${PACKHOME}/soft/periodic/GIZMO', '/obs/periodic/GIZMO']}})
        self.launcher = self.root / 'frozen-launcher.sh'

    def test_only_active_known_solver_is_replaced_resources_unchanged(self):
        original = ('#!/bin/bash\n#SBATCH --nodes=4\n#SBATCH --mem-per-cpu=3800M\n'
                    '#mpirun /obs/periodic/GIZMO debug.PARAM\n'
                    'mpirun "${PACKHOME}/soft/periodic/GIZMO" GZWJL.PARAM ${rerun}\n')
        self.launcher.write_text(original)
        selected = models.select_production_solver(self.root, self.launcher)
        self.assertEqual(self.launcher.read_text(), original.replace(
            '"${PACKHOME}/soft/periodic/GIZMO"', '"' + str(self.binary) + '"'))
        self.assertEqual(selected['alignment_bytes'], 16)

    def test_failed_validation_never_activates_or_modifies_launcher(self):
        self.launcher.write_text('mpirun /obs/periodic/GIZMO GZWJL.PARAM\n')
        runtime.write_json(self.root / 'validation.json', {'state': 'failed'})
        with self.assertRaises(ValueError):
            models.select_production_solver(self.root, self.launcher)
        self.assertEqual(self.launcher.read_text(), 'mpirun /obs/periodic/GIZMO GZWJL.PARAM\n')

    def test_unrelated_isolation_or_custom_solver_is_not_replaced(self):
        self.launcher.write_text('mpirun /obs/isolation/GIZMO GZWJL.PARAM\n')
        self.assertIsNone(models.select_production_solver(self.root, self.launcher))
        self.assertEqual(self.launcher.read_text(), 'mpirun /obs/isolation/GIZMO GZWJL.PARAM\n')

    def test_similar_custom_binary_name_is_not_a_path_match(self):
        self.launcher.write_text('mpirun /obs/periodic/GIZMO-custom GZWJL.PARAM\n')
        self.assertIsNone(models.select_production_solver(self.root, self.launcher))
        self.assertEqual(self.launcher.read_text(), 'mpirun /obs/periodic/GIZMO-custom GZWJL.PARAM\n')
