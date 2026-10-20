"""Exercise export, real package imports and prepare manifest enforcement."""
from pathlib import Path
import contextlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

WORKFLOW = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(WORKFLOW))
if os.name == 'nt' and 'fcntl' not in sys.modules:
    # These tests exercise prepare and file verification, not POSIX locking.
    sys.modules['fcntl'] = SimpleNamespace(flock=lambda *args: None, LOCK_EX=2, LOCK_NB=4)
import workflow_runtime as runtime
import workflow_github as github


class OrbitPackageFreezingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='dsph-orbit-freezing-',
                                                dir=os.environ.get('DSPH_TEST_SCRATCH'))
        self.root = Path(self.temp.name)
        self.scratch = self.root.parent.resolve()
        self.addCleanup(self.cleanup_fixture)
        self.model = self.root / 'Fornax7004'
        (self.model / 'output').mkdir(parents=True)
        for name, text in {
            'Fornax.ini': '1 0 0 0\n', 'GZWJL.PARAM': 'synthetic parameters\n',
            'mergegalaxy.par': 'MW.ini Fornax.ini\n', 'mergecube.par': 'cube.ini\n',
        }.items():
            (self.model / name).write_text(text, encoding='utf-8')
        (self.model / 'output/snapshot_001.hdf5').write_bytes(b'synthetic snapshot fixture')
        self.code = self.root / 'source'
        self.code.mkdir()
        self.top_names = {'GetInfo.py', 'PlotFig.py', 'data_processing.py', 'canonical_observer.py'}
        for name in self.top_names:
            text = ("CANONICAL_OBSERVER_ID = 'astropy_native_direct_v1'\n"
                    if name == 'canonical_observer.py' else '# unchanged synthetic entry\n')
            (self.code / name).write_text(text, encoding='utf-8')
        cohort = self.model / 'output/analysis/astropy_native_direct_v1/cohort_ids.npz'
        cohort.parent.mkdir(parents=True)
        np.savez_compressed(cohort, particle_id=np.array([1], dtype=np.uint64))
        self.cohort = cohort
        (self.root / 'workflow.json').write_text(json.dumps({'analysis_code': 'source'}), encoding='utf-8')
        self.args = SimpleNamespace(root=self.root, model=self.model.name, name='package-test',
                                    snapshots='output', code=None, dry_run=False)

    def cleanup_fixture(self):
        self.assertEqual(self.root.resolve().parent, self.scratch)
        self.assertTrue(self.root.name.startswith('dsph-orbit-freezing-'))
        self.temp.cleanup()

    def make_package(self):
        package = self.code / 'orbit_calibration_v1'
        files = {
            '__init__.py': 'from .nested.solver import solve\n',
            '__main__.py': 'from . import solve\nprint(solve(21))\n',
            'nested/__init__.py': '# nested import fixture\n',
            'nested/solver.py': 'def solve(value):\n    return value * 2\n',
            'README.md': '# Synthetic canonical package\n',
            'config.example.json': '{"observer_id":"astropy_native_direct_v1"}\n',
            'nested/config_schema.json': '{"type":"object"}\n',
        }
        for name, text in files.items():
            path = package / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text, encoding='utf-8')
        return package, files

    def prepare(self):
        with contextlib.ExitStack() as stack:
            stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
            stack.enter_context(patch.dict(os.environ, {'DSPH_PREEXISTING_ID_PATH': str(self.cohort)}))
            if os.name == 'nt':
                # Only model input links are emulated. Package files are copied
                # by the real prepare implementation and imported below.
                def fixture_link(link, target, target_is_directory=False):
                    source = (link.parent / target).resolve()
                    self.assertTrue(source.is_relative_to(self.root.resolve()))
                    shutil.copy2(source, link)
                stack.enter_context(patch.object(Path, 'symlink_to', fixture_link))
            self.assertEqual(runtime.prepare(self.args), 0)
        directory = runtime.run_dir(self.args)
        return directory, json.loads((directory / 'manifest.json').read_text(encoding='utf-8'))

    def make_link(self, link, target):
        try:
            link.symlink_to(target, target_is_directory=target.is_dir())
        except OSError as error:
            self.skipTest('Host cannot create a real source link: ' + str(error))

    def test_export_and_prepare_exclude_data_caches_hidden_and_output_code(self):
        package, allowed = self.make_package()
        excluded = {
            'private.json', 'measurements.csv', 'stars.npz', 'image.png', 'REPORT.md',
            '.secret.py', '.git/private.py', '__pycache__/cached.py', 'cache/cached.py',
            'output/accidental.py', 'outputs/config_leaked.json', 'results/report.py',
            'data/private.py', 'figures/plot.py', 'logs/trace.py',
        }
        for name in excluded:
            path = package / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('not package source', encoding='utf-8')
        (self.code / 'unlisted_package').mkdir()
        (self.code / 'unlisted_package/private.py').write_text('secret = 1\n', encoding='utf-8')
        selected = {path.relative_to(self.code).as_posix() for path in github.source_files(self.code)}
        expected = self.top_names | {'orbit_calibration_v1/' + name for name in allowed}
        self.assertEqual(selected, expected)
        directory, manifest = self.prepare()
        frozen = directory / 'code/orbit_calibration_v1'
        self.assertEqual({p.relative_to(frozen).as_posix() for p in frozen.rglob('*') if p.is_file()}, set(allowed))
        for name in allowed:
            key = str(Path('code/orbit_calibration_v1') / name)
            self.assertEqual(manifest['code_sha256'][key], runtime.sha(frozen / name))
        for name in self.top_names:
            self.assertEqual((directory / 'code' / name).read_bytes(), (self.code / name).read_bytes())
        self.assertFalse((directory / 'code/unlisted_package').exists())

    def test_nested_import_is_frozen_and_verify_inputs_rejects_tampering(self):
        package, _ = self.make_package()
        directory, manifest = self.prepare()
        # The editable implementation changes after prepare; the frozen import
        # must still execute its original nested dependency.
        (package / 'nested/solver.py').write_text('def solve(value):\n    return -999\n', encoding='utf-8')
        result = subprocess.run(
            [sys.executable, '-I', '-B', '-c',
             'import sys;sys.path.insert(0,sys.argv[1]);from orbit_calibration_v1 import solve;print(solve(21))',
             str(directory / 'code')], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.assertEqual(result.stdout.strip(), '42')
        runtime.verify_inputs(directory, manifest)
        (directory / 'code/orbit_calibration_v1/nested/solver.py').write_text('def solve(value):\n    return 43\n', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'Pinned code changed'):
            runtime.verify_inputs(directory, manifest)

    def test_symlinked_module_is_refused(self):
        package, _ = self.make_package()
        outside = self.root / 'outside.py'
        outside.write_text('private = 1\n', encoding='utf-8')
        self.make_link(package / 'borrowed.py', outside)
        with self.assertRaisesRegex(ValueError, 'Linked analysis package'):
            github.source_files(self.code)
        with self.assertRaisesRegex(ValueError, 'Linked analysis package'):
            self.prepare()
        self.assertFalse((self.root / 'workspace').exists())

    def test_symlinked_package_is_refused(self):
        outside = self.root / 'outside_package'
        outside.mkdir()
        (outside / '__init__.py').write_text('private = 1\n', encoding='utf-8')
        self.make_link(self.code / 'orbit_calibration_v1', outside)
        with self.assertRaisesRegex(ValueError, 'Linked analysis package'):
            github.source_files(self.code)
        with self.assertRaisesRegex(ValueError, 'Linked analysis package'):
            self.prepare()
        self.assertFalse((self.root / 'workspace').exists())

    def test_invalid_or_large_sources_fail_before_creating_a_run(self):
        package, _ = self.make_package()
        oversized = package / 'oversized.py'
        oversized.write_bytes(b'#' + b'x' * (2 * 1024**2))
        with self.assertRaisesRegex(ValueError, 'Unexpected large source file'):
            github.source_files(self.code)
        with self.assertRaisesRegex(ValueError, 'Unexpected large source file'):
            self.prepare()
        self.assertFalse((self.root / 'workspace').exists())
        oversized.unlink()
        (package / 'broken.py').write_text('def broken(:\n', encoding='utf-8')
        with self.assertRaises(SyntaxError):
            self.prepare()
        self.assertFalse((self.root / 'workspace').exists())

    def test_invalid_config_json_is_refused(self):
        package, _ = self.make_package()
        (package / 'config.invalid.json').write_text('{invalid', encoding='utf-8')
        with self.assertRaises(json.JSONDecodeError):
            github.source_files(self.code)
        with self.assertRaises(json.JSONDecodeError):
            self.prepare()
        self.assertFalse((self.root / 'workspace').exists())

    def test_absent_optional_package_keeps_existing_prepare_paths(self):
        directory, manifest = self.prepare()
        self.assertEqual(set(manifest['code_sha256']), {str(Path('code') / name) for name in self.top_names})
        self.assertEqual(manifest['observation']['analysis_relative_dir'], 'output/analysis/astropy_native_direct_v1')
        self.assertTrue((directory / self.model.name / 'output/snapshot_001.hdf5').is_file())
        self.assertFalse((directory / 'code/orbit_calibration_v1').exists())
        runtime.verify_inputs(directory, manifest)


if __name__ == '__main__':
    unittest.main()
