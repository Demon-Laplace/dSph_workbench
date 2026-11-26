"""Thread policy across IC entry points, without running an orbit solver."""
from contextlib import redirect_stderr, redirect_stdout
import importlib.util
import io
import json
import os
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import workflow as cli
import workflow_models as models

REPO = Path(__file__).resolve().parents[2]
GENERATOR = REPO / 'local_tools/ic_generation'


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


with patch.dict(sys.modules, catalog=load('thread_policy_catalog', GENERATOR / 'catalog.py')):
    builder = load('thread_policy_builder', GENERATOR / 'build_ic.py')
tool = load('thread_policy_model_tool', REPO / 'model_tool.py')


class StopWorker(BaseException):
    """Stop at the solver boundary, leaving the captured manifest in a temp dir."""


class ICThreadPolicyTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.ini = self.root / 'IC_Fornax9001.ini'
        self.ini.write_text('[Density_gas]\nmass=0.01\n[Component_gas]\nnbody=2\n')

    def test_workflow_cli_preserves_default_and_explicit_threads(self):
        for command, owner, method, arguments in [
                ('build-ic', cli.models, 'build_ic', ['--model', 'Fornax9001']),
                ('isolation', cli.isolation, 'pipeline', ['Fornax9001'])]:
            for extra, expected in [([], None), (['--threads', '2'], 2)]:
                with self.subTest(command=command, extra=extra), \
                     patch.object(sys, 'argv', ['workflow.py', '--root', str(self.root),
                                              command, *arguments, '--name', 'v1', *extra]), \
                     patch.object(owner, method, return_value=0) as handler:
                    self.assertEqual(cli.main(), 0)
                    self.assertEqual(handler.call_args.args[0].threads, expected)

    def test_workflow_forwards_only_explicit_thread_limits(self):
        (self.root / 'workflow.json').write_text(json.dumps({'ic': {'root': '.'}}))
        for threads in (None, 2):
            args = SimpleNamespace(root=self.root, repo=REPO, ini=str(self.ini), model=None,
                                   name='v1', python='python3', threads=threads, seed=42)
            with self.subTest(threads=threads), \
                 patch.object(models.subprocess, 'run', return_value=SimpleNamespace(returncode=0)) as run:
                self.assertEqual(models.build_ic(args), 0)
                self.assert_thread_argument(run.call_args.args[0], threads)

    def assert_thread_argument(self, command, threads):
        if threads is None:
            self.assertNotIn('--threads', command)
        else:
            self.assertEqual(command[command.index('--threads') + 1], str(threads))

    def test_model_tool_forwards_only_explicit_thread_limits(self):
        for threads in (None, 2):
            with self.subTest(threads=threads), \
                 patch.object(tool.subprocess, 'run', return_value=SimpleNamespace(returncode=0)) as run:
                result = tool.build_one(('Fornax9001', self.ini, threads, 42,
                                         self.root, self.root / 'workflow.py'))
                self.assertTrue(result[1])
                self.assert_thread_argument(run.call_args.args[0], threads)

    def test_batch_without_thread_limit_is_sequential(self):
        second = self.root / 'IC_Fornax9002.ini'
        second.write_text(self.ini.read_text())
        args = tool.make_parser().parse_args(['build', str(self.ini), str(second)])
        with patch.object(tool, 'sync_generated_models') as sync, \
             patch.object(tool, 'build_workbench', return_value=(self.root, self.root / 'workflow.py')), \
             patch.object(tool, 'build_one', side_effect=lambda task: (task[0], True, '')) as build, \
             patch.object(tool, 'ProcessPoolExecutor') as pool, redirect_stdout(io.StringIO()):
            self.assertEqual(tool.command_build(args), 0)
            pool.assert_not_called()
            self.assertEqual(build.call_count, 2)
            self.assertTrue(all(call.args[0][2] is None for call in build.call_args_list))
            sync.assert_called_once_with([self.ini, second])

    def test_worker_inherits_environment_unless_limit_is_explicit(self):
        for inherited, extra, expected in [
                ({}, [], {}),
                ({'OMP_NUM_THREADS': '7', 'OPENBLAS_NUM_THREADS': '3'}, [],
                 {'OMP_NUM_THREADS': '7', 'OPENBLAS_NUM_THREADS': '3'}),
                ({'OMP_NUM_THREADS': '7', 'OPENBLAS_NUM_THREADS': '3'}, ['--threads', '2'],
                 {'OMP_NUM_THREADS': '2', 'OPENBLAS_NUM_THREADS': '1'})]:
            with self.subTest(extra=extra), tempfile.TemporaryDirectory() as temporary:
                output = Path(temporary) / 'build'
                captured = {}

                def stop(command, **kwargs):
                    captured.update(kwargs['env'])
                    raise StopWorker()

                environment = {key: value for key, value in os.environ.items()
                               if key not in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS')}
                environment.update(inherited)
                argv = ['build_ic.py', '--ini', str(self.ini), '--output', str(output), *extra]
                with patch.dict(os.environ, environment, clear=True), patch.object(sys, 'argv', argv), \
                     patch.object(builder.importlib.util, 'find_spec', return_value=True), \
                     patch.object(builder.subprocess, 'check_output', return_value='test-revision'), \
                     patch.object(builder.subprocess, 'run', side_effect=stop):
                    with self.assertRaises(StopWorker):
                        builder.main()
                actual = {key: captured[key] for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS')
                          if key in captured}
                self.assertEqual(actual, expected)
                manifest = json.loads((output / 'build.json').read_text())
                self.assertEqual(manifest['threads'], 2 if extra else None)
                self.assertEqual(manifest['thread_policy'], 'explicit' if extra else 'inherited')

    def test_local_id_requires_no_catalog_and_keeps_default_plots(self):
        output = self.root / 'local-build'
        argv = ['build_ic.py', '--ini', str(self.ini), '--output', str(output), '--model', 'X001']
        with patch.object(sys, 'argv', argv), redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as rejected:
                builder.main()
            self.assertEqual(rejected.exception.code, 2)
        self.assertFalse(output.exists())

        def worker(command, **kwargs):
            particles = np.zeros((2, 7))
            particles[:, 6] = 0.005
            np.savetxt(Path(kwargs['cwd']) / 'model_component_gas.nbody', particles)

        with patch.object(sys, 'argv', argv + ['--no-catalog']), \
             patch.object(builder.importlib.util, 'find_spec', return_value=True), \
             patch.object(builder.subprocess, 'check_output', return_value='test-revision'), \
             patch.object(builder.subprocess, 'run', side_effect=worker), \
             patch.object(builder, 'generate_default_plots',
                          return_value={'state': 'complete', 'directory': str(output)}) as plots, \
             patch.object(builder.catalog, 'register') as register, redirect_stdout(io.StringIO()):
            self.assertEqual(builder.main(), 0)
            register.assert_not_called()
            plots.assert_called_once()
        self.assertTrue((output / 'X001.ini').is_file())


if __name__ == '__main__':
    unittest.main()
