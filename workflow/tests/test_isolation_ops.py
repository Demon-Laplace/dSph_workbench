import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import workflow_isolation as isolation
import workflow_isolation_ops as ops
import workflow_runtime as runtime


class IsolationOpsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        runtime.write_json(self.root / 'workflow.json', {'isolation_check': {
            'script': 'density.py', 'setup': ':', 'processes': 2,
            'sbatch': {'clusters': 'test', 'ntasks': 1}}})
        (self.root / 'density.py').write_text('fixture')
        self.old = self.make_run('failed', 'failed', '100')
        self.keep = self.make_run('success', 'complete', '101')

    def make_run(self, name, state, job):
        run = isolation.case_dir(self.root, name, 'Fornax7002')
        run.mkdir(parents=True)
        (run / 'output').mkdir()
        for f in ('Fornax.ini', 'parameters.ini'): (run / f).write_text('same IC')
        files = {f: isolation.sha(run / f) for f in ('Fornax.ini', 'parameters.ini')}
        runtime.write_json(run / 'request.json', {'name': name, 'model': 'Fornax7002', 'files': files})
        runtime.write_json(run / 'prepared.json', files)
        runtime.write_json(run / 'status.json', {'state': state, 'job_id': job})
        runtime.write_json(run / 'profile.json', {'sbatch': {'clusters': 'test'}})
        (run / 'solver.log').write_text('EXIT CODE: 11\n')
        return run

    def clean(self, accounting=None, active=None, dry_run=False):
        record = accounting or {'state': 'FAILED', 'elapsed_seconds': 53, 'exit_code': '1:0'}
        with patch.object(ops, 'accounting', return_value=record), patch.object(ops, 'active_directories', return_value=active or []):
            return ops.cleanup(self.root, 'Fornax7002', 'success', dry_run)

    def test_remove_short_failure_keep_replacement_and_small_tombstone(self):
        self.assertIn(str(self.old), self.clean(dry_run=True)['removed'])
        self.assertTrue(self.old.exists())
        self.assertIn(str(self.old), self.clean()['removed'])
        self.assertFalse(self.old.exists())
        self.assertTrue(self.keep.exists())
        receipt = self.root / 'workspace/receipts/isolation-cleaned/Fornax7002/failed.json'
        self.assertEqual(json.loads(receipt.read_text())['accounting']['elapsed_seconds'], 53)
        self.assertEqual(self.clean()['removed'], [])

    def test_replacement_must_pass_startup_before_cleanup(self):
        runtime.write_json(self.keep / 'status.json', {'state': 'running', 'job_id': '101'})
        for elapsed, log in [(100, 'Sync-Point 2, Time: 0.001'), (150, 'initializing')]:
            (self.keep / 'solver.log').write_text(log)
            report = self.clean(accounting={'state': 'RUNNING', 'elapsed_seconds': elapsed, 'exit_code': '0:0'})
            self.assertIn('deferred', report)
            self.assertTrue(self.old.exists())
        (self.keep / 'solver.log').write_text('Sync-Point 2, Time: 0.001')
        with patch.object(ops, 'active_directories', return_value=[self.keep]), patch.object(ops, 'accounting', side_effect=[
            {'state': 'RUNNING', 'elapsed_seconds': 150, 'exit_code': '0:0'},
            {'state': 'FAILED', 'elapsed_seconds': 53, 'exit_code': '1:0'}]):
            ops.cleanup(self.root, 'Fornax7002', 'success')
        self.assertFalse(self.old.exists())

    def test_solver_monitor_cleans_after_threshold_and_propagates_failure(self):
        from unittest.mock import MagicMock
        proc = MagicMock()
        proc.__enter__.return_value = proc
        proc.wait.side_effect = [subprocess.TimeoutExpired('solver', 140), 0]
        with patch.object(isolation.subprocess, 'Popen', return_value=proc), \
             patch.object(isolation, 'cleanup_started_replacement', return_value=True) as clean:
            isolation.monitored_solver(['solver'], self.keep, None, self.root, 'Fornax7002')
        self.assertEqual(proc.wait.call_args_list[0].kwargs['timeout'], 140)
        self.assertEqual(proc.wait.call_args_list[1].kwargs['timeout'], None)
        clean.assert_called_once()
        proc.wait.side_effect = [11]
        with patch.object(isolation.subprocess, 'Popen', return_value=proc):
            with self.assertRaises(subprocess.CalledProcessError):
                isolation.monitored_solver(['solver'], self.keep, None, self.root, 'Fornax7002')

    def test_verified_transfers_removed_but_other_attempts_retained(self):
        import shutil
        transfers = self.root / 'workspace/transfers/isolation'
        transfers.mkdir(parents=True)
        shutil.copytree(self.old, transfers / 'failed-upload')
        shutil.copytree(self.keep, transfers / 'success-upload')
        pending = self.make_run('pending', 'prepared', '102')
        shutil.copytree(pending, transfers / 'pending-upload')
        preview = self.clean(dry_run=True)
        actual = self.clean()
        self.assertEqual(preview['removed'], actual['removed'])
        self.assertEqual(preview['bytes'], actual['bytes'])
        self.assertFalse((transfers / 'failed-upload').exists())
        self.assertFalse((transfers / 'success-upload').exists())
        self.assertTrue((transfers / 'pending-upload').exists())
        self.assertTrue(pending.exists())

    def test_model_label_comes_from_receipt_not_run_date(self):
        import ast
        source = Path(isolation.__file__).parent / 'isolation_python/check_density_profile.py'
        module = ast.parse(source.read_text())
        functions = [n for n in module.body if isinstance(n, ast.FunctionDef) and
                     n.name in ('workflow_model_name', 'detect_model_name', 'detect_model_number')]
        namespace = {'json': json, 'Path': Path, 're': __import__('re')}
        exec(compile(ast.Module(body=functions, type_ignores=[]), str(source), 'exec'), namespace)
        self.assertEqual(namespace['detect_model_name'](self.keep / 'output'), 'Fornax7002')
        self.assertEqual(namespace['detect_model_number'](self.keep / 'output'), '7002')

    def test_active_long_running_or_unknown_failures_are_preserved(self):
        self.clean(active=[self.old / 'output'])
        self.assertTrue(self.old.exists())
        for state, elapsed in [('RUNNING', 50), ('FAILED', 141), ('COMPLETED', 50)]:
            self.clean(accounting={'state': state, 'elapsed_seconds': elapsed, 'exit_code': '1:0'})
            self.assertTrue(self.old.exists())
        (self.old / 'solver.log').write_text('Unknown solver failure')
        self.clean()
        self.assertTrue(self.old.exists())

    def test_other_ics_completed_cases_and_changed_replacement_preserved(self):
        request = json.loads((self.old / 'request.json').read_text())
        request['files']['Fornax.ini'] = 'different'
        runtime.write_json(self.old / 'request.json', request)
        self.clean()
        self.assertTrue(self.old.exists())
        (self.keep / 'Fornax.ini').write_text('changed')
        with self.assertRaises(ValueError): self.clean()
        self.assertTrue(self.old.exists())

    def test_scheduler_failure_aborts_deletion(self):
        with patch.object(ops, 'active_directories', side_effect=subprocess.CalledProcessError(1, 'squeue')):
            with self.assertRaises(subprocess.CalledProcessError):
                ops.cleanup(self.root, 'Fornax7002', 'success')
        self.assertTrue(self.old.exists())

    def test_check_selects_three_snapshots_and_deduplicates_submission(self):
        for n in (0, 300, 600): (self.keep / 'output' / f'snapshot_{n:03d}.hdf5').touch()
        with patch.object(isolation, 'verify_completion'), patch.object(ops, 'active_directories', return_value=[]), \
             patch.object(ops.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, '123;test\n', '')) as submit, \
             patch.object(ops, 'accounting', return_value={'state': 'RUNNING'}):
            first = ops.check(self.root, 'Fornax7002', 'success', 'quick')
            second = ops.check(self.root, 'Fornax7002', 'success', 'quick')
        self.assertEqual(first['selection'], '0,300,600')
        self.assertEqual(second['job_id'], '123')
        submit.assert_called_once()
        self.assertIn('--no-video', submit.call_args.args[0][-1])

    def test_default_check_uses_all_snapshots_and_enables_video(self):
        for n in (0, 300, 600): (self.keep / 'output' / f'snapshot_{n:03d}.hdf5').touch()
        with patch.object(isolation, 'verify_completion'), patch.object(ops, 'active_directories', return_value=[]), \
             patch.object(ops.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, '124;test\n', '')) as submit:
            result = ops.check(self.root, 'Fornax7002', 'success')
        self.assertEqual(result['selection'], 'all')
        self.assertNotIn('--range', submit.call_args.args[0][-1])
        self.assertNotIn('--no-video', submit.call_args.args[0][-1])

    def test_check_refuses_active_directory(self):
        (self.keep / 'output/snapshot_000.hdf5').touch()
        with patch.object(isolation, 'verify_completion'), patch.object(ops, 'active_directories', return_value=[self.keep]):
            with self.assertRaisesRegex(ValueError, 'already using'):
                ops.check(self.root, 'Fornax7002', 'success')


if __name__ == '__main__':
    unittest.main()
