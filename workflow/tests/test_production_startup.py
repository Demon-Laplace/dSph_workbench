import contextlib
import io
import json
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import workflow_cleanup as cleanup
import workflow_isolation_ops as ops
import workflow_models as models
import workflow_runtime as runtime


class ProductionStartupTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.model = self.root / 'Fornax7004'
        self.output = self.model / 'output'
        self.output.mkdir(parents=True)
        self.parameter = self.model / 'GZWJL.PARAM'
        self.parameter.write_text('InitCondFile mergecube.IC_000\nOutputDir output\n')
        (self.model / 'mergecube.IC_000').write_bytes(b'fixture')
        (self.model / 'subastrofat.sh').write_text('#!/bin/sh\n')
        self.job = self.root / 'workspace/jobs/orbit'
        self.job.mkdir(parents=True)
        runtime.write_json(self.job / 'submission.json', dict(submitted=True, model=str(self.model),
            stdout='19012;astro_fat\n', parameter_sha256=runtime.sha(self.parameter)))
        self.args = SimpleNamespace(root=self.root, model='Fornax7004', name='orbit', timeout=60)
        self.info = self.output / 'info.txt'
        self.info.write_text('Sync-Point 8, Time: 0.001, Nf = 10\n')

    def monitor(self):
        with contextlib.redirect_stdout(io.StringIO()):
            return models.simulation_startup(self.args)

    def test_waits_for_140_running_seconds_then_exits_without_submission(self):
        rows = [dict(state='RUNNING', elapsed_seconds=n, exit_code='0:0') for n in [139, 140]]
        with patch.object(ops, 'accounting', side_effect=rows) as accounting, \
             patch.object(models.time, 'sleep') as sleep, patch.object(models.subprocess, 'run') as process:
            self.assertEqual(self.monitor(), 0)
        self.assertEqual(accounting.call_count, 2)
        accounting.assert_called_with('19012', 'astro_fat')
        sleep.assert_called_once()
        process.assert_not_called()
        result = json.loads((self.job / 'startup.json').read_text())
        self.assertEqual(result['state'], 'startup_verified')
        self.assertEqual(result['progress']['sync_point'], 8)

    def test_initial_snapshot_alone_does_not_pass_and_timeout_does_not_cancel(self):
        self.info.write_text('Sync-Point 0, Time: 0, Nf = 100\n')
        with patch.object(ops, 'accounting', return_value=dict(state='RUNNING', elapsed_seconds=900, exit_code='0:0')), \
             patch.object(models.time, 'monotonic', side_effect=[0, 61]), \
             patch.object(models.time, 'sleep') as sleep, patch.object(models.subprocess, 'run') as process:
            with self.assertRaisesRegex(ValueError, 'not cancelled'):
                self.monitor()
        sleep.assert_not_called()
        process.assert_not_called()
        self.assertEqual(json.loads((self.job / 'startup.json').read_text())['state'], 'timed_out')

    def test_failed_or_short_completed_job_cannot_be_accepted(self):
        for state, code, elapsed in [('FAILED', '1:0', 150), ('OUT_OF_MEMORY', '0:125', 150),
                                     ('CANCELLED', '0:0', 150), ('COMPLETED', '0:0', 139)]:
            with self.subTest(state=state), patch.object(ops, 'accounting',
                    return_value=dict(state=state, elapsed_seconds=elapsed, exit_code=code)), \
                 patch.object(models.time, 'sleep') as sleep:
                with self.assertRaises(ValueError): self.monitor()
                sleep.assert_not_called()
                self.assertEqual(json.loads((self.job / 'startup.json').read_text())['state'], 'failed')

    def test_missing_accounting_is_retried_but_changed_parameters_are_rejected(self):
        with patch.object(ops, 'accounting', side_effect=[
                ValueError('Slurm accounting is unavailable or ambiguous; retry later'),
                dict(state='RUNNING', elapsed_seconds=140, exit_code='0:0')]), \
             patch.object(models.time, 'sleep') as sleep:
            self.assertEqual(self.monitor(), 0)
            sleep.assert_called_once()
        self.parameter.write_text(self.parameter.read_text() + 'PartAllocFactor 6.5\n')
        with patch.object(ops, 'accounting') as accounting:
            with self.assertRaisesRegex(ValueError, 'changed since submission'): self.monitor()
            accounting.assert_not_called()

    def test_pending_job_with_old_progress_does_not_pass(self):
        with patch.object(ops, 'accounting', return_value=dict(state='PENDING', elapsed_seconds=900, exit_code='0:0')), \
             patch.object(models.time, 'monotonic', side_effect=[0, 61]):
            with self.assertRaisesRegex(ValueError, 'timed out'): self.monitor()

    def test_successful_submit_checks_startup_by_default_and_can_return_immediately(self):
        for no_wait in [False, True]:
            args = SimpleNamespace(root=self.root, model='Fornax7004', name='submit-' + str(no_wait),
                launcher='subastrofat.sh', submit=True, no_wait_startup=no_wait)
            with patch.object(cleanup, 'active_jobs', return_value=[]), \
                 patch.object(models.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, '19012;astro_fat\n', '')), \
                 patch.object(models, 'simulation_startup', return_value=0) as monitor, \
                 contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(models.simulate(args), 0)
                self.assertEqual(monitor.call_count, 0 if no_wait else 1)

    def test_cli_can_monitor_existing_job_without_resubmission(self):
        import workflow
        with patch('sys.argv', ['workflow.py', '--root', str(self.root), 'simulation-startup',
                               'Fornax7004', '--name', 'orbit']), \
             patch.object(models, 'simulation_startup', return_value=0) as monitor:
            self.assertEqual(workflow.main(), 0)
            self.assertEqual(monitor.call_args.args[0].name, 'orbit')


if __name__ == '__main__':
    unittest.main()
