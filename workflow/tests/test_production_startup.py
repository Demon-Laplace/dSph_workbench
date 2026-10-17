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
import workflow_management as management
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
        self.details_patch = patch.object(models, 'simulation_details', return_value={'model': 'Fornax7004'})
        self.details_patch.start()
        self.addCleanup(self.details_patch.stop)

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
        self.assertEqual(result['details_sha256'], runtime.sha(self.job / 'simulation-details.json'))

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

    @unittest.skipUnless(__import__('importlib.util').util.find_spec('numpy'), 'numpy required for mass resolution')
    def test_details_use_actual_component_mass_resolution_and_orbit(self):
        import struct
        self.details_patch.stop()
        for name, header in [('cube.ini', '2 0 0 0'), ('MW.ini', '1 1 1 0'), ('dwarf.ini', '1 0 2 0')]:
            (self.model / name).write_text(header + '\n')
        (self.model / 'mergegalaxy.par').write_text('MW.ini dwarf.ini\n1 2 3 4 5 6 90 0 0 103.42 35.84 -1.8\n')
        (self.model / 'mergecube.par').write_text('cube.ini mergegalaxy.dat\n500 500 500 0 0 0 0 0 0 0 0 0\n')
        runtime.write_json(self.model / 'case.json', {'orbit_reference': 'Fornax5028', 'orientation_reference': 'Fornax5028'})
        header = bytearray(256)
        struct.pack_into('<6I', header, 0, 4, 1, 3, 0, 0, 0)
        blocks = [header, bytes(12 * 8), bytes(12 * 8), bytes(4 * 8),
                  struct.pack('<8f', *[n * 1e-6 for n in [2, 2, 3, 4, 5, 6, 7, 9]])]
        with (self.model / 'mergecube.IC_000').open('wb') as stream:
            for block in blocks:
                stream.write(struct.pack('<I', len(block))); stream.write(block); stream.write(struct.pack('<I', len(block)))
        details = models.simulation_details(self.root, self.model, self.job, {'job_id': '19012'},
            {'InitCondFile': 'mergecube.IC_000', 'UnitMass_in_g': '1.989e43'})
        self.assertEqual(details['total_particles'], 8)
        self.assertEqual(details['components'][1]['model'], 'MW')
        dwarf = details['components'][2]
        self.assertEqual(dwarf['total_particles'], 3)
        self.assertAlmostEqual(dwarf['resolution']['0']['particle_mass_mean_msun'], 40000, places=1)
        self.assertAlmostEqual(dwarf['resolution']['2']['particle_mass_min_msun'], 70000, places=1)
        self.assertAlmostEqual(dwarf['resolution']['2']['particle_mass_max_msun'], 90000, places=1)
        self.assertAlmostEqual(dwarf['resolution']['2']['particle_mass_mean_msun'], 80000, places=1)
        self.assertEqual(details['orbit']['dwarf_rotation_deg'], [103.42, 35.84, -1.8])
        self.assertEqual(details['orbit_reference'], 'Fornax5028')

    def test_remote_monitor_downloads_only_checksum_verified_report(self):
        payload = b'{"model":"Fornax7004","total_particles":8}\n'
        import hashlib
        digest = hashlib.sha256(payload).hexdigest()
        args = SimpleNamespace(root=self.root, arguments=['simulation-startup', 'Fornax7004', '--name', 'orbit'])
        metadata = {'model': 'Fornax7004', 'name': 'orbit', 'sha256': digest}
        transfers = []
        def run(command, **kwargs):
            if command[0] == 'ssh':
                return subprocess.CompletedProcess(command, 0, 'DSPH_SIMULATION_DETAILS ' + json.dumps(metadata) + '\n', '')
            transfers.append(command)
            Path(command[-1]).write_bytes(payload)
            return subprocess.CompletedProcess(command, 0, '', '')
        with patch.object(management, 'remote_command', return_value=['ssh', 'fixture']), \
             patch.object(management, 'ssh_config', return_value=({'root': '/data', 'host': 'server'}, ['ssh'])), \
             patch.object(management.subprocess, 'run', side_effect=run), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(management.remote(args), 0)
            metadata['sha256'] = '0' * 64
            with self.assertRaisesRegex(ValueError, 'checksum mismatch'): management.remote(args)
        self.assertEqual(transfers[0][-2], 'server:/data/workspace/jobs/orbit/simulation-details.json')
        saved = self.root / 'workspace/receipts/production/Fornax7004/orbit/simulation-details.json'
        self.assertEqual(saved.read_bytes(), payload)

    def test_in_job_monitor_waits_without_querying_or_allocating_slurm_resources(self):
        import os
        launcher = self.job / 'subastrofat.sh'
        launcher.write_text((self.model / 'subastrofat.sh').read_text())
        models.prepare_startup_monitor(self.args, self.model, self.job, launcher)
        with patch.dict(os.environ, SLURM_JOB_ID='19012', SLURM_CLUSTER_NAME='astro_fat'), \
             patch.object(models.time, 'monotonic', side_effect=[0, 139, 140]), \
             patch.object(models.time, 'sleep') as sleep, patch.object(ops, 'accounting') as accounting, \
             patch.object(models.subprocess, 'run') as process, contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(models.in_job_startup(self.args), 0)
        sleep.assert_called_once_with(1)
        accounting.assert_not_called()
        process.assert_not_called()
        saved = json.loads((self.job / 'startup.json').read_text())
        self.assertEqual(saved['monitor_location'], 'inside_batch_allocation')
        self.assertEqual(saved['scheduler']['elapsed_seconds'], 140)

    def test_batch_report_is_returned_without_polling_after_monitor_passes(self):
        with patch.object(ops, 'accounting', return_value=dict(state='RUNNING', elapsed_seconds=140, exit_code='0:0')):
            self.monitor()
        runtime.write_json(self.job / 'startup-request.json', {})
        with patch.object(ops, 'accounting') as accounting:
            self.assertEqual(self.monitor(), 0)
            accounting.assert_not_called()

    def test_batch_wrapper_preserves_launcher_and_simulation_exit_status(self):
        import sys
        marker = self.job / 'monitor-complete'
        done = self.job / 'solver-complete'
        launcher = self.job / 'original.sh'
        launcher.write_text('#!/bin/sh\n#SBATCH --nodes=1\n#SBATCH --ntasks=96\n'
            'for n in 1 2 3 4 5 6 7 8 9 10; do\n'
            '  test -f "' + str(marker) + '" && break\n  sleep 0.1\ndone\n'
            'test -f "' + str(marker) + '" || exit 8\n'
            'echo done > "' + str(done) + '"\nexit 7\n')
        original = launcher.read_bytes()
        monitor = self.job / 'monitor.py'
        monitor.write_text('from pathlib import Path\nPath(' + repr(str(marker)) + ').write_text("done")\n')
        wrapper = models.write_monitored_launcher(self.job, launcher, [sys.executable, str(monitor)])
        completed = subprocess.run(['bash', str(wrapper)], timeout=5, capture_output=True)
        self.assertEqual(completed.returncode, 7)
        self.assertTrue(done.exists())
        self.assertEqual((self.job / 'solver-exit-status').read_text().strip(), '7')
        self.assertEqual(launcher.read_bytes(), original)
        self.assertIn('#SBATCH --ntasks=96', wrapper.read_text())
        self.assertNotIn('srun ', wrapper.read_text())
        self.assertNotIn('sbatch ', wrapper.read_text())

    def test_batch_wrapper_stops_sleeping_monitor_when_solver_exits(self):
        import sys
        launcher = self.job / 'original.sh'
        launcher.write_text('#!/bin/sh\nexit 9\n')
        monitor = self.job / 'monitor.py'
        monitor.write_text('import time\ntime.sleep(140)\n')
        wrapper = models.write_monitored_launcher(self.job, launcher, [sys.executable, str(monitor)])
        completed = subprocess.run(['bash', str(wrapper)], timeout=5, capture_output=True)
        self.assertEqual(completed.returncode, 9)


if __name__ == '__main__':
    unittest.main()
