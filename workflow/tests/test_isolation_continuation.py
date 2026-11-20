import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import workflow_isolation as isolation
import workflow_isolation_ops as ops
import workflow_runtime as runtime
import workflow_versions as versions


class ContinuationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.model, self.name = 'Fornax8003', 'disk'
        self.run = isolation.case_dir(self.root, self.name, self.model)
        self.run.mkdir(parents=True)
        (self.run / 'output/restartfiles').mkdir(parents=True)
        for rank in range(2):
            (self.run / ('output/restartfiles/restart.' + str(rank))).write_bytes(b'checkpoint')
        runtime.write_json(self.run / 'request.json',
                           {'model': self.model, 'name': self.name, 'profile': 'disk', 'build_id': 'same-IC'})
        runtime.write_json(self.run / 'status.json', {'state': 'complete', 'job_id': '10'})
        runtime.write_json(self.run / 'profile.json',
                           {'solver': ['mpirun', '/readonly/GIZMO', '{parameters}'],
                            'converter': ['must-not-run'], 'sbatch': {'ntasks': 2}})
        (self.run / 'GZWJL.PARAM').write_text(
            'OutputDir ' + str(self.run / 'output') + '\nTimeMax 0.5\nRestartFile restart\nSofteningGas .01\n')
        (self.run / 'parameters.ini').write_text('[Density_stellar]\nmass=.002\n')
        (self.root / 'density.py').write_text('# existing script')
        runtime.write_json(self.root / 'workflow.json', {'isolation_check': {
            'script': 'density.py', 'setup': ':', 'processes': 16,
            'sbatch': {'clusters': 'test', 'ntasks': 1, 'cpus-per-task': 16}}})

    def test_preview_does_not_query_submit_or_modify_parameters(self):
        original = (self.run / 'GZWJL.PARAM').read_text()
        with patch.object(isolation.subprocess, 'run') as process:
            result = isolation.resume(self.root, self.name, self.model, 3.)
        process.assert_not_called()
        self.assertFalse(result['submitted'])
        self.assertEqual((self.run / 'GZWJL.PARAM').read_text(), original)

    def test_restart_preserves_output_skips_converter_and_versions_attempt(self):
        snapshot = self.run / 'output/snapshot_000.hdf5'
        snapshot.write_bytes(b'old snapshot')
        inode = snapshot.stat().st_ino
        answers = [subprocess.CompletedProcess([], 0, '', ''),
                   subprocess.CompletedProcess([], 0, '20;test\n', '')]
        with patch.object(isolation.subprocess, 'run', side_effect=answers):
            result = isolation.resume(self.root, self.name, self.model, 3., submit=True)
        with patch.object(isolation.subprocess, 'run') as process:
            duplicate = isolation.resume(self.root, self.name, self.model, 3., submit=True)
        process.assert_not_called()
        self.assertEqual(duplicate['job_id'], result['job_id'])
        self.assertEqual(result['directory'], str(self.run))
        self.assertEqual(snapshot.stat().st_ino, inode)
        attempt = self.run / '.resume-001'
        self.assertEqual(versions.parameter_values((attempt / 'GZWJL.PARAM').read_text()),
                         {'TimeMax': '3.0', 'RestartFile': 'restart', 'SofteningGas': '.01'})
        with patch.dict(os.environ, SLURM_JOB_ID='20'), \
             patch.object(isolation.subprocess, 'run',
                          return_value=subprocess.CompletedProcess([], 0)) as solver:
            isolation.run_job(self.root, self.name, self.model, attempt)
        self.assertEqual(solver.call_count, 1)
        self.assertEqual(solver.call_args.args[0][-1], '1')
        self.assertFalse((self.run / 'converter.log').exists())
        self.assertEqual((self.run / 'output/restartfiles/restart.0').read_bytes(), b'checkpoint')
        ledger = versions.lookup(isolation.version_directory(self.root, self.model), self.model, self.name)
        self.assertEqual(len(ledger['continuations']), 1)

    def test_active_solver_is_not_resubmitted(self):
        runtime.write_json(self.run / 'status.json', {'state': 'running', 'job_id': '10'})
        with patch.object(isolation.subprocess, 'run') as process, self.assertRaises(ValueError):
            isolation.resume(self.root, self.name, self.model, 3., submit=True)
        process.assert_not_called()

    def test_partial_analysis_freezes_inputs_and_ignores_only_newest_file(self):
        runtime.write_json(self.run / 'status.json', {'state': 'running', 'job_id': '10'})
        for n in range(4):
            (self.run / ('output/snapshot_%03d.hdf5' % n)).touch()
        with patch.object(ops, 'active_directories', return_value=[self.run]), \
             patch.object(ops.subprocess, 'run',
                          return_value=subprocess.CompletedProcess([], 0, '30;test\n', '')), \
             patch.object(isolation, 'sha', side_effect=AssertionError('no routine hashing')):
            result = ops.check(self.root, self.model, self.name)
        self.assertEqual(result['snapshots'], [0, 1, 2])
        self.assertTrue(result['incomplete_simulation'])
        self.assertNotEqual(result['output_directory'], str(self.run / 'density_output'))
        output = Path(result['output_directory'])
        output.mkdir(parents=True, exist_ok=True)
        (output / 'snapshot_000_density_profile.png').touch()
        (output / 'video.mp4').touch()
        (self.run / 'output/snapshot_004.hdf5').touch()
        with patch.object(ops, 'accounting', return_value={'state': 'COMPLETED', 'exit_code': '0:0'}), \
             patch.object(ops, 'active_directories', return_value=[self.run]), \
             patch.object(isolation, 'sha', side_effect=AssertionError('no routine hashing')):
            fetched = ops.results(self.root, self.model, self.name)
        self.assertEqual(fetched['videos'][0]['frames'], 3)
        self.assertTrue(fetched['ready'])
        later = ops.density_plan(self.root, self.model, self.name, 'all')
        self.assertEqual(later['identity']['snapshots'], [0, 1, 2, 3])
        self.assertNotEqual(later['identity']['output_directory'], result['output_directory'])
