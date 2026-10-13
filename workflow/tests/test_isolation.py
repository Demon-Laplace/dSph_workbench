import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import workflow_isolation as isolation
import workflow_runtime as runtime


class IsolationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.stage = self.root / 'transfer'
        self.stage.mkdir()
        template = self.root / 'template'
        template.mkdir()
        self.parameters = 'InitCondFile old\nOutputDir old-output\nResubmitOn 1\nTimeMax 1\nSofteningGas 0.1\n'
        (template / 'GZWJL.PARAM').write_text(self.parameters)
        config = {'isolation_profiles': {'default': {'template': 'template',
                  'converter': [sys.executable, '-c', 'from pathlib import Path; Path("{binary}").write_bytes(b"binary")'],
                  'solver': [sys.executable, '-c', 'print("solver fixture")'],
                  'sbatch': {'ntasks': 2, 'time': '00:10:00'}}}}
        runtime.write_json(self.root / 'workflow.json', config)
        for name in ('Fornax.ini', 'parameters.ini'):
            (self.stage / name).write_text('fixture ' + name)
        for name in ('workflow_isolation.py', 'workflow_runtime.py'):
            shutil.copy2(Path(isolation.__file__).with_name(name), self.stage / name)
        self.request = {'name': 'Fornax7003-test', 'model': 'Fornax7003', 'profile': 'default',
                        'build_id': 'fixture', 'files': {n: isolation.sha(self.stage / n)
                        for n in ('Fornax.ini', 'parameters.ini')}}
        runtime.write_json(self.stage / 'request.json', self.request)
        self.run = isolation.case_dir(self.root, self.request['name'], self.request['model'])

    def receive(self, submit=False):
        with patch.object(isolation.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, '12345\n', '')) as sbatch:
            result = isolation.receive(self.root, self.stage, submit)
            return result, sbatch

    def test_preflight_preserves_physics_and_does_not_submit(self):
        result, sbatch = self.receive()
        self.assertEqual(result['state'], 'prepared')
        self.assertIn('--test-only', sbatch.call_args.args[0])
        params = (self.run / 'GZWJL.PARAM').read_text()
        self.assertIn('TimeMax 1\nSofteningGas 0.1', params)
        self.assertIn('ResubmitOn 0', params)
        self.assertEqual((self.root / 'template/GZWJL.PARAM').read_text(), self.parameters)

    def test_repeat_submission_returns_original_job_without_sbatch(self):
        result, _ = self.receive(submit=True)
        self.assertEqual(result['job_id'], '12345')
        again, sbatch = self.receive(submit=True)
        self.assertEqual(again['job_id'], '12345')
        sbatch.assert_not_called()

    def test_changed_inputs_or_settings_cannot_reuse_run(self):
        self.receive()
        self.request['profile'] = 'other'
        runtime.write_json(self.stage / 'request.json', self.request)
        with self.assertRaisesRegex(ValueError, 'different ICs'):
            isolation.receive(self.root, self.stage)
        self.request['profile'] = 'default'
        runtime.write_json(self.stage / 'request.json', self.request)
        (self.stage / 'Fornax.ini').write_text('bad upload')
        with self.assertRaisesRegex(ValueError, 'checksum'):
            isolation.receive(self.root, self.stage)

    def test_prepared_parameter_edit_and_existing_data_refused(self):
        self.receive()
        (self.run / 'GZWJL.PARAM').write_text('changed')
        with self.assertRaisesRegex(ValueError, 'Prepared case changed'):
            isolation.receive(self.root, self.stage, True)

    def test_uncertain_submission_never_blindly_retries(self):
        with patch.object(isolation.subprocess, 'run', side_effect=OSError('connection lost')):
            with self.assertRaises(OSError):
                isolation.receive(self.root, self.stage, True)
        with self.assertRaisesRegex(ValueError, 'uncertain'):
            isolation.receive(self.root, self.stage, True)

    def test_conversion_and_solver_run_once_in_case_directory(self):
        self.receive(submit=True)
        with patch.dict(os.environ, SLURM_JOB_ID='12345'), patch.object(isolation, 'verify_completion', return_value={'final_snapshot_time': 1}):
            isolation.run_job(self.root, self.run.name)
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                isolation.run_job(self.root, self.run.name)
        self.assertEqual((self.run / 'Fornax.IC_000').read_bytes(), b'binary')
        self.assertEqual(json.loads((self.run / 'status.json').read_text())['state'], 'complete')

    def test_solver_without_final_snapshot_is_failed(self):
        self.receive(submit=True)
        with patch.dict(os.environ, SLURM_JOB_ID='12345'), patch.object(isolation, 'verify_completion', side_effect=ValueError('TimeMax not reached')):
            with self.assertRaisesRegex(ValueError, 'TimeMax'):
                isolation.run_job(self.root, self.run.name)
        self.assertEqual(json.loads((self.run / 'status.json').read_text())['state'], 'failed')

    def test_symlink_escape_cannot_create_server_case(self):
        outside = self.root / 'outside'
        outside.mkdir()
        data = self.root / 'data'
        data.mkdir()
        (data / 'Fornax7003').symlink_to(outside, target_is_directory=True)
        with self.assertRaises(ValueError):
            isolation.case_dir(data, self.request['name'], self.request['model'])
        self.assertFalse((outside / 'models').exists())

    def test_local_transfer_and_remote_receiver_end_to_end_without_network(self):
        local = self.root / 'local'
        local.mkdir()
        ic = local / 'generated/existing'
        ic.mkdir(parents=True)
        for name in ('parameters.ini', 'Fornax.ini'):
            (ic / name).write_text('transport fixture')
        runtime.write_json(ic / 'build.json', {'state': 'complete', 'source_ini': 'IC_Fornax7003.ini',
                           'seed': 42, 'build_id': 'local-fixture', 'ini_sha256': isolation.sha(ic / 'parameters.ini'),
                           'particle_sha256': isolation.sha(ic / 'Fornax.ini'), 'catalog_sync': {'state': 'pushed'}})
        init = self.root / 'modules.sh'
        init.write_text('module() { :; }\n')
        runtime.write_json(local / 'workflow.json', {'ic': {'root': '.'}, 'remote': {'host': 'fixture-host',
                           'root': str(self.root), 'module_init': str(init), 'modules': []}})
        binaries = self.root / 'bin'
        binaries.mkdir()
        scripts = {
            'ssh': 'import subprocess,sys; sys.exit(subprocess.run(["bash","-c",sys.argv[-1]]).returncode)',
            'rsync': 'import shutil,sys; shutil.copytree(sys.argv[-2],sys.argv[-1].split(":",1)[1],dirs_exist_ok=True)',
            'sbatch': 'from pathlib import Path; p=Path(' + repr(str(self.root / 'sbatch-count')) + '); p.write_text(p.read_text()+"1" if p.exists() else "1"); print("45678")',
        }
        for name, script in scripts.items():
            path = binaries / name
            path.write_text('#!' + sys.executable + '\n' + script + '\n')
            path.chmod(0o755)
        args = SimpleNamespace(root=local, repo=Path(isolation.__file__).parents[1],
                               model='Fornax7003', name='transport', ic_run='existing',
                               profile='default', threads=2, seed=42, python=None, submit=True, dry_run=False)
        with patch.dict(os.environ, PATH=str(binaries) + os.pathsep + os.environ['PATH']):
            self.assertEqual(isolation.pipeline(args), 0)
            self.assertEqual(isolation.pipeline(args), 0)
        self.assertEqual((self.root / 'sbatch-count').read_text(), '1')
        self.assertEqual(len(list((self.root / 'workspace/transfers/isolation').iterdir())), 1)
        receipt = json.loads((local / 'Fornax7003/isolation/transport/client/server.json').read_text())
        self.assertEqual(receipt['job_id'], '45678')
        self.assertEqual((self.root / 'Fornax7003/isolation/transport/Fornax.ini').read_bytes(), (ic / 'Fornax.ini').read_bytes())


if __name__ == '__main__':
    unittest.main()
