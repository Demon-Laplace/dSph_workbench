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
                  'sbatch': {'clusters': 'test_cluster', 'qos': 'test_qos', 'ntasks': 2, 'time': '00:10:00'}}}}
        runtime.write_json(self.root / 'workflow.json', config)
        for name in ('Fornax.ini', 'parameters.ini'):
            (self.stage / name).write_text('fixture ' + name)
        for name in ('workflow_isolation.py', 'workflow_runtime.py'):
            shutil.copy2(Path(isolation.__file__).with_name(name), self.stage / name)
            shutil.copy2(Path(isolation.__file__).with_name(name), self.root / name)
        self.request = {'name': 'Fornax7003-test', 'model': 'Fornax7003', 'profile': 'default',
                        'build_id': 'fixture', 'files': {n: isolation.sha(self.stage / n)
                        for n in ('Fornax.ini', 'parameters.ini')}}
        runtime.write_json(self.stage / 'request.json', self.request)
        self.run = isolation.case_dir(self.root, self.request['name'], self.request['model'])

    def receive(self, submit=False):
        with patch.object(isolation.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, '12345\n', '')) as sbatch:
            result = isolation.receive(self.root, self.stage, submit)
            return result, sbatch

    def test_server_layout_is_separate_and_old_frozen_run_is_reused(self):
        self.assertEqual(self.run, self.root / 'isoFornax/Fornax7003/runs/Fornax7003-test')
        old = self.root / 'Fornax7003/isolation/Fornax7003-test'
        old.mkdir(parents=True)
        self.assertEqual(isolation.case_dir(self.root, self.run.name, 'Fornax7003'), old)
        self.run.mkdir(parents=True)
        with self.assertRaisesRegex(ValueError, 'Ambiguous'):
            isolation.case_dir(self.root, self.run.name, 'Fornax7003')

    def test_preflight_preserves_physics_and_does_not_submit(self):
        result, sbatch = self.receive()
        self.assertEqual(result['state'], 'prepared')
        self.assertIn('--test-only', sbatch.call_args.args[0])
        self.assertIn('--clusters', sbatch.call_args.args[0])
        self.assertIn('--qos', sbatch.call_args.args[0])
        self.assertIn('export TMPDIR="$PWD/.tmp"', (self.run / 'run.sh').read_text())
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
        (data / 'isoFornax').mkdir()
        (data / 'isoFornax/Fornax7003').symlink_to(outside, target_is_directory=True)
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
        receipt = json.loads((local / 'workspace/receipts/isolation/Fornax7003/transport/server.json').read_text())
        self.assertEqual(receipt['job_id'], '45678')
        self.assertEqual(list((local / 'Fornax7003/isolation').iterdir()), [])
        self.assertEqual((self.root / 'isoFornax/Fornax7003/runs/transport/Fornax.ini').read_bytes(), (ic / 'Fornax.ini').read_bytes())

    def test_cooling_is_materialized_and_required_table_cannot_disappear(self):
        table = self.root / 'table'
        table.write_text('cooling data')
        (self.root / 'template/TREECOOL').symlink_to(table)
        self.receive()
        self.assertFalse((self.run / 'TREECOOL').is_symlink())
        self.assertEqual((self.run / 'TREECOOL').read_text(), 'cooling data')
        isolation.verify_case(self.run)
        shutil.rmtree(self.run)
        table.unlink()
        config = json.loads((self.root / 'workflow.json').read_text())
        config['isolation_profiles']['default']['cooling_required'] = True
        runtime.write_json(self.root / 'workflow.json', config)
        with self.assertRaisesRegex(ValueError, 'TREECOOL'):
            self.receive()
        self.assertFalse(self.run.exists())

    def manual_model(self, name='Fornax9001'):
        model = self.root / 'isoFornax' / name
        model.mkdir(parents=True)
        (model / (name + '.ini')).write_text('numbered particle fixture ' + name)
        (model / 'GZWJL.PARAM').write_text('InitCondFile ./Fornax.IC_000\nOutputDir ./output\n')
        (self.root / 'template/TREECOOL').write_text('cooling table')
        return model

    def test_manual_numbered_ic_conversion_and_cooling_are_idempotent(self):
        model = self.manual_model()
        (model / 'TREECOOL').symlink_to('broken-old-link')
        result = isolation.prepare_legacy(self.root, model)
        binary = model / 'Fornax.IC_000'
        original = binary.stat().st_mtime_ns
        self.assertEqual(result['state'], 'prepared')
        self.assertEqual((model / 'Fornax.ini').resolve(), model / 'Fornax9001.ini')
        self.assertFalse((model / 'TREECOOL').is_symlink())
        isolation.prepare_legacy(self.root, model)
        self.assertEqual(binary.stat().st_mtime_ns, original)
        (model / 'Fornax9001.ini').write_text('changed input')
        with self.assertRaisesRegex(ValueError, 'changed'):
            isolation.prepare_legacy(self.root, model)
        self.assertEqual(binary.stat().st_mtime_ns, original)

    def test_manual_wrong_model_and_existing_binary_are_refused(self):
        model = self.manual_model()
        (model / 'Fornax9001.ini').rename(model / 'Fornax9002.ini')
        with self.assertRaisesRegex(ValueError, 'numbered'):
            isolation.prepare_legacy(self.root, model)
        self.assertFalse((model / 'Fornax.IC_000').exists())
        (model / 'Fornax9002.ini').rename(model / 'Fornax9001.ini')
        (model / 'Fornax.IC_000').write_bytes(b'keep original')
        with self.assertRaisesRegex(ValueError, 'refusing to overwrite'):
            isolation.prepare_legacy(self.root, model)
        self.assertEqual((model / 'Fornax.IC_000').read_bytes(), b'keep original')

    def test_both_local_workspaces_share_receipts(self):
        analysis = self.root / 'analysis'
        ic = self.root / 'ic'
        analysis.mkdir(); ic.mkdir()
        runtime.write_json(analysis / 'workflow.json', {'models_root': '.'})
        runtime.write_json(ic / 'workflow.json', {'models_root': '../analysis'})
        expected = analysis / 'workspace/receipts/isolation/Fornax7003/run'
        self.assertEqual(isolation.client_directory(analysis, 'Fornax7003', 'run'), expected)
        self.assertEqual(isolation.client_directory(ic, 'Fornax7003', 'run'), expected)

    def test_batch_keeps_model_identity_and_stops_at_first_failure(self):
        args = SimpleNamespace(model=['Fornax9001', 'Fornax9002', 'Fornax9003'])
        with patch.object(isolation, 'pipeline_one', side_effect=[0, 3]) as one:
            self.assertEqual(isolation.pipeline(args), 3)
            self.assertEqual([c.args[0].model for c in one.call_args_list], ['Fornax9001', 'Fornax9002'])
        args.model = ['Fornax9001', 'Fornax9001']
        with patch.object(isolation, 'pipeline_one') as one:
            with self.assertRaisesRegex(ValueError, 'Duplicate'):
                isolation.pipeline(args)
            one.assert_not_called()


if __name__ == '__main__':
    unittest.main()
