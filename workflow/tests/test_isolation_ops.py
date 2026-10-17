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

    def test_batch_submits_once_runs_sequentially_and_preserves_model_outputs(self):
        import os
        import shutil
        second = isolation.case_dir(self.root, 'other-success', 'Fornax7003')
        shutil.copytree(self.keep, second)
        for run in [self.keep, second]:
            for n in (0, 300, 600): (run / 'output' / f'snapshot_{n:03d}.hdf5').touch()
        config = json.loads((self.root / 'workflow.json').read_text())
        config['isolation_check']['sbatch'].update(time='01:00:00', nodes=2, ntasks=2)
        runtime.write_json(self.root / 'workflow.json', config)
        (self.root / 'density.py').write_text(
            'import argparse\nfrom pathlib import Path\n'
            'p=argparse.ArgumentParser();p.add_argument("--model-dir");p.add_argument("--output-dir");'
            'p.add_argument("--processes");a=p.parse_args();out=Path(a.output_dir);out.mkdir(exist_ok=True);'
            '(out/"snapshot_000_density_profile.png").write_text(str(Path.cwd()));'
            '(out/"density.mp4").write_bytes(b"video fixture");print(Path.cwd())\n')
        with patch.object(isolation, 'verify_completion'), patch.object(ops, 'active_directories', return_value=[]), \
             patch.object(ops.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, '125;test\n', '')) as submit, \
             patch.object(ops, 'accounting', return_value={'state': 'RUNNING'}):
            overrides = ['Fornax7003=other-success']
            first = ops.check_many(self.root, ['Fornax7002', 'Fornax7003'], 'success', run_names=overrides)
            again = ops.check_many(self.root, ['Fornax7002', 'Fornax7003'], 'success', run_names=overrides)
            with self.assertRaisesRegex(ValueError, 'active batch'):
                ops.check(self.root, 'Fornax7003', 'other-success', 'quick')
        submit.assert_called_once()
        command = submit.call_args.args[0]
        self.assertEqual(command[command.index('--time') + 1], '02:00:00')
        self.assertEqual(command[command.index('--nodes') + 1], '1')
        self.assertEqual(command[command.index('--ntasks') + 1], '1')
        self.assertEqual(first['job_id'], '125')
        self.assertEqual([r['job_id'] for r in again['existing']], ['125', '125'])
        import sys
        probe = self.root / 'ffprobe'
        probe.write_text('#!' + sys.executable + '\nprint(\'{"streams":[{"nb_frames":"3"}]}\')\n')
        probe.chmod(0o755)
        environment = {**os.environ, 'SLURM_JOB_ID': '125', 'PATH': str(self.root) + os.pathsep + os.environ['PATH']}
        subprocess.run(['bash', '-c', command[-1]], env=environment, check=True)
        for run in [self.keep, second]:
            self.assertEqual((run / 'density_output/snapshot_000_density_profile.png').read_text(), str(run))
            self.assertTrue((run / 'density_output/density.mp4').exists())
            self.assertIn(str(run), (run / 'density-125.log').read_text())
        probe.write_text('#!' + sys.executable + '\nprint(\'{"streams":[{"nb_frames":"2"}]}\')\n')
        (second / 'density_output/snapshot_000_density_profile.png').unlink()
        with self.assertRaises(subprocess.CalledProcessError):
            subprocess.run(['bash', '-c', command[-1]], env=environment, check=True)
        self.assertFalse((second / 'density_output/snapshot_000_density_profile.png').exists())

    def test_run_names_are_explicit_complete_and_validated_before_connection(self):
        from types import SimpleNamespace
        models = ['Fornax7002', 'Fornax7003']
        self.assertEqual(ops.density_cases(models, 'success'), [(m, 'success') for m in models])
        self.assertEqual(ops.density_cases(models, run_names=['Fornax7002=one', 'Fornax7003=two']),
                         [('Fornax7002', 'one'), ('Fornax7003', 'two')])
        invalid = [[], ['Fornax7002=one'], ['Fornax7004=one'], ['Fornax7002'],
                   ['Fornax7002='], ['Fornax7002=../escape'],
                   ['Fornax7002=one', 'Fornax7002=two']]
        for names in invalid:
            with self.subTest(names=names), patch.object(isolation, 'remote_exec') as remote:
                args = SimpleNamespace(model=models, name=None, run_names=names)
                with self.assertRaises(ValueError): ops.client(args)
                remote.assert_not_called()

    def test_cli_accepts_per_model_runs_without_shared_name(self):
        import workflow
        overrides = ['Fornax7002=one', 'Fornax7003=two']
        for operation, flag in [('isolation-check', '--fetch'), ('isolation-fetch', '--wait')]:
            arguments = ['workflow.py', '--root', str(self.root), operation, 'Fornax7002', 'Fornax7003',
                         '--run', overrides[0], '--run', overrides[1], flag]
            with self.subTest(operation=operation), patch('sys.argv', arguments), \
                 patch.object(ops, 'client', return_value=0) as client:
                self.assertEqual(workflow.main(), 0)
                parsed = client.call_args.args[0]
                self.assertIsNone(parsed.name)
                self.assertEqual(parsed.run_names, overrides)

    def test_client_sends_one_batch_request_from_either_local_workspace(self):
        from types import SimpleNamespace
        import workflow_management as management
        args = SimpleNamespace(root=self.root, model=['Fornax7002', 'Fornax7003'],
                               name='success', command='isolation-check', range='all')
        # An IC checkout can link model data outside its own root; client validation
        # must not treat it as a server filesystem layout.
        (self.root / 'Fornax7002').symlink_to(self.root.parent)
        response = subprocess.CompletedProcess([], 0, json.dumps({'job_id': '125',
            'submitted_models': args.model, 'existing': []}), '')
        with patch.object(management, 'ssh_config', return_value=({'root': '/data'}, [])), \
             patch.object(isolation, 'remote_exec', return_value=response) as remote:
            self.assertEqual(ops.client(args), 0)
        remote.assert_called_once()
        command = remote.call_args.args[1]
        self.assertIn('check-batch', command)
        self.assertEqual(command[command.index('--models') + 1:command.index('--name')], args.model)

    @unittest.skipUnless(__import__('shutil').which('ffmpeg') and __import__('shutil').which('ffprobe'), 'ffmpeg required')
    def test_video_sequence_keeps_every_nonconsecutive_snapshot(self):
        import ast
        import re
        import shutil
        import struct
        import zlib
        source = Path(isolation.__file__).parent / 'isolation_python/check_density_profile.py'
        functions = [n for n in ast.parse(source.read_text()).body if isinstance(n, ast.FunctionDef) and
                     n.name in ('create_density_video', 'get_output_paths', 'build_video_output_path',
                                'sanitize_filename_part', 'detect_model_number', 'workflow_model_name')]
        namespace = dict(Path=Path, json=json, re=re, shutil=shutil, subprocess=subprocess,
                         tempfile=tempfile, CODE_VERSION='test')
        exec(compile(ast.Module(body=functions, type_ignores=[]), str(source), 'exec'), namespace)
        output = self.keep / 'density_output'
        output.mkdir()
        def chunk(kind, data):
            return struct.pack('>I', len(data)) + kind + data + struct.pack('>I', zlib.crc32(kind + data) & 0xffffffff)
        for number, color in zip([0, 301, 601], [(255, 0, 0), (0, 255, 0), (0, 0, 255)]):
            raw = (b'\0' + bytes(color) * 16) * 16
            png = b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('>IIBBBBB', 16, 16, 8, 2, 0, 0, 0))
            png += chunk(b'IDAT', zlib.compress(raw)) + chunk(b'IEND', b'')
            (output / f'snapshot_{number:03d}_density_profile.png').write_bytes(png)
        video = namespace['create_density_video'](output, self.keep / 'output', 'Fornax', [601, 0, 301], '16:16', 18)
        metadata = json.loads(subprocess.check_output(['ffprobe', '-v', 'error', '-select_streams', 'v:0',
                              '-show_entries', 'stream=nb_frames', '-of', 'json', str(video)]))
        self.assertEqual(metadata['streams'][0]['nb_frames'], '3')
        self.assertEqual(list(output.glob('.density-video-*')), [])

    def density_receipt(self):
        directory = self.root / 'workspace/receipts/isolation-check/Fornax7002/success'
        directory.mkdir(parents=True)
        runtime.write_json(directory / 'check.json', {'job_id': '125', 'selection': 'all'})
        return directory

    def test_pending_results_can_be_polled_but_failures_abort(self):
        self.density_receipt()
        with patch.object(ops, 'accounting', return_value={'state': 'RUNNING', 'exit_code': '0:0'}):
            self.assertFalse(ops.results(self.root, 'Fornax7002', 'success', True)['ready'])
            with self.assertRaises(ValueError): ops.results(self.root, 'Fornax7002', 'success')
        for state, code in [('FAILED', '1:0'), ('COMPLETED', '1:0')]:
            with patch.object(ops, 'accounting', return_value={'state': state, 'exit_code': code}):
                with self.assertRaises(ValueError): ops.results(self.root, 'Fornax7002', 'success', True)

    def test_results_verify_full_video_frames_and_hash(self):
        self.density_receipt()
        out = self.keep / 'density_output'; out.mkdir()
        (out/'snapshot_000_density_profile.png').touch()
        (out/'video.mp4').write_bytes(b'video')
        for n in (0, 1, 2): (self.keep/'output'/f'snapshot_{n:03d}.hdf5').touch()
        with patch.object(ops, 'accounting', return_value={'state':'COMPLETED', 'exit_code':'0:0'}), \
             patch.object(ops, 'active_directories', return_value=[]), \
             patch.object(ops.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0,
                 json.dumps({'streams':[{'nb_frames':'3', 'duration':'0.12'}]}), '')) as probe:
            result = ops.results(self.root, 'Fornax7002', 'success')
            self.assertEqual(result['videos'][0]['frames'], 3)
            ops.verify_download(out, result)
            (out/'video.mp4').write_bytes(b'corrupt')
            with self.assertRaisesRegex(ValueError, 'checksum'): ops.verify_download(out, result)
            probe.return_value.stdout = json.dumps({'streams':[{'nb_frames':'2', 'duration':'0.08'}]})
            with self.assertRaisesRegex(ValueError, 'frame count'): ops.results(self.root, 'Fornax7002', 'success')

    def test_full_download_requires_server_video_verification(self):
        with self.assertRaisesRegex(ValueError, 'Missing verified video metadata'):
            ops.verify_download(self.root, {'analysis': {'selection': 'all'}})

    def test_wait_download_polls_until_ready_and_obeys_timeout(self):
        from types import SimpleNamespace
        args = SimpleNamespace(timeout=60)
        pending = {'ready':False, 'analysis':{'job_id':'125', 'scheduler':{'state':'RUNNING'}}}
        ready = {'ready':True}
        replies = [subprocess.CompletedProcess([], 0, json.dumps(x), '') for x in (pending, ready)]
        with patch.object(isolation, 'remote_exec', side_effect=replies) as remote, patch.object(ops.time, 'sleep') as sleep:
            self.assertEqual(ops.await_results(args, ['results'], True), ready)
            self.assertEqual(remote.call_count, 2)
            self.assertIn('--allow-pending', remote.call_args.args[1])
            sleep.assert_called_once()
        with patch.object(isolation, 'remote_exec', return_value=replies[0]), \
             patch.object(ops.time, 'monotonic', side_effect=[0, 61]), patch.object(ops.time, 'sleep') as sleep:
            with self.assertRaisesRegex(ValueError, 'Timed out'): ops.await_results(args, ['results'], True)
            sleep.assert_not_called()

    def test_check_fetch_downloads_both_models_after_single_batch_request(self):
        from types import SimpleNamespace
        import workflow_management as management
        args = SimpleNamespace(root=self.root, model=['Fornax7002', 'Fornax7003'], name='success',
                               run_names=['Fornax7003=other-success'],
                               command='isolation-check', range='all', fetch=True, timeout=60)
        batch = {'job_id':'125', 'submitted_models':args.model, 'existing':[]}
        replies = [batch] + [{'ready':True, 'directory':'/data/'+model+'/density_output',
                              'analysis':{'job_id':'125'}, 'videos':[]} for model in args.model]
        with patch.object(management, 'ssh_config', return_value=({'root':'/data', 'host':'server'}, ['ssh'])), \
             patch.object(isolation, 'remote_exec', side_effect=[subprocess.CompletedProcess([], 0, json.dumps(r), '') for r in replies]) as remote, \
             patch.object(ops.subprocess, 'run') as rsync:
            self.assertEqual(ops.client(args), 0)
        self.assertEqual(remote.call_count, 3)
        self.assertIn('check-batch', remote.call_args_list[0].args[1])
        command = remote.call_args_list[0].args[1]
        self.assertEqual(command[command.index('--run') + 1], 'Fornax7003=other-success')
        self.assertEqual(rsync.call_count, 2)
        for model in args.model:
            self.assertTrue((self.root/model/'isolation').is_dir())
            self.assertFalse((self.root/model/'isolation/density_output').exists())
            name = 'success' if model == 'Fornax7002' else 'other-success'
            request = remote.call_args_list[1 if model == 'Fornax7002' else 2].args[1]
            self.assertEqual(request[request.index('--name') + 1], name)
            saved = isolation.client_directory(self.root, model, name)/'density-result.json'
            self.assertIn('download_verified_utc', json.loads(saved.read_text()))
            source = saved.parent.parent / 'density-source.json'
            self.assertEqual(json.loads(source.read_text())['run'], name)
        # The same mapping also works when resuming downloads without submission.
        args.command = 'isolation-fetch'
        args.fetch = False
        args.wait = True
        with patch.object(management, 'ssh_config', return_value=({'root':'/data', 'host':'server'}, ['ssh'])), \
             patch.object(isolation, 'remote_exec', side_effect=[subprocess.CompletedProcess([], 0, json.dumps(r), '') for r in replies[1:]]) as remote, \
             patch.object(ops.subprocess, 'run'):
            self.assertEqual(ops.client(args), 0)
        self.assertEqual(remote.call_count, 2)
        self.assertTrue(all('check-batch' not in c.args[1] for c in remote.call_args_list))

    @unittest.skipUnless(__import__('shutil').which('rsync'), 'rsync required')
    def test_fetch_only_video_and_summary_preserves_existing_ic(self):
        from types import SimpleNamespace
        import workflow_management as management
        source = self.root/'remote-output'; source.mkdir()
        (source/'video.mp4').write_bytes(b'verified video')
        (source/'density_evolution.csv').write_text('snapshot,density\n0,1\n')
        (source/'snapshot_000_density_profile.png').write_bytes(b'image')
        (source/'Fornax.ini').write_text('must not replace IC')
        (source/'other.csv').write_text('unrelated table')
        (source/'nested').mkdir(); (source/'nested/other.mp4').write_bytes(b'nested video')
        destination = self.root/'Fornax7002/isolation'; destination.mkdir(parents=True)
        (destination/'Fornax.ini').write_text('original IC')
        (destination/'components.png').write_bytes(b'original diagnostic')
        result = {'ready': True, 'directory':'/remote/output', 'analysis':{'selection':'all'},
                  'videos':[{'name':'video.mp4', 'sha256':isolation.sha(source/'video.mp4'), 'bytes':14}]}
        args = SimpleNamespace(root=self.root, model=['Fornax7002'], name='success', command='isolation-fetch')
        real_run = subprocess.run
        def transfer(command, **kwargs):
            command = list(command); command[-2] = str(source)+'/'
            return real_run(command, **kwargs)
        with patch.object(management, 'ssh_config', return_value=({'root':'/data', 'host':'server'}, ['ssh'])), \
             patch.object(isolation, 'remote_exec', return_value=subprocess.CompletedProcess([], 0, json.dumps(result), '')), \
             patch.object(ops.subprocess, 'run', side_effect=transfer):
            self.assertEqual(ops.client(args), 0)
            self.assertEqual(set(p.name for p in destination.iterdir()),
                             {'video.mp4','density_evolution.csv','Fornax.ini','components.png'})
            self.assertEqual((destination/'Fornax.ini').read_text(),'original IC')
            source_receipt = isolation.client_directory(self.root, 'Fornax7002', 'success').parent/'density-source.json'
            source_receipt.unlink()
            with self.assertRaisesRegex(ValueError, 'unknown provenance'): ops.client(args)
            self.assertEqual((destination/'video.mp4').read_bytes(), b'verified video')

    def test_check_refuses_active_directory(self):
        (self.keep / 'output/snapshot_000.hdf5').touch()
        with patch.object(isolation, 'verify_completion'), patch.object(ops, 'active_directories', return_value=[self.keep]):
            with self.assertRaisesRegex(ValueError, 'already using'):
                ops.check(self.root, 'Fornax7002', 'success')


if __name__ == '__main__':
    unittest.main()
