from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import workflow_isolation as isolation
import workflow_runtime as runtime
import workflow_versions as versions


class RunVersionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.template = self.root / 'template'
        self.template.mkdir()
        self.template_text = ('InitCondFile unused\nOutputDir unused\nResubmitOn 1\n'
                              'TimeMax 6\nSofteningGas 0.1\n')
        (self.template / 'GZWJL.PARAM').write_text(self.template_text)
        self.settings = {
            'template': 'template',
            'converter': ['converter', '{ascii}', '{binary}'],
            'solver': ['solver', '{parameters}'],
            'sbatch': {'ntasks': 2, 'time': '01:00:00'},
        }
        runtime.write_json(self.root / 'workflow.json', {
            'isolation_profiles': {'default': self.settings}})
        self.stage = self.root / 'stage'
        self.stage.mkdir()
        for name in ('Fornax.ini', 'parameters.ini'):
            (self.stage / name).write_text('particle input placeholder')
        self.request = {
            'model': 'Fornax9001', 'name': 'isolate', 'profile': 'default',
            'build_id': 'build-a', 'source_directory': '/client/ic/build-a',
            'files': {'Fornax.ini': 'Fornax.ini', 'parameters.ini': 'parameters.ini'}}
        self.directory = self.root / 'isoFornax/Fornax9001'

    def receive(self, request=None, submit=False):
        runtime.write_json(self.stage / 'request.json', request or self.request)
        return isolation.receive(self.root, self.stage, submit=submit)

    def ledger(self):
        return json.loads((self.directory / versions.LEDGER).read_text())

    def test_first_version_has_json_and_generated_text_and_no_submission(self):
        with patch.object(isolation.subprocess, 'run') as process, \
                patch.object(isolation, 'sha', side_effect=AssertionError('hashing forbidden')):
            result = self.receive()
        process.assert_not_called()
        self.assertEqual((result['ic_version'], result['parameter_version']), ('IC001', 'P001'))
        self.assertEqual(result['name'], 'isolate')
        self.assertFalse(result['reused'])
        self.assertEqual(self.ledger()['runs']['isolate']['state'], 'prepared')
        text = (self.directory / versions.TEXT).read_text()
        self.assertIn('ic_version: IC001', text)
        self.assertIn('state: prepared', text)

    def test_changed_build_gets_a_new_run_without_touching_old_inputs(self):
        old = self.receive()
        path = Path(old['directory']) / 'Fornax.ini'
        original = path.read_text()
        (self.stage / 'Fornax.ini').write_text('new particles')
        new = dict(self.request, build_id='build-b', source_directory='/client/ic/build-b')
        result = self.receive(new)
        self.assertEqual(result['name'], 'isolate-IC002-P001-01')
        self.assertEqual((result['ic_version'], result['parameter_version']), ('IC002', 'P001'))
        self.assertEqual(path.read_text(), original)
        self.assertEqual((Path(result['directory']) / 'Fornax.ini').read_text(), 'new particles')
        self.assertEqual(self.ledger()['runs']['isolate']['request']['build_id'], 'build-a')

    def test_changed_effective_parameters_get_new_version_and_stable_retry(self):
        first = self.receive()
        cfg = runtime.config(self.root)
        cfg['isolation_profiles']['default']['parameter_overrides'] = {'TimeMax': '7'}
        runtime.write_json(self.root / 'workflow.json', cfg)
        second = self.receive()
        third = self.receive()
        self.assertEqual(second['ic_version'], first['ic_version'])
        self.assertEqual(second['parameter_version'], 'P002')
        self.assertEqual(second['name'], third['name'])
        self.assertTrue(third['reused'])
        self.assertIn('TimeMax 7', (Path(second['directory']) / 'GZWJL.PARAM').read_text())
        self.assertIn('TimeMax 6', (Path(first['directory']) / 'GZWJL.PARAM').read_text())

    def test_template_parameters_are_snapshotted_at_resolution(self):
        selected = isolation.resolve_request(self.root, self.request)
        (self.template / 'GZWJL.PARAM').write_text(self.template_text.replace('TimeMax 6', 'TimeMax 9'))
        original = self.receive(selected['request'])
        self.assertIn('TimeMax 6', (Path(original['directory']) / 'GZWJL.PARAM').read_text())
        next_version = self.receive()
        self.assertEqual(next_version['parameter_version'], 'P002')
        self.assertIn('TimeMax 9', (Path(next_version['directory']) / 'GZWJL.PARAM').read_text())

    def test_same_build_can_move_without_fabricating_a_new_ic_version(self):
        self.receive()
        result = self.receive(dict(self.request, source_directory='/different/client/build-a'))
        self.assertTrue(result['reused'])
        self.assertEqual(len(self.ledger()['ics']), 1)

    def test_explicit_import_version_changes_even_at_same_source_directory(self):
        a = dict(self.request, source_version='import-v1')
        b = dict(self.request, source_version='import-v2')
        self.receive(a)
        result = self.receive(b)
        self.assertEqual(result['ic_version'], 'IC002')
        self.assertEqual(len(self.ledger()['runs']), 2)

    def test_source_alias_never_hides_a_changed_build_id(self):
        self.receive(dict(self.request, source_version='label'))
        changed = dict(self.request, source_version='label', build_id='build-b')
        self.assertEqual(self.receive(changed)['ic_version'], 'IC002')

    def test_changed_build_without_source_directory_never_reuses_old_binding(self):
        old = {key: value for key, value in self.request.items() if key != 'source_directory'}
        self.receive(old)
        result = self.receive(dict(old, build_id='build-b'))
        self.assertFalse(result['reused'])
        self.assertEqual(result['ic_version'], 'IC002')
        self.assertFalse(isolation.request_matches(old, dict(old, build_id='build-b')))

    def test_repeating_submitted_version_does_not_submit_again(self):
        with patch.object(isolation.subprocess, 'run',
                          return_value=subprocess.CompletedProcess([], 0, '12345\n', '')) as process:
            first = self.receive(submit=True)
            second = self.receive(submit=True)
        self.assertEqual(process.call_count, 1)
        self.assertEqual(second['job_id'], first['job_id'])
        self.assertTrue(second['reused'])
        entry = self.ledger()['runs']['isolate']
        self.assertEqual((entry['state'], entry['job_id']), ('submitted', '12345'))
        self.assertIn('job_id: 12345', (self.directory / versions.TEXT).read_text())

    def test_actual_generated_name_also_reuses_its_existing_binding(self):
        self.receive()
        changed = dict(self.request, build_id='build-b')
        second = self.receive(changed)
        direct = self.receive(dict(changed, name=second['name']))
        self.assertEqual(direct['name'], second['name'])
        self.assertTrue(direct['reused'])

    def test_concurrent_resolutions_reserve_one_run(self):
        def resolve(_):
            return isolation.resolve_request(self.root, self.request)
        with ThreadPoolExecutor(max_workers=8) as workers:
            records = list(workers.map(resolve, range(16)))
        self.assertEqual({item['request']['name'] for item in records}, {'isolate'})
        self.assertEqual(len(self.ledger()['runs']), 1)

    def test_concurrent_receivers_submit_once(self):
        selected = isolation.resolve_request(self.root, self.request)
        stages = []
        for number in range(2):
            stage = self.root / ('upload-' + str(number))
            stage.mkdir()
            runtime.write_json(stage / 'request.json', selected['request'])
            for name in ('Fornax.ini', 'parameters.ini'):
                (stage / name).write_text('same declared build')
            stages.append(stage)
        with patch.object(isolation.subprocess, 'run',
                          return_value=subprocess.CompletedProcess([], 0, '12345\n', '')) as process:
            with ThreadPoolExecutor(max_workers=2) as workers:
                results = list(workers.map(lambda p: isolation.receive(self.root, p, True), stages))
        self.assertEqual(process.call_count, 1)
        self.assertEqual({item['job_id'] for item in results}, {'12345'})

    def legacy(self, request):
        run = isolation.case_dir(self.root, request['name'], request['model'])
        run.mkdir(parents=True)
        runtime.write_json(run / 'request.json', request)
        runtime.write_json(run / 'status.json', {'state': 'submitted', 'job_id': 'old-job'})
        (run / 'Fornax.ini').write_text('historical particles')
        return run

    def test_known_legacy_build_reuses_as_legacy_without_rewriting_it(self):
        run = self.legacy(self.request)
        original = (run / 'request.json').read_bytes()
        selected = isolation.resolve_request(self.root, self.request)
        self.assertEqual(selected['request']['version_contract'], 'legacy')
        self.assertNotIn('ic_version', selected['request'])
        result = isolation.receive(self.root, run, submit=True)
        self.assertEqual(result['job_id'], 'old-job')
        self.assertEqual((run / 'request.json').read_bytes(), original)

    def test_unidentified_legacy_run_is_preserved_and_new_input_gets_new_run(self):
        legacy = {key: value for key, value in self.request.items()
                  if key not in ('build_id', 'source_directory')}
        run = self.legacy(legacy)
        result = self.receive()
        self.assertNotEqual(result['name'], 'isolate')
        self.assertEqual((run / 'Fornax.ini').read_text(), 'historical particles')
        self.assertEqual(self.ledger()['runs']['isolate']['request']['version_contract'], 'legacy')
        self.assertEqual(result['ic_version'], 'IC001')

    def test_legacy_direct_resume_without_source_fields_remains_allowed(self):
        legacy = {key: value for key, value in self.request.items()
                  if key not in ('build_id', 'source_directory')}
        run = self.legacy(legacy)
        self.assertEqual(isolation.receive(self.root, run, True)['job_id'], 'old-job')

    def test_wrong_version_metadata_cannot_write_into_an_existing_run(self):
        selected = isolation.resolve_request(self.root, self.request)
        with self.assertRaisesRegex(ValueError, 'binding'):
            self.receive(dict(selected['request'], build_id='different-build'))
        self.assertFalse(Path(selected['directory']).exists())

    def test_failed_or_uncertain_submission_does_not_create_duplicate_jobs(self):
        with patch.object(isolation.subprocess, 'run', side_effect=OSError('lost response')):
            with self.assertRaises(OSError):
                self.receive(submit=True)
        with patch.object(isolation.subprocess, 'run') as process:
            with self.assertRaisesRegex(ValueError, 'uncertain'):
                self.receive(submit=True)
        process.assert_not_called()
        self.assertEqual(self.ledger()['runs']['isolate']['state'], 'submitting')

    def test_definitive_scheduler_rejection_keeps_prepared_binding(self):
        with patch.object(isolation.subprocess, 'run',
                          return_value=subprocess.CompletedProcess([], 1, '', 'resource limit')):
            with self.assertRaisesRegex(ValueError, 'rejected'):
                self.receive(submit=True)
        self.assertEqual(self.ledger()['runs']['isolate']['state'], 'prepared')
        self.assertNotIn('job_id: 12345', (self.directory / versions.TEXT).read_text())
        with patch.object(isolation.subprocess, 'run',
                          return_value=subprocess.CompletedProcess([], 0, '12345\n', '')):
            result = self.receive(submit=True)
        self.assertEqual(result['job_id'], '12345')
        self.assertEqual(len(self.ledger()['runs']), 1)

    def test_reservation_is_not_reported_as_a_real_job(self):
        isolation.resolve_request(self.root, self.request)
        entry = self.ledger()['runs']['isolate']
        self.assertEqual(entry['state'], 'reserved')
        self.assertNotIn('job_id', entry)
        self.assertIn('job_id: -', (self.directory / versions.TEXT).read_text())

    def test_failed_upload_does_not_rebind_an_existing_run(self):
        first = self.receive()
        old_request = (Path(first['directory']) / 'request.json').read_text()
        newer = dict(self.request, build_id='build-b')
        (self.stage / 'parameters.ini').unlink()
        with self.assertRaises(FileNotFoundError):
            self.receive(newer)
        self.assertEqual((Path(first['directory']) / 'request.json').read_text(), old_request)
        reserved = self.ledger()['runs']['isolate-IC002-P001-01']
        self.assertEqual(reserved['state'], 'reserved')
        self.assertFalse(Path(reserved['directory']).exists())

    def test_complete_state_is_process_completion_not_scientific_validation(self):
        result = self.receive()
        run = Path(result['directory'])
        runtime.write_json(run / 'status.json', {'state': 'submitted', 'job_id': '12345'})
        def process(command, **kwargs):
            if command[0] == 'converter':
                (run / 'Fornax.IC_000').write_text('converted fixture')
            return subprocess.CompletedProcess(command, 0)
        with patch.dict(os.environ, SLURM_JOB_ID='12345'), \
                patch.object(isolation.subprocess, 'run', side_effect=process):
            isolation.run_job(self.root, 'isolate', 'Fornax9001')
        self.assertEqual(self.ledger()['runs']['isolate']['state'], 'complete')


if __name__ == '__main__':
    unittest.main()
