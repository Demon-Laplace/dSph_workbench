import contextlib
import io
import json
from pathlib import Path
import tempfile
import subprocess
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import workflow_cleanup as cleanup
import workflow_runtime as runtime
import workflow_models as models
import workflow_release as release
import workflow_github as github
import workflow_archive as archive
import workflow_management as management


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        model = self.root / 'Fornax1'
        (model / 'output').mkdir(parents=True)
        for name, text in {'Fornax.ini': '1 0 0 0\n', 'GZWJL.PARAM': 'test',
                           'mergegalaxy.par': 'MW.ini Fornax.ini\n', 'mergecube.par': 'cube.ini\n'}.items():
            (model / name).write_text(text)
        (model / 'output/snapshot_001.hdf5').write_bytes(b'snapshot fixture')
        code = self.root / 'code'
        code.mkdir()
        for name in ('GetInfo.py', 'PlotFig.py', 'data_processing.py'):
            (code / name).write_text('# fixture\n')
        (self.root / 'workflow.json').write_text(json.dumps({'analysis_code': 'code'}))
        self.args = SimpleNamespace(root=self.root, model='Fornax1', name='test', snapshots=None, dry_run=False)

    def prepare(self):
        with contextlib.redirect_stdout(io.StringIO()):
            runtime.prepare(self.args)
        return self.root / 'workspace/analysis/test'

    def test_dry_run_creates_nothing(self):
        self.args.dry_run = True
        self.prepare()
        self.assertFalse((self.root / 'workspace').exists())

    def test_pins_code_and_links_data_without_importing_legacy_summary(self):
        (self.root / 'Fornax1/elinfo_Fornax1.csv').write_text('old schema')
        run = self.prepare()
        self.assertTrue((run / 'Fornax1/output/snapshot_001.hdf5').is_symlink())
        self.assertFalse((run / 'Fornax1/elinfo_Fornax1.csv').exists())
        self.assertFalse((run / 'code/GetInfo.py').is_symlink())
        manifest = json.loads((run / 'manifest.json').read_text())
        runtime.verify_inputs(run, manifest)
        (self.root / 'Fornax1/output/snapshot_001.hdf5').write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'Source changed'):
            runtime.verify_inputs(run, manifest)

    def test_existing_run_and_escape_rejected(self):
        self.prepare()
        with self.assertRaisesRegex(ValueError, 'already exists'):
            self.prepare()
        self.args.name = '../../outside'
        with self.assertRaises(ValueError):
            self.prepare()

    def test_ambiguous_snapshot_locations_rejected(self):
        (self.root / 'Fornax1/snapshot_002.hdf5').write_bytes(b'test')
        with self.assertRaisesRegex(ValueError, 'ambiguous'):
            self.prepare()

    def test_summary_validation_catches_upstream_false_success(self):
        path = self.root / 'summary.csv'
        path.write_text('numsp,star_mass,distance,sigma\n1,2,3,4\n')
        runtime.verify_summary(path, [1])
        with self.assertRaisesRegex(ValueError, 'missing or duplicate'):
            runtime.verify_summary(path, [1, 2])

    def test_cleanup_only_recognizes_scheduler_logs(self):
        (self.root / 'slurm-123.out').write_text('log')
        (self.root / 'density_001.out').write_text('science data')
        (self.root / 'custom_123.err').write_text('error')
        (self.root / 'job.slurm').write_text('#SBATCH --error=custom_%j.err\n')
        (self.root / 'slurm-456.out').symlink_to(self.root / 'density_001.out')
        names = {p.name for p, _, _ in cleanup.candidates(self.root)}
        self.assertEqual(names, {'slurm-123.out', 'custom_123.err'})

    def test_new_case_preserves_template_and_rejects_existing_case(self):
        template = self.root / 'Fornax1'
        for name in ('MW.ini', 'cube.ini'):
            (template / name).write_text('particle fixture')
        (template / 'subastrothin.sh').write_text('#!/bin/bash\n')
        (template / 'subgetinfo.sh').write_text('obsolete relative paths')
        args = SimpleNamespace(root=self.root, model='Fornax2', template='Fornax1',
                               dwarf_ic='Fornax1/Fornax.ini', mw_ic=None)
        with contextlib.redirect_stdout(io.StringIO()):
            models.new_case(args)
        case = self.root / 'workspace/models/production/Fornax2'
        self.assertEqual((self.root / 'Fornax2').resolve(), case)
        self.assertEqual((template / 'mergegalaxy.par').read_text(), 'MW.ini Fornax.ini\n')
        self.assertFalse((case / 'subgetinfo.sh').exists())
        self.assertFalse((case / 'output').exists())
        with contextlib.redirect_stdout(io.StringIO()):
            management.organize(SimpleNamespace(root=self.root, repo=self.root / 'missing-checkout'))
        self.assertTrue(case.is_dir())
        self.assertFalse(case.is_symlink())
        with self.assertRaisesRegex(ValueError, 'already exists'):
            models.new_case(args)

    def test_simulation_refuses_existing_snapshots(self):
        model = self.root / 'Fornax1'
        (model / 'GZWJL.PARAM').write_text('InitCondFile Fornax.ini\nOutputDir output\n')
        (model / 'subastrothin.sh').write_text('#!/bin/bash\n')
        args = SimpleNamespace(root=self.root, model='Fornax1', launcher='subastrothin.sh', name='sim', submit=False)
        with self.assertRaisesRegex(ValueError, 'already contains data'):
            models.simulate(args)

    def test_github_source_allowlist_excludes_model_data(self):
        source = self.root / 'code'
        for name in ('snapshot_001.hdf5', 'Fornax.ini', 'slurm-123.out', 'config.json', 'figure.png'):
            (source / name).write_text('not source')
        (source / 'link.py').symlink_to(self.root / 'workflow.json')
        self.assertEqual({p.name for p in github.source_files(source)}, {'GetInfo.py', 'PlotFig.py', 'data_processing.py'})
        (source / 'GetInfo.py').write_text('invalid python !')
        with self.assertRaises(SyntaxError):
            github.source_files(source)

    def test_tampered_release_cannot_change_active_configuration(self):
        directory = self.root / 'workspace/releases/test/repo'
        directory.mkdir(parents=True)
        code = directory / 'GetInfo.py'
        code.write_text('# original\n')
        runtime.write_json(directory.parent / 'release.json', {'git_commit': 'fixture',
                           'files': {'repo/GetInfo.py': runtime.sha(code)}})
        settings = (self.root / 'workflow.json').read_bytes()
        code.write_text('# changed\n')
        with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
            release.activate(SimpleNamespace(root=self.root, name='test'))
        self.assertEqual(settings, (self.root / 'workflow.json').read_bytes())
        code.write_text('# original\n')
        (directory / 'untracked.py').write_text('# unexpected')
        with self.assertRaisesRegex(ValueError, 'Unexpected files'):
            release.verify(directory.parent)

    def test_archive_roundtrip_and_no_overwrite(self):
        source = self.root / 'old-check'
        source.mkdir(); (source / 'figure.png').write_bytes(b'original')
        plan = self.root / 'plan.json'
        runtime.write_json(plan, {'moves': [{'source': 'old-check', 'destination': 'archive/old-check'}]})
        with contextlib.redirect_stdout(io.StringIO()):
            archive.execute(self.root, plan)
            archive.execute(self.root, plan, restore=True)
        self.assertEqual((source / 'figure.png').read_bytes(), b'original')
        (self.root / 'archive/old-check').mkdir()
        with self.assertRaisesRegex(ValueError, 'Destination exists'):
            archive.execute(self.root, plan)

    def test_github_sync_refuses_stale_source_over_new_remote_code(self):
        remote = self.root / 'remote.git'
        checkout = self.root / 'checkout'
        def git_at(path, *arguments):
            return subprocess.run(['git', '-C', str(path), *arguments], check=True,
                                  text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
        subprocess.run(['git', 'init', '--bare', str(remote)], check=True, capture_output=True)
        checkout.mkdir()
        git_at(checkout, 'init', '-b', 'main')
        git_at(checkout, 'config', 'user.name', 'Workflow test')
        git_at(checkout, 'config', 'user.email', 'workflow-test@example.invalid')
        git_at(checkout, 'remote', 'add', 'origin', str(remote))
        (checkout / 'repo').mkdir(); (checkout / 'workflow').mkdir()
        code = checkout / 'repo/GetInfo.py'
        code.write_text('# baseline\n')
        runtime.write_json(checkout / 'workflow/source_manifest.json', {'files': {'GetInfo.py': runtime.sha(code)}})
        git_at(checkout, 'add', '.')
        git_at(checkout, 'commit', '-m', 'baseline')
        git_at(checkout, 'push', '-u', 'origin', 'main')
        source = self.root / 'source'; source.mkdir()
        (source / 'GetInfo.py').write_text('# baseline\r\n')
        (source / 'secret-config.json').write_text('must stay outside Git')
        settings = {'github': {'checkout': 'checkout', 'source': 'source', 'branch': 'main', 'url': str(remote)}}
        runtime.write_json(self.root / 'workflow.json', settings)
        args = SimpleNamespace(root=self.root, message='test source sync')
        with contextlib.redirect_stdout(io.StringIO()):
            github.sync(args)
        self.assertFalse((checkout / 'repo/secret-config.json').exists())
        peer = self.root / 'peer'
        subprocess.run(['git', 'clone', '-b', 'main', str(remote), str(peer)], check=True, capture_output=True)
        git_at(peer, 'config', 'user.name', 'Workflow test')
        git_at(peer, 'config', 'user.email', 'workflow-test@example.invalid')
        (peer / 'repo/GetInfo.py').write_text('# newer remote change\n')
        runtime.write_json(peer / 'workflow/source_manifest.json', {'files': {'GetInfo.py': runtime.sha(peer / 'repo/GetInfo.py')}})
        git_at(peer, 'add', '.')
        git_at(peer, 'commit', '-m', 'peer change')
        git_at(peer, 'push', 'origin', 'main')
        # Even a previously updated checkout must not make stale server code look new.
        git_at(checkout, 'fetch', 'origin')
        git_at(checkout, 'merge', '--ff-only', 'origin/main')
        with self.assertRaisesRegex(ValueError, 'both changed'):
            github.sync(args)
        self.assertEqual(code.read_text(), '# newer remote change\n')
        self.assertEqual((source / 'GetInfo.py').read_text(), '# baseline\n')


if __name__ == '__main__':
    unittest.main()
