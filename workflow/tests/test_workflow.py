import contextlib
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import workflow_cleanup as cleanup
import workflow_runtime as runtime
import workflow_models as models
import workflow_release as release
import workflow_github as github
import workflow_archive as archive


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
        directory = self.root / 'workspace/releases/test/production_analysis'
        directory.mkdir(parents=True)
        code = directory / 'GetInfo.py'
        code.write_text('# original\n')
        runtime.write_json(directory.parent / 'release.json', {'git_commit': 'fixture',
                           'files': {'production_analysis/GetInfo.py': runtime.sha(code)}})
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


if __name__ == '__main__':
    unittest.main()
