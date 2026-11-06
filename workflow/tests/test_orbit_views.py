import contextlib
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import workflow_management as management


class OrbitViewTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.model = self.root / 'Fornax1'
        self.source = self.model / 'historical/Fornax1'
        (self.source / 'output').mkdir(parents=True)
        (self.source / 'GZWJL.PARAM').write_text('frozen physics\n')
        (self.source / 'subastrofat.sh').write_text('protected launcher\n')
        (self.source / 'output/snapshot_000.hdf5').write_bytes(b'frozen snapshot')
        self.run = self.root / 'workspace/analysis-contexts/old/workspace/analysis/run'
        (self.run / 'display').mkdir(parents=True)
        (self.run / 'display/Fornax1.mp4').write_bytes(b'existing video')
        (self.run / 'manifest.json').write_text(json.dumps({'model': 'Fornax1', 'source': str(self.source)}))
        self.spec = {'source': str(self.source.relative_to(self.root)), 'entries': {
            'inputs/GZWJL.PARAM': str((self.source / 'GZWJL.PARAM').relative_to(self.root)),
            'output': str((self.source / 'output').relative_to(self.root)),
            'analysis/old-run': str(self.run.relative_to(self.root)),
            'videos/Fornax1.mp4': str((self.run / 'display/Fornax1.mp4').relative_to(self.root)),
        }}
        self.settings = {'orbit_views': {'Fornax1': {'trial': self.spec}}}
        self.args = SimpleNamespace(root=self.root, model='Fornax1', name=None, refresh=True)

    def run_command(self):
        (self.root / 'workflow.json').write_text(json.dumps(self.settings))
        with contextlib.redirect_stdout(io.StringIO()):
            management.orbit_views(self.args)

    def test_flat_entries_preserve_frozen_files_and_are_idempotent(self):
        before = {p: (p.stat().st_ino, p.read_bytes()) for p in self.source.rglob('*') if p.is_file()}
        self.run_command()
        self.run_command()
        view = self.model / 'orbit_trial'
        self.assertEqual((view / 'output/snapshot_000.hdf5').read_bytes(), b'frozen snapshot')
        self.assertEqual((view / 'videos/Fornax1.mp4').read_bytes(), b'existing video')
        self.assertEqual((view / 'analysis/old-run').resolve(), self.run)
        self.assertEqual(before, {p: (p.stat().st_ino, p.read_bytes()) for p in before})

    def test_listing_does_not_create_navigation(self):
        self.args.refresh = False
        self.run_command()
        self.assertFalse((self.model / 'orbit_trial').exists())

    def test_other_orbit_analysis_is_rejected_before_any_writes(self):
        (self.run / 'manifest.json').write_text(json.dumps({'model': 'Fornax1', 'source': str(self.model)}))
        with self.assertRaisesRegex(ValueError, 'different orbit'):
            self.run_command()
        self.assertFalse((self.model / 'orbit_trial').exists())

    def test_conflicting_link_is_preserved(self):
        self.run_command()
        link = self.model / 'orbit_trial/videos/Fornax1.mp4'
        link.unlink()
        other = self.model / 'other.mp4'
        other.write_bytes(b'other orbit video')
        link.symlink_to(other)
        with self.assertRaisesRegex(ValueError, 'points elsewhere'):
            self.run_command()
        self.assertEqual(link.resolve(), other)

    def test_external_sources_and_other_orbit_output_are_rejected(self):
        self.spec['entries']['videos/outside.mp4'] = str(self.root.parent)
        with self.assertRaisesRegex(ValueError, 'escape'):
            self.run_command()
        del self.spec['entries']['videos/outside.mp4']
        other = self.model / 'other/output'
        other.mkdir(parents=True)
        self.spec['entries']['output'] = str(other.relative_to(self.root))
        with self.assertRaisesRegex(ValueError, 'different orbit'):
            self.run_command()

    def test_removed_entries_cannot_leave_stale_orbit_links(self):
        self.run_command()
        del self.spec['entries']['videos/Fornax1.mp4']
        with self.assertRaisesRegex(ValueError, 'identity or links changed'):
            self.run_command()

    def test_local_video_only_view_does_not_require_remote_snapshots(self):
        self.spec.pop('source')
        self.spec['remote_source'] = '/remote/data/Fornax1'
        self.spec['entries'] = {'videos/Fornax1.mp4': str((self.run / 'display/Fornax1.mp4').relative_to(self.root))}
        self.run_command()
        view = self.model / 'orbit_trial'
        self.assertTrue((view / 'videos/Fornax1.mp4').is_file())
        self.assertFalse((view / 'output').exists())

    def test_unmanaged_directory_is_not_adopted(self):
        (self.model / 'orbit_trial').mkdir()
        with self.assertRaisesRegex(ValueError, 'unmanaged'):
            self.run_command()


if __name__ == '__main__':
    unittest.main()
