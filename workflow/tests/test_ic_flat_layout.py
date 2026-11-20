import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SOURCE = Path(__file__).resolve().parents[2] / 'local_tools/ic_generation/flat_layout.py'
SPEC = importlib.util.spec_from_file_location('ic_flat_layout_test', SOURCE)
LAYOUT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(LAYOUT)


class FlatICLayoutTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.model = self.root / 'Fornax8003/isolation'

    def build(self, name):
        directory = self.root / 'workspace/ic-builds/Fornax8003' / name
        directory.mkdir(parents=True)
        (directory / 'build.json').write_text(json.dumps(
            {'state': 'complete', 'build_id': name, 'particle_ic': 'Fornax8003.ini'}))
        (directory / 'Fornax8003.ini').write_text('particle fixture ' + name)
        (directory / 'Fornax.ini').symlink_to('Fornax8003.ini')
        (directory / 'parameters.ini').write_text('parameter fixture')
        return directory

    def test_files_are_moved_to_flat_isolation_and_fixed_names_still_resolve(self):
        directory = self.build('v1')
        inode = (directory / 'Fornax8003.ini').stat().st_ino
        result = LAYOUT.publish(directory, self.model)
        self.assertEqual(result['prefix'], '')
        self.assertEqual((self.model / 'Fornax8003.ini').stat().st_ino, inode)
        self.assertTrue((directory / 'Fornax8003.ini').is_symlink())
        self.assertEqual((self.model / 'Fornax.ini').resolve(), self.model / 'Fornax8003.ini')
        self.assertEqual((directory / 'Fornax.ini').resolve(), self.model / 'Fornax8003.ini')
        self.assertEqual(json.loads((directory / 'build.json').read_text())['build_id'], 'v1')

    def test_second_version_uses_filename_prefix_without_overwriting_first(self):
        LAYOUT.publish(self.build('v1'), self.model)
        directory = self.build('v2')
        result = LAYOUT.publish(directory, self.model)
        self.assertEqual(result['prefix'], 'v2__')
        self.assertEqual((self.model / 'Fornax8003.ini').read_text(), 'particle fixture v1')
        self.assertEqual((self.model / 'v2__Fornax.ini').resolve(),
                         self.model / 'v2__Fornax8003.ini')

    def test_repeat_publication_moves_only_new_catalogue_file(self):
        directory = self.build('v1')
        LAYOUT.publish(directory, self.model)
        self.assertEqual(LAYOUT.publish(directory, self.model)['moved_files'], 0)
        (directory / 'catalog-entry.json').write_text('{}')
        self.assertEqual(LAYOUT.publish(directory, self.model)['moved_files'], 1)
        self.assertTrue((self.model / 'catalog-entry.json').is_file())

    def test_conflicting_explicit_prefix_preserves_existing_files(self):
        LAYOUT.publish(self.build('v1'), self.model)
        with self.assertRaises(FileExistsError):
            LAYOUT.publish(self.build('v2'), self.model, prefix='')
        self.assertEqual((self.model / 'Fornax8003.ini').read_text(), 'particle fixture v1')
