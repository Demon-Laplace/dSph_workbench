import importlib.util
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location('centralize', Path(__file__).parents[1] / 'centralize_isolation_python.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class CentralizeTests(unittest.TestCase):
    def test_preserves_shell_and_source_paths_and_resumes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model = root / 'isoFornax/Fornax7001'
            model.mkdir(parents=True)
            source = root / 'source'
            source.mkdir()
            code = b'from pathlib import Path\nREPO_ROOT = Path(__file__).resolve().parents[1]\n'
            (source / 'check.py').write_bytes(code)
            (model / 'check.py').write_bytes(code)
            shell = model / 'job.sh'
            shell.write_text('python3 check.py\n')
            shell_stat = shell.stat()
            files = {p.name: {'sha256': module.digest(p), 'mtime_ns': p.stat().st_mtime_ns} for p in model.iterdir()}
            rows = [{'model': model.name, 'files': files}]
            archive = root / 'observer_review_delete_20261006/test'
            for verify in (False, False, True):
                module.centralize(root, rows, source, archive, verify)
            self.assertTrue((model / 'check.py').is_symlink())
            self.assertEqual((model / 'check.py').resolve().parents[1], model.parent)
            self.assertEqual((archive / model.name / 'check.py').read_bytes(), code)
            self.assertEqual(shell.stat(), shell_stat)
            # A changed shell or nonidentical shared version is never silently accepted.
            shell.write_text('changed\n')
            with self.assertRaisesRegex(ValueError, 'Shell changed'):
                module.centralize(root, rows, source, archive)

    def test_version_conflict_preserves_original(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model = root / 'isoFornax/Fornax7001'
            model.mkdir(parents=True)
            (model / 'check.py').write_text('original\n')
            source = root / 'source'
            source.mkdir()
            (source / 'check.py').write_text('different\n')
            rows = [{'model': model.name, 'files': {'check.py': {'sha256': module.digest(model / 'check.py')}}}]
            with self.assertRaisesRegex(ValueError, 'Different Python version'):
                module.centralize(root, rows, source, root / 'observer_review_delete_20261006/test')
            self.assertFalse((model / 'check.py').is_symlink())
            self.assertEqual((model / 'check.py').read_text(), 'original\n')
