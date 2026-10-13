import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import flatten_local_isolation as flatten

spec = importlib.util.spec_from_file_location('model_paths', Path(__file__).parents[2] / 'local_tools/isolation_analysis/model_paths.py')
paths = importlib.util.module_from_spec(spec)
spec.loader.exec_module(paths)


class FlatIsolationTests(unittest.TestCase):
    def test_model_identity_survives_flat_layout(self):
        for directory in ('Fornax7003', 'Fornax7003/isolation', 'Fornax7003/isolation/output',
                          'Fornax7003/isolation/legacy/Fornax7003'):
            self.assertEqual(paths.model_id(directory), 'Fornax7003')
        self.assertEqual(paths.model_id('/tmp/Fornax'), 'Fornax')

    def test_flatten_preserves_particle_and_alias_and_refuses_collisions(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ic = root / 'old'; ic.mkdir()
            (root / 'workflow.json').write_text(json.dumps({'ic': {'root': 'old'}}))
            old = root / 'Fornax7003/isolation/legacy/Fornax7003'
            old.mkdir(parents=True)
            particle = old / 'Fornax.ini'; particle.write_bytes(b'scientific data')
            dest = old.parents[1] / particle.name
            alias = ic / 'Fornax7003'; alias.symlink_to(old)
            plan = {'root': str(root), 'source_root': str(ic), 'models': 1, 'locks': [],
                    'moves': [{'source': str(particle), 'target': str(dest), 'identity': flatten.identity(particle)}],
                    'links': [{'path': str(alias), 'old': str(old), 'new': str(dest.parent)}],
                    'metadata_updates': [], 'remove_links': [], 'remove_dirs': [str(old), str(old.parent)]}
            dest.write_text('keep')
            with self.assertRaisesRegex(ValueError, 'occupied'):
                flatten.apply(plan, root / 'journal')
            self.assertTrue(particle.exists())
            dest.unlink()
            result = flatten.apply(plan, root / 'journal')
            self.assertEqual(result['empty_directories_removed'], 2)
            self.assertEqual(flatten.identity(dest), plan['moves'][0]['identity'])
            self.assertEqual((alias / 'Fornax.ini').read_bytes(), b'scientific data')
            self.assertFalse(old.parent.exists())
