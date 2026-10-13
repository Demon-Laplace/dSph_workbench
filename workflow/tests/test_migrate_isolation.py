import json
from pathlib import Path
import tempfile
import unittest
import migrate_isolation_layout as migration


class MigrationTests(unittest.TestCase):
    def test_rename_preserves_data_and_relative_links_and_is_idempotent(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); old = root/'old'; source = old/'Fornax7001'
            source.mkdir(parents=True); (source/'output').mkdir()
            particle = source/'output/snapshot'; particle.write_bytes(b'unchanged science')
            inode = particle.stat().st_ino
            (old/'shared').write_text('shared'); (source/'external').symlink_to('../shared')
            (source/'internal').symlink_to('output/snapshot')
            plan = dict(root=str(root), source_root=str(old), models=['Fornax7001'], receipt=str(root/'receipt.jsonl'))
            migration.migrate(plan); migration.migrate(plan)
            target = root/'Fornax7001/isolation/legacy/Fornax7001'
            self.assertTrue(source.is_symlink())
            self.assertEqual((target/'output/snapshot').stat().st_ino, inode)
            self.assertEqual((target/'external').read_text(), 'shared')
            self.assertEqual((target/'internal').read_bytes(), b'unchanged science')
            self.assertEqual(len((root/'receipt.jsonl').read_text().splitlines()), 2)

    def test_server_and_occupied_target_refuse_move(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); source=root/'old/Fornax7001';source.mkdir(parents=True)
            plan=dict(root=str(root),source_root=str(source.parent),models=['Fornax7001'],receipt=str(root/'receipt.jsonl'))
            with self.assertRaisesRegex(ValueError,'local-only'):
                migration.migrate(plan, server=True)
            target=root/'Fornax7001/isolation/legacy/Fornax7001';target.mkdir(parents=True)
            with self.assertRaisesRegex(ValueError,'occupied'):migration.migrate(plan)
            self.assertFalse(source.is_symlink())
