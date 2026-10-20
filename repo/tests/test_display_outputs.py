"""Regression checks preventing reuse or export of unbound display media."""
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import display_outputs as outputs


class DisplayOutputTests(unittest.TestCase):
    def test_frame_requires_matching_identity_and_current_bytes(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'frame.png'
            path.write_bytes(b'original image')
            identity = {'display_contract_id': 'TAN', 'observer_config_sha256': 'physical'}
            self.assertFalse(outputs.frame_is_current(path, identity))
            outputs.record_frame(path, identity, 285)
            self.assertTrue(outputs.frame_is_current(path, identity))
            self.assertFalse(outputs.frame_is_current(path, dict(identity, observer_config_sha256='other')))
            path.write_bytes(b'changed image')
            self.assertFalse(outputs.frame_is_current(path, identity))

    def test_manifest_exports_only_verified_selected_media(self):
        import json
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            image = root / 'output_img/check_nature_285.png'
            image.parent.mkdir()
            image.write_bytes(b'frame')
            identity = {'display_contract_id': 'TAN'}
            outputs.record_frame(image, identity, 285)
            (root / 'unrelated.mp4').write_bytes(b'old video')
            outputs.record_display_manifest(root, identity, [285])
            manifest = json.loads((root / 'display_manifest.json').read_text())
            self.assertEqual(manifest['videos'], [])
            self.assertEqual([row['numsp'] for row in manifest['frames']], [285])
            image.write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'Missing or different-display'):
                outputs.record_display_manifest(root, identity, [285])


if __name__ == '__main__':
    unittest.main()
