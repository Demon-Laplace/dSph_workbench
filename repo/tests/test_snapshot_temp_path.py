"""Runtime markers stay in the configured analysis scratch directory."""
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from snapshot_io import SnapshotLoaderMixin


class SnapshotTempPathTests(unittest.TestCase):
    def test_auto_boxsize_marker_respects_analysis_tmpdir(self):
        with tempfile.TemporaryDirectory() as folder:
            with patch.dict(os.environ, TMPDIR=folder, DSPH_RUN_ID='tan-path-regression'):
                SnapshotLoaderMixin._announce_auto_boxsize_once(100.)
                SnapshotLoaderMixin._announce_auto_boxsize_once(100.)
            self.assertEqual([p.name for p in Path(folder).iterdir()],
                             ['dsph_boxsize_tan-path-regression.flag'])


if __name__ == '__main__':
    unittest.main()
