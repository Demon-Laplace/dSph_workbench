import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import workflow_models as models

SOURCE = Path(__file__).resolve().parents[2] / 'local_tools/ic_generation/catalog.py'
spec = importlib.util.spec_from_file_location('ic_catalog', SOURCE)
catalog = importlib.util.module_from_spec(spec)
spec.loader.exec_module(catalog)

INI = '''[Density_gas]
mass=0.013
scaleRadius=0.9
scaleHeight=-0.45
[Component_gas]
nbody=2
ickappa=0.85
[Density_stellar]
mass=0.006
scaleRadius=0.52
axisRatioZ=0.55
[Component_stellar]
nbody=3
beta=-0.1
'''


class ICCatalogTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.repo = self.root / 'repo'
        self.remote = self.root / 'remote.git'
        self.repo.mkdir()
        self.git(self.repo, 'init', '-b', 'main')
        self.git(self.repo, 'config', 'user.name', 'IC tests')
        self.git(self.repo, 'config', 'user.email', 'ic@example.invalid')
        self.git(self.root, 'init', '--bare', str(self.remote))
        self.git(self.repo, 'remote', 'add', 'origin', str(self.remote))
        (self.repo / 'README').write_text('fixture')
        self.git(self.repo, 'add', 'README')
        self.git(self.repo, 'commit', '-m', 'baseline')
        self.git(self.repo, 'push', '-u', 'origin', 'main')

    def git(self, root, *args):
        return subprocess.check_output(['git', '-C', str(root), *args], stderr=subprocess.PIPE, text=True)

    def build(self, name='run1'):
        run = self.root / name
        run.mkdir()
        (run / 'parameters.ini').write_text(INI)
        (run / 'Fornax.ini').write_text('2 0 3 0\nparticle fixture\n')
        build = dict(state='complete', source_ini='IC_Fornax7003.ini', seed=42,
                     build_id=name, generated_utc='2026-10-12T00:00:00+00:00',
                     particle_ic='Fornax.ini', ini_sha256=catalog.sha(run / 'parameters.ini'),
                     particle_sha256=catalog.sha(run / 'Fornax.ini'))
        (run / 'build.json').write_text(json.dumps(build))
        return run

    def rows(self):
        with (self.repo / catalog.TABLE).open(newline='') as stream:
            return list(csv.DictReader(stream))

    def test_units_signed_height_and_multiple_components(self):
        ini = self.root / 'params.ini'
        ini.write_text(INI + '[Density_gas_2]\nmass=0.001\nscaleRadius=1\n'
                       '[Component_gas_2]\nnbody=4\n')
        row = catalog.parameters(ini, 'Fornax7003')
        self.assertEqual(row['gas_mass_Msun'], '140000000')
        self.assertEqual(row['gas_particles'], '6')
        self.assertEqual(row['gas_scaleHeight_kpc'], '-0.45')
        self.assertEqual(row['gas_scaleRadius_kpc'], 'gas=0.9; gas_2=1')
        self.assertEqual(row['stellar_mass_Msun'], '60000000')

    def test_historical_missing_mismatched_and_frozen_parameters(self):
        (self.root / 'IniFile').mkdir()
        (self.root / 'IniFile/IC_Fornax1.ini').write_text(INI)
        self.assertEqual(catalog.historical_rows(self.root)[0]['status'], 'parameters_only')
        model = self.root / 'Fornax1'
        model.mkdir()
        (model / 'IC_Fornax1.ini').write_text(INI.replace('0.013', '0.01'))
        (model / 'Fornax.ini').write_text('1 0 3 0\n')
        rows = catalog.historical_rows(self.root)
        self.assertEqual([r['status'] for r in rows], ['existing_ic_check_needed', 'parameters_only'])
        self.assertEqual(rows[0]['gas_mass_Msun'], '100000000')

    def test_registration_pushes_once_and_preserves_peer_update(self):
        run = self.build()
        catalog.register(run, self.repo)
        first = self.git(self.repo, 'rev-parse', 'HEAD')
        catalog.register(run, self.repo)
        self.assertEqual(first, self.git(self.repo, 'rev-parse', 'HEAD'))
        self.assertEqual(len(self.rows()), 1)
        peer = self.root / 'peer'
        self.git(self.root, 'clone', '-b', 'main', str(self.remote), str(peer))
        self.git(peer, 'config', 'user.name', 'IC tests')
        self.git(peer, 'config', 'user.email', 'ic@example.invalid')
        catalog.update(peer, [catalog.generated_row(self.build('peer-run'))])
        catalog.register(self.build('run2'), self.repo)
        self.assertEqual(len(self.rows()), 3)
        self.assertEqual(self.git(self.repo, 'rev-parse', 'HEAD'), self.git(self.remote, 'rev-parse', 'main'))
        self.assertEqual(self.git(self.repo, 'status', '--porcelain'), '')
        self.assertEqual(json.loads((run / 'build.json').read_text())['catalog_sync']['state'], 'pushed')

    def test_push_failure_preserves_completed_ic_and_retries(self):
        hook = self.remote / 'hooks/pre-receive'
        hook.write_text('#!/bin/sh\nexit 1\n')
        hook.chmod(0o755)
        run = self.build()
        before = catalog.sha(run / 'Fornax.ini')
        with self.assertRaises(subprocess.CalledProcessError):
            catalog.register(run, self.repo)
        receipt = json.loads((run / 'build.json').read_text())
        self.assertEqual(receipt['state'], 'complete')
        self.assertEqual(receipt['catalog_sync']['state'], 'pending')
        self.assertEqual(before, catalog.sha(run / 'Fornax.ini'))
        hook.unlink()
        catalog.register(run, self.repo)
        self.assertEqual(len(self.rows()), 1)

    def test_failed_or_changed_data_cannot_enter_catalogue(self):
        run = self.build()
        (run / 'Fornax.ini').write_text('changed')
        with self.assertRaisesRegex(ValueError, 'changed'):
            catalog.register(run, self.repo)
        receipt = json.loads((run / 'build.json').read_text())
        receipt['state'] = 'failed'
        (run / 'build.json').write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError, 'successful'):
            catalog.register(run, self.repo)
        self.assertFalse((self.repo / catalog.TABLE).exists())

    def test_workflow_uses_configured_sibling_and_refuses_escape(self):
        analysis = self.root / 'analysis'
        analysis.mkdir()
        ic = self.root / 'ic'
        (ic / 'IniFile').mkdir(parents=True)
        (ic / 'IniFile/IC_Fornax7003.ini').write_text(INI)
        (analysis / 'workflow.json').write_text(json.dumps({'ic': {'root': '../ic'}}))
        args = SimpleNamespace(root=analysis, repo=SOURCE.parents[2], ini=None,
                               model='Fornax7003', name='new-run', python='python3', threads=2, seed=42)
        with patch.object(models.subprocess, 'run', return_value=SimpleNamespace(returncode=0)) as run:
            self.assertEqual(models.build_ic(args), 0)
            command = run.call_args.args[0]
            self.assertEqual(command[command.index('--output') + 1], str(ic / 'generated/new-run'))
        (ic / 'generated').symlink_to(analysis, target_is_directory=True)
        with self.assertRaises(ValueError):
            models.build_ic(args)


if __name__ == '__main__':
    unittest.main()
