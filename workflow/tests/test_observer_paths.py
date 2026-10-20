"""Synthetic workflow integration only; no simulation or server changes."""
from pathlib import Path
import base64
import contextlib
import csv
import hashlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

WORKFLOW = Path(__file__).resolve().parent.parent
REPOSITORY = WORKFLOW.parent
PRODUCTION = REPOSITORY / 'repo'
SCRATCH = Path(tempfile.gettempdir())
if os.name == 'nt':
    # Windows has no POSIX file locking. The server implementation remains
    # unchanged; this shim lets us test path/schema behavior on this host.
    sys.modules['fcntl'] = SimpleNamespace(flock=lambda *args: None, LOCK_EX=2, LOCK_NB=4)

sys.path.insert(0, str(WORKFLOW))
import workflow_runtime as runtime
import workflow_management as management
import workflow_isolation as isolation_stub


PNG = base64.b64decode('iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAusB9Wl6GaAAAAAASUVORK5CYII=')
MODEL = 'Fornax5028'
ID = runtime.CANONICAL_OBSERVER_ID

class ObserverWorkflowTests(unittest.TestCase):
    def setUp(self):
        SCRATCH.mkdir(exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(prefix='test_',dir=SCRATCH)
        self.root = Path(self.temp.name)
        self.model = self.root / MODEL
        (self.model / 'output').mkdir(parents=True)
        for name in ['mergegalaxy.par', 'mergecube.par', 'GZWJL.PARAM', 'Fornax.ini']:
            (self.model / name).write_text('synthetic fixture\n', encoding='utf-8')
        (self.model / 'output' / 'snapshot_300.hdf5').write_bytes(b'synthetic snapshot placeholder')
        (self.model / ('elinfo_' + MODEL + '.csv')).write_text('LEGACY_SENTINEL\n', encoding='utf-8')
        self.code = self.root / 'candidate_code'
        self.code.mkdir()
        shutil.copy2(PRODUCTION / 'canonical_observer.py', self.code / 'canonical_observer.py')
        (self.code / 'data_processing.py').write_text('# synthetic dependency\n', encoding='utf-8')
        columns = sorted(runtime.CANONICAL_SUMMARY_COLUMNS)
        getinfo = ("import csv,os\nfrom pathlib import Path\nfrom canonical_observer import canonical_csv_comment\n"
                   "root=Path(os.environ['DSPH_CANONICAL_OUTPUT_ROOT']); root.mkdir(parents=True,exist_ok=True)\n"
                   f"columns={columns!r}\n"
                   f"with (root/'elinfo_{MODEL}.csv').open('w',encoding='utf-8',newline='') as out:\n"
                   " out.write(canonical_csv_comment()); writer=csv.DictWriter(out,fieldnames=columns);writer.writeheader()\n"
                   " writer.writerow(dict.fromkeys(columns,1)|{'numsp':300})\n")
        plot = ("import os\nfrom pathlib import Path\n"
                "root=Path(os.environ['DSPH_CANONICAL_OUTPUT_ROOT'])/'output_img';root.mkdir(parents=True,exist_ok=True)\n"
                f"(root/'check_nature_300.png').write_bytes({PNG!r})\n")
        (self.code / 'GetInfo.py').write_text(getinfo, encoding='utf-8')
        (self.code / 'PlotFig.py').write_text(plot, encoding='utf-8')
        self.cohort = self.model / 'output' / 'analysis' / ID / 'cohort_ids.npz'
        self.cohort.parent.mkdir(parents=True)
        np.savez_compressed(self.cohort, particle_id=np.array([1,2,3], dtype=np.uint64))
        (self.root / 'workflow.json').write_text(json.dumps({'analysis_code': 'candidate_code'}), encoding='utf-8')
        self.args = SimpleNamespace(root=self.root, model=MODEL, name='canonical_fixture',
                                    snapshots='output', code=str(self.code), dry_run=False,
                                    range='300,300', processes=1, stage='all')

    def tearDown(self):
        self.assertTrue(self.root.resolve().is_relative_to(SCRATCH.resolve()))
        self.temp.cleanup()

    def prepare(self):
        with patch.dict(os.environ, {'DSPH_PREEXISTING_ID_PATH': str(self.cohort)}):
            if os.name == 'nt':
                def fixture_link(link, target, target_is_directory=False):
                    source=(link.parent/target).resolve()
                    self.assertTrue(source.is_relative_to(self.root.resolve()))
                    shutil.copy2(source,link)
                with patch.object(Path,'symlink_to',fixture_link):
                    self.assertEqual(runtime.prepare(self.args), 0)
            else:
                self.assertEqual(runtime.prepare(self.args), 0)
        directory = runtime.run_dir(self.args)
        manifest = json.loads((directory / 'manifest.json').read_text())
        return directory, manifest

    def run_canonical(self):
        directory, manifest = self.prepare()
        try:
            self.assertEqual(runtime.run(self.args), 0)
        except Exception:
            log=directory/'logs'/'info.log'
            if log.exists():
                print(log.read_text(encoding='utf-8',errors='replace'),file=sys.stderr)
            raise
        return directory, manifest

    def enable_tan(self):
        (self.code / 'canonical_display.py').write_text(
            "TAN_DISPLAY_CONTRACT_ID = 'heliocentric_tan_xi_eta_2p1_v1'\n"
            "DISPLAY_RELATIVE_DIR = 'display/' + TAN_DISPLAY_CONTRACT_ID\n", encoding='utf-8')
        (self.code / 'display_outputs.py').write_text('# synthetic renderer identity helper\n', encoding='utf-8')
        plot = ("import os,json,hashlib\nfrom pathlib import Path\n"
                "from canonical_display import TAN_DISPLAY_CONTRACT_ID,DISPLAY_RELATIVE_DIR\n"
                "canonical=Path(os.environ['DSPH_CANONICAL_OUTPUT_ROOT'])\n"
                "root=canonical/DISPLAY_RELATIVE_DIR;images=root/'output_img';images.mkdir(parents=True,exist_ok=True)\n"
                f"image=images/'check_nature_300.png';image.write_bytes({PNG!r})\n"
                "sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()\n"
                "module_sha=sha(Path(__file__).with_name('canonical_display.py'))\n"
                f"lines=(canonical/'elinfo_{MODEL}.csv').read_text().splitlines()\n"
                "identity=dict(item.split('=',1) for item in next(l for l in lines if l.startswith('# canonical_observer_id=')).split()[1:])\n"
                "identity.update(display_contract_id=TAN_DISPLAY_CONTRACT_ID,display_module_sha256=module_sha,relative_directory=DISPLAY_RELATIVE_DIR,render_source_sha256={p.name:sha(p) for p in sorted(Path(__file__).parent.glob('*.py'))})\n"
                "sidecar=dict(identity,png_sha256=sha(image),numsp=300,bytes=image.stat().st_size)\n"
                "image.with_name(image.name+'.display.json').write_text(json.dumps(sidecar))\n"
                "manifest=dict(identity,frames=[dict(numsp=300,bytes=image.stat().st_size,sha256=sha(image))],videos=[])\n"
                "(root/'display_manifest.json').write_text(json.dumps(manifest))\n")
        (self.code / 'PlotFig.py').write_text(plot, encoding='utf-8')

    def test_tan_run_pins_display_and_rejects_stale_frame_sidecar(self):
        self.enable_tan()
        directory, manifest = self.run_canonical()
        canonical = directory / MODEL / manifest['observation']['analysis_relative_dir']
        display = manifest['display']
        root = runtime.display_output_directory(canonical, display)
        self.assertEqual(display['display_module_sha256'], runtime.sha(directory/'code/canonical_display.py'))
        self.assertTrue((root/'output_img/check_nature_300.png').is_file())
        self.assertFalse((canonical/'output_img').exists())
        kwargs = dict(summary_path=canonical/f'elinfo_{MODEL}.csv', code=directory/'code')
        self.assertEqual(runtime.verify_display_outputs(canonical, display, [300], **kwargs)['videos'], [])
        sidecar = root/'output_img/check_nature_300.png.display.json'
        original=json.loads(sidecar.read_text())
        for changed in [dict(original,display_contract_id='raw'),dict(original,numsp=301),
                        dict(original,bytes=0),dict(original,observer_config_sha256='0'*64),
                        dict(original,render_source_sha256={})]:
            sidecar.write_text(json.dumps(changed))
            with self.assertRaisesRegex(ValueError, 'sidecar mismatch'):
                runtime.verify_display_outputs(canonical, display, [300], **kwargs)
        sidecar.write_text(json.dumps(original))
        with self.assertRaisesRegex(ValueError, 'missing, extra or duplicate'):
            runtime.verify_display_outputs(canonical, display, [301], **kwargs)

    def test_tan_results_export_only_manifest_videos_and_keep_flat_local_media(self):
        self.enable_tan()
        directory, manifest = self.run_canonical()
        canonical = directory / MODEL / manifest['observation']['analysis_relative_dir']
        root = runtime.display_output_directory(canonical, manifest['display'])
        video = root/(MODEL+'_v6_tan_xi_eta_2p1_v1.mp4')
        video.write_bytes(b'verified synthetic video')
        (canonical/'old_raw.mp4').write_bytes(b'old raw video')
        (root/'unlisted_tan_xi_eta_2p1_v1.mp4').write_bytes(b'unlisted video')
        (canonical/'observer_snapshot_300.json').write_text('{}')
        sidecar = root/'display_manifest.json'
        metadata = json.loads(sidecar.read_text())
        metadata['videos'] = [dict(name=video.name, bytes=video.stat().st_size, sha256=runtime.sha(video))]
        sidecar.write_text(json.dumps(metadata))
        result = management.production_results(self.root, name=self.args.name)
        names = {record['name'] for record in result['files']}
        self.assertEqual(names, {video.name, 'display_manifest.json', f'elinfo_{MODEL}.csv', 'observer_snapshot_300.json'})
        video_record = next(record for record in result['files'] if record['name']==video.name)
        self.assertEqual(video_record['source_relative_path'], manifest['display']['relative_directory']+'/'+video.name)
        def transfer(command, **kwargs):
            source=Path(command[-2].split(':',1)[1]); target=Path(command[-1])
            for option in command:
                if option.startswith('--include=/'):
                    filename=option.removeprefix('--include=/')
                    shutil.copy2(source/filename,target/filename)
        args=SimpleNamespace(root=self.root,name=self.args.name,dry_run=False)
        with patch.object(management,'ssh_config',return_value=({'root':'/remote','host':'fixture'},['ssh'])), \
             patch.object(isolation_stub,'remote_exec',return_value=SimpleNamespace(stdout=json.dumps(result)),create=True), \
             patch.object(management.subprocess,'run',side_effect=transfer), \
             patch.object(management.subprocess,'check_output',return_value=json.dumps({'streams':[{'nb_frames':'1','duration':'1'}]})):
            self.assertEqual(management.fetch_production(args,name=self.args.name),0)
        receipt=self.root/'workspace/receipts/production'/MODEL/self.args.name
        self.assertEqual((self.model/video.name).read_bytes(),video.read_bytes())
        self.assertTrue((self.model/f'elinfo_{MODEL}.csv').is_file())
        self.assertTrue((receipt/'display/display_manifest.json').is_file())
        self.assertTrue((receipt/'observer/observer_snapshot_300.json').is_file())
        self.assertFalse((self.model/'display_manifest.json').exists())

    def test_tan_contract_and_transfer_paths_reject_overrides(self):
        self.enable_tan()
        directory, manifest=self.run_canonical()
        canonical=directory/MODEL/manifest['observation']['analysis_relative_dir']
        display=dict(manifest['display'],relative_directory='../../raw')
        with self.assertRaisesRegex(ValueError,'inconsistent display'):
            runtime.display_output_directory(canonical,display)
        result=management.production_results(self.root,name=self.args.name)
        result['files'][0]['source_relative_path']='../../display_manifest.json'
        args=SimpleNamespace(root=self.root,name=self.args.name,dry_run=True)
        with patch.object(management,'ssh_config',return_value=({'root':'/remote','host':'fixture'},['ssh'])), \
             patch.object(isolation_stub,'remote_exec',return_value=SimpleNamespace(stdout=json.dumps(result)),create=True), \
             patch.object(management.subprocess,'run') as transfer:
            with self.assertRaisesRegex(ValueError,'source-relative'):
                management.fetch_production(args,name=self.args.name)
            transfer.assert_not_called()

    def test_tan_manifest_rejects_changed_physical_identity(self):
        self.enable_tan()
        directory, manifest=self.run_canonical()
        canonical=directory/MODEL/manifest['observation']['analysis_relative_dir']
        root=runtime.display_output_directory(canonical,manifest['display'])
        sidecar=root/'display_manifest.json'; metadata=json.loads(sidecar.read_text())
        metadata['observer_config_sha256']='0'*64
        sidecar.write_text(json.dumps(metadata))
        with self.assertRaisesRegex(ValueError,'identity mismatch'):
            management.production_results(self.root,name=self.args.name)
        metadata['observer_config_sha256']=runtime.canonical_csv_identity(canonical/f'elinfo_{MODEL}.csv')['observer_config_sha256']
        sidecar.write_text(json.dumps(metadata))
        for name in ['PlotFig.py','data_processing.py']:
            module=directory/'code'/name; saved=module.read_bytes()
            module.write_text('# changed rendering source\n')
            with self.assertRaisesRegex(ValueError,'identity mismatch: render_source_sha256'):
                management.production_results(self.root,name=self.args.name)
            module.write_bytes(saved)

    def test_tan_info_only_fetch_requires_explicit_complete_named_run(self):
        self.enable_tan()
        self.args.stage='info'
        directory, manifest=self.run_canonical()
        result=management.production_results(self.root,name=self.args.name)
        self.assertEqual([entry['name'] for entry in result['files']],[f'elinfo_{MODEL}.csv'])
        self.assertEqual(result['display'],manifest['display'])
        status_path=directory/'status.json'; original=json.loads(status_path.read_text())
        for state in [dict(original,status='running'),dict(original,status='failed'),
                      dict(original,stages=[]),dict(original,stages=[dict(stage='plot',status='complete')]),
                      dict(original,stages=[dict(stage='info',status='failed')])]:
            status_path.write_text(json.dumps(state))
            with self.assertRaisesRegex(ValueError,'display manifest is missing'):
                management.production_results(self.root,name=self.args.name)
        status_path.write_text(json.dumps(original))
        canonical=directory/MODEL/manifest['observation']['analysis_relative_dir']
        model_output=self.model/'output/analysis'/ID
        shutil.copy2(canonical/f'elinfo_{MODEL}.csv',model_output/f'elinfo_{MODEL}.csv')
        with self.assertRaisesRegex(ValueError,'display manifest is missing'):
            management.production_results(self.root,model=MODEL)

    def test_tan_named_fetch_rejects_manifest_from_other_selected_range(self):
        self.enable_tan()
        directory, manifest=self.run_canonical()
        manifest['snapshots'][MODEL+'/output/snapshot_301.hdf5']={}
        (directory/'manifest.json').write_text(json.dumps(manifest))
        status=json.loads((directory/'status.json').read_text())
        status['range']='301,301'
        (directory/'status.json').write_text(json.dumps(status))
        with self.assertRaisesRegex(ValueError,'missing, extra or duplicate'):
            management.production_results(self.root,name=self.args.name)

    def test_prepare_and_run_new_paths_and_pinned_cohort(self):
        directory, manifest = self.run_canonical()
        observation = manifest['observation']
        self.assertEqual(observation['canonical_observer_id'], ID)
        self.assertEqual(observation['cohort_definition'], 'selected_then_frozen')
        self.assertEqual(runtime.sha(directory / observation['cohort_relative_path']), runtime.sha(self.cohort))
        output = directory / MODEL / 'output' / 'analysis' / ID
        self.assertTrue((output / ('elinfo_' + MODEL + '.csv')).is_file())
        self.assertTrue((output / 'output_img' / 'check_nature_300.png').is_file())
        self.assertFalse((directory / MODEL / ('elinfo_' + MODEL + '.csv')).exists())
        self.assertEqual((self.model / ('elinfo_' + MODEL + '.csv')).read_text(), 'LEGACY_SENTINEL\n')
        state = json.loads((directory / 'status.json').read_text())
        self.assertEqual(state['canonical_observer_id'], ID)

    def test_wrong_observer_schema_and_cohort_are_rejected(self):
        directory, manifest = self.run_canonical()
        observation = manifest['observation']
        summary = directory / MODEL / observation['analysis_relative_dir'] / ('elinfo_' + MODEL + '.csv')
        kwargs = dict(observation=observation,code=directory/'code',cohort_path=directory/observation['cohort_relative_path'])
        original = summary.read_text()
        identity = runtime.canonical_csv_identity(summary)
        for replacement in [original.replace(ID, 'wrong_observer'),
                            original.replace('observer_config_sha256=', 'missing_config='),
                            original.replace(identity['observer_config_sha256'], '0'*64),
                            original.replace('walker_R5_nstar', 'wrong_schema_column')]:
            summary.write_text(replacement)
            with self.assertRaises((ValueError, subprocess.CalledProcessError)):
                runtime.verify_summary(summary, [300], **kwargs)
        summary.write_text(original)
        kwargs['cohort_path'].write_bytes(b'changed frozen fixture')
        with self.assertRaises(ValueError): runtime.verify_summary(summary, [300], **kwargs)

    def test_new_prepare_does_not_allow_legacy_code_or_missing_cohort(self):
        (self.code/'canonical_observer.py').unlink()
        with self.assertRaises(ValueError): runtime.prepare(self.args)
        self.assertFalse(runtime.run_dir(self.args).exists())
        shutil.copy2(PRODUCTION/'canonical_observer.py',self.code/'canonical_observer.py')
        self.cohort.unlink()
        with self.assertRaises(ValueError): runtime.prepare(self.args)
        self.assertFalse(runtime.run_dir(self.args).exists())

    def test_model_fetch_does_not_fallback_to_old_csv(self):
        with self.assertRaises(ValueError): management.production_results(self.root, model=MODEL)

    def test_default_model_fetch_accepts_only_new_contract_and_sidecars(self):
        directory, manifest = self.run_canonical()
        output = self.model / 'output' / 'analysis' / ID
        source = directory / MODEL / manifest['observation']['analysis_relative_dir']
        shutil.copy2(source/('elinfo_'+MODEL+'.csv'), output/('elinfo_'+MODEL+'.csv'))
        (output/'observer_snapshot_300.json').write_text('{}')
        result = management.production_results(self.root, model=MODEL)
        self.assertEqual(result['canonical_observer_id'], ID)
        self.assertEqual({r['name'] for r in result['files']}, {'elinfo_'+MODEL+'.csv','observer_snapshot_300.json'})

    def test_server_results_stay_versioned_and_local_media_use_model_root(self):
        directory, manifest = self.run_canonical()
        result = management.production_results(self.root, name=self.args.name)
        self.assertEqual(result['canonical_observer_id'], ID)
        self.assertEqual(result['local_relative_directory'], 'output/analysis/'+ID)
        self.assertEqual(Path(result['directory']).parts[-3:],('output','analysis',ID))
        self.assertEqual([f['name'] for f in result['files']], ['elinfo_'+MODEL+'.csv'])
        args = SimpleNamespace(root=self.root, name=self.args.name, dry_run=True)
        with patch.object(management, 'ssh_config', return_value=({'root':'/remote','host':'fixture'},['ssh'])), \
             patch.object(isolation_stub, 'remote_exec', return_value=SimpleNamespace(stdout=json.dumps(result)), create=True), \
             patch.object(management.subprocess, 'run') as rsync:
            self.assertEqual(management.fetch_production(args,name=self.args.name),0)
        self.assertEqual(rsync.call_args.args[0][-1],str(self.model)+'/')

    def test_fetch_keeps_observer_sidecars_in_receipts(self):
        directory, manifest = self.run_canonical()
        output = directory / MODEL / manifest['observation']['analysis_relative_dir']
        (output / 'observer_snapshot_300.json').write_text('{}')
        result = management.production_results(self.root, name=self.args.name)
        args = SimpleNamespace(root=self.root, name=self.args.name, dry_run=False)
        def transfer(command, **kwargs):
            target = Path(command[-1])
            for option in command:
                if option.startswith('--include=/'):
                    filename = option.removeprefix('--include=/')
                    shutil.copy2(output / filename, target / filename)
        with patch.object(management, 'ssh_config', return_value=({'root':'/remote','host':'fixture'},['ssh'])), \
             patch.object(isolation_stub, 'remote_exec', return_value=SimpleNamespace(stdout=json.dumps(result)), create=True), \
             patch.object(management.subprocess, 'run', side_effect=transfer):
            self.assertEqual(management.fetch_production(args,name=self.args.name),0)
        receipt_root = self.root / 'workspace/receipts/production' / MODEL
        metadata = receipt_root / self.args.name / 'observer'
        self.assertTrue((self.model / ('elinfo_'+MODEL+'.csv')).is_file())
        self.assertTrue((metadata / 'observer_snapshot_300.json').is_file())
        self.assertFalse((self.model / 'observer_snapshot_300.json').exists())
        receipt = json.loads((receipt_root / (self.args.name+'.json')).read_text())
        self.assertEqual(receipt['local_relative_directory'], '.')
        self.assertEqual(receipt['local_directory'], str(self.model))
        self.assertEqual(receipt['local_metadata_directory'], str(metadata))

    def test_named_frozen_legacy_run_is_explicit_and_manifest_unchanged(self):
        run = self.root/'workspace'/'analysis'/'historical'
        (run/MODEL).mkdir(parents=True)
        old_manifest={'model':MODEL,'snapshots':{},'code_sha256':{},'inputs':{}}
        (run/'manifest.json').write_text(json.dumps(old_manifest))
        (run/MODEL/('elinfo_'+MODEL+'.csv')).write_text('numsp,star_mass,distance,sigma\n300,1,1,1\n')
        before=(run/'manifest.json').read_bytes()
        result=management.production_results(self.root,name='historical')
        self.assertEqual(result['canonical_observer_id'],'legacy_unversioned')
        self.assertEqual(result['local_relative_directory'],'.')
        self.assertEqual((run/'manifest.json').read_bytes(),before)

    def test_manifest_path_or_schema_override_is_rejected(self):
        for field,value in [('analysis_relative_dir','../../legacy'),('schema_id','legacy_elinfo'),('canonical_observer_id','legacy')]:
            observation={'canonical_observer_id':ID,'schema_id':runtime.CANONICAL_SCHEMA_ID,'analysis_relative_dir':'output/analysis/'+ID}
            observation[field]=value
            with self.assertRaises(ValueError): runtime.observation_output_directory(self.model,observation)

if __name__=='__main__':
    with contextlib.redirect_stdout(io.StringIO()):
        unittest.main(verbosity=2)
