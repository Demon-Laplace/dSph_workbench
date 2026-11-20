import json
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import workflow_isolation as isolation
import workflow_local_isolation as local
import workflow_management as management
import workflow_runtime as runtime


class NativeTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.ic=self.root/'Fornax7004/ic/build';self.ic.mkdir(parents=True)
        (self.ic/'Fornax.ini').write_text('1 0 0 0\n0 0 0 0 0 0 0 1 0 0\n')
        (self.ic/'parameters.ini').write_text('[Potential]\n')
        runtime.write_json(self.ic/'build.json',{'state':'complete','source_ini':'IC_Fornax7004.ini',
            'particle_sha256':isolation.sha(self.ic/'Fornax.ini'),'ini_sha256':isolation.sha(self.ic/'parameters.ini')})
        template=self.root/'template';template.mkdir();(template/'TREECOOL').write_text('cooling')
        (template/'GZWJL.PARAM').write_text('InitCondFile old\nOutputDir output\nResubmitOn 1\nNumFilesWrittenInParallel 4\nTimeMax 6\n')
        (self.root/'solver').write_text('fixture');(self.root/'converter').write_text('fixture')
        runtime.write_json(self.root/'workflow.json',{'isolation_local':{'template':'template','solver':'solver',
            'converter':'converter','python':sys.executable,'command':['solver','{parameters}'],'processes':2}})
        self.args=SimpleNamespace(root=self.root,repo=self.root,model='Fornax7004',ic_run='build',name='test',
                                  processes=2,time_max=None,dry_run=False)

    def test_freezes_relative_solver_paths_and_never_duplicates(self):
        with patch.object(local.subprocess,'Popen') as popen:
            popen.return_value.pid=123
            local.start(self.args);local.start(self.args)
            popen.assert_called_once()
        run=local.directory(self.root,'Fornax7004','test')
        params=(run/'GZWJL.PARAM').read_text()
        self.assertIn('InitCondFile Fornax.IC_000',params)
        self.assertIn('OutputDir output',params)
        self.assertTrue((run/'Fornax.ini').is_file())
        (self.ic/'Fornax.ini').write_text('changed outside workflow')
        with patch.object(local.subprocess,'Popen') as again:
            self.assertEqual(local.start(self.args),0)
            again.assert_not_called()

    def test_dry_run_does_not_create_case_and_invalid_limit_is_rejected(self):
        self.args.dry_run=True;local.start(self.args)
        self.assertFalse((self.root/'workspace').exists())
        self.args.time_max=-1
        with self.assertRaises(ValueError):local.start(self.args)
        self.assertFalse((self.root/'workspace').exists())

    def test_declared_new_build_is_supported_without_digests(self):
        receipt=self.ic/'build.json'
        build=json.loads(receipt.read_text())
        build.update(build_id='explicit-build-version',particle_sha256=None,ini_sha256=None)
        runtime.write_json(receipt,build)
        with patch.object(local.subprocess,'Popen') as popen, \
             patch.object(isolation,'sha',side_effect=AssertionError('no hashing')):
            popen.return_value.pid=123
            local.start(self.args);local.start(self.args)
            popen.assert_called_once()
        identity=json.loads((local.directory(self.root,'Fornax7004','test')/'identity.json').read_text())
        self.assertEqual(identity['build_id'],'explicit-build-version')

    def test_legacy_run_keeps_its_identity_without_rereading_data(self):
        with patch.object(local.subprocess,'Popen') as popen:
            popen.return_value.pid=123;local.start(self.args)
        run=local.directory(self.root,'Fornax7004','test')
        identity=json.loads((run/'identity.json').read_text())
        identity.pop('build_id')
        identity['particle_sha256']=json.loads((self.ic/'build.json').read_text())['particle_sha256']
        runtime.write_json(run/'identity.json',identity)
        before=(run/'identity.json').read_text()
        with patch.object(local.subprocess,'Popen') as popen, \
             patch.object(isolation,'sha',side_effect=AssertionError('no hashing')):
            local.start(self.args);popen.assert_not_called()
        self.assertEqual((run/'identity.json').read_text(),before)

    def test_normal_execution_does_not_recheck_frozen_file_contents(self):
        with patch.object(local.subprocess,'Popen') as popen:
            popen.return_value.pid=123;local.start(self.args)
        run=local.directory(self.root,'Fornax7004','test')
        (run/'Fornax.ini').write_text('corrupted after freeze')
        process = unittest.mock.MagicMock()
        process.__enter__.return_value.wait.return_value=0
        with patch.object(local.subprocess,'run'), \
             patch.object(local.subprocess,'Popen',return_value=process), \
             patch.object(local.os,'nice'), \
             patch.object(isolation,'sha',side_effect=AssertionError('no hashing')):
            local.execute(self.root,run)
        self.assertEqual(json.loads((run/'status.json').read_text())['state'],'complete')

    def test_eta_uses_measured_progress_after_startup(self):
        run=self.root/'run';run.mkdir();(run/'solver.log').write_text('Sync-Point 12, Time: 0.2\n')
        samples=[{'elapsed_seconds':100,'simulation_time':0.1}]
        with patch.object(local.time,'time',return_value=200):local.progress(run,{'started_epoch':0,'target_time':6},samples)
        value=json.loads((run/'progress.json').read_text())
        self.assertAlmostEqual(value['estimated_remaining_hours'],5.8/0.001/3600)


class ProductionFetchTests(unittest.TestCase):
    def test_video_and_table_go_to_model_root_without_images_or_inputs(self):
        import workflow_isolation as isolation
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'remote';source.mkdir();target=root/'Fornax1';target.mkdir()
            (target/'Fornax.ini').write_text('original IC')
            for name,data in [('movie.mp4',b'video'),('elinfo_Fornax1.csv',b'summary'),
                              ('Fornax.ini',b'wrong IC'),('frame.png',b'frame')]: (source/name).write_bytes(data)
            (source/'output_img').mkdir();(source/'output_img/frame.png').write_bytes(b'frame')
            records=[{'name':p.name,'bytes':p.stat().st_size,'sha256':isolation.sha(p)}
                     for p in [source/'movie.mp4',source/'elinfo_Fornax1.csv']]
            result={'model':'Fornax1','directory':'/data/model','files':records,'expected_frames':3}
            args=SimpleNamespace(root=root,model='Fornax1',dry_run=False)
            real_run=subprocess.run
            def transfer(command,**kwargs):
                command=list(command);command[-2]=str(source)+'/'
                return real_run(command,**kwargs)
            with patch.object(management,'ssh_config',return_value=({'root':'/data','host':'server'},['ssh'])), \
                 patch.object(isolation,'remote_exec',return_value=subprocess.CompletedProcess([],0,json.dumps(result),'')), \
                 patch.object(management.subprocess,'run',side_effect=transfer), \
                 patch.object(management.subprocess,'check_output',return_value=json.dumps({'streams':[{'nb_frames':'3'}]})):
                management.fetch_model(args)
            self.assertEqual(set(p.name for p in target.iterdir()),{'movie.mp4','elinfo_Fornax1.csv','Fornax.ini'})
            self.assertEqual((target/'Fornax.ini').read_text(),'original IC')
            self.assertTrue((root/'workspace/receipts/production/Fornax1/model.json').exists())
