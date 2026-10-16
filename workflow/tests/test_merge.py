import json
import os
from pathlib import Path
import struct
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import workflow_merge as merge


def fake_binary(ascii_path, binary):
    with Path(ascii_path).open() as f:
        counts = [int(x) for x in next(f).split()]
        data = np.loadtxt(f, ndmin=2)
    header = bytearray(256)
    for offset in (0, 96): struct.pack_into('<6I', header, offset, *counts, 0, 0)
    blocks = [header, data[:, 1:4].astype('<f4').tobytes(), data[:, 4:7].astype('<f4').tobytes(),
              data[:, 0].astype('<i4').tobytes(), data[:, 7].astype('<f4').tobytes()]
    if counts[0]: blocks += [data[:counts[0], 8].astype('<f4').tobytes(), bytes(4*counts[0])]
    with Path(binary).open('wb') as f:
        for block in blocks:
            f.write(struct.pack('<I', len(block))); f.write(block); f.write(struct.pack('<I', len(block)))


class MergeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.build = tempfile.TemporaryDirectory()
        cls.executable = Path(cls.build.name) / 'ic_merge'
        subprocess.run(['gcc', '-O2', '-std=c99', '-ffp-contract=off', '-fno-fast-math', '-Wall', '-Wextra',
                        str(Path(merge.__file__).with_name('ic_merge.c')), '-lm', '-o', str(cls.executable)], check=True)

    @classmethod
    def tearDownClass(cls): cls.build.cleanup()

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.model = self.root / 'Fornax9999'; self.model.mkdir()
        self.work = self.root / 'work'; self.work.mkdir()
        self.fixture()

    def fixture(self, seed=1, galaxy=None, cube=None, gas=True):
        rng = np.random.default_rng(seed)
        for name, n in [('mw.ini', 101), ('dwarf.ini', 79), ('cube.ini', 33)]:
            types = rng.integers(0 if gas else 1, 4, n)
            values = rng.normal(size=(n, 8)) * 7
            values[:, 6] = rng.uniform(1e-7, 1e-4, n); values[:, 7] = rng.uniform(0, 0.1, n)
            with (self.model / name).open('w') as f:
                f.write(' '.join(map(str, np.bincount(types, minlength=4)))+'\n')
                for i in rng.permutation(n):
                    f.write(f'{i} '+' '.join(format(v, '.17g') for v in values[i])+f' {types[i]}\n')
        g = galaxy if galaxy is not None else [-192.3313,164.1997,-183.42,38.24,-89.24,3.28,90,0,0,151.354,37.547,-46.683]
        c = cube if cube is not None else [500,500,500,0,0,0,0,0,0,0,0,0]
        (self.model/'mergegalaxy.par').write_text('mw.ini dwarf.ini\n'+' '.join(map(str,g))+'\n\nobj=mergegalaxy\n')
        (self.model/'mergecube.par').write_text('cube.ini mergegalaxy.dat\n'+' '.join(map(str,c))+'\n')
        (self.model/'GZWJL.PARAM').write_text('InitCondFile mergecube.IC_000\nOutputDir output\n')
        (self.root/'workflow.json').write_text(json.dumps({'ic_merge': {'converter':['fake-converter','{ascii}','{binary}']}}))

    def onepass(self, name='new.dat'):
        return subprocess.run([str(self.executable),'mergegalaxy.par','mergecube.par',str(self.work/name)],
                              cwd=self.model,capture_output=True,text=True)

    def test_counts_stable_order_ids_and_cube_unchanged(self):
        self.assertEqual(self.onepass().returncode,0)
        new = np.loadtxt(self.work/'new.dat',skiprows=1)
        self.assertEqual(len(np.unique(new[:,0])),213)
        self.assertTrue(np.all(np.diff(new[:,9])>=0))
        for row in np.loadtxt(self.model/'cube.ini',skiprows=1):
            result = new[new[:,0]==row[0]][0]
            expected = np.array([float(format(float(v),'.7g')) for v in row[1:9].astype(np.float32)])
            np.testing.assert_array_equal(result[1:9],expected)

    def test_duplicate_ids_and_truncated_input_refused(self):
        p = self.model/'cube.ini'; lines=p.read_text().splitlines()
        fields=lines[1].split(); fields[0]=lines[2].split()[0]; lines[1]=' '.join(fields)
        p.write_text('\n'.join(lines)+'\n')
        self.assertNotEqual(self.onepass().returncode,0); self.assertFalse((self.work/'new.dat').exists())
        p.write_text('\n'.join(lines[:-1])+'\n')
        self.assertNotEqual(self.onepass().returncode,0)

    def test_refuses_unrelated_input_or_existing_ic(self):
        p=self.model/'mergecube.par'; p.write_text(p.read_text().replace('mergegalaxy.dat','other.dat'))
        with self.assertRaises(ValueError): merge.plan(self.root,'Fornax9999')
        self.fixture(); (self.model/'mergecube.IC_000').write_bytes(b'keep')
        with self.assertRaises(ValueError): merge.plan(self.root,'Fornax9999')
        self.assertEqual((self.model/'mergecube.IC_000').read_bytes(),b'keep')

    def test_binary_all_fields_and_truncation(self):
        self.assertEqual(self.onepass().returncode,0)
        ascii_path,binary=self.work/'new.dat',self.work/'new.IC_000'
        fake_binary(ascii_path,binary)
        self.assertEqual(merge.validate_binary(ascii_path,binary)['particles'],213)
        with binary.open('r+b') as f: f.seek(268); f.write(struct.pack('<f',12345))
        with self.assertRaises(ValueError): merge.validate_binary(ascii_path,binary)
        fake_binary(ascii_path,binary)
        with binary.open('r+b') as f: f.truncate(binary.stat().st_size-1)
        with self.assertRaises(ValueError): merge.validate_binary(ascii_path,binary)

    def test_cleanup_only_after_validated_publication(self):
        real_run=subprocess.run
        def run(cmd,**kw):
            if cmd[0]=='fake-converter':
                fake_binary(cmd[1],cmd[2]); return subprocess.CompletedProcess(cmd,0)
            return real_run(cmd,**kw)
        legacy=self.model/'mergegalaxy.dat'; legacy.write_text('unrelated existing data')
        with patch.object(merge.subprocess,'run',side_effect=run), patch.object(merge.shutil,'which',return_value='/bin/true'):
            receipt=merge.prepare(self.root,'Fornax9999',self.work)
        self.assertTrue(receipt['intermediate_deleted']); self.assertFalse((self.work/'merged.dat').exists())
        self.assertTrue((self.model/'mergecube.IC_000').is_file()); self.assertTrue(legacy.is_file())

    def test_failed_converter_preserves_intermediate(self):
        (self.root/'workflow.json').write_text(json.dumps({'ic_merge': {'converter':['/bin/false','{ascii}','{binary}']}}))
        with self.assertRaises(subprocess.CalledProcessError): merge.prepare(self.root,'Fornax9999',self.work)
        self.assertTrue((self.work/'merged.dat').is_file()); self.assertFalse((self.model/'mergecube.IC_000').exists())

    def test_queued_input_change_refused(self):
        identity=merge.plan(self.root,'Fornax9999')[-1]
        with (self.model/'dwarf.ini').open('a') as f: f.write('\n')
        with self.assertRaises(ValueError): merge.prepare(self.root,'Fornax9999',self.work,identity)
        self.assertFalse((self.work/'merged.dat').exists())

    def test_submit_once_and_freeze_settings(self):
        from types import SimpleNamespace
        cfg = {'ic_merge': {'converter':['/bin/true','{ascii}','{binary}'], 'sbatch':{'ntasks':1}, 'setup':'module load python'}}
        (self.root/'workflow.json').write_text(json.dumps(cfg))
        args = SimpleNamespace(root=self.root,model='Fornax9999',name='v1',submit=True)
        with patch('workflow_cleanup.active_jobs',return_value=[]), patch.object(merge.subprocess,'run',return_value=subprocess.CompletedProcess([],0,'123;cluster\n','')) as run:
            self.assertEqual(merge.command(args),0)
            self.assertEqual(merge.command(args),0)
            self.assertEqual(run.call_count,1)
        work=self.root/'workspace/merge-ic/Fornax9999/v1'
        self.assertEqual(json.loads((work/'request.json').read_text())['settings'],cfg['ic_merge'])
        self.assertTrue((work/'ic_merge.c').is_file())
        self.assertFalse((self.model/'mergecube.IC_000').exists())

    def test_active_model_blocks_preparation(self):
        from types import SimpleNamespace
        args=SimpleNamespace(root=self.root,model='Fornax9999',name='v1',submit=True)
        with patch('workflow_cleanup.active_jobs',return_value=[('99',self.model)]):
            with self.assertRaises(ValueError): merge.command(args)

    @unittest.skipUnless(os.environ.get('IC_LEGACY_MERGER') and os.environ.get('IC_LEGACY_CONVERTER'),'Set reference binaries')
    def test_reference_equivalence(self):
        merger=os.environ['IC_LEGACY_MERGER']; converter=os.environ['IC_LEGACY_CONVERTER']
        rng=np.random.default_rng(20261016)
        cases=[(None,None),([-12.7,220.1,-94.9,-23.7,-67.2,-68.4,90,0,0,83,45,0],[1000,1000,1000,0,0,0,0,0,0,0,0,0]),([0]*12,[0]*12)]
        cases += [(rng.uniform(-360,360,12).tolist(),rng.uniform(-360,360,12).tolist()) for _ in range(17)]
        for index,(g,c) in enumerate(cases):
            with self.subTest(index=index):
                self.fixture(index+1,g,c,gas=index!=19)
                for par,out in [('mergegalaxy.par','mergegalaxy.dat'),('mergecube.par','mergecube.dat')]:
                    subprocess.run([merger,par,out],cwd=self.model,check=True,capture_output=True)
                name=f'new-{index}.dat'; self.assertEqual(self.onepass(name).returncode,0)
                self.assertEqual((self.model/'mergecube.dat').read_bytes(),(self.work/name).read_bytes())
                for ap,bp in [(self.model/'mergecube.dat',self.work/'old.IC_000'),(self.work/name,self.work/'new.IC_000')]:
                    subprocess.run([converter,str(ap),str(bp)],check=True,capture_output=True)
                    merge.validate_binary(ap,bp)
                self.assertEqual((self.work/'old.IC_000').read_bytes(),(self.work/'new.IC_000').read_bytes())


if __name__=='__main__': unittest.main()
