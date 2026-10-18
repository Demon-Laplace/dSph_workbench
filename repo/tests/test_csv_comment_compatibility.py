"""Regress actual GetInfo comment generation against the actual legacy reader.

Portable: defaults to the repo containing tests/, with --source-dir override.
No particle snapshot reads and no production writes.
"""
from pathlib import Path
import argparse, ast, copy, hashlib, importlib.util, json, os, sys, tempfile, unittest
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument('--source-dir',type=Path,default=Path(__file__).resolve().parents[1])
parser.add_argument('--output-json',type=Path)
ARGS,_=parser.parse_known_args();SOURCE=ARGS.source_dir.resolve()
spec=importlib.util.spec_from_file_location('csv_compat_canonical_observer',SOURCE/'canonical_observer.py')
co=importlib.util.module_from_spec(spec);spec.loader.exec_module(co)

def real_getinfo_comments():
    tree=ast.parse((SOURCE/'GetInfo.py').read_text(encoding='utf-8'))
    assignments=[n for n in ast.walk(tree) if isinstance(n,ast.Assign) and isinstance(n.value,ast.List)
        and any(isinstance(t,ast.Name) and t.id=='output_comments' for t in n.targets)]
    if len(assignments)!=1:raise AssertionError('Expected one actual GetInfo output_comments assignment')
    expression=ast.fix_missing_locations(ast.Expression(copy.deepcopy(assignments[0].value)))
    orbit=[1.25,-2.5,3.75,-40.5,50.25,-60.125];fd=16.6
    comments=eval(compile(expression,str(SOURCE/'GetInfo.py'),'eval'),
        {'orbit_data':orbit,'fd_value':fd,'canonical_csv_comment':co.canonical_csv_comment})
    return comments,tuple(orbit+[fd])

def real_legacy_reader():
    path=SOURCE/'data_processing.py'
    if not path.is_file():path=SOURCE.parent/'original/data_processing.py'
    tree=ast.parse(path.read_text(encoding='utf-8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='DataProcessor')
    method=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='get_ini_data')
    # Execute the actual reader method byte-for-byte as AST; unrelated imports
    # and model-dependent module globals are intentionally not executed.
    shell=ast.ClassDef(name='ActualDataProcessor',bases=[],keywords=[],body=[copy.deepcopy(method)],decorator_list=[])
    module=ast.fix_missing_locations(ast.Module(body=[shell],type_ignores=[]));ns={}
    exec(compile(module,str(path),'exec'),ns)
    return ns['ActualDataProcessor'].get_ini_data

class CSVCompatibilityTests(unittest.TestCase):
    def setUp(self):
        scratch=Path(os.environ.get('TMPDIR',str(Path.cwd()/'test_scratch'))).resolve()
        scratch.mkdir(parents=True,exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(prefix='csv-comment-compat-',dir=scratch);self.root=Path(self.temp.name)
        self.previous=os.environ.get('DSPH_PREEXISTING_ID_PATH')
        fixture=self.root/'cohort_ids.npz';np.savez_compressed(fixture,particle_id=np.array([2**60+1,2**60+2],dtype=np.uint64))
        os.environ['DSPH_PREEXISTING_ID_PATH']=str(fixture)
    def tearDown(self):
        if self.previous is None:os.environ.pop('DSPH_PREEXISTING_ID_PATH',None)
        else:os.environ['DSPH_PREEXISTING_ID_PATH']=self.previous
        self.temp.cleanup()
    def test_actual_output_comments_compatible_with_actual_plot_reader(self):
        comments,expected=real_getinfo_comments()
        self.assertTrue(comments[0].startswith('# Orbit Data:'))
        self.assertTrue(comments[1].startswith('# Timescale_fd:'))
        self.assertEqual(comments[2],co.canonical_csv_comment())
        path=self.root/'elinfo_synthetic.csv';path.write_text(''.join(comments)+'numsp,star_mass\n300,1\n',encoding='utf-8')
        self.assertEqual(real_legacy_reader()(path),expected)
        co.require_canonical_csv(path)
    def test_contract_still_rejects_changed_cohort_without_modifying_csv(self):
        comments,_=real_getinfo_comments();path=self.root/'elinfo_synthetic.csv'
        path.write_text(''.join(comments)+'numsp\n300\n',encoding='utf-8');before=path.read_bytes()
        different=self.root/'different.npz';np.savez_compressed(different,particle_id=np.array([2**60+3],dtype=np.uint64))
        os.environ['DSPH_PREEXISTING_ID_PATH']=str(different)
        with self.assertRaises(ValueError):co.require_canonical_csv(path)
        self.assertEqual(path.read_bytes(),before)

if __name__=='__main__':
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CSVCompatibilityTests))
    receipt=dict(state='PASS' if result.wasSuccessful() else 'FAIL',tests_run=result.testsRun,source_dir=str(SOURCE),
        GetInfo_sha256=hashlib.sha256((SOURCE/'GetInfo.py').read_bytes()).hexdigest(),
        test_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),snapshot_reads=False)
    if ARGS.output_json:
        with ARGS.output_json.open('x',encoding='utf-8') as f:json.dump(receipt,f,indent=2)
    print(json.dumps(receipt,indent=2));raise SystemExit(not result.wasSuccessful())
