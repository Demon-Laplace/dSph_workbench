import contextlib
import io
import json
from pathlib import Path
import re
import sys
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np

CODE = Path(__file__).resolve().parents[2]/'local_tools/analysis'
sys.path.insert(0, str(CODE))
import fornax_viewer_gui as desktop
import visualize_simulation as viewer


def disk():
    theta = np.linspace(0,2*np.pi,80,endpoint=False)
    x = np.column_stack((np.cos(theta),np.sin(theta),np.zeros_like(theta)))
    v = np.column_stack((-np.sin(theta),np.cos(theta),np.zeros_like(theta))) * 7
    return x,v,np.ones(len(x))


def ascii_dataset(path):
    x,v,m=disk()
    with path.open('w') as f:
        f.write(f'0 0 {len(x)} 0\n')
        for a,b,c in zip(x,v,m):
            # Verify that the reader never parses an ID.
            f.write('ignored '+' '.join(map(str, [*a,*b,c,0,2]))+'\n')


class ViewerTests(unittest.TestCase):
    def setUp(self):
        td = tempfile.TemporaryDirectory()
        self.addCleanup(td.cleanup)
        self.root=Path(td.name)

    def test_view_basis_and_full_aperture_stats(self):
        x,v,m=disk()
        center,bulk,basis,shape=viewer.define_frame((x+20,v+1,m),5,2,'short')
        np.testing.assert_allclose(center,[20,20,20],atol=1e-12)
        np.testing.assert_allclose(bulk,[1,1,1],atol=1e-12)
        np.testing.assert_allclose(basis[:,2],[0,0,1],atol=1e-12)
        self.assertAlmostEqual(np.linalg.det(basis),1)
        self.assertAlmostEqual(shape[1],0)
        g=viewer.group_display(x,v,m,5,100,8,42,15)
        self.assertEqual(g['n'],80)
        self.assertAlmostEqual(g['mean'],7)
        self.assertAlmostEqual(g['r50'],1)
        self.assertEqual(g['prograde'],1)
        self.assertEqual(len(g['points']['x']),8)
        self.assertEqual(len(g['light']['xyz']),15)
        self.assertAlmostEqual(sum(g['light']['mass']),8000)

    def test_offline_html_no_overwrite(self):
        path=self.root/'Fornax8004.ini'
        ascii_dataset(path)
        output=self.root/'viewer.html'
        args=[str(path),'--output',str(output),'--units','fornax','--points','7','--profile-points','18']
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(viewer.main(args),0)
        page=output.read_text()
        data=json.loads(re.search(r'<script id="payload" type="application/json">(.*?)</script>',page).group(1))
        g=data['frames'][0]['groups']['2']
        self.assertEqual(g['n'],80)
        self.assertAlmostEqual(g['mass'],80e10)
        self.assertEqual(len(g['points']['x']),7)
        self.assertEqual(len(g['light']['xyz']),18)
        self.assertIn('plotly_relayout',page)
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            viewer.main(args)

    def test_split_hdf5_masstable_and_desktop_discovery(self):
        import h5py
        d=self.root/'Fornax8004'/'orbit_example'
        d.mkdir(parents=True)
        x,v,_=disk()
        for i in (0,1):
            with h5py.File(d/f'snapshot_002.{i}.hdf5','w') as f:
                h=f.create_group('Header')
                h.attrs['Time']=.2
                h.attrs['NumFilesPerSnapshot']=2
                h.attrs['NumPart_Total']=[0,0,80,0,0,0]
                h.attrs['MassTable']=[0,0,.3,0,0,0]
                g=f.create_group('PartType2')
                g['Coordinates']=x[i*40:(i+1)*40]
                g['Velocities']=v[i*40:(i+1)*40]
        candidates=desktop.candidate_inputs(d.parent)
        self.assertEqual(len(candidates),1)
        self.assertIn('Old stars: 80',desktop.source_info(candidates[0]))
        key,paths=viewer.discover(candidates,8)[0]
        self.assertEqual(len(paths),2)
        args=SimpleNamespace(center=None,bulk=None,center_type=2,types=[0,2],radius=5,axis_radius=2,
                             axis='short',points=10,light_types=[2,4],profile_points=20,seed=42)
        frame=viewer.build_frame(key,paths,args,1.,1.,1.)
        self.assertAlmostEqual(frame['time'],.2)
        self.assertAlmostEqual(frame['groups']['2']['mass'],24,places=5)
        self.assertEqual(frame['groups']['2']['n'],80)
        with self.assertRaisesRegex(ValueError,'Incomplete'):
            viewer.hdf_header(paths[:1])

    def test_ignores_parameter_ini_and_selects_latest_frame(self):
        d=self.root/'Fornax8004'
        d.mkdir()
        ini=self.root/'IniFile'
        ini.mkdir()
        (ini/'IC_Fornax8004.ini').write_text('[Density_stellar]\nmass=.004\n[Component_stellar]\nnbody=400000\n')
        self.assertIn('400,000',desktop.model_info(d))
        (d/'IC_Fornax8004.ini').write_text('[Density_stellar]\nmass=.004\n')
        p=d/'Fornax8004.ini'; ascii_dataset(p)
        self.assertEqual(desktop.candidate_inputs(d),[p])
        self.assertIn('80',desktop.source_info(p))
        for i in (2,3,10):
            (d/f'snapshot_{i:03d}.hdf5').touch()
        latest=viewer.discover([d],1,skip_latest=True)
        self.assertEqual(latest[0][0].name,'snapshot_003')


if __name__=='__main__':
    unittest.main()
