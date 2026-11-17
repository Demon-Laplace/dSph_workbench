"""Property-only loader and centre-seed tests; no real snapshot or particle IDs."""
from pathlib import Path
import tempfile
import unittest
import sys
import types
import numpy as np
try:
    import h5py
    HAS_H5PY=True
except ImportError:
    sys.modules['h5py']=types.ModuleType('h5py')
    HAS_H5PY=False
import spatial_context as S


class PropertyContext(unittest.TestCase):
    def test_merger_seed_translation_and_legacy_rotation(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)
            (p/'mergegalaxy.par').write_text('mw.ini dwarf.ini\n1 2 3 0 0 0 0 0 0 19 22 25\n')
            (p/'mergecube.par').write_text('cube.ini mergegalaxy.dat\n500 500 500 0 0 0 0 0 0 0 0 0\n')
            seed,_=S.initial_seed_from_mergers(p)
            np.testing.assert_allclose(seed,[501,502,503])
            (p/'mergecube.par').write_text('cube.ini mergegalaxy.dat\n500 500 500 0 0 0 0 0 0 0 0 90\n')
            seed,_=S.initial_seed_from_mergers(p)
            theta=90*.017453
            np.testing.assert_allclose(seed,[500+np.cos(theta)+2*np.sin(theta),
                 500-np.sin(theta)+2*np.cos(theta),503],atol=1e-12)

    def test_local_center_does_not_jump_to_more_massive_host(self):
        rng=np.random.default_rng(711)
        dwarf=rng.normal(0,.2,(2000,3))+[10,0,0]
        host=rng.normal(0,2,(20000,3))+[100,0,0]
        phase=np.concatenate([dwarf,host]);mass=np.ones(len(phase))
        center=S.shrinking_center(phase,mass,[10.1,0,0],1.)
        np.testing.assert_allclose(center,[10,0,0],atol=.05)
        translated=S.shrinking_center(phase+500,mass,np.array([10.1,0,0])+500,1.)
        np.testing.assert_allclose(translated,center+500,atol=1e-11)

    def test_sparse_local_center_is_explicit_failure(self):
        with self.assertRaisesRegex(ValueError,'Too few local'):
            S.shrinking_center(np.zeros((10,3)),np.ones(10),[0,0,0],1.)

    @unittest.skipUnless(HAS_H5PY,'h5py unavailable in this local Python')
    def test_loader_reads_hdf5_without_any_ParticleIDs_dataset(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'snapshot_000.hdf5'
            with h5py.File(path,'w') as f:
                head=f.create_group('Header')
                head.attrs['Time']=2.;head.attrs['MassTable']=[0,0,3,0,0,0]
                g=f.create_group('PartType2')
                g['Coordinates']=np.arange(12,dtype='f4').reshape(4,3)
                g['Velocities']=np.arange(12,dtype='f4').reshape(4,3)*2
            st=path.stat()
            ctx=dict(snapshot_paths={'0':[str(path)]},units=dict(time_to_gyr=.5,mass_to_msun=10.),
               input_snapshot_fingerprints={'0':[dict(path=str(path),size_bytes=st.st_size,mtime_ns=st.st_mtime_ns)]})
            result=S.read_stellar_properties(ctx,0)
            self.assertEqual(set(result),{'p','v','m','tp','birth','age_Gyr','age_source','time','source_stats','source_fields_sha256'})
            self.assertEqual(result['time'],1.)
            np.testing.assert_array_equal(result['m'],[30.]*4)
            np.testing.assert_array_equal(result['tp'],[2]*4)
            self.assertTrue(np.isnan(result['birth']).all())
            self.assertTrue(np.isnan(result['age_Gyr']).all())
            self.assertEqual(result['p'].dtype,np.dtype('f4'))
            self.assertNotIn('ids',result)

    @unittest.skipUnless(HAS_H5PY,'h5py unavailable in this local Python')
    def test_loader_rejects_changed_snapshot(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'snapshot_000.hdf5';path.write_bytes(b'abc')
            ctx=dict(snapshot_paths={'0':[str(path)]},input_snapshot_fingerprints={'0':[
              dict(path=str(path),size_bytes=100,mtime_ns=0)]})
            with self.assertRaisesRegex(RuntimeError,'Source snapshot changed'):
                S.read_stellar_properties(ctx,0)

    @unittest.skipUnless(HAS_H5PY,'h5py unavailable in this local Python')
    def test_recorded_formation_time_zero_is_not_discarded(self):
        with tempfile.TemporaryDirectory() as directory:
            with h5py.File(str(Path(directory)/'age.hdf5'),'w') as f:
                g=f.create_group('PartType4')
                g['StellarFormationTime']=[0.,1.,2.]
                age,birth,source=S.stellar_age_properties(g,1.,dict(time_to_gyr=.5),3)
                np.testing.assert_allclose(age,[1.,.5,0.])
                np.testing.assert_allclose(birth,[0.,.5,1.])
                self.assertEqual(source['field'],'StellarFormationTime')
                self.assertIn('including_zero',source['method'])

    @unittest.skipUnless(HAS_H5PY,'h5py unavailable in this local Python')
    def test_explicit_age_attribute_units_take_precedence(self):
        with tempfile.TemporaryDirectory() as directory:
            with h5py.File(str(Path(directory)/'age.hdf5'),'w') as f:
                g=f.create_group('PartType2')
                g['Age']=[100.,500.,-1.];g['Age'].attrs['units']='Myr'
                g['StellarFormationTime']=[0.,0.,0.]
                age,birth,source=S.stellar_age_properties(g,3.,dict(time_to_gyr=.5),3)
                np.testing.assert_allclose(age[:2],[.1,.5]);self.assertTrue(np.isnan(age[2]))
                self.assertEqual(source['unit_source'],'dataset_attribute')
                self.assertEqual(source['method'],'recorded_age_field')
                self.assertEqual(source['invalid_count'],1)

    @unittest.skipUnless(HAS_H5PY,'h5py unavailable in this local Python')
    def test_unannotated_explicit_age_declares_code_time_default(self):
        with tempfile.TemporaryDirectory() as directory:
            with h5py.File(str(Path(directory)/'age.hdf5'),'w') as f:
                g=f.create_group('PartType3');g['StellarAge']=[2.,4.]
                age,birth,source=S.stellar_age_properties(g,3.,dict(time_to_gyr=.5),2)
                np.testing.assert_allclose(age,[1.,2.])
                self.assertEqual(source['unit'],'code_time')
                self.assertEqual(source['unit_source'],'explicit_loader_default')


if __name__=='__main__':unittest.main()
