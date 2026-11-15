import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
import h5py
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import cgm_shell_archive as archive


class ShellArchiveTest(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.path=self.root/'snapshot_000.hdf5'
        self.parameters=dict(UnitLength_in_cm='3.085678e21',UnitMass_in_g='1.989e43',
                             UnitVelocity_in_cm_per_s='1e5',ComovingIntegrationOn='0')
        self.frame=dict(host_center_input_code=[50.,50.,50.],host_bulk_input_code=[5.,0.,0.])
        with h5py.File(self.path,'w') as f:
            h=f.create_group('Header');h.attrs['NumFilesPerSnapshot']=1
            h.attrs['NumPart_Total']=[4,1,0,0,0,0];h.attrs['MassTable']=[0.,2.,0.,0.,0.,0.]
            h.attrs['BoxSize']=100.;h.attrs['Time']=.5;h.attrs['Flag_Cooling']=1
            gas=f.create_group('PartType0')
            gas['Coordinates']=np.array([[51,50,50],[52,50,50],[49,50,50],[151,50,50]],np.float32)
            gas['Velocities']=np.array([[6,0,0],[3,0,0],[7,0,0],[6,1,0]],np.float32)
            gas['ParticleIDs']=np.array([0,1,2,3],np.uint64)
            gas['Masses']=np.array([1.,2.,3.,4.],np.float32)
            gas['Density']=np.array([1.,2.,3.,4.],np.float32)
            gas['InternalEnergy']=np.array([1.,1.e4,100.,1.e5],np.float32)
            gas['ElectronAbundance']=np.ones(4,np.float32)
            gas['NeutralHydrogenAbundance']=np.array([1.,0.,.5,0.],np.float32)
            gas['Metallicity']=np.array([[.1],[.2],[.3],[.4]],np.float32)
            dm=f.create_group('PartType1');dm['Coordinates']=np.array([[50,50,50]],np.float32)
            dm['Velocities']=np.array([[5,0,0]],np.float32);dm['ParticleIDs']=np.array([4],np.uint64)

    def tearDown(self):self.tmp.cleanup()

    def test_periodic_mass_phases_moments_and_chunk_independence(self):
        records=[]
        for chunk in (1,3):
            path=self.root/('archive-%d.npz'%chunk)
            archive.reduce_snapshot(self.path,path,self.parameters,chunk_size=chunk,frame=self.frame)
            with np.load(path,allow_pickle=False) as f:records.append({k:f[k].copy() for k in f.files})
        data=records[0]
        np.testing.assert_allclose(data['mass_by_type_Msun'].sum(axis=1),[1.e11,2.e10,0,0,0,0])
        self.assertEqual(int(data['gas_count'][0].sum()),4)
        self.assertAlmostEqual(float(data['gas_sums'][0,:,15].sum()),.76*(1+3*.5)*1.e10)
        self.assertAlmostEqual(float(data['gas_sums'][0,:,7].sum()),(1-4-6+4)*1.e10)
        np.testing.assert_allclose(data['mass_second_moments_Msun_kpc2'][0,:,0].sum(),
                                   (1+2*4+3+4)*1.e10*archive.units(self.parameters)['length_to_kpc']**2)
        for key in records[0]:
            if key=='metadata_json':continue
            np.testing.assert_allclose(records[0][key],records[1][key],rtol=1.e-12,atol=1.e-10)
        self.assertEqual(json.loads(str(records[0]['metadata_json']))['field_sha256'],
                         json.loads(str(records[1]['metadata_json']))['field_sha256'])
        self.assertTrue(self.path.exists())

    def test_unsupported_units_and_nonfinite_source_fail(self):
        with self.assertRaises(ValueError):archive.units(dict(self.parameters,ComovingIntegrationOn='1'))
        with self.assertRaises(ValueError):archive.units(dict(self.parameters,UnitMass_in_g='1e44'))
        with h5py.File(self.path,'a') as f:f['PartType0/Velocities'][0,0]=np.nan
        with self.assertRaises(ValueError):archive.reduce_snapshot(self.path,self.root/'bad.npz',self.parameters,frame=self.frame)
        self.assertFalse((self.root/'bad.npz').exists())

    def test_histogram_edge_split_is_retained(self):
        self.assertIn(2.e4,archive.TEMP_EDGES)
        self.assertIn(10.,archive.RADIAL_EDGES)
        self.assertIn(150.,archive.RADIAL_EDGES)

    def test_usedvalues_six_digit_serialization_keeps_original_units(self):
        rounded=archive.units(dict(self.parameters,UnitLength_in_cm='3.08568e21'))
        exact=archive.units(self.parameters)
        self.assertEqual(rounded['native_units'],exact['native_units'])
        self.assertEqual(rounded['density_to_cgs'],exact['density_to_cgs'])
        self.assertEqual(rounded['recorded_unit_values']['UnitLength_in_cm'],3.08568e21)
        with self.assertRaises(ValueError):archive.units(dict(self.parameters,UnitLength_in_cm='3.085679e21'))

    def test_stellar_reference_across_periodic_faces(self):
        positions=np.array([[.2,.1,99.8],[99.8,99.9,.2]])
        unwrapped,reference,concentration=archive.unwrap_stellar_positions(
            positions,np.ones(2),np.ones(2,dtype=bool),100.)
        center=unwrapped.mean(axis=0)
        delta=unwrapped-center
        self.assertTrue(np.all(np.linalg.norm(delta,axis=1)<1.))
        np.testing.assert_allclose(np.mod(center+50,100)-50,0.,atol=1.e-12)
        self.assertTrue(np.all(concentration>.99))


if __name__=='__main__':unittest.main()
