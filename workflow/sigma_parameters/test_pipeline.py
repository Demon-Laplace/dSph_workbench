"""Synthetic end-to-end physical catalogue; never reads a real snapshot."""
from pathlib import Path
import tempfile,unittest
import numpy as np
import h5py
import stellar_parameters as S

class Pipeline(unittest.TestCase):
 def test_property_catalogue_with_missing_or_repeated_IDs(self):
  rng=np.random.default_rng(712)
  with tempfile.TemporaryDirectory(dir=str(Path(__file__).parent)) as folder:
   p=Path(folder)/'snapshot_000.hdf5'
   host=rng.normal(0,1,(1600,3));dwarf=rng.normal(0,.3,(600,3))+[150,0,0]
   positions=np.vstack([host,dwarf]);vel=rng.normal(0,5,(2200,3));vel[1600:]+=[0,90,0]
   with h5py.File(p,'w') as f:
    hd=f.create_group('Header');hd.attrs['Time']=1.;hd.attrs['MassTable']=np.zeros(6)
    g=f.create_group('PartType2');g['Coordinates']=positions;g['Velocities']=vel;g['Masses']=np.ones(2200);g['StellarFormationTime']=np.full(2200,.2)
   ctx=dict(snapshot_paths={'0':[str(p)]},units=dict(length_to_kpc=1.,velocity_to_kms=1.,mass_to_msun=1.,time_to_gyr=1.),
       center_seeds={'0':[150,0,0]},methods_path=str(Path(__file__).with_name('observation.py')))
   row,center,bulk=S.snapshot(ctx,0)
   self.assertFalse(row['particle_IDs_read']);self.assertEqual(row['stats']['all']['rlt5.0']['N'],600)
   self.assertEqual(row['Walker']['N'],600)
   self.assertTrue(np.isfinite(row['stats']['all']['Rlt5.0']['sigma_los_Walker_kms']))
   with h5py.File(p,'a') as f:f['PartType2']['ParticleIDs']=np.ones(2200,dtype='u8')
   changed,_,_=S.snapshot(ctx,0)
   self.assertEqual(row['source_stellar_fields_sha256'],changed['source_stellar_fields_sha256'])
   self.assertEqual(row['stats'],changed['stats'])
if __name__=='__main__':unittest.main()
