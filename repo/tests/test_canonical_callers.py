"""Focused synthetic caller integration, with snapshot I/O/classification replaced.

Real centering, per-component six-D adapters, context geometry, motion and Walker
are exercised. No real snapshot/file loader or publication output is executed.
"""
from __future__ import annotations
import argparse, ast, importlib.util, json, os, sys, types, unittest, hashlib, tempfile
from pathlib import Path
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
PARSER=argparse.ArgumentParser()
PARSER.add_argument('--source-dir',type=Path,default=HERE.parent)
PARSER.add_argument('--output-json',type=Path)
ARGS,_=PARSER.parse_known_args();SOURCE=ARGS.source_dir.resolve()
sys.path.insert(0,str(SOURCE));os.environ['DSPH_SKIP_CUBE_PATH']='1'
import numpy as np
import pandas as pd
_fixture_dir=tempfile.TemporaryDirectory(prefix='caller_fixture_',)
_fixture_cohort=Path(_fixture_dir.name)/'cohort_ids.npz'
np.savez(_fixture_cohort,particle_id=np.arange(2**60+132,2**60+140,dtype=np.uint64))
import astropy.units as un
from astropy.coordinates import SkyCoord
import canonical_observer as co
import canonical_display as display
from particle_selection import ParticleSelectionMixin
from reference_observer import expected_sixd

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def fixture(include_gas=True):
    H=np.array([480.,500.,525.]);B=np.array([-17.,9.,22.])
    core=np.tile(np.vstack([np.eye(3),-np.eye(3)])*.2,(20,1))
    bulge=np.tile(np.vstack([np.eye(3),-np.eye(3)])*.15,(2,1))
    center=np.array([-15.,-80.,-113.]);xy=np.array([[-1.,0],[1.,0],[0.,-1],[0.,1],[-.7,-.7],[.7,.7],[-.3,.8],[.3,-.8]])
    dwarf=center+np.column_stack([xy,np.zeros(len(xy))])
    grf=np.column_stack([96.+xy[:,0],-56.+xy[:,1]*2,70.+xy[:,0]*.4])
    p=np.vstack([core,bulge,dwarf/co.LENGTH_CODE_TO_KPC,(center+[.3,.2,.1])/co.LENGTH_CODE_TO_KPC])
    v=np.vstack([np.zeros_like(core),np.tile([2.,-4.,6.],(len(bulge),1)),grf/co.VELOCITY_CODE_TO_KMS,[150.,10.,90.]])+B
    masses=np.r_[np.ones(len(core)),np.ones(len(bulge))*2,np.ones(len(dwarf)),7.]
    tp=np.r_[np.tile([2,4],60),np.full(len(bulge),3),np.full(len(dwarf),2),4]
    temp=np.zeros(len(p))
    if include_gas:
        p=np.vstack([p,dwarf[0]/co.LENGTH_CODE_TO_KPC,dwarf[0]/co.LENGTH_CODE_TO_KPC])
        v=np.vstack([v,grf[0]/co.VELOCITY_CODE_TO_KMS+B,grf[0]/co.VELOCITY_CODE_TO_KMS+B])
        masses=np.r_[masses,3.,4.];tp=np.r_[tp,0,0];temp=np.r_[temp,100.,1e6]
    p=p+H
    birth=np.zeros(len(p));birth[140]=1.
    df=pd.DataFrame(dict(id=np.arange(2**60,2**60+len(p),dtype=np.uint64),
        x=p[:,0],y=p[:,1],z=p[:,2],vx=v[:,0],vy=v[:,1],vz=v[:,2],
        m=masses,tp=tp,temp=temp,nh=np.ones(len(p)),birth=birth))
    # Positions/velocities are independently expected after core H and R5 B.
    expectedB=B+np.array([2.,-4.,6.])*(24./144.)
    return df,H,expectedB

class SyntheticSimulation(ParticleSelectionMixin):
    include_gas=True
    def __init__(self,folder_path=None,snapshot_num=None):self.df,self.H,self.B=fixture(self.include_gas)
    def load_snapshot(self,**kwargs):self.loader_kwargs=kwargs;return self.df,2.7
    def classify_mw_dwarf(self,df=None,**kwargs):
        df=self.df if df is None else df
        stars=df['tp'].to_numpy()>=2
        dwarf=stars & (df['r'].to_numpy()>20.)
        return stars&~dwarf,dwarf

class MinimalAnalysis:
    @staticmethod
    def find_center_3d(x,y,z,mass):return tuple(np.average(np.column_stack([x,y,z]),axis=0,weights=mass))

def context_module():
    # Only remove the loader-owned basefunc dependency: all methods under test
    # below come directly from the production files.
    old=sys.modules.get('basefunc')
    stub=types.ModuleType('basefunc');stub.Analysis=MinimalAnalysis;stub.GalaxySimulation=SyntheticSimulation
    sys.modules['basefunc']=stub
    try:
        spec=importlib.util.spec_from_file_location('caller_review_context',SOURCE/'snapshot_context.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    finally:
        if old is None:sys.modules.pop('basefunc',None)
        else:sys.modules['basefunc']=old
    return module

class CallerTests(unittest.TestCase):
    def setUp(self):
        environment=patch.dict(os.environ,DSPH_PREEXISTING_ID_PATH=str(_fixture_cohort),
            DSPH_CANONICAL_OUTPUT_ROOT=str(Path(_fixture_dir.name)/'outputs'))
        environment.start();self.addCleanup(environment.stop)

    def test_MW_center_units_bulk_once_and_tp3(self):
        simulation=SyntheticSimulation();raw=simulation.df.copy(deep=True)
        result=simulation.center_on_MW()
        state=simulation.observer_host_state
        np.testing.assert_allclose(state['host_center_input_code'],simulation.H,rtol=0,atol=2e-12)
        np.testing.assert_allclose(state['host_bulk_input_code'],simulation.B,rtol=0,atol=2e-12)
        self.assertEqual(state['bulk_particle_count'],132)
        p=(raw[['x','y','z']].to_numpy()-simulation.H)*co.LENGTH_CODE_TO_KPC
        v=(raw[['vx','vy','vz']].to_numpy()-simulation.B)*co.VELOCITY_CODE_TO_KMS
        np.testing.assert_allclose(result[['x','y','z']],p,rtol=0,atol=2e-12)
        np.testing.assert_allclose(result[['vx','vy','vz']],v,rtol=0,atol=2e-12)
        before=result.copy(deep=True)
        simulation.center_on_MW()
        pd.testing.assert_frame_equal(simulation.df,before,check_exact=True)

    def test_star_cold_hot_same_sixd_and_ID_integrity(self):
        simulation=SyntheticSimulation();simulation.center_on_MW()
        first_dwarf=np.flatnonzero((simulation.df['tp'].to_numpy()==2)&(simulation.df['r'].to_numpy()>20.))[0]
        cold=len(simulation.df)-2;hot=len(simulation.df)-1
        outputs=[]
        for index,with_velocity in [(first_dwarf,True),(cold,False),(hot,False)]:
            mask=np.zeros(len(simulation.df),dtype=bool);mask[index]=True
            out=simulation.convert_coordinates_for_mask(mask,with_velocity=with_velocity)
            self.assertEqual(out['id'][0],simulation.df['id'].iloc[index])
            self.assertEqual(out['id'].dtype,np.dtype('uint64'))
            outputs.append(out)
        for key in ('ra','dec','rh','vlos','pmra','pmdec','vxh','vyh','vzh'):
            np.testing.assert_allclose([out[key][0] for out in outputs],outputs[0][key][0],rtol=0,atol=2e-11)
        row=simulation.df.iloc[first_dwarf]
        expected=expected_sixd(row[['x','y','z']].to_numpy(dtype=float)[None,:],row[['vx','vy','vz']].to_numpy(dtype=float)[None,:])
        for key in ('ra','dec','rh','vlos','pmra','pmdec'):
            np.testing.assert_allclose(outputs[0][key],expected[key],rtol=0,atol=2e-11)

    def test_conversion_before_center_rejected(self):
        simulation=SyntheticSimulation()
        with self.assertRaises(ValueError):simulation.convert_coordinates_for_mask(np.ones(len(simulation.df),dtype=bool))

    def test_context_real_adapters_share_geometry_and_R5(self):
        SyntheticSimulation.include_gas=True
        context=context_module().prepare_snapshot_context('UNOPENED_SYNTHETIC_ONLY',285,core_radius=5.)
        self.assertEqual(context['simulation'].loader_kwargs['box'],False)
        self.assertEqual(context['canonical_observer_id'],co.CANONICAL_OBSERVER_ID)
        self.assertEqual(context['walker_R5']['N_fit'],9)
        self.assertEqual(len(context['star_coords']['id']),9)
        self.assertEqual(int(np.sum(context['fixed_old_local_mask'])),8)
        self.assertEqual(context['motion_R5']['HRF_EN_masyr'].shape,(2,))
        for key in ('star_coords','cold_gas_coords','hot_gas_coords'):
            self.assertEqual(context[key]['canonical_observer_id'],co.CANONICAL_OBSERVER_ID)
            self.assertAlmostEqual(context[key]['_sky']['center_RA_deg'],context['dw_cra'],places=12)
            self.assertAlmostEqual(context[key]['_sky']['center_Dec_deg'],context['dw_cdec'],places=12)
        self.assertAlmostEqual(context['cold_gas_coords']['vlos'][0],context['hot_gas_coords']['vlos'][0],places=12)
        old_only=co.walker61_zero_error(context['star_coords']['_sky'],
            context['star_coords']['vlos'],
            fit_mask=context['fit_R5_mask'] & context['fixed_old_local_mask'])
        self.assertEqual(old_only['N_fit'],8)
        self.assertFalse(np.isclose(context['walker_R5']['sigma_raw_mle_kms'],
            old_only['sigma_raw_mle_kms']))
        np.testing.assert_array_equal(context['walker_R5']['fit_particle_id'],
            context['star_coords']['id'][context['fit_R5_mask']])

    def test_context_empty_gases(self):
        SyntheticSimulation.include_gas=False
        try:
            context=context_module().prepare_snapshot_context('UNOPENED_SYNTHETIC_ONLY',285,core_radius=5.)
        finally:SyntheticSimulation.include_gas=True
        self.assertEqual(len(context['cold_gas_coords']['ra']),0)
        self.assertEqual(len(context['hot_gas_coords']['pmra']),0)
        self.assertEqual(context['walker_R5']['N_fit'],9)

    def test_geodesic_arrow_origin_and_field_edge(self):
        direction=np.array([.6,-.8]);length=.65
        for ra,dec in [(41.,-34.),(38.9,-36.1),(43.1,-31.9),(359.8,-70.)]:
            result=display.center_motion_arrow({'center_RA_deg':ra,'center_Dec_deg':dec},direction,length_deg=length)
            center=SkyCoord(ra=ra*un.deg,dec=dec*un.deg,frame='icrs')
            tip=center.directional_offset_by(np.arctan2(direction[0],direction[1])*un.rad,length*un.deg)
            expected=np.array([(tip.ra.deg-ra+180)%360-180,tip.dec.deg-dec])
            np.testing.assert_allclose(result['origin_delta_RA_Dec_deg'],0,rtol=0,atol=0)
            np.testing.assert_allclose(result['tip_delta_RA_Dec_deg'],expected,rtol=0,atol=2e-11)
            self.assertEqual(result['length_represents'],'direction_only')

    def test_entry_schema_and_guard_calls(self):
        getinfo=ast.parse((SOURCE/'GetInfo.py').read_text(encoding='utf-8'))
        columns=next(ast.literal_eval(node.value) for node in getinfo.body if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='ELINFO_COLUMNS' for t in node.targets))
        for key in ('walker_R5_k_kms_per_deg','walker_R5_v0_kms','walker_R5_sigma_raw_mle_kms','walker_R5_sigma_after_mle_kms','walker_R5_nstar'):
            self.assertIn(key,columns)
        for name in ('GetInfo.py','PlotFig.py'):
            tree=ast.parse((SOURCE/name).read_text(encoding='utf-8'))
            called=[n.func.id for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name)]
            self.assertIn('require_canonical_csv',called)
            self.assertIn('DSPH_CANONICAL_OUTPUT_ROOT',(SOURCE/name).read_text(encoding='utf-8'))
        metrics=(SOURCE/'snapshot_metrics.py').read_text(encoding='utf-8')
        self.assertIn("snapshot['motion_R5']['HRF_EN_masyr']",metrics)
        plot=(SOURCE/'PlotFig.py').read_text(encoding='utf-8')
        self.assertIn("snapshot['motion_R5']['HRF_EN_masyr']",plot)
        self.assertIn("label='Heliocentric proper motion'",plot)

    def test_inverse_gnomonic_raw_offsets(self):
        for ra0,dec0 in [(41.,-34.),(359.8,-70.)]:
            coords=dict(center_RA_deg=ra0,center_Dec_deg=dec0,center_distance_kpc=140.)
            ra=np.deg2rad(ra0);dec=np.deg2rad(dec0)
            e=np.array([-np.sin(ra),np.cos(ra),0.])
            n=np.array([-np.sin(dec)*np.cos(ra),-np.sin(dec)*np.sin(ra),np.cos(dec)])
            ray=np.array([np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)])
            x,y=np.meshgrid(np.array([-5.,-.1,0.,2.,5.]),np.array([-5.,0.,5.]))
            dra,ddec=co.tangent_kpc_to_raw_offsets(coords,x,y)
            # Independent Cartesian SkyCoord round-trip to ICRS angles.
            q=140.*ray+x[...,None]*e+y[...,None]*n
            sky=SkyCoord(x=q[...,0]*un.kpc,y=q[...,1]*un.kpc,z=q[...,2]*un.kpc,
                         representation_type='cartesian',frame='icrs').spherical
            expected_dra=(sky.lon.deg-ra0+180)%360-180
            np.testing.assert_allclose(dra,expected_dra,rtol=0,atol=2e-11)
            np.testing.assert_allclose(ddec,sky.lat.deg-dec0,rtol=0,atol=2e-11)

    def test_raw_sky_pixel_area_conservation(self):
        for dec0 in (-70.,-34.,0.,40.):
            n=120;half=2.1
            row_area=co.raw_sky_pixel_area_arcsec2(dec0,n,half_width_deg=half)
            self.assertEqual(row_area.shape,(n,));self.assertTrue(np.all(row_area>0))
            # Exact integral dOmega=dRA*d(sin Dec), independently integrated
            # over the complete rectangular RA/Dec patch.
            solid=np.deg2rad(2*half)*(np.sin(np.deg2rad(dec0+half))-np.sin(np.deg2rad(dec0-half)))
            expected=solid*(180/np.pi*3600)**2
            self.assertAlmostEqual(float(row_area.sum()*n)/expected,1.,places=12)

    def test_frozen_cohort_missing_ID_fails_closed(self):
        current=np.arange(2**60+132,2**60+140,dtype=np.uint64)
        mask,state=co.frozen_cohort_mask(current,'UNOPENED_SYNTHETIC_ONLY')
        self.assertTrue(mask.all());self.assertEqual(state['cohort_N'],8)
        self.assertEqual(state['cohort_definition'],'selected_then_frozen')
        with self.assertRaises(ValueError):co.frozen_cohort_mask(current[:-1],'UNOPENED_SYNTHETIC_ONLY')

if __name__=='__main__':
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CallerTests))
    names=['canonical_observer.py','canonical_display.py','coordinate_transform.py','particle_selection.py',
           'snapshot_context.py','snapshot_metrics.py','GetInfo.py','PlotFig.py','snapshot_io.py']
    document=dict(state='PASS' if result.wasSuccessful() else 'FAIL',tests_run=result.testsRun,
        source_dir=str(SOURCE),source_sha256={name:sha(SOURCE/name) for name in names},
        test_source_sha256=sha(Path(__file__)),errors=[str(t) for t,_ in result.errors],
        failures=[str(t) for t,_ in result.failures],snapshot_IO='mocked; no real snapshot read',
        classification='synthetic fixed selection; no claim to validate density classification',
        real_algorithms='production MW center/bulk, full 6D adapters, common context, Walker and PM',production_mutations=False)
    if ARGS.output_json:
        if ARGS.output_json.exists():raise FileExistsError(ARGS.output_json)
        ARGS.output_json.parent.mkdir(parents=True,exist_ok=True)
        ARGS.output_json.write_text(json.dumps(document,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(document,ensure_ascii=False));_fixture_dir.cleanup();raise SystemExit(not result.wasSuccessful())
