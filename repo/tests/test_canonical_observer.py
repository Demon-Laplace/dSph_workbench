"""Independent, snapshot-free numerical contract checks for canonical observing implementation.

Run with --source-dir <source-directory> --output-json <new-result.json>.
No production data or outputs are opened by these tests.
"""
from __future__ import annotations
import argparse, importlib.util, json, hashlib, os, sys, tempfile, unittest
from pathlib import Path
from unittest.mock import patch

HERE=Path(__file__).resolve().parent
import numpy as np
_fixture_dir=tempfile.TemporaryDirectory(prefix='pure_contract_fixture_',)
_fixture_cohort=Path(_fixture_dir.name)/'cohort_ids.npz'
np.savez(_fixture_cohort,particle_id=np.arange(1000,1010,dtype=np.uint64))
from reference_observer import (expected_sixd, true_spherical_offsets,
    independent_walker, center_pm_from_h_and_u, S_SUN, V_SUN, OBSERVER_ID)

PARSER=argparse.ArgumentParser()
PARSER.add_argument('--source-dir',type=Path,default=HERE.parent)
PARSER.add_argument('--output-json',type=Path)
ARGS,REMAINING=PARSER.parse_known_args()
SOURCE=ARGS.source_dir.resolve()
spec=importlib.util.spec_from_file_location('reviewed_canonical_observer',SOURCE/'canonical_observer.py')
co=importlib.util.module_from_spec(spec);spec.loader.exec_module(co)

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

class ContractTests(unittest.TestCase):
    def setUp(self):
        environment=patch.dict(os.environ,DSPH_PREEXISTING_ID_PATH=str(_fixture_cohort))
        environment.start();self.addCleanup(environment.stop)
        self.p=np.array([[-14.2,-77.5,-110.3],[-15.,-80.,-113.],[-17.9,-83.7,-108.4],
                         [55.,2.1,30.],[.1,-90.,10.],[8.,2.,-5.]])
        self.v=np.array([[80.,-30.,71.],[96.,-56.,70.],[70.,-77.,65.],[-120.,33.,-8.],
                         [1.,234.,77.],[-91.,14.,-82.]])
        self.ids=np.arange(2**60,2**60+len(self.p),dtype=np.uint64)

    def assertAngle(self,a,b,atol=2e-10):
        np.testing.assert_allclose((np.asarray(a)-np.asarray(b)+180)%360-180,0,rtol=0,atol=atol)

    def test_complete_astropy_sixd_and_no_R_LOS(self):
        result=co.production_observables(self.p,self.v,self.ids)
        expected=expected_sixd(self.p,self.v)
        self.assertEqual(result['canonical_observer_id'],OBSERVER_ID)
        self.assertAngle(result['ra'],expected['ra'])
        for key in ('dec','rh','vlos','pmra','pmdec'):
            np.testing.assert_allclose(result[key],expected[key],rtol=0,atol=2e-10)
        np.testing.assert_allclose(result['vlos'],expected['direct_vlos'],rtol=0,atol=2e-10)
        np.testing.assert_array_equal(result['id'],self.ids)
        self.assertEqual(result['id'].dtype,np.dtype('uint64'))

    def test_comoving_sun_and_solar_reflex_probes(self):
        result=co.production_observables(self.p,np.tile(V_SUN,(len(self.p),1)),self.ids)
        for key in ('vlos','pmra','pmdec'):
            np.testing.assert_allclose(result[key],0,rtol=0,atol=2e-10)
        at_rest=co.production_observables(self.p,np.zeros_like(self.p),self.ids)
        h=self.p-S_SUN
        np.testing.assert_allclose(at_rest['vlos'],-(h@V_SUN)/np.linalg.norm(h,axis=1),rtol=0,atol=2e-10)

    def test_units_host_center_bulk_once(self):
        H=np.array([400.,500.,600.]);B=np.array([-17.,9.,22.]);L=1.7;V=.37
        raw_p=self.p/L+H;raw_v=self.v/V+B
        result=co.astropy_native_direct_observables(raw_p,raw_v,self.ids,
            host_center_raw=H,host_bulk_raw=B,length_to_kpc=L,velocity_to_kms=V,**co.OBSERVER_FRAME)
        direct=co.production_observables(self.p,self.v,self.ids)
        for key,newkey in [('distance_kpc','rh'),('vlos_native_kms','vlos'),
                           ('mu_alpha_cosdec_masyr','pmra'),('mu_delta_masyr','pmdec')]:
            np.testing.assert_allclose(result[key],direct[newkey],rtol=0,atol=3e-10)
        np.testing.assert_allclose(result['host_relative_position_kpc'],self.p,rtol=0,atol=2e-13)
        np.testing.assert_allclose(result['grf_velocity_native_kms'],self.v,rtol=0,atol=2e-13)
        # A common inertial boost cancels with the same boosted host bulk.
        boost=np.array([91.,-133.,4.])
        boosted=co.astropy_native_direct_observables(raw_p,raw_v+boost,self.ids,
            host_center_raw=H,host_bulk_raw=B+boost,length_to_kpc=L,velocity_to_kms=V,**co.OBSERVER_FRAME)
        np.testing.assert_allclose(boosted['vlos_native_kms'],result['vlos_native_kms'],rtol=0,atol=3e-10)

    def test_input_array_immutability_and_invalid_IDs(self):
        p=self.p.copy();v=self.v.copy();ids=self.ids.copy()
        co.production_observables(p,v,ids)
        np.testing.assert_array_equal(p,self.p);np.testing.assert_array_equal(v,self.v)
        np.testing.assert_array_equal(ids,self.ids)
        duplicate=self.ids.copy();duplicate[-1]=duplicate[0]
        for invalid in [duplicate,self.ids.astype(float),np.arange(len(self.ids),dtype=np.int64)]:
            with self.assertRaises(ValueError):co.production_observables(self.p,self.v,invalid)
        bad=self.p.copy();bad[0,0]=np.nan
        with self.assertRaises(ValueError):co.production_observables(bad,self.v,self.ids)

    def test_empty_component_supported(self):
        result=co.production_observables(np.empty((0,3)),np.empty((0,3)),np.empty(0,dtype=np.uint64))
        for key in ('id','ra','dec','rh','vlos','pmra','pmdec'):
            self.assertEqual(len(result[key]),0)

    def fixture(self):
        # Physical gnomonic square corners are display members, never R5 members.
        base=co.production_observables(self.p,self.v,self.ids)
        h0=np.array([-15.,-80.,-113.])
        sky0=co.production_sky(base,h0)
        east,north=sky0['east_native'],sky0['north_native']
        xy=np.array([[-4.,0.],[4.,0.],[0.,-4.],[0.,4.],[-2.,-2.],[2.,2.],
                     [-1.,2.],[1.,-2.],[4.4,4.4],[-4.4,-4.4]])
        h=h0+xy[:,0,None]*east+xy[:,1,None]*north
        p=h+S_SUN
        ids=np.arange(1000,1000+len(h),dtype=np.uint64)
        # Start a physically varied velocity field; no hand-assigned source LOS.
        v=np.column_stack([30+xy[:,0]*2,-60+xy[:,1]*3,80+xy[:,0]*.7])
        coords=co.production_observables(p,v,ids)
        sky=co.production_sky(coords,h0)
        return coords,sky,p,v,ids

    def test_spherical_offsets_and_Walker61_reference(self):
        coords,sky,*_=self.fixture()
        expected=true_spherical_offsets(coords['ra'],coords['dec'],sky['center_RA_deg'],
            sky['center_Dec_deg'],sky['center_distance_kpc'])
        for key in ('east_deg','north_deg','rho_deg','projected_radius_kpc'):
            np.testing.assert_allclose(sky[key],expected[key],rtol=0,atol=2e-10)
        np.testing.assert_array_equal(sky['fit_R5_mask'],expected['fit_R5_mask'])
        self.assertFalse(sky['fit_R5_mask'][-1]);self.assertFalse(sky['fit_R5_mask'][-2])
        result=co.walker61_zero_error(sky,coords['vlos'],fit_mask=sky['fit_R5_mask'])
        expectedfit=independent_walker(coords['ra'],coords['dec'],coords['vlos'],
            ra0=sky['center_RA_deg'],dec0=sky['center_Dec_deg'],d0=sky['center_distance_kpc'])
        self.assertEqual(len(result['theta0_scan_deg']),61)
        for key,ex in [('k_scan_kms_per_deg','k_scan'),('v0_scan_kms','v0_scan'),('sigma_mle_scan_kms','sigma_scan')]:
            np.testing.assert_allclose(result[key],expectedfit[ex],rtol=0,atol=2e-9)
        self.assertEqual(result['theta_index'],expectedfit['index'])
        # Extreme corner LOS values remain outside the immutable R5 circle.
        contaminated=coords['vlos'].copy();contaminated[-2:]=[5000.,-4000.]
        result2=co.walker61_zero_error(sky,contaminated,fit_mask=sky['fit_R5_mask'])
        np.testing.assert_array_equal(result2['k_scan_kms_per_deg'],result['k_scan_kms_per_deg'])

    def test_HRF_motion_at_center_and_field_edge(self):
        coords,sky,*_=self.fixture();mass=np.arange(1,len(coords['id'])+1,dtype=float)
        selected=sky['fit_R5_mask']
        summary=co.production_motion(coords,sky,mass,selected)
        u=np.average(coords['_native']['u_native_kms'][selected],axis=0,weights=mass[selected])
        expected=np.array([u@sky['east_native'],u@sky['north_native']])/(4.740470463533349*sky['center_distance_kpc'])
        np.testing.assert_allclose(summary['HRF_EN_masyr'],expected,rtol=0,atol=2e-12)
        M=coords['_native']['physical_to_icrs']
        ra0,dec0=sky['center_RA_deg'],sky['center_Dec_deg'];D=sky['center_distance_kpc']
        for dra,ddec in [(0.,0.),(-2.1,-2.1),(2.1,2.1),(-2.1,2.1),(2.1,-2.1)]:
            ra,dec=np.deg2rad([ra0+dra,dec0+ddec])
            h=M.T@(D*np.array([np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)]))
            at_point=co.proper_motion_components(h[None,:],u[None,:],physical_to_icrs=M)
            expected_point=center_pm_from_h_and_u(h,u)
            np.testing.assert_allclose([at_point['mu_alpha_cosdec_masyr'][0],at_point['mu_delta_masyr'][0]],expected_point,rtol=0,atol=2e-11)

    def test_raw_delta_grid_40_and_2p1_extent(self):
        x=np.array([-2.1,-.01,.01,2.1,2.1001]);y=np.zeros(5)
        v=np.array([1.,2.,4.,8.,9999.])
        grid,extent=co.angle_velocity_grid(x,y,v,1.)
        self.assertEqual(grid.shape,(40,40))
        self.assertEqual(extent,(-2.1,2.1,-2.1,2.1))
        self.assertEqual(np.count_nonzero(np.isfinite(grid)),4)
        np.testing.assert_array_equal(np.sort(grid[np.isfinite(grid)]),[0.,1.,3.,7.])

    def test_CSV_observer_guard_preserves_legacy(self):
        with tempfile.TemporaryDirectory(prefix='observer_review_',) as tmp:
            path=Path(tmp)/'elinfo.csv'
            co.require_canonical_csv(path)
            old=b'# legacy observer\nnumsp,pmra,pmdec\n285,0,0\n';path.write_bytes(old)
            with self.assertRaises(ValueError):co.require_canonical_csv(path)
            self.assertEqual(path.read_bytes(),old)
            path.write_text(co.canonical_csv_comment()+'numsp,pmra,pmdec\n285,0,0\n',encoding='utf-8')
            co.require_canonical_csv(path)
            modified=path.read_text().replace(co.observer_contract_sha256(),'0'*64)
            path.write_text(modified,encoding='utf-8')
            with self.assertRaises(ValueError):co.require_canonical_csv(path)

if __name__=='__main__':
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(ContractTests)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    document=dict(state='PASS' if result.wasSuccessful() else 'FAIL',
        source_dir=str(SOURCE),canonical_source_sha256=sha(SOURCE/'canonical_observer.py'),
        test_source_sha256=sha(Path(__file__)),reference_source_sha256=sha(HERE/'reference_observer.py'),
        tests_run=result.testsRun,errors=[str(t) for t,_ in result.errors],
        failures=[str(t) for t,_ in result.failures],snapshot_reads=False,production_mutations=False)
    if ARGS.output_json:
        if ARGS.output_json.exists():raise FileExistsError(ARGS.output_json)
        ARGS.output_json.parent.mkdir(parents=True,exist_ok=True)
        ARGS.output_json.write_text(json.dumps(document,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(document,ensure_ascii=False))
    _fixture_dir.cleanup()
    raise SystemExit(not result.wasSuccessful())
