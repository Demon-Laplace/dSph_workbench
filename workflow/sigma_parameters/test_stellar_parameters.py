"""Focused invariants for the exported aperture moments (synthetic inputs)."""
from pathlib import Path
import importlib.util
import sys
import types
import unittest
import numpy as np

# The moment routines do not read snapshots. Avoid importing production I/O
# just to exercise their mathematical definitions in the local bundled Python.
stub = types.ModuleType('physical_helpers')
sys.modules.setdefault('physical_helpers', stub)
spec = importlib.util.spec_from_file_location('stellar_parameters_under_test',
          Path(__file__).with_name('stellar_parameters.py'))
S = importlib.util.module_from_spec(spec)
spec.loader.exec_module(S)


class ApertureMoments(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(7719)
        self.x = self.rng.normal(size=(120,3))
        self.v = self.rng.normal(size=(120,3))*[2.,3.,4.]
        self.row_labels = np.arange(120)
        self.m = np.linspace(1.,2.,120)
        self.axis = np.array([1.,2.,3.]);self.axis /= np.linalg.norm(self.axis)
        self.q = np.ones(120,dtype=bool)
        self.los = self.v@self.axis
        self.resid = self.los-2.

    def calculate(self, **kw):
        args = dict(row_labels=self.row_labels,x=self.x,v=self.v,m=self.m,los=self.los,
                    resid=self.resid,q=self.q,axis=self.axis)
        args.update(kw)
        return S.moments(**args)

    def test_known_three_dimensional_definition(self):
        got = self.calculate()
        expected = np.var(self.v,axis=0,ddof=0).sum()
        self.assertAlmostEqual(got['sigma_3d_kms']**2,expected,places=11)
        self.assertAlmostEqual(got['sigma_3d_kms'],np.sqrt(3)*got['sigma_1d_equiv_kms'],places=12)
        self.assertAlmostEqual(got['sigma_common_LOS_tensor_kms'],np.std(self.los),places=12)
        self.assertAlmostEqual(got['rms_Walker_residual_kms']**2,
             got['sigma_los_Walker_kms']**2+got['mean_Walker_residual_kms']**2,places=11)

    def test_rigid_rotation_preserves_trace_and_shared_LOS(self):
        Q,_ = np.linalg.qr(self.rng.normal(size=(3,3)))
        if np.linalg.det(Q)<0:Q[:,0]*=-1
        original = self.calculate()
        rotated = self.calculate(x=self.x@Q.T,v=self.v@Q.T,axis=Q@self.axis)
        for key in ('sigma_3d_kms','sigma_1d_equiv_kms','sigma_common_LOS_tensor_kms',
                    'rms_axis_c_over_a','los_depth_std_kpc','sigma_los_raw_kms'):
            self.assertAlmostEqual(original[key],rotated[key],places=11,msg=key)
        np.testing.assert_allclose(rotated['velocity_covariance_kms2'],
                    Q@np.asarray(original['velocity_covariance_kms2'])@Q.T,rtol=1e-12,atol=1e-12)

    def test_velocity_zero_point_preserves_dispersion(self):
        shift=np.array([35.,-20.,70.])
        original = self.calculate()
        translated = self.calculate(v=self.v+shift,los=self.los+shift@self.axis)
        for key in ('sigma_3d_kms','sigma_1d_equiv_kms','sigma_common_LOS_tensor_kms',
                    'sigma_los_raw_kms','sigma_3d_mass_kms'):
            self.assertAlmostEqual(original[key],translated[key],places=11,msg=key)
        np.testing.assert_allclose(translated['velocity_mean_relative_core_kms'],
                      np.asarray(original['velocity_mean_relative_core_kms'])+shift,atol=1e-12)

    def test_sphere_and_true_ray_projection_select_different_rows(self):
        x=np.array([[.1,0.,.1],[.1,.1,-.1],[.1,.1,4.],[.2,0.,-3.]])
        h=x+[0.,0.,10.]
        R=10.*np.linalg.norm(h[:,:2],axis=1)/h[:,2]
        r=np.linalg.norm(x,axis=1)
        rows=np.arange(4)
        v=np.arange(12,dtype=float).reshape(4,3)
        a=S.moments(rows,x,v,np.ones(4),v[:,2],v[:,2],r<.5,np.array([0.,0.,1.]))
        b=S.moments(rows,x,v,np.ones(4),v[:,2],v[:,2],R<.5,np.array([0.,0.,1.]))
        self.assertEqual(a['N'],2);self.assertEqual(b['N'],4)
        self.assertFalse(np.array_equal(r<.5,R<.5))

    def test_constant_velocity_zero_variance(self):
        x=np.tile(np.eye(3),(2,1));v=np.zeros((6,3));ids=np.arange(6,dtype=np.uint64)
        got=S.moments(ids,x,v,np.ones(6),np.ones(6)*4,np.ones(6)*2,
                        np.ones(6,dtype=bool),np.array([0.,0.,1.]))
        self.assertEqual(got['sigma_3d_kms'],0.)
        self.assertEqual(got['sigma_los_raw_kms'],0.)
        self.assertEqual(got['mean_Walker_residual_kms'],2.)
        self.assertEqual(got['rms_Walker_residual_kms'],2.)
        self.assertIsNone(got['beta_sigma'])
        self.assertAlmostEqual(got['N_eff_mass'],6.)

    def test_empty_and_singleton_are_explicitly_unsupported(self):
        for n in (0,1):
            q=np.zeros(120,dtype=bool);q[:n]=True
            got=self.calculate(q=q)
            self.assertEqual(got['N'],n)
            self.assertFalse(got['supported'])
            self.assertNotIn('sigma_3d_kms',got)

    def test_synthetic_row_labels_do_not_affect_property_statistics(self):
        original=self.calculate()
        relabeled=self.calculate(row_labels=self.row_labels[::-1])
        self.assertEqual(original,relabeled)


if __name__=='__main__':unittest.main()
