"""Small synthetic-HDF5 tests for gas statistics and missing-field boundaries."""
import copy
import importlib.util
from pathlib import Path
import tempfile
import unittest

import h5py
import numpy as np

SPEC = importlib.util.spec_from_file_location("gas_parameters_under_test", Path(__file__).with_name("gas_parameters.py"))
G = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(G)


class GasParametersTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="gas-test-", dir=str(Path(__file__).parent))
        self.root = Path(self.temporary.name)
        self.path = self.root / "snapshot_000.hdf5"
        self.units = dict(length_to_kpc=1.0, velocity_to_kms=1.0, mass_to_msun=1.0,
                          UnitLength_in_cm=G.KPC_CM, UnitMass_in_g=G.MSUN_G,
                          UnitVelocity_in_cm_per_s=1e5, time_to_gyr=1.0)
        self.ctx = dict(snapshot_paths={"0": [str(self.path)]}, initial_source=dict(parent_bits=25),
                        origin_capabilities=dict(gas_parent=True), config={})
        self.roster = dict(gas_parent=np.array([1001, 1002, 1003, 1004], np.uint64))
        self.position = np.zeros((132, 3), np.float64)
        self.position[:4, 0] = [0.1, 0.2, 0.3, 0.4]
        self.position[4:, 0] = 35.0
        self.position[4:, 1] = np.linspace(-3, 3, 128)
        self.nH = np.r_[[0.05, 0.2, 0.8, 2.0], np.full(128, 0.01)]
        self.gas_fields = dict(
            ParticleIDs=np.r_[np.arange(1001, 1005), np.arange(2000, 2128)].astype(np.uint64),
            Coordinates=self.position, Velocities=np.zeros_like(self.position), Masses=np.ones(132),
            Density=self.nH / (G.MSUN_G / G.KPC_CM ** 3 * G.XH / G.MP_G),
            InternalEnergy=np.full(132, 100.0), ElectronAbundance=np.zeros(132),
            NeutralHydrogenAbundance=np.full(132, 0.5), StarFormationRate=np.r_[[0.1, 0.2, 0.0, 0.0], np.zeros(128)],
            default=np.full(132, 10.0))

    def tearDown(self):
        self.temporary.cleanup()

    def write(self, omit=(), no_gas=False):
        with h5py.File(self.path, "w") as snapshot:
            header = snapshot.create_group("Header")
            header.attrs["Time"] = 0.125
            header.attrs["MassTable"] = np.zeros(6)
            if not no_gas:
                group = snapshot.create_group("PartType0")
                for name, value in self.gas_fields.items():
                    if name not in omit:
                        group.create_dataset(name, data=value)

    def run_gas(self, ctx=None, bulk=(10, 0, 0)):
        return G.analyze_gas(ctx or self.ctx, 0, [0, 0, 0], bulk, self.roster, self.units)

    def test_no_gas(self):
        self.write(no_gas=True)
        result = self.run_gas()
        self.assertEqual(result["gas"]["5.0"]["all"]["Mgas_Msun"], 0)
        self.assertIsNone(result["gas"]["5.0"]["all"]["Mdense"])
        self.assertFalse(result["environment"]["supported"])
        self.assertTrue(result["source_unchanged"])

    def test_fixed_physical_density_proxies(self):
        self.write()
        result = self.run_gas()
        own = result["gas"]["5.0"]["dwarf_spatial"]
        self.assertEqual(own["Mgas_Msun"], 4)
        for key, expected in (("Mdensity_proxy_nHgt0p1_Msun", 3), ("Mdensity_proxy_nHgt0p5_Msun", 2),
                              ("Mdensity_proxy_nHgt1_Msun", 1), ("Mdensity_proxy_nHgt10_Msun", 0)):
            self.assertEqual(own[key], expected)
        self.assertIsNone(own["Mdense"])
        self.assertIsNone(result["dwarf_spatial_gas_r20"]["Mdense"])
        self.assertAlmostEqual(own["MHI"], 1.52, places=6)
        self.assertAlmostEqual(own["SFR"], 0.3, places=6)

    def test_missing_optional_fields_remain_independent(self):
        self.write(omit=("ElectronAbundance", "NeutralHydrogenAbundance", "StarFormationRate", "default"))
        result = self.run_gas()
        own = result["gas"]["5.0"]["dwarf_spatial"]
        self.assertEqual(own["Mgas_Msun"], 4)
        for key in ("MHI", "Mcold", "SFR", "mFd"):
            self.assertIsNone(own[key])
        self.assertEqual(own["Mdensity_proxy_nHgt1_Msun"], 1)
        self.assertTrue(result["environment"]["upstream20_40"]["P_ram_bulk_dyn_cm2"] > 0)
        self.assertTrue(result["environment"]["upstream20_40"]["P_thermal_dyn_cm2"] > 0)
        self.assertIsNone(result["environment"]["upstream20_40"]["P_feedback_stored_proxy_dyn_cm2"])

    def test_unknown_origin_does_not_change_spatial_membership(self):
        self.write()
        ctx = copy.deepcopy(self.ctx)
        ctx["origin_capabilities"]["gas_parent"] = False
        result = self.run_gas(ctx)
        self.assertEqual(result["gas"]["5.0"]["all"]["Mgas_Msun"], 4)
        self.assertEqual(result["gas"]["5.0"]["dwarf_spatial"]["Mgas_Msun"], 4)
        self.assertTrue(result["environment"]["supported"])
        self.assertEqual(result["status"], "COMPLETE")
        self.assertFalse(result["membership_status"]["particle_IDs_read"])

    def test_known_effective_threshold(self):
        self.write()
        ctx = copy.deepcopy(self.ctx)
        ctx["config"]["effective_sf_nH_threshold_cm3"] = 0.5
        ctx["config"]["effective_sf_threshold_evidence"] = "synthetic fixture explicit threshold"
        result = self.run_gas(ctx)
        self.assertEqual(result["gas"]["5.0"]["dwarf_spatial"]["Mdense"], 2)
        self.assertAlmostEqual(result["effective_sf_nH_threshold_cm3"], 0.5)

    def test_velocity_translation_invariance(self):
        self.write()
        original = self.run_gas()
        translation = np.array([100.0, -40.0, 80.0])
        self.gas_fields["Velocities"] = self.gas_fields["Velocities"] + translation
        self.write()
        translated = self.run_gas(bulk=np.array([10.0, 0.0, 0.0]) + translation)
        a = original["environment"]["upstream20_40"]
        b = translated["environment"]["upstream20_40"]
        self.assertAlmostEqual(a["vrel_kms"], b["vrel_kms"], places=12)
        self.assertAlmostEqual(a["P_ram_bulk_dyn_cm2"] / b["P_ram_bulk_dyn_cm2"], 1.0, places=12)
        np.testing.assert_allclose(original["environment"]["wind_direction_native"], translated["environment"]["wind_direction_native"], atol=1e-12)

    def test_missing_velocity_still_retains_thermal_pressure(self):
        self.write(omit=("Velocities",))
        result = self.run_gas()
        seed = result["environment"]["seed_shell30_60"]
        self.assertIsNone(seed["P_ram_bulk_dyn_cm2"])
        self.assertTrue(seed["P_thermal_dyn_cm2"] > 0)
        self.assertTrue(seed["nH_geom_cm3"] > 0)

    def test_missing_file_is_recorded(self):
        result = self.run_gas()
        self.assertEqual(result["status"], "PARTIAL")
        self.assertFalse(result["input_snapshot_stats"][0]["stat_available"])
        self.assertIn("error", result["input_snapshot_stats"][0])
        self.assertTrue(result["errors"])
        self.assertFalse(result["source_unchanged"])

    def test_changed_frozen_source_is_not_silently_consumed(self):
        self.write()
        ctx = copy.deepcopy(self.ctx)
        st = self.path.stat()
        ctx["input_snapshot_fingerprints"] = {"0": [dict(path=str(self.path), size_bytes=st.st_size, mtime_ns=st.st_mtime_ns - 1)]}
        result = self.run_gas(ctx)
        self.assertEqual(result["status"], "INVALID_SOURCE")
        self.assertFalse(result["source_unchanged"])
        self.assertTrue(result["input_snapshot_stats"][0]["changed_before_read"])
        self.assertEqual(result["counts"]["gas_particles_read"], 0)

    def test_host_enclosed_gas_excludes_spatial_dwarf_core(self):
        self.write()
        ctx = copy.deepcopy(self.ctx)
        ctx.update(snapshot_host_center_raw=[0.0, 0.0, 0.0], snapshot_Rgc_kpc=40.0)
        result = self.run_gas(ctx)
        self.assertEqual(result["host_enclosed_gas_mass_Msun"], 128.0)
        self.assertEqual(result["host_enclosed_gas_N"], 128)
        self.assertEqual(result["host_enclosed_gas_mass_status"], "COMPLETE")
        ctx["origin_capabilities"]["gas_parent"] = False
        result = self.run_gas(ctx)
        self.assertEqual(result["host_enclosed_gas_mass_Msun"], 128.0)
        self.assertEqual(result["host_enclosed_gas_mass_status"], "COMPLETE")

    def test_duplicate_or_absent_IDs_do_not_change_any_physical_result(self):
        self.write()
        original = self.run_gas()
        self.gas_fields["ParticleIDs"][:] = 1
        self.write()
        duplicated = self.run_gas()
        self.assertEqual(original["gas"], duplicated["gas"])
        self.assertEqual(original["environment"], duplicated["environment"])
        self.assertEqual(original["gas_flux"], duplicated["gas_flux"])
        self.assertEqual(original["source_gas_fields_sha256"], duplicated["source_gas_fields_sha256"])
        self.write(omit=("ParticleIDs",))
        absent = self.run_gas()
        self.assertEqual(original["gas"], absent["gas"])
        self.assertEqual(original["environment"], absent["environment"])
        self.assertEqual(original["source_gas_fields_sha256"], absent["source_gas_fields_sha256"])
        self.assertNotIn("ParticleIDs", absent["field_availability"])

    def test_strict_loose_membership_and_annulus_are_physical_apertures(self):
        self.gas_fields["Coordinates"][2:4, 0] = [6.0, 14.0]
        self.write()
        ctx = copy.deepcopy(self.ctx)
        ctx.update(snapshot_host_center_raw=[0.0, 0.0, 0.0], snapshot_Rgc_kpc=40.0)
        result = self.run_gas(ctx)
        self.assertEqual(result["membership_strict_r5"]["Mgas_Msun"], 2)
        self.assertEqual(result["membership_loose_r15"]["Mgas_Msun"], 4)
        self.assertEqual(result["annulus5_15"]["Mgas_Msun"], 2)
        self.assertEqual(result["annulus5_15"]["Mdensity_proxy_nHgt0p5_Msun"], 2)
        self.assertAlmostEqual(result["annulus5_15"]["nH_massmean_cm3"], 1.4, places=5)
        self.assertEqual(result["host_enclosed_gas_mass_Msun"], 128)
        self.assertIn("15.0", result["gas_flux"])
        self.assertIn("upstream15_25", result["environment"])
        self.assertNotIn("upstream10_20", result["environment"])


if __name__ == "__main__":
    unittest.main()
