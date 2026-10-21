import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch
import io
import ast
import math

import h5py
import numpy as np

from snapshot_metrics import hi_particle_masses, hi_fraction, _sum_hi_mass, require_hi_fraction_csv
from snapshot_io import SnapshotLoaderMixin


class HIMassTests(unittest.TestCase):
    def test_neutral_number_fraction_is_converted_to_hydrogen_mass(self):
        result=hi_particle_masses([100,100,100],[1,0.5,0])
        np.testing.assert_allclose(result,[76,38,0],rtol=1e-7)
        self.assertEqual(result.dtype,np.float32)
        self.assertAlmostEqual(hi_fraction(38,62),0.38)
        self.assertTrue(np.isnan(hi_fraction(0,0)))

    def test_hi_sums_both_temperature_phases_and_respects_aperture(self):
        # Cold HI=0.76*10*.5, warm HI=0.76*100*.01; a massive outside particle is excluded.
        cold=_sum_hi_mass(np.array([10,10000]),np.array([.5,1]),np.array([True,False]))
        warm=_sum_hi_mass(np.array([100]),np.array([.01]),np.array([True]))
        self.assertAlmostEqual(cold+warm,4.56,places=6)
        self.assertAlmostEqual(hi_fraction(cold+warm,10),4.56/14.56,places=7)

    def test_fraction_limits_and_invalid_masses(self):
        self.assertEqual(hi_fraction(0, 62), 0)
        self.assertEqual(hi_fraction(38, 0), 1)
        for hi, stars in [(0, 0), (-1, 2), (1, -2), (np.nan, 2), (1, np.inf)]:
            self.assertTrue(np.isnan(hi_fraction(hi, stars)))

    def test_csv_definition_rejects_legacy_and_missing_headers(self):
        good = "# hi_fraction=M_HI/(M_star+M_HI) in the same aperture; gas_mass includes all gas\n"
        with patch("builtins.open", return_value=io.StringIO(good + "numsp\n285\n")):
            require_hi_fraction_csv("unused.csv")
        for text in ["# hi_fraction=M_HI/M_star in the same aperture;\n", "numsp\n285\n"]:
            with patch("builtins.open", return_value=io.StringIO(text)):
                with self.assertRaises(ValueError):
                    require_hi_fraction_csv("unused.csv")

    def test_ic_report_uses_the_same_bounded_fraction(self):
        source = Path(__file__).resolve().parents[1] / "build_production_ic" / "cli.py"
        tree = ast.parse(source.read_text(encoding="utf-8"))
        assignments = {target.id: node.value for node in ast.walk(tree)
                       if isinstance(node, ast.Assign) for target in node.targets
                       if isinstance(target, ast.Name) and target.id in ("hi_denominator", "gas_fraction")}
        for hi, stars in [(38, 62), (0, 62), (38, 0), (0, 0), (-1, 2), (1, np.inf)]:
            ns = dict(math=math, hi_mass=hi, star_mass=stars)
            for name in ("hi_denominator", "gas_fraction"):
                ns[name] = eval(compile(ast.Expression(assignments[name]), str(source), "eval"), ns)
            expected = hi_fraction(hi, stars)
            if np.isnan(expected):
                self.assertTrue(np.isnan(ns["gas_fraction"]))
            else:
                self.assertAlmostEqual(ns["gas_fraction"], expected)

    def test_invalid_abundances_are_not_silently_treated_as_hi(self):
        for neutral in ([-0.1],[1.1],[np.nan]):
            with self.assertRaises(ValueError):hi_particle_masses([1],neutral)
        with self.assertRaises(ValueError):hi_particle_masses([1,2],[0.5])

    def test_hi_loading_does_not_depend_on_electron_field(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'snapshot_000.hdf5'
            with h5py.File(path,'w') as f:
                header=f.create_group('Header')
                for name,value in {'Time':0,'NumPart_Total':[2,0,0,0,0,0],
                    'NumPart_ThisFile':[2,0,0,0,0,0],'MassTable':[0]*6,'NumFilesPerSnapshot':1,'Flag_Cooling':1}.items():header.attrs[name]=value
                gas=f.create_group('PartType0')
                gas['Coordinates']=np.zeros((2,3));gas['Velocities']=np.zeros((2,3))
                gas['ParticleIDs']=[1,2];gas['Masses']=[1,2];gas['NeutralHydrogenAbundance']=[.25,.75]
            class Loader(SnapshotLoaderMixin):
                boxsize=0
                def _resolve_snapshot_hdf5(self):return str(path),str(path.with_suffix('')),'.hdf5'
            frame,_=Loader()._load_snapshot_direct_hdf5(dftype=2,box=False,gas_fields=('nh',),particle_types=(0,))
            np.testing.assert_allclose(frame.nh,[.25,.75])
            with h5py.File(path,'a') as f:del f['PartType0/NeutralHydrogenAbundance']
            with self.assertRaisesRegex(ValueError,'NeutralHydrogenAbundance'):
                Loader()._load_snapshot_direct_hdf5(dftype=2,box=False,gas_fields=('nh',),particle_types=(0,))


if __name__=='__main__':unittest.main()
