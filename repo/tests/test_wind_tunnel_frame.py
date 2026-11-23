"""Regression for reconstructing a drifting stellar centre from box snapshots."""

import json
from pathlib import Path
import sys

import h5py
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from basefunc import GalaxySimulation
from wind_tunnel_frame import KPC_CM, MSUN_G, MYR_S, require_summary_frame


def test_wind_snapshot_keeps_stellar_drift_and_native_age(tmp_path):
    output = tmp_path / 'output'
    output.mkdir()
    support = tmp_path / '_support'
    support.mkdir()
    np.savetxt(support / 'orbit.dat', [
        [0, 0, 0, 0, 10, 20, 30, 0, 0, 0, -10, -20, -30, 100, 1e-15],
        [100, 2, 0, 0, 12, 24, 36, 0, 0, 0, -12, -24, -36, 200, 3e-15],
    ])
    np.savetxt(support / 'mass.dat', [[0, 0], [2, 8e8], [4, 64e8]])
    (support / 'wind_tunnel_analysis.json').write_text(json.dumps(dict(
        model_name='Fornax1', run_id='wind_case', orbit_table='orbit.dat',
        host_mass_table='mass.dat', time_begin_code=.5, unit_length_cm=2 * KPC_CM,
        unit_velocity_cm_s=2e5, unit_mass_g=MSUN_G, feedback_timescale_myr=2,
    )), encoding='utf-8')
    time = .5 + 50 * MYR_S / (KPC_CM / 1.e5)
    offsets = np.stack(np.meshgrid(*([np.arange(-3, 4) * .04] * 3)), axis=-1).reshape(-1, 3)
    core = np.array([.3, -.1, .2])
    positions = np.vstack([offsets + core, [10., 0., 0.]]) + 5.
    count = len(positions)
    with h5py.File(output / 'snapshot_007.hdf5', 'w') as h5:
        header = h5.create_group('Header')
        particle_counts = [0, 0, count, 0, 0, 0]
        header.attrs.update(Time=time, BoxSize=10., NumPart_Total=particle_counts,
                            NumPart_ThisFile=particle_counts, NumFilesPerSnapshot=1,
                            MassTable=[0., 0., 2., 0., 0., 0.], Flag_Cooling=0)
        stars = h5.create_group('PartType2')
        stars['Coordinates'] = positions
        stars['Velocities'] = np.tile([.25, -.5, .75], (count, 1))
        stars['StellarFormationTime'] = np.full(count, .4)

    simulation = GalaxySimulation(str(output), 7)
    df, tsnap = simulation.load_snapshot(particle_types=(2,))
    assert tsnap == time
    assert df.m.iloc[0] == 2.
    assert df.stellar_age_gyr.iloc[0] == pytest.approx((time - .4) * (KPC_CM / 1.e5) / (1000 * MYR_S))
    centered = simulation.center_on_MW()
    assert simulation.center_on_MW() is centered
    np.testing.assert_allclose(centered[['vx', 'vy', 'vz']].iloc[0], [11.5, 21., 34.5])
    np.testing.assert_allclose(centered[['x', 'y', 'z']].iloc[-1], [1., 0., 0.])
    host = simulation.observer_host_state
    assert host['elapsed_since_insertion_Myr'] == pytest.approx(50.)
    assert host['bulk_particle_count'] == 0
    assert simulation.wind_tunnel.enclosed_host_mass(3.) == pytest.approx(27e8)
    mw, dwarf = simulation.classify_mw_dwarf(dwarf_radius=.9)
    np.testing.assert_allclose(simulation.dwarf_center_grf_kpc, np.array([1., 0., 0.]) + 2 * core, atol=.05)
    assert not mw.any()
    assert dwarf[:-1].all()
    assert dwarf[-1]

    require_summary_frame([simulation.wind_tunnel.csv_comment()], simulation.wind_tunnel)
    with pytest.raises(ValueError, match='different simulation frame'):
        require_summary_frame([], simulation.wind_tunnel)
    with pytest.raises(ValueError, match='outside'):
        simulation.wind_tunnel.at_time(time - 60 / simulation.wind_tunnel.time_to_myr)
