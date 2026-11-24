"""Regression for reconstructing a drifting stellar centre from box snapshots."""

import json
from pathlib import Path
import sys
from types import SimpleNamespace

import h5py
import numpy as np
import pandas as pd
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


def test_wind_sfh_uses_member_birth_times_and_partial_final_bin():
    from PlotFig import _wind_sfr_series, SFH_SFR_UNIT_MSUN_PER_YEAR

    frame = SimpleNamespace(time_to_myr=1000., time_begin_code=.5)
    snapshot = dict(
        simulation=SimpleNamespace(wind_tunnel=frame), tsnap=.623,
        total_dw_star_mask=np.array([True, True, True, True, True, False, True]),
        df=pd.DataFrame(dict(
            tp=[2, 4, 4, 4, 4, 4, 4],
            birth=[.42, .475, .515, .615, np.nan, .60, .64],
            stellar_age_gyr=[.203, .148, .108, .008, .063, .023, -1.],
            m=[900., 50., 75., 23., 10., 130., 300.],
        )),
    )
    edges, sfr = _wind_sfr_series(snapshot, .123)
    np.testing.assert_allclose(edges, [-.05, 0., .05, .10, .123])
    recovered_mass = sfr * np.diff(edges) * 1.e9 * SFH_SFR_UNIT_MSUN_PER_YEAR
    np.testing.assert_allclose(recovered_mass, [50., 75., 10., 23.])


def test_wind_sfh_observations_use_the_explicit_lookback_reference():
    from PlotFig import _observed_sfh_on_wind_clock

    edges, sfr = _observed_sfh_on_wind_clock(pd.DataFrame(dict(age=[.25, 1.25], sfr=[1., 3.])), .6)
    np.testing.assert_allclose(edges, [-.9, -.4, .1, .6])
    np.testing.assert_allclose(sfr, [3., np.nan, 1.])


def test_wind_scene_projects_gas_across_the_whole_box():
    from matplotlib.collections import QuadMesh
    from matplotlib.figure import Figure
    from PlotFig import WIND_SURFACE_DENSITY_LIMITS, _plot_scene_panel

    # The two exterior gas cells lie far beyond the stellar-centred aperture.
    # The non-gas row lies inside the box and carries a larger mass.
    snapshot = dict(
        wind_tunnel={}, cy=11., cz=-21.,
        observer_host_state=dict(box_center_grf_kpc=[0., 10., -20.],
                                 box_size_code=50., length_to_kpc=2.),
        membership_state=dict(gas_temperature_split_K=20000.),
        df=pd.DataFrame(dict(
            tp=[0, 0, 0, 2],
            y=[11., 59., -39.8, 11.], z=[-21., 29.8, -69., -21.],
            m=[2., 6., 9., 1000.], temp=[100., 1.e6, 1.e6, 0.],
        )),
    )
    fig = Figure()
    ax = fig.subplots()
    _plot_scene_panel(fig, ax, snapshot, np.array([11.]), np.array([-21.]))
    mesh = next(artist for artist in ax.collections if isinstance(artist, QuadMesh))
    assert mesh.get_array().sum() * 2.**2 == pytest.approx(17.)
    assert ax.get_xlim() == ax.get_ylim() == (-50., 50.)
    assert mesh.get_clim() == WIND_SURFACE_DENSITY_LIMITS


def test_wind_hi_panel_keeps_its_scale_when_hi_fades():
    from matplotlib.figure import Figure
    from PlotFig import WIND_HI_COLOR_MAX, WIND_HI_HALF_WIDTH_KPC, _plot_hi_sky_panel

    fig = Figure()
    panels = fig.subplots(1, 2)
    x, y = np.meshgrid(np.linspace(-4., 4., 3), np.linspace(-4., 4., 3))
    for ax, peak in zip(panels, (8.e20, 0.)):
        column = np.array([[0., 0., 0.], [0., peak, 0.], [0., 0., 0.]])
        _plot_hi_sky_panel(fig, ax, x, y, column, WIND_HI_HALF_WIDTH_KPC, wind_tunnel=True)
        assert ax.get_xlim() == (WIND_HI_HALF_WIDTH_KPC, -WIND_HI_HALF_WIDTH_KPC)
        assert ax.get_ylim() == (-WIND_HI_HALF_WIDTH_KPC, WIND_HI_HALF_WIDTH_KPC)
        assert ax.collections[0].get_clim() == (0., WIND_HI_COLOR_MAX)
    assert len(fig.axes) == 4  # Both panels reserve the same space for a colour bar.


def test_wind_sfh_moves_the_time_marker_on_one_reference_history():
    from matplotlib.figure import Figure
    from PlotFig import _plot_wind_sfh_panel

    fig = Figure()
    early, late = fig.subplots(1, 2)
    series = np.array([-.05, 0., .05, .1]), np.array([0., 50., 20.])
    _plot_wind_sfh_panel(early, series, .01, (-.5, 3.))
    _plot_wind_sfh_panel(late, series, 2.8, (-.5, 3.))
    assert early.get_xlim() == late.get_xlim() == (-.5, 3.)
    assert early.get_ylim() == late.get_ylim()
    np.testing.assert_allclose(early.patches[0].get_data().values, late.patches[0].get_data().values)
    assert early.lines[-1].get_xdata()[0] == .01
    assert late.lines[-1].get_xdata()[0] == 2.8


def test_wind_fixed_view_uses_circular_physical_annuli():
    from PlotFig import _fixed_view_light_profile
    from basefunc import Analysis

    center = np.array([4., -5., 2.])
    positions = center + np.array([[.2, 0., 5.], [0., .3, -6.], [.6, .8, 3.]])
    mass = np.array([1.e5, 2.e5, 4.e5])
    edges = [0., .25, .5, 1.2]
    fixed = _fixed_view_light_profile(positions, mass, center, np.array([0., 0., 1.]), 100., r_edges=edges)
    _, expected = Analysis.radial_magnitude_profile([.2, .3, 1.], [0., 0., 0.], mass, 300., ep=0., r_edges=edges)
    np.testing.assert_allclose(fixed, expected)


def test_wind_media_identity_tracks_the_summary_and_sfh_reference(tmp_path, monkeypatch):
    from canonical_observer import canonical_csv_comment
    from display_outputs import WIND_RENDER_REVISION, display_identity

    monkeypatch.setenv('DSPH_MEDIA_REUSE', 'scientific')
    summary = tmp_path / 'elinfo.csv'
    header = canonical_csv_comment() + '# wind_tunnel_run_id=wind_case\n'
    summary.write_text(header + 'numsp\n1\n', encoding='utf-8')
    before = display_identity(summary, wind_sfh_snapshot=1)
    assert before['render_revision'] == WIND_RENDER_REVISION
    assert before['wind_sfh_reference_snapshot'] == 1
    summary.write_text(header + 'numsp\n1\n2\n', encoding='utf-8')
    after = display_identity(summary, wind_sfh_snapshot=2)
    assert after['wind_summary_sha256'] != before['wind_summary_sha256']
