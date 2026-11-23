"""The fixed-axis, translating GIZMO wind-tunnel snapshot frame."""

from functools import lru_cache
import json
from pathlib import Path

import numpy as np


KPC_CM = 3.0856775814913673e21
MSUN_G = 1.989e33
MYR_S = 3.15576e13
WIND_CSV_COLUMNS = (
    'wind_elapsed_Gyr',
    'wind_prescribed_Pram_dyn_cm2',
    'wind_upstream12_20_Pram_dyn_cm2',
    'wind_upstream20_40_Pram_dyn_cm2',
    'wind_dwarf_gas_r2_msun',
    'wind_dwarf_gas_r5_msun',
    'wind_ambient_gas_r5_msun',
)


def _table(path, columns):
    rows = np.loadtxt(path, ndmin=2)
    if (rows.shape[0] < 2 or rows.shape[1] != columns or not np.isfinite(rows).all()
            or rows[0, 0] != 0 or np.any(np.diff(rows[:, 0]) <= 0)):
        raise ValueError(f'Expected {columns} finite columns with an increasing first column starting at zero: {path}')
    return rows


class WindTunnelFrame:
    """Run inputs for an axis-aligned periodic box with a prescribed orbit."""

    def __init__(self, path):
        self.path = path
        config = json.loads(path.read_text(encoding='utf-8-sig'))
        self.model_name = config['model_name']
        self.run_id = config['run_id']
        if not isinstance(self.run_id, str) or not self.run_id or any(c.isspace() for c in self.run_id):
            raise ValueError(f'Expected a non-empty wind tunnel run_id without whitespace in {path}')
        self.time_begin_code = float(config['time_begin_code'])
        length_cm = float(config['unit_length_cm'])
        velocity_cm_s = float(config['unit_velocity_cm_s'])
        mass_g = float(config['unit_mass_g'])
        if not all(np.isfinite(unit) and unit > 0 for unit in (length_cm, velocity_cm_s, mass_g)):
            raise ValueError(f'Wind tunnel length, velocity and mass units must be positive: {path}')
        self.length_to_kpc = length_cm / KPC_CM
        self.velocity_to_kms = velocity_cm_s / 1.e5
        self.mass_to_msun = mass_g / MSUN_G
        self.time_to_myr = length_cm / velocity_cm_s / MYR_S
        self.feedback_timescale_myr = float(config['feedback_timescale_myr'])
        self.orbit_path = (path.parent / config['orbit_table']).resolve()
        self.orbit = _table(self.orbit_path, 15)
        self.host_mass_path = self.host_mass = None
        if config.get('host_mass_table'):
            self.host_mass_path = (path.parent / config['host_mass_table']).resolve()
            self.host_mass = _table(self.host_mass_path, 2)
        self.gas_origin = config.get('gas_origin')
        self.free_flow_half_width_fraction = config.get('free_flow_half_width_fraction')

    def at_time(self, time_code):
        elapsed = (float(time_code) - self.time_begin_code) * self.time_to_myr
        times = self.orbit[:, 0]
        # Match the solver's tolerance for the rounded first and final epochs.
        if not np.isfinite(elapsed) or not times[0] - .001 <= elapsed <= times[-1] + .001:
            raise ValueError(f'Snapshot time {time_code} lies outside the wind tunnel orbit in {self.orbit_path}')
        elapsed = float(np.clip(elapsed, times[0], times[-1]))
        values = [np.interp(elapsed, times, column) for column in self.orbit[:, 1:].T]
        return dict(
            elapsed_since_insertion_Myr=elapsed,
            box_center_grf_kpc=np.asarray(values[:3]),
            box_velocity_grf_kms=np.asarray(values[3:6]),
            prescribed_wind_box_kms=np.asarray(values[9:12]),
            prescribed_density_msun_kpc3=values[12],
            prescribed_thermal_pressure_dyn_cm2=values[13],
        )

    def enclosed_host_mass(self, galactocentric_radius_kpc):
        if self.host_mass is None:
            return np.nan
        # The moving-box force interpolates the tabulated mass in r^3.
        radius, mass = self.host_mass.T
        return float(np.interp(galactocentric_radius_kpc ** 3, radius ** 3, mass, right=np.nan))

    def csv_comment(self):
        return f'# wind_tunnel_run_id={self.run_id}\n'


@lru_cache(maxsize=8)
def _load_frame(path):
    return WindTunnelFrame(path)


def wind_tunnel_frame(folder_path):
    path = Path(folder_path).resolve().parent / '_support' / 'wind_tunnel_analysis.json'
    return _load_frame(path) if path.is_file() else None


def require_summary_frame(comments, frame):
    found = [line.strip() for line in comments if line.startswith('# wind_tunnel_run_id=')]
    expected = [frame.csv_comment().strip()] if frame else []
    if found != expected:
        raise ValueError('The elinfo CSV belongs to a different simulation frame or wind tunnel run. Use a fresh analysis output path.')


def box_displacement(raw_coordinate, box_size_code):
    displacement = raw_coordinate - .5 * box_size_code
    return displacement - box_size_code * np.rint(displacement / box_size_code)


def measure_wind_snapshot(simulation, radius_squared, stellar_rows):
    """Measure the local upstream flow and gas origins about the measured star centre."""
    frame = simulation.wind_tunnel
    host = simulation.observer_host_state
    box_center = host['box_center_grf_kpc']
    box_velocity = host['box_velocity_grf_kms']
    wind = host['prescribed_wind_box_kms']
    speed = np.linalg.norm(wind)
    df = simulation.df
    inner_stars = stellar_rows & (radius_squared < 25.)
    velocity_columns = ['vx', 'vy', 'vz']
    star_bulk = np.average(df.loc[inner_stars, velocity_columns], axis=0, weights=df.loc[inner_stars, 'm'])
    center = simulation.dwarf_center_grf_kpc
    record = dict(
        stellar_center_grf_kpc=center,
        stellar_center_box_kpc=center - box_center,
        stellar_bulk_r5_grf_kms=star_bulk,
        stellar_bulk_r5_box_kms=star_bulk - box_velocity,
        prescribed_Pram_dyn_cm2=host['prescribed_density_msun_kpc3'] * (MSUN_G / KPC_CM**3) * (speed * 1.e5)**2,
        host_mass_table=str(frame.host_mass_path) if frame.host_mass_path else None,
    )
    origin = frame.gas_origin
    if origin is None:
        return record

    gas = df['tp'].to_numpy() == 0
    covered = np.zeros(len(df), dtype=bool)
    ambient = np.zeros(len(df), dtype=bool)
    mask = np.uint64((1 << int(origin['parent_id_bits'])) - 1)
    first_ambient_id = int(origin['ambient_first_id'])
    for block in simulation.source_particle_id_blocks:
        labels = block['values']
        if labels.dtype.kind in 'iu':
            block_rows = slice(block['start'], block['stop'])
            ambient[block_rows] = (labels.astype(np.uint64, copy=False) & mask) >= first_ambient_id
            covered[block_rows] = True
    if not covered[gas].all():
        record['gas_origin_status'] = 'source labels unavailable for some gas rows'
        return record
    ambient &= gas
    dwarf = gas & ~ambient
    mass = df['m'].to_numpy()
    record.update(gas_origin_status='classified by the configured parent-ID convention',
                  gas_origin=origin, origin_apertures={})
    for radius in (0.5, 1, 2, 5, 15):
        inside = radius_squared < radius ** 2
        record['origin_apertures'][str(radius)] = dict(
            dwarf_gas_msun=mass[dwarf & inside].sum(), ambient_gas_msun=mass[ambient & inside].sum())

    if speed == 0:
        return record
    # This is the same geometric cone estimator used in the live-MW comparison.
    # The snapshot velocities and star bulk are in the same fixed-axis frame.
    from gas_parameters import _ambient

    extent = np.max(np.abs(center - box_center))
    box_width_kpc = simulation.boxsize * frame.length_to_kpc
    probe_rows = ambient & (radius_squared < 40.**2)
    selected = df.loc[probe_rows]
    distance = np.sqrt(radius_squared[probe_rows])
    gas_data = dict(
        r=distance, x=selected[['x', 'y', 'z']].to_numpy() - center,
        m=selected['m'].to_numpy(), v=selected[velocity_columns].to_numpy(),
        bulk=star_bulk, u=selected['ugas'].to_numpy(), T=selected['temp'].to_numpy(),
        fd=np.broadcast_to(np.nan, distance.shape), UV2=(frame.velocity_to_kms * 1.e5)**2,
    )
    for inner, outer in ((12, 20), (20, 40)):
        label = f'upstream{inner}_{outer}'
        if extent + outer > .5 * box_width_kpc:
            record[label] = dict(supported=False, reason='The bounding sphere extends outside the physical box')
            continue
        probe = _ambient(gas_data, inner, outer, axis=-wind / speed)
        relative = probe.get('relative_velocity_native_kms')
        if relative is not None:
            probe['gas_minus_stars_kms'] = [-component for component in relative]
        if frame.free_flow_half_width_fraction is not None:
            probe['bounding_sphere_within_free_flow'] = (
                extent + outer <= frame.free_flow_half_width_fraction * box_width_kpc)
        record[label] = probe
    return record
