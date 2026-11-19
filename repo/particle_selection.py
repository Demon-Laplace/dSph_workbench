import importlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.neighbors import KDTree

from analysis_core import Analysis as AnalysisCore
from coordinate_transform import CoordinateTransform
from canonical_observer import LENGTH_CODE_TO_KPC, VELOCITY_CODE_TO_KMS, CANONICAL_OBSERVER_ID, MEMBER_SELECTION_ID, observer_contract_sha256
from variable import fornax_core_radius


core_radius = fornax_core_radius
_SPATIAL_MODULES = None
_SPATIAL_CONTEXTS = {}


def _shared_spatial_helpers():
    """Load the maintained spatial implementation, or its job package copy."""
    global _SPATIAL_MODULES
    if _SPATIAL_MODULES is None:
        configured = os.environ.get('DSPH_SPATIAL_HELPERS_DIRECTORY')
        if configured:
            directory = Path(configured)
        else:
            code = next((p for p in Path(__file__).resolve().parents if p.name == 'Code'), None)
            directory = code / 'simulation' / 'sigma_parameters' if code else Path(__file__).resolve().parent
        if str(directory) not in sys.path:
            sys.path.insert(0, str(directory))
        spatial = importlib.import_module('spatial_context')
        physical = importlib.import_module('physical_helpers')
        _SPATIAL_MODULES = spatial, physical
    return _SPATIAL_MODULES


def _vector(value):
    try:
        result = np.asarray(value, dtype=float)
        return result if result.shape == (3,) and np.isfinite(result).all() else None
    except (TypeError, ValueError):
        return None


def _snapshot_seed(table, snapshot):
    """Reuse geometric reference points only as local spatial search seeds."""
    direct = _vector(table.get(str(int(snapshot))))
    if direct is not None:
        return direct, 'saved_snapshot_geometric_seed'
    samples = sorted((int(key), vector) for key, value in table.items()
                     if (vector := _vector(value)) is not None)
    if not samples:
        return None, None
    if len(samples) == 1:
        return samples[0][1].copy(), 'nearest_saved_geometric_seed'
    indices = np.asarray([item[0] for item in samples])
    values = np.asarray([item[1] for item in samples])
    n = int(snapshot)
    if n < indices[0]:
        seed = values[0] + (n - indices[0]) * (values[1] - values[0]) / (indices[1] - indices[0])
    elif n > indices[-1]:
        seed = values[-1] + (n - indices[-1]) * (values[-1] - values[-2]) / (indices[-1] - indices[-2])
    else:
        seed = np.array([np.interp(n, indices, values[:, axis]) for axis in range(3)])
    return seed, 'interpolated_or_endpoint_extrapolated_geometric_seed'


class ParticleSelectionMixin:
    def _spatial_seed_context(self):
        if hasattr(self, '_spatial_context'):
            return self._spatial_context
        folder = Path(self.folder_path).resolve()
        run = Path(os.environ.get('DSPH_SOURCE_RUN_DIRECTORY', str(folder.parent if folder.name == 'output' else folder)))
        explicit = os.environ.get('DSPH_SPATIAL_CONTEXT')
        candidates = [Path(explicit)] if explicit else sorted((run / 'output' / 'analysis' / 'sigma_spatial_parameters').glob('*_sigma_parameters_context.json'))
        context = {}
        if candidates:
            path = candidates[0]
            key = str(path)
            if key not in _SPATIAL_CONTEXTS:
                _SPATIAL_CONTEXTS[key] = json.loads(path.read_text(encoding='utf-8-sig'))
            context = dict(_SPATIAL_CONTEXTS[key])
            context['seed_context_path'] = str(path)
        context.setdefault('run_dir', str(run))
        # A new run's merger origin is an initial search seed, not a late-time
        # dwarf position. Per-snapshot trajectories are supplied by the batch.
        if _vector(context.get('initial_host_seed_raw')) is None:
            spatial, _ = _shared_spatial_helpers()
            _, transform = spatial.merger_parameters(Path(context['run_dir']) / 'mergecube.par')
            context['initial_host_seed_raw'] = transform[:3].tolist()
        self._spatial_context = context
        return context

    def find_MW_particles(self, dw_particles_num, cube_file):
        """Historical ID-count utility; never used for current spatial membership."""
        if cube_file is None:
            cube_count = 0
        else:
            with open(cube_file, 'r') as f:
                cube_count = sum(1 for _ in f) - 1

        df_sorted = self.df.sort_values(by='id')
        front_ids = df_sorted.head(cube_count)['id']
        back_ids = df_sorted.tail(dw_particles_num)['id']
        return ~self.df['id'].isin(front_ids) & ~self.df['id'].isin(back_ids)

    def find_dwarf_particles(self, dw_particles_num):
        """Historical ID-count utility; never used for current spatial membership."""
        if dw_particles_num <= 0:
            return pd.Series([False] * len(self.df), index=self.df.index)
        dw_all = self.df.nlargest(dw_particles_num, 'id')
        return self.df['id'].isin(dw_all['id'])

    def classify_mw_dwarf(self, r_exclude=5, dwarf_radius=5, k_density=32, df=None):
        df = self.df if df is None else df
        if not getattr(self, '_canonical_MW_centered', False):
            raise ValueError('Locate the physical host before selecting dwarf members')
        spatial, _ = _shared_spatial_helpers()
        context = self._spatial_seed_context()
        tp = df['tp'].to_numpy()
        rows = np.flatnonzero(np.isin(tp, [2, 3, 4]))
        position = df.iloc[rows][['x', 'y', 'z']].to_numpy()
        mass = df['m'].to_numpy()[rows]
        valid = np.isfinite(position).all(axis=1) & np.isfinite(mass) & (mass > 0)
        rows, p, m = rows[valid], position[valid], mass[valid]
        stellar = np.zeros(len(df), dtype=bool)
        stellar[rows] = True
        if len(p) < 2:
            raise ValueError('Fewer than two finite physical stellar rows')
        raw_seed, method = _snapshot_seed(context.get('center_seeds', {}), self.snapshot_num)
        if raw_seed is None and int(self.snapshot_num) == 0:
            raw_seed = _vector(context.get('initial_seed_raw', context.get('initial_center_raw')))
            if raw_seed is None:
                raw_seed, _ = spatial.initial_seed_from_mergers(context['run_dir'])
            method = 'initial_merger_search_seed'
        if raw_seed is not None:
            host = np.asarray(self.observer_host_state['host_center_input_code'])
            seed = (raw_seed - host) * self.observer_host_state['length_to_kpc']
            search_radius = 30. if int(self.snapshot_num) == 0 else 5.
        else:
            # This explicitly labelled property-only fallback is retained for
            # genuinely new runs without a reference trajectory. It is never
            # an ID/cohort recovery and never uses a global stellar COM.
            radius = np.linalg.norm(p, axis=1)
            outer = radius > float(r_exclude)
            if not np.any(outer):
                raise ValueError('No spatial dwarf seed outside the host exclusion sphere')
            density = AnalysisCore.estimate_local_density(p[:, 0], p[:, 1], p[:, 2], k=min(int(k_density), len(p)))
            peak = np.argmax(np.where(outer & ~np.isnan(density), density, -np.inf))
            seed = p[peak]
            method = 'untracked_current_stellar_density_peak_outside_host_exclusion'
            search_radius = max(5., float(dwarf_radius))
        center = spatial.shrinking_center(p, m, seed, 1., initial_radius_kpc=search_radius,
                                          min_radius_kpc=.3, min_count=128)
        self.dwarf_center_grf_kpc = np.asarray(center)
        self.dwarf_center_raw = (self.dwarf_center_grf_kpc / self.observer_host_state['length_to_kpc']
                                 + np.asarray(self.observer_host_state['host_center_input_code']))
        self.dwarf_center_method = dict(seed_method=method, seed_grf_kpc=np.asarray(seed).tolist(),
            search_radius_kpc=search_radius, min_radius_kpc=.3, minimum_local_stars=128,
            method='shared spatial_context.shrinking_center on current TP2/3/4 properties',
            seed_context_path=context.get('seed_context_path'))
        dwarf = np.zeros(len(df), dtype=bool)
        dwarf[rows] = np.linalg.norm(p - center, axis=1) < float(dwarf_radius)
        return stellar & ~dwarf, dwarf

    def center_on_MW(self, r_init=300.0, r_min=5.0, shrink=0.8):
        """Use the shared host-seeded centre and the unchanged strict r5 bulk."""
        df = self.df
        if getattr(self, '_canonical_MW_centered', False):
            return df
        _, physical = _shared_spatial_helpers()
        context = self._spatial_seed_context()
        host_seed, seed_method = _snapshot_seed(context.get('host_center_seeds', {}), self.snapshot_num)
        if host_seed is None:
            host_seed = _vector(context.get('initial_host_seed_raw'))
            seed_method = 'initial_host_origin_after_cube_merger'
        dwarf_seed, _ = _snapshot_seed(context.get('center_seeds', {}), self.snapshot_num)
        if dwarf_seed is None and int(self.snapshot_num) == 0:
            dwarf_seed = _vector(context.get('initial_seed_raw', context.get('initial_center_raw')))
        stars = np.isin(df['tp'].to_numpy(), [2, 3, 4])
        current = df.loc[stars, ['x', 'y', 'z', 'vx', 'vy', 'vz', 'm', 'tp']]
        properties = dict(p=current[['x', 'y', 'z']].to_numpy(),
                          v=current[['vx', 'vy', 'vz']].to_numpy(),
                          m=current['m'].to_numpy(), tp=current['tp'].to_numpy())
        units = context.get('units', {})
        length = float(units.get('length_to_kpc', context.get('length_to_kpc', LENGTH_CODE_TO_KPC)))
        velocity = float(units.get('velocity_to_kms', context.get('velocity_to_kms', VELOCITY_CODE_TO_KMS)))
        if not (np.isfinite(length) and length > 0 and np.isfinite(velocity) and velocity > 0):
            raise ValueError('Explicit positive physical length/velocity units required')
        host = physical._mw_center_bulk(properties, length, seed=host_seed,
                                        dwarf_center=dwarf_seed, return_info=True)
        if host['center_raw'] is None or host['bulk_raw'] is None:
            raise ValueError('Canonical observer host reference unavailable: ' + str(host.get('reason')))
        H, B = np.asarray(host['center_raw']), np.asarray(host['bulk_raw'])
        for axis, center in zip(('x', 'y', 'z'), H):
            df[axis] = (df[axis].to_numpy(dtype=float) - center) * length
        df['r'] = np.sqrt(df['x'] ** 2 + df['y'] ** 2 + df['z'] ** 2)
        for axis, bulk in zip(('vx', 'vy', 'vz'), B):
            df[axis] = (df[axis].to_numpy(dtype=float) - bulk) * velocity
        self.observer_host_state = dict(host_center_input_code=H.tolist(), host_bulk_input_code=B.tolist(),
            length_to_kpc=length, velocity_to_kms=velocity,
            bulk_particle_count=host['N'], host_reference_status=host['status'],
            low_bulk_sample=host.get('low_bulk_sample', False), bulk_radius_kpc=5.,
            minimum_bulk_stars=2, low_bulk_sample_reference_N=100,
            host_center_method=host['center_method'], host_seed_method=seed_method,
            host_seed_input_code=host_seed.tolist(), seed_context_path=context.get('seed_context_path'),
            canonical_observer_id=CANONICAL_OBSERVER_ID, member_selection_id=MEMBER_SELECTION_ID,
            observer_config_sha256=observer_contract_sha256())
        self._canonical_MW_centered = True
        return df
    def convert_coordinates(self, df=None, copy_df=True, include_galactic=True):
        df=self.df if df is None else df
        if copy_df:df=df.copy()
        result=self.convert_coordinates_for_mask(np.ones(len(df),dtype=bool),df=df,
            with_velocity=True,include_galactic=include_galactic)
        for key in ('xh','yh','zh','vxh','vyh','vzh','rh','ra','dec','vlos','pmra','pmdec'):
            df[key]=result[key]
        if include_galactic:df['l'],df['b']=result['l'],result['b']
        return df

    def convert_coordinates_for_mask(self, mask, df=None, with_velocity=True, include_galactic=False):
        df=self.df if df is None else df
        if not getattr(self,'_canonical_MW_centered',False):
            raise ValueError('Call center_on_MW exactly once before observation conversion')
        mask=np.asarray(mask,dtype=bool)
        selected=df.loc[mask]
        # All temperature/particle components use one complete six-dimensional contract.
        result=CoordinateTransform.observe_grf(*(selected[a].to_numpy() for a in ('x','y','z','vx','vy','vz')),
            particle_id=selected['analysis_row'].to_numpy() if 'analysis_row' in selected else np.flatnonzero(mask).astype(np.uint64))
        result['catalogue_row_index'] = np.flatnonzero(mask)
        result['particle_id_definition'] = 'analysis_row_index; source IDs are optional separate metadata'
        result['_native']['particle_id_definition'] = result['particle_id_definition']
        result['source_particle_id'] = (selected['source_particle_id'].to_numpy()
            if 'source_particle_id' in selected else selected['id'].to_numpy()
            if 'id' in selected else None)
        if result['source_particle_id'] is None:
            result['source_particle_id_blocks'] = []
            for block in getattr(self, 'source_particle_id_blocks', []):
                use = mask[block['start']:block['stop']]
                if np.any(use):
                    result['source_particle_id_blocks'].append(dict(
                        analysis_row=np.arange(block['start'], block['stop'], dtype=np.uint64)[use],
                        source_particle_id=block['values'][use]))
        result['source_id_available'] = (selected['source_id_available'].to_numpy()
            if 'source_id_available' in selected else np.full(len(selected), 'id' in selected, dtype=bool))
        if include_galactic:
            result['l'],result['b']=CoordinateTransform.to_galactic(result['x'],result['y'],result['z'])
        return result

    def classify_new_dwarf_stars(self, dw_star_mask, r_cut=core_radius):
        try:
            new_stars_mask = self.df.tp == 4

            if not new_stars_mask.any() or not dw_star_mask.any():
                empty_mask = pd.Series(False, index=self.df.index, dtype=bool)
                return empty_mask, empty_mask

            new_pos = self.df.loc[new_stars_mask, ['x', 'y', 'z']].values
            dw_pos = self.df.loc[dw_star_mask, ['x', 'y', 'z']].values

            dw_tree = KDTree(dw_pos)
            dw_dist, _ = dw_tree.query(new_pos, k=1)
            is_dw = dw_dist.squeeze() < r_cut

            new_dw_stars = pd.Series(False, index=self.df.index, dtype=bool)
            new_mw_stars = pd.Series(False, index=self.df.index, dtype=bool)

            new_dw_stars.loc[new_stars_mask] = is_dw.astype(bool)
            new_mw_stars.loc[new_stars_mask] = (~is_dw).astype(bool)

            return new_dw_stars, new_mw_stars

        except Exception as e:
            raise RuntimeError(f"Error classifying new stars: {str(e)}")

    def find_dwarf_gas(self, dw_star_mask, radius=25.0):
        df = self.df

        x = df["x"].to_numpy()
        y = df["y"].to_numpy()
        z = df["z"].to_numpy()
        m = df["m"].to_numpy()

        if dw_star_mask.sum() == 0:
            return np.zeros(len(df), dtype=bool)

        xc = np.average(x[dw_star_mask], weights=m[dw_star_mask])
        yc = np.average(y[dw_star_mask], weights=m[dw_star_mask])
        zc = np.average(z[dw_star_mask], weights=m[dw_star_mask])

        gas_sel = df["tp"].to_numpy() == 0
        dx = x - xc
        dy = y - yc
        dz = z - zc
        r2 = dx * dx + dy * dy + dz * dz

        dw_gas_mask = gas_sel & (r2 < radius ** 2)
        return dw_gas_mask

    def find_MW_gas(self, radius=200.0):
        df = self.df

        x = df["x"].to_numpy()
        y = df["y"].to_numpy()
        z = df["z"].to_numpy()
        tp = df["tp"].to_numpy()

        gas_sel = tp == 0
        r2 = x * x + y * y + z * z

        mw_gas_mask = gas_sel & (r2 < radius ** 2)
        return mw_gas_mask

    @staticmethod
    def get_dw_num(file_path):
        try:
            with open(file_path, 'r') as file:
                first_line = file.readline().strip()
                numbers = [
                    float(num)
                    for num in first_line.replace(',', ' ').split()
                    if num.replace('.', '', 1).isdigit()
                ]
                total_sum = int(sum(numbers))

                if total_sum <= 0:
                    raise ValueError("No valid particle count found")

                return total_sum

        except FileNotFoundError:
            print(f"Error: File not found - {file_path}")
            return 0
        except (ValueError, TypeError) as e:
            print(f"Error parsing dwarf particle count: {e}")
            return 0
        except Exception as e:
            print(f"Unexpected error reading dwarf config: {e}")
            return 0
