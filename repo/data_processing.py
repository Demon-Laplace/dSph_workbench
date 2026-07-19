import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

from variable import gizmo_config, orbit_file

ELINFO_DTYPE_MAP = {
    'star_mass': np.float32,
    'star_half_mass': np.float32,
    'hotgas_mass': np.float32,
    'hotgas_half_mass': np.float32,
    'coldgas_mass': np.float32,
    'coldgas_half_mass': np.float32,
    'mw_mass_r': np.float32,
    'eps': np.float32,
    'pa': np.float32,
    'gas_density': np.float32,
    'vr': np.float32,
    'vtheta': np.float32,
    'vphi': np.float32,
    'distance_gal': np.float32,
    'rhalf': np.float32,
    'rhalf_circularized': np.float32,
    'shape_center_x_kpc': np.float32,
    'shape_center_y_kpc': np.float32,
    'distance': np.float32,
    'age': np.float32,
    'sigma': np.float32,
    'tsigma': np.float32,
    'tsigma_xyz': np.float32,
    'sigma_x': np.float32,
    'sigma_y': np.float32,
    'sigma_z': np.float32,
    'pmra': np.float32,
    'pmdec': np.float32,
    'cold_gas_center_ra': np.float32,
    'cold_gas_center_dec': np.float32,
    'numsp': np.int32,
}


class DataProcessor:
    MODEL_NAME_PATTERN = re.compile(r'([A-Za-z]+)(\d+)(?:_(.+))?$')
    ELINFO_DTYPE_MAP = ELINFO_DTYPE_MAP
    IC_SNAPSHOT_INI_TEMPLATE = r"IC_snapshot(?P<snapshot>\d+)_{dwarf_name}(?P<model>\d+)\.ini"

    @staticmethod
    def read_feedback_value(file_path=gizmo_config):
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith('%'):
                    continue
                if line.startswith("timescale_fd"):
                    parts = line.split()
                    if len(parts) > 1:
                        return float(parts[1].split(";")[0])
        return None

    @staticmethod
    def read_sfreff_value(file_path=gizmo_config):
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith('%'):
                    continue
                if line.startswith("sfr_effiency"):
                    parts = line.split()
                    if len(parts) > 1:
                        return float(parts[1].split(";")[0])
        return None

    @staticmethod
    def get_ini_data(elinfo_path):
        with open(elinfo_path, 'r') as f:
            line1 = f.readline().strip().removeprefix("# Orbit Data:").strip()
            x, y, z, vx, vy, vz = map(float, line1.split(','))

            line2 = f.readline().strip().removeprefix("# Timescale_fd:").strip()
            fd_value = float(line2)

        return x, y, z, vx, vy, vz, fd_value

    @staticmethod
    def read_orbit_parameters(orbit_file=orbit_file):
        if not os.path.exists(orbit_file):
            print(f"Orbit file not found: {orbit_file}")
            return None

        with open(orbit_file, 'r') as f:
            f.readline()
            second_line = f.readline().strip()
            columns = second_line.split()
            if len(columns) < 6:
                print("Orbit file does not have enough columns.")
                return None
            return [float(val) for val in columns[:6]]

    @staticmethod
    def calsfh(sfr_csv, Tage):
        sfr_idw = sfr_csv.copy()
        sfr_idw['sfr'] = sfr_idw['sfr'] * 1e-4
        sfr_idw['sfr0'] = sfr_idw['sfr'] / sfr_idw['sfr'].sum()

        sfr_idw['dage'] = sfr_idw['age'].shift(-1) - sfr_idw['age']
        sfr_idw.iloc[-1, sfr_idw.columns.get_loc('dage')] = sfr_idw.iloc[-2]['dage']
        sfr_idw['mt'] = sfr_idw['sfr'] * sfr_idw['dage']

        max_age = sfr_idw['age'].max()
        sfr_idw['cosmic_time'] = max_age - sfr_idw['age']
        sfr_idw['cosmic_time_shifted'] = sfr_idw['cosmic_time'] - sfr_idw['cosmic_time'].max() + Tage

        sfr_idw.sort_values('cosmic_time_shifted', inplace=True)
        sfr_idw['cumulative_mass'] = sfr_idw['mt'].cumsum()

        mask_positive = sfr_idw['cosmic_time_shifted'] > 0
        if mask_positive.any():
            min_positive_time = sfr_idw.loc[mask_positive, 'cosmic_time_shifted'].min()
            mt_at_min = sfr_idw.loc[
                sfr_idw['cosmic_time_shifted'] == min_positive_time,
                'cumulative_mass',
            ].values[0]
        else:
            mt_at_min = 0.0

        sfr_idw['pure_sfh'] = sfr_idw['cumulative_mass'] - mt_at_min
        return sfr_idw

    @staticmethod
    def GetSFR(sfr_csv, bins=None):
        sfr = sfr_csv.copy()

        sfr['dage'] = sfr['age'].shift(-1) - sfr['age']
        sfr.iloc[-1, sfr.columns.get_loc('dage')] = sfr.iloc[-2]['dage']

        lookback = sfr['age'].to_numpy()
        sfr_values = sfr['sfr'].to_numpy()

        if bins is None:
            bin_width = 0.5
            bins = np.arange(0 - bin_width, np.max(lookback), bin_width)

        sfr_in_bin, edges = np.histogram(lookback, bins=bins, weights=sfr_values)
        count_in_bin, _ = np.histogram(lookback, bins=bins)
        avg_sfr = np.divide(
            sfr_in_bin,
            count_in_bin,
            out=np.zeros_like(sfr_in_bin, dtype=float),
            where=count_in_bin > 0,
        )

        return edges[:-1], avg_sfr

    @staticmethod
    def parse_model_name(path=None):
        model_info = DataProcessor.parse_model_name_details(path=path)
        return model_info['modelname'], model_info['dwarf_name'], model_info['model_num']

    @staticmethod
    def parse_model_name_details(path=None):
        modelname = (Path.cwd() if path is None else Path(path).resolve()).name
        match = DataProcessor.MODEL_NAME_PATTERN.fullmatch(modelname)
        if not match:
            raise ValueError(
                f"Invalid model name format: '{modelname}'. "
                "Expected <dwarf_name><number> or <dwarf_name><number>_<description>."
            )
        dwarf_name, num_str, description = match.groups()
        return {
            'modelname': modelname,
            'base_modelname': f'{dwarf_name}{num_str}',
            'dwarf_name': dwarf_name,
            'model_num': int(num_str),
            'description': description,
        }

    @staticmethod
    def resolve_dwarf_ini_path(dwarf_name, directory=None):
        search_dir = Path.cwd() if directory is None else Path(directory)
        default_ini = search_dir / f'{dwarf_name}.ini'
        if default_ini.exists():
            return default_ini

        ic_pattern = re.compile(
            DataProcessor.IC_SNAPSHOT_INI_TEMPLATE.format(dwarf_name=re.escape(dwarf_name))
        )
        matched_ic_files = sorted(
            path for path in search_dir.iterdir()
            if path.is_file() and ic_pattern.fullmatch(path.name)
        )

        if len(matched_ic_files) == 1:
            return matched_ic_files[0]
        if len(matched_ic_files) > 1:
            matched_names = ', '.join(path.name for path in matched_ic_files)
            raise ValueError(
                f"Multiple IC snapshot INI files found for {dwarf_name} under {search_dir}: {matched_names}"
            )

        raise FileNotFoundError(
            f"INI file not found: expected {default_ini.name} or exactly one "
            f"{DataProcessor.IC_SNAPSHOT_INI_TEMPLATE.format(dwarf_name=dwarf_name)} under {search_dir}"
        )

    @staticmethod
    def read_csv_with_comments(file_path, dtype=None):
        return pd.read_csv(file_path, comment='#', dtype=dtype)

    @staticmethod
    def format_float_for_csv(value, precision=8):
        value = float(value)
        if np.isnan(value) or np.isinf(value):
            return str(value)
        return format(value, f'.{precision}g')

    @staticmethod
    def format_row_for_csv(row, float_precision=8):
        formatted = []
        for value in row:
            if isinstance(value, (int, np.integer)):
                formatted.append(int(value))
            elif isinstance(value, (float, np.floating)):
                formatted.append(DataProcessor.format_float_for_csv(value, precision=float_precision))
            else:
                formatted.append(value)
        return formatted

    @staticmethod
    def list_snapshot_numbers(folder):
        folder_path = Path(folder)
        numbers = []
        for file_path in folder_path.glob("snapshot_*.hdf5"):
            match = re.match(r"snapshot_(\d+)\.hdf5", file_path.name)
            if match:
                numbers.append(int(match.group(1)))
        return sorted(numbers)

    @staticmethod
    def max_snapshot_number(folder):
        numbers = DataProcessor.list_snapshot_numbers(folder)
        return numbers[-1] if numbers else None

    @staticmethod
    def find_snapshot_range(folder):
        return DataProcessor.max_snapshot_number(folder)

    @staticmethod
    def sci_notation(value, precision=2):
        if value == 0:
            return f"{0:.{precision}f} × 10⁰"
        exponent = int(np.floor(np.log10(abs(value))))
        base = value / (10 ** exponent)

        superscript_map = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
        exponent_str = str(exponent).translate(superscript_map)

        return f"{base:.{precision}f} × 10{exponent_str}"
