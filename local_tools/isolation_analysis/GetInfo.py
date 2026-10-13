import numpy as np
# import pandas as pd
import os
from multiprocessing import get_context
import argparse
from basefunc import DataProcessor, GalaxySimulation, Analysis, CoordinateTransform
from variable import *
from pathlib import Path
from model_paths import model_id
import sys
import csv
import re
import gc

class GetInfo:

    def __init__(self):
        # Get model name from current directory
        self.modelname = model_id()
        match = re.match(r'([A-Za-z]+)(\d+)', self.modelname)
        if not match:
            raise ValueError(f"Invalid model name format: '{self.modelname}'. "
                           f"Expected format: <dwarf_name><number>")

        # Extract dwarf name and number
        self.dwarf_name, num_str = match.groups()
        self.num = int(num_str)
        
        # Set up file paths using Path objects
        self.inifile = Path(f'./{self.dwarf_name}.ini')
        if not self.inifile.exists():
            raise FileNotFoundError(f"INI file not found: {self.inifile}")
        
        # Get dwarf parameters
        self.dw_particles = GalaxySimulation.get_dw_num(self.inifile)
        if self.dw_particles <= 0:
            raise ValueError(f"Invalid number of dwarf particles: {self.dw_particles}")
            
        self.core_radius = core_radius_dict.get(self.dwarf_name)
        if self.core_radius is None:
            raise ValueError(f"Core radius not found for dwarf: {self.dwarf_name}")
            
        # Set output file path
        self.output_csv_file = Path(f'./elinfo_{self.modelname}.csv')

    def get_info(self, numsp):

        simulation = GalaxySimulation(folder_path=folder_path, snapshot_num=numsp)
        df, tsnap = simulation.load_snapshot()
        mw_mask = simulation.find_MW_particles(dw_particles_num=self.dw_particles, cube_file=cube_path)
        dw_mask = simulation.find_dwarf_particles(dw_particles_num=self.dw_particles)

        try:
            # Only attempt centering if we have MW particles
            if mw_mask.any():
                df = simulation.center_on_MW(mw_mask, iterations=2, radius=100)
            else:
                print(f"Warning: No MW particles found in snapshot {numsp}, skipping centering")
        except Exception as e:
            print(f"Warning: Failed to center on MW in snapshot {numsp}: {str(e)}")

        df = simulation.convert_coordinates(df=df)

        mw_star_mask = mw_mask & ((df.tp == 2) | (df.tp == 3))
        dw_star_mask = dw_mask & ((df.tp == 2) | (df.tp == 3))
        dw_gas_mask = dw_mask & (df.tp == 0)
        mw_gas_mask = mw_mask & (df.tp == 0)
        # new_star_mask = df.tp == 4
        dw_hot_gas_mask = dw_gas_mask & (df.temp >= 20000)
        dw_cold_gas_mask = dw_gas_mask & (df.temp < 20000)
        new_dw_star_mask, new_mw_star_mask = simulation.classify_new_dwarf_stars(dw_star_mask)
        total_dw_star_mask = dw_star_mask | new_dw_star_mask
        total_mw_star_mask = mw_star_mask | new_mw_star_mask 
        
        dw_cra, dw_cdec = Analysis.find_center(df[total_dw_star_mask].ra, df[total_dw_star_mask].dec)
        rotra_all, rotdec_all = Analysis.rotate_to_sky(df.ra, df.dec, center=[dw_cra, dw_cdec])

        d_mean = np.mean(df[total_dw_star_mask].rh)
        d_mean_gal = np.mean(df[total_dw_star_mask].r)

        mw_star_r = total_mw_star_mask & (df['r'] <= d_mean_gal)
        mw_gas_r = mw_gas_mask & (df['r'] <= d_mean_gal)
        mass_mw_star_r = np.sum(df.loc[mw_star_r, 'm'])
        mass_mw_gas_r = np.sum(df.loc[mw_gas_r, 'm'])

        mw_mass_r = mass_mw_star_r + mass_mw_gas_r

        rotra_dw_star, rotdec_dw_star = rotra_all[total_dw_star_mask], rotdec_all[total_dw_star_mask]
        rotra_dw_hot_gas, rotdec_dw_hot_gas = rotra_all[dw_hot_gas_mask], rotdec_all[dw_hot_gas_mask]
        rotra_dw_cold_gas, rotdec_dw_cold_gas = rotra_all[dw_cold_gas_mask], rotdec_all[dw_cold_gas_mask]

        x_kpc = np.asarray(np.radians(rotra_dw_star) * d_mean)
        y_kpc = np.asarray(np.radians(rotdec_dw_star) * d_mean)
        hot_gas_x_kpc = np.asarray(np.radians(rotra_dw_hot_gas) * d_mean)
        hot_gas_y_kpc = np.asarray(np.radians(rotdec_dw_hot_gas) * d_mean)
        cold_gas_x_kpc = np.asarray(np.radians(rotra_dw_cold_gas) * d_mean)
        cold_gas_y_kpc = np.asarray(np.radians(rotdec_dw_cold_gas) * d_mean)

        #--- Get total mass in dwarf galaxy for different components ---#

        star_points = np.column_stack([x_kpc, y_kpc])
        cold_points = np.column_stack([cold_gas_x_kpc, cold_gas_y_kpc])
        hot_points = np.column_stack([hot_gas_x_kpc, hot_gas_y_kpc])

        boundary_path = Analysis.estimate_galaxy_boundary(x_kpc, y_kpc)

        star_in = boundary_path.contains_points(star_points)
        cold_in = boundary_path.contains_points(cold_points)
        hot_in = boundary_path.contains_points(hot_points)

        M_star = df.loc[total_dw_star_mask, 'm'][star_in].sum()
        M_cold_gas = df.loc[dw_cold_gas_mask, 'm'][cold_in].sum()
        M_hot_gas = df.loc[dw_hot_gas_mask, 'm'][hot_in].sum()

        #--- Get rhalf mass in dwarf galaxy for different components ---#

        ellipticity, pa = Analysis.calculate_ellipticity(x_kpc, y_kpc, mass=df[total_dw_star_mask].m.values, n_neighbors=30)

        Rotate = np.array([[np.cos(-pa), -np.sin(-pa)],
                   [np.sin(-pa),  np.cos(-pa)]])
        coords = np.vstack((x_kpc, y_kpc))  # shape: (2, N)
        rotated_coords = Rotate @ coords        # shape: (2, N)

        x_rot, y_rot = rotated_coords[0], rotated_coords[1]
        r_half = Analysis.half_light_radius(x_rot, y_rot, df[total_dw_star_mask].m, ep=ellipticity)

        dw_star_rhalf_mask = Analysis.get_radial_mask(x_kpc, y_kpc, r_half)
        dw_hot_gas_rhalf_mask = Analysis.get_radial_mask(hot_gas_x_kpc, hot_gas_y_kpc, r_half)
        dw_cold_gas_rhalf_mask = Analysis.get_radial_mask(cold_gas_x_kpc, cold_gas_y_kpc, r_half)

        M_star_rhalf = df.loc[total_dw_star_mask, :].loc[dw_star_rhalf_mask, 'm'].sum()
        M_hot_gas_rhalf = df.loc[dw_hot_gas_mask, :].loc[dw_hot_gas_rhalf_mask, 'm'].sum()
        M_cold_gas_rhalf = df.loc[dw_cold_gas_mask, :].loc[dw_cold_gas_rhalf_mask, 'm'].sum()

        vlos_array = Analysis.cal_vlos(df[total_dw_star_mask].xh, df[total_dw_star_mask].yh, df[total_dw_star_mask].zh,
                                   df[total_dw_star_mask].vxh, df[total_dw_star_mask].vyh, df[total_dw_star_mask].vzh)
        bin_centers, vlos_dispersion = Analysis.get_vlos_dispersion(x_kpc, y_kpc, vlos_array)

        rhalf_mask = bin_centers < r_half
        mean_vlos_dispersion = np.mean(vlos_dispersion[rhalf_mask])

        v_r, v_theta, v_phi = Analysis.cartesian_to_spherical(
            df.loc[total_dw_star_mask, 'x'][star_in].to_numpy(), df.loc[total_dw_star_mask, 'y'][star_in].to_numpy(), 
            df.loc[total_dw_star_mask, 'z'][star_in].to_numpy(), 
            df.loc[total_dw_star_mask, 'vx'][star_in].to_numpy(), df.loc[total_dw_star_mask, 'vy'][star_in].to_numpy(), 
            df.loc[total_dw_star_mask, 'vz'][star_in].to_numpy()
            )
        
        mean_vr, mean_vtheta, mean_vphi = np.mean(v_r), np.mean(v_theta), np.mean(v_phi)
        theoretical_sigma = np.sqrt((M_star_rhalf + M_hot_gas_rhalf + M_cold_gas_rhalf) / 930 / (r_half * 1000))

        masses = {
            'star': M_star,
            'star_half': M_star_rhalf,
            'hotgas': M_hot_gas,
            'hotgas_half': M_hot_gas_rhalf,
            'coldgas': M_cold_gas,
            'coldgas_half': M_cold_gas_rhalf,
            'mw_mass_r': mw_mass_r
        }
        
        return [masses['star'], masses['star_half'], 
                masses['hotgas'], masses['hotgas_half'], 
                masses['coldgas'], masses['coldgas_half'], 
                masses['mw_mass_r'],
                mean_vr, mean_vtheta, mean_vphi, d_mean_gal, r_half, 
                d_mean, tsnap, mean_vlos_dispersion, theoretical_sigma, 
                numsp]

    def process_and_stream(self, pool, numsp_list):
    
        columns = ['star_mass', 'star_half_mass', 
                   'hotgas_mass', 'hotgas_half_mass', 
                   'coldgas_mass', 'coldgas_half_mass', 
                   'mw_mass_r',
                   'vr', 'vtheta', 'vphi', 'distance_gal', 'rhalf', 
                   'distance', 'age', 'sigma', 'tsigma', 'numsp']

        results = []
        for result in pool.imap_unordered(self.get_info, numsp_list):
            results.append(result)
            gc.collect()

        results.sort(key=lambda x: x[-1])
    
        # Write sorted results to file
        with open(self.output_csv_file, 'a', newline='', buffering=16384) as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(columns)
            writer.writerows(results)

    def main(self):

        orbit_data = DataProcessor.read_orbit_parameters(orbit_file)
        if orbit_data is None:
            print("Failed to read orbit data.")
            sys.exit(1)

        fd_value = DataProcessor.read_feedback_value(gizmo_config)
        if fd_value is None:
            print("Failed to extract timescale value.")
            sys.exit(1)

        snapshot_files = [f for f in os.listdir(folder_path) if f.startswith('snapshot')]
        max_snapshot_num = len(snapshot_files) - 1

        parser = argparse.ArgumentParser(description='Process some integers.')
        default_range = f'0,{max_snapshot_num}'
        parser.add_argument('--range', type=str, default=default_range, 
                            help='the range of snapshots, e.g., "1,20"')
        parser.add_argument('--processes', type=int, default=1, 
                            help='number of processes to use')
        args = parser.parse_args()

        start, end = map(int, args.range.split(','))
        numsp_list = list(range(start, end + 1))

        with open(self.output_csv_file, 'w') as f:
            f.write(f"# Orbit Data: {','.join(map(str, orbit_data))}\n")
            f.write(f"# Timescale_fd: {fd_value}\n")

        with get_context("spawn").Pool(processes=args.processes) as pool:
            self.process_and_stream(pool, numsp_list)

        print(f"\nData saved to {self.output_csv_file}")

if __name__ == "__main__":
    processor = GetInfo()
    processor.main()