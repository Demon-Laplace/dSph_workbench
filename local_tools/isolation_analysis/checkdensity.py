from basefunc import *
from variable import *
from pathlib import Path
from model_paths import model_id
import re
from astropy.coordinates import CartesianRepresentation, CartesianDifferential, SkyCoord
# from scipy.interpolate import UnivariateSpline
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib import gridspec
from astropy.coordinates import SkyCoord
import warnings
import gc  # 内存回收

import argparse
from multiprocessing import get_context

warnings.filterwarnings("ignore")
plt.switch_backend('Agg')

def remove_outliers(group):
    q1 = group['rho'].quantile(0.30)
    q3 = group['rho'].quantile(0.70)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    return group[(group['rho'] >= lower_bound) & (group['rho'] <= upper_bound)]

modelname = model_id()
match = re.match(r'([A-Za-z]+)(\d+)', modelname)
if not match:
    raise ValueError(f"can't get '{modelname}' to dwarf name and number")

output_folder = './output_density_img'
output_folder_csv = 'output_density_csv'

dwarf_name, num_str = match.groups()
num = int(num_str)

modelfolder = './output'
elinfo_path = f'./elinfo_{modelname}.csv'

# core_radius = core_radius_dict.get(dwarf_name, None)
# if core_radius is None:
#     raise ValueError(f"please input '{dwarf_name}''s radius")

# dw_elinfo = pd.read_csv(elinfo_path, skiprows=2)

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    print(f"Created output folder: {output_folder}")

# ini_x, ini_y, ini_z, ini_vx, ini_vy, ini_vz, fd_value = DataProcessor.get_ini_data(elinfo_path=elinfo_path)
# fd_factor_raw = fd_value / 8.3
# if fd_factor_raw >= 0.90:
#     fd_factor_str = f"{round(fd_factor_raw):.0f}"
# else:
#     fd_factor_str = f"{fd_factor_raw:.1f}"

core_radius = fornax_core_radius

try:
    num = int(modelname.replace('Fornax', ''))
except:
    num = None

try:
    inifile = f'./{dwarf_name}.ini'
    dw_particles = GalaxySimulation.get_dw_num(inifile)
except:
    dw_particles = 0

sample_rate = 0.7

try:
    dw_elinfo = pd.read_csv(elinfo_path, skiprows=2)
except:
    pass

def plot(numsp, dw_particles=dw_particles):
    try:
        simulation = GalaxySimulation(folder_path=folder_path, snapshot_num=numsp)
        df, tsnap = simulation.load_snapshot()
        mw_mask = simulation.find_MW_particles(dw_particles_num=dw_particles, cube_file=cube_path)
        dw_mask = simulation.find_dwarf_particles(dw_particles_num=dw_particles)
        df = simulation.center_on_MW(mw_mask, iterations=2, radius=100) 
        df = simulation.convert_coordinates(df=df)

        istar = df.tp == 2
        igas = df.tp == 0

        o_dw_star = (dw_mask & istar)
        dw_gas = dw_mask & igas

        # hgas = df[df.temp >= 20000]
        num_bins = 1000
        hot_gas_mask = df.temp >= 20000

        if dw_particles == 0:

            if 'r' in df.columns:
                df = df.dropna(subset=['r'])
                
                if df[hot_gas_mask]['r'].nunique() > 1:
                    # Create bins for all particles
                    df.loc[hot_gas_mask, 'bin'] = pd.cut(df.loc[hot_gas_mask, 'r'], num_bins, duplicates='drop')
                    
                    # Calculate density profile
                    df_hot = df[hot_gas_mask].groupby('bin').apply(remove_outliers).reset_index(drop=True)
                    bin_means = df_hot.groupby('bin')['rho'].median()
                    bin_centers = df_hot.groupby('bin')['r'].mean()
                    
                    # Calculate cumulative mass within each radius
                    if 'm' in df.columns:
                        df_sorted = df.sort_values('r')  # Sort by radius
                        cumulative_mass = df_sorted['m'].cumsum()  # Cumulative sum of masses
                        # Interpolate cumulative mass at bin centers
                        mass = np.interp(bin_centers, df_sorted['r'], cumulative_mass)
                        del df_sorted, cumulative_mass  # 及时释放临时 DataFrame
                    else:
                        print("Warning: No mass column found")
                        mass = np.zeros(len(bin_centers))
                else:
                    print("Warning: 'r' has only one unique value, cannot create bins")
                    
        else:
            # Handle case where dw_particles != 0: build hot-gas table from df
            hgas = df[hot_gas_mask].copy()
            # Ensure we have radius column and enough points to bin
            if 'r' in hgas.columns and hgas['r'].nunique() > 1:
                # assign bins on the hot-gas subset
                hgas.loc[:, 'bin'] = pd.cut(hgas['r'], num_bins, duplicates='drop')
                df_hot = hgas.groupby('bin').apply(remove_outliers).reset_index(drop=True)
                bin_means = df_hot.groupby('bin')['rho'].median()
                bin_centers = df_hot.groupby('bin')['r'].mean()

                # Calculate cumulative mass within each radius using full df if mass exists
                if 'm' in df.columns:
                    df_sorted = df.sort_values('r')  # Sort by radius
                    cumulative_mass = df_sorted['m'].cumsum()  # Cumulative sum of masses
                    mass = np.interp(bin_centers, df_sorted['r'], cumulative_mass)
                    del df_sorted, cumulative_mass
                else:
                    mass = np.zeros(len(bin_centers))
            else:
                print("Warning: insufficient hot-gas data for binning when dwarf particles > 0; skipping snapshot")
                return

        bin_centers = bin_centers.to_numpy()
        bin_means = bin_means.to_numpy()

        indices = np.random.choice(len(df_hot.r), int(len(df_hot.r) * sample_rate), replace=False)

        Tage = round(tsnap, 2)

        # if 'rho_hgas' in dw_elinfo.columns:
        #     idx = np.abs(dw_elinfo['age'] - Tage).idxmin()
        #     red_x = dw_elinfo['distance'].iloc[idx]
        #     red_y = dw_elinfo['rho_hgas'].iloc[idx]

        fig = plt.figure(figsize=(9, 8)) 

        gs = gridspec.GridSpec(2, 3, height_ratios=[2, 1], width_ratios=[1, 1, 1]) 

        ax1 = fig.add_subplot(gs[0, :])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[1, 1])
        ax4 = fig.add_subplot(gs[1, 2])

        try:

            if 'rho_hgas' in dw_elinfo.columns:
                traj_mask = dw_elinfo.index <= numsp
                traj_x = dw_elinfo.loc[traj_mask, 'distance']
                traj_y = dw_elinfo.loc[traj_mask, 'rho_hgas']
                # Plot the trajectory line
                ax1.plot(traj_x.to_numpy(), traj_y.to_numpy(), color='skyblue', lw=3, alpha=0.8, label=f'Fornax{num}', zorder=20)
                # Highlight the current snapshot as a larger dot
                ax1.scatter(traj_x.iloc[-1], traj_y.iloc[-1], color='skyblue', s=18, zorder=10)
                del traj_mask, traj_x, traj_y  # 释放临时变量
        except:
            pass

        ax1.scatter(df_hot.r.iloc[indices], df_hot.rho.iloc[indices], s=0.3, alpha=0.6, linewidth=0, c='black')
        ax1.plot(bin_centers, bin_means, c='red', lw=2)
        # ax1.scatter(red_x, red_y, color='skyblue', s=30, zorder=10, label=f'Fornax{num}')

        ref_x = np.array([30, 50, 100, 200, 300, 450])
        ref_y = np.array([1.0e-3, 6.0e-4, 2.0e-4, 6.0e-5, 2.5e-5, 1.0e-5])

        # Plot the reference line
        ax1.plot(ref_x, ref_y, ls='--', color='blue', label='WJL MW1', zorder=5)
        ax1.set_xlabel('R (kpc)')
        ax1.set_ylabel('rho')
        ax1.set_xlim(40,500)
        ax1.set_ylim(1e-6,1e-2)
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        ax1.yaxis.set_ticks_position('both')
        ax1.xaxis.set_ticks_position('both')
        ax1.grid(which='both', linestyle='--', linewidth=0.5, alpha=0.5, zorder=0, color='gray')
        ax1.tick_params(direction='in', which='both', top=True, right=True)
        ax1.set_title(f'{Tage:.2f} Gyr')
        ax1.legend()

        ax2.scatter(df[o_dw_star].x, df[o_dw_star].y, s=4, linewidth=0, c='black', alpha=0.6)
        ax2.scatter(df[dw_gas].x, df[dw_gas].y, s=0.2, lw=0, alpha=0.2, c='orange')
        ax2.scatter(0, 0, s=60, marker='x', c='red')
        ax2.set_xlabel('X (kpc)')
        ax2.set_ylabel('Y (kpc)')
        ax2.set_title('X-Y', fontsize=10)
        ax2.tick_params(direction='in')
        ax2.set_xlim(-300,300)
        ax2.set_ylim(-300,300)

        ax3.scatter(df[o_dw_star].x, df[o_dw_star].z, c='black', s=4, linewidth=0, alpha=0.6)
        ax3.scatter(df[dw_gas].x, df[dw_gas].z, s=0.2, lw=0, alpha=0.2, c='orange')
        ax3.scatter(0, 0, s=60, marker='x', c='red')
        ax3.set_xlabel('X (kpc)')
        ax3.set_ylabel('Z (kpc)')
        ax3.set_title('X-Z', fontsize=10)
        ax3.tick_params(direction='in')
        ax3.set_xlim(-300,300)
        ax3.set_ylim(-300,300)

        ax4.scatter(df[o_dw_star].y, df[o_dw_star].z, c='black', s=4, linewidth=0, alpha=0.6)
        ax4.scatter(df[dw_gas].y, df[dw_gas].z, s=0.2, lw=0, alpha=0.2, c='orange')
        ax4.scatter(0, 0, s=60, marker='x', c='red')
        ax4.set_xlabel('Y (kpc)')
        ax4.set_ylabel('Z (kpc)')
        ax4.set_title('Y-Z', fontsize=10)
        ax4.tick_params(direction='in')
        ax4.set_xlim(-300,300)
        ax4.set_ylim(-300,300)

        plt.tight_layout()

        output_file = os.path.join(output_folder, f'density_{numsp:03d}.png')
        os.makedirs(output_folder, exist_ok=True)
        
        plt.savefig(output_file, dpi=200)
        plt.close('all')  # 关闭所有图表

        bin_df = pd.DataFrame({
            "rho": bin_means,
            "radius": bin_centers,
            # "rho2": bin_means,
            "mass": mass
        })

        os.makedirs(output_folder_csv, exist_ok=True)
        bin_csv_file = os.path.join(output_folder_csv, f'density_profile_{numsp:03d}.csv')
        bin_df.to_csv(bin_csv_file, index=False)
        
    finally:
        # 无条件执行内存清理，即使出现异常也会运行
        try:
            plt.close('all')
        except:
            pass
        
        # 删除大对象引用
        for var_name in ['df', 'df_hot', 'df_sorted', 'fig', 'ax1', 'ax2', 'ax3', 'ax4', 
                         'indices', 'bin_centers', 'bin_means', 'bin_df', 'simulation']:
            if var_name in locals():
                try:
                    del locals()[var_name]
                except:
                    pass
        
        # 强制垃圾回收
        gc.collect()


def main():

    snapshot_files = [f for f in os.listdir(folder_path) if f.startswith('snapshot')]
    max_snapshot_num = len(snapshot_files) - 1

    parser = argparse.ArgumentParser(description='Process some integers.')
    
    default_range = f'1,{max_snapshot_num}'
    
    parser.add_argument('--range', type=str, default=default_range, help='the range of snapshots, e.g., "1,20"')
    parser.add_argument('--processes', type=int, default=1, help='number of processes to use')
    # parser.add_argument('--plot', action='store_true', help='generate plots for snapshots')
    args = parser.parse_args()

    start, end = map(int, args.range.split(','))
    processes = args.processes
    # plotfig = args.plot

    numsp_list = range(start, end + 1)

    with get_context("spawn").Pool(processes=processes) as pool:
        pool.starmap(plot, [(numsp, ) for numsp in numsp_list])

if __name__ == "__main__":
    main()