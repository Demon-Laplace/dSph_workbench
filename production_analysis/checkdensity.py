from commonfunc import *
import matplotlib.pyplot as plt
import os
from multiprocessing import get_context
from variable import *
import warnings
from pathlib import Path
from matplotlib import gridspec
import argparse

warnings.filterwarnings("ignore")
plt.switch_backend('Agg')

output_folder = './output_density_img'
output_folder_csv = 'output_density_csv'
modelfolder = './output'
modelname = Path().resolve().name

try:
    elinfo_path = f'./elinfo_{modelname}.csv'   
except:
    print("elinfo file not found")

core_radius = fornax_core_radius

try:
    num = int(modelname.replace('Fornax', ''))
except:
    num = None

try:
    inifile = f'./Fornax.ini'
    dw_particles = cal_dwsum(inifile)
except:
    dw_particles = 0

sample_rate = 1

try:
    dw_elinfo = pd.read_csv(elinfo_path, skiprows=2)
except:
    pass

def plot(numsp, dw_particles=dw_particles):     

    df, tsnap = create_dataframe(numsp=numsp, path=modelfolder, dftype=1)

    if num is not None and num <= 55:

        x_arr = df['x'].to_numpy()
        y_arr = df['y'].to_numpy()
        z_arr = df['z'].to_numpy()

        vx_arr = df['vx'].to_numpy()
        vy_arr = df['vy'].to_numpy()
        vz_arr = df['vz'].to_numpy()

        df['x'] = z_arr
        df['y'] = x_arr
        df['z'] = y_arr

        df['vx'] = vz_arr
        df['vy'] = vx_arr
        df['vz'] = vy_arr

    idw = find_dw(df, dw_particles)
    df, iMW = modelbias(df, dw_particles)
    df = convertCoor(df)

    istar = df.tp == 2
    igas = df.tp == 0

    o_dw_star = (idw & istar)
    dw_gas = idw & igas

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
                else:
                    print("Warning: No mass column found")
                    mass = np.zeros(len(bin_centers))
            else:
                print("Warning: 'r' has only one unique value, cannot create bins")
                
    else:
        print("please confirmed you have removed the dwarf galaxy")

    # hgas = hgas.groupby('bin').apply(remove_outliers).reset_index(drop=True)
    # bin_means = hgas.groupby('bin')['rho'].median()
    # bin_centers = hgas.groupby('bin')['rh'].mean()

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
    plt.close("all")

    bin_df = pd.DataFrame({
        "rho": bin_means,
        "radius": bin_centers,
        "rho2": bin_means,
        "mass": mass
    })

    os.makedirs(output_folder_csv, exist_ok=True)
    bin_csv_file = os.path.join(output_folder_csv, f'density_profile_{numsp:03d}.csv')
    bin_df.to_csv(bin_csv_file, index=False)

    del df
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