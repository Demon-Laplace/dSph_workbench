from basefunc import *
from variable import *
from pathlib import Path
import re
# from astropy.coordinates import CartesianRepresentation, CartesianDifferential, SkyCoord
# from scipy.interpolate import UnivariateSpline
# import astropy.units as u
import matplotlib.pyplot as plt
# from matplotlib import gridspec
# from astropy.coordinates import SkyCoord
import warnings
import gc

import argparse
from multiprocessing import get_context

warnings.filterwarnings("ignore")
plt.switch_backend('Agg')

modelname = Path().resolve().name
match = re.match(r'([A-Za-z]+)(\d+)', modelname)
if not match:
    raise ValueError(f"can't get '{modelname}' to dwarf name and number")

output_folder = './output_density_img'
output_folder_csv = 'output_density_csv'

dwarf_name, num_str = match.groups()
num = int(num_str)

modelfolder = './output'
elinfo_path = f'./elinfo_{modelname}.csv'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    print(f"Created output folder: {output_folder}")

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

try:
    dw_elinfo = pd.read_csv(elinfo_path, skiprows=2)
except:
    pass

ref_x = np.array([30, 50, 100, 200, 300, 450])
ref_y = np.array([1.0e-3, 6.0e-4, 2.0e-4, 6.0e-5, 2.5e-5, 1.0e-5])
ref2_y = np.array([5.0e-4, 3.0e-4, 2.0e-4, 1.2e-4, 5.0e-5, 3.0e-5])

def plot(numsp):

    simulation = GalaxySimulation(folder_path=folder_path, snapshot_num=numsp)
    df, tsnap = simulation.load_snapshot()
    df = simulation.center_on_MW()

    mw_star_mask, dw_star_mask = simulation.classify_mw_dwarf(df=df, r_exclude=5, dwarf_radius=core_radius, k_density=32)

    xc, yc = Analysis.find_center(df[dw_star_mask].x, df[dw_star_mask].y, units='kpc')
    yc, zc = Analysis.find_center(df[dw_star_mask].y, df[dw_star_mask].z, units='kpc')
    dw_r = np.sqrt(xc ** 2 + yc ** 2 + zc ** 2)

    hot_gas = df.loc[(df['tp'] == 0) & (df['temp'] > 20000), ['r', 'rho', 'm']]

    nbins = 600
    sample_rate = 0.05

    bins = np.linspace(hot_gas['r'].min(), 600, nbins+1)
    bins = np.sort(bins)
    bins = np.unique(bins)
    bins = bins[~np.isnan(bins)]

    if len(bins) < 2:
        print("Invalid bins:", bins)
        return

    hot_gas['r_bin'] = pd.cut(hot_gas['r'], bins=bins)

    rho_median = hot_gas.groupby('r_bin')['rho'].median()
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    idx = np.abs(bin_centers - dw_r).argmin()
    closest_bin = bin_centers[idx]
    closest_rho = rho_median[idx]

    plt.figure(figsize=(9,5))

    plt.plot(bin_centers, rho_median.values, c='red')

    plt.plot(ref_x, ref_y, c = 'blue', ls='--', lw=3, label='WJL MW1')
    plt.plot(ref_x, ref2_y, c = 'green', ls='--', lw=3, label='WJL MW2')

    mask = np.random.rand(len(hot_gas.r)) < sample_rate
    plt.scatter(hot_gas.r[mask], hot_gas.rho[mask], s=2, c='black', lw=0)
    plt.scatter(dw_r, closest_rho, s=60, c='skyblue', zorder=5, alpha=0.8)
    # plt.scatter(hot_gas.r, hot_gas.rho, s=1, c='gray', lw=0)

    plt.xlabel('Radius [kpc]')
    plt.ylabel('Median Gas Density')
    plt.title(f'{tsnap:.2f} Gyr')
    plt.grid(which='both', linestyle='--', linewidth=0.5, alpha=0.5, zorder=0, color='gray')
    plt.minorticks_on()
    plt.legend(loc='upper right')
    # plt.yaxis.set_ticks_position('both')
    # plt.xaxis.set_ticks_position('both')
    plt.tick_params(direction='in', which='both', top=True, right=True)
    plt.xlim(30, 500)
    plt.ylim(1e-6,1e-2)
    plt.xscale('log')
    plt.yscale('log')
    
    plt.tight_layout()
    output_file = os.path.join(output_folder, f'density_{numsp:03d}.png')
    os.makedirs(output_folder, exist_ok=True)
    
    plt.savefig(output_file, dpi=200)
    plt.close("all")

    del df
    del simulation
    del hot_gas
    
    gc.collect()

def main():

    snapshot_files = [f for f in os.listdir(folder_path) if f.startswith('snapshot')]
    max_snapshot_num = len(snapshot_files) - 1

    parser = argparse.ArgumentParser(description='Process some integers.')
    
    default_range = f'0,{max_snapshot_num}'
    
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
