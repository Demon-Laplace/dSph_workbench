from basefunc import *
from variable import *
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import sys
import argparse
sys.stdout.reconfigure(encoding='utf-8')

warnings.filterwarnings("ignore")

# ======================================================
# --- User Configuration ---
# ======================================================
modelname = Path().resolve().name
match = re.match(r'([A-Za-z]+)(\d+)', modelname)
if not match:
    raise ValueError(f"can't get '{modelname}' to dwarf name and number")

dwarf_name, num_str = match.groups()
num = int(num_str)
inifile = f'./{dwarf_name}.ini'

output_dir = Path("./density_output")
output_dir.mkdir(exist_ok=True)

max_snap = 400           # process snapshots 0–400
num_threads = 8          # number of parallel workers

# Observational comparison data
wlm_r = np.load("../radiusbin.npy")
wlm_f = np.load("../flux.npy")

# ======================================================
# --- Helper Function ---
# ======================================================
def compute_center(df, mask):
    """Return (x,y,z) center of mass for the given mask."""
    x = df.loc[mask, "x"].values
    y = df.loc[mask, "y"].values
    z = df.loc[mask, "z"].values
    m = df.loc[mask, "m"].values
    msum = m.sum()
    if msum == 0:
        return 0.0, 0.0, 0.0
    return np.sum(x*m)/msum, np.sum(y*m)/msum, np.sum(z*m)/msum


def process_snapshot(numsp):
    """Compute and save density profile figure for one snapshot."""
    try:
        folder_path = f"./output/"
        # inifile = folder_path + f"{dwarf_name}.ini"

        # --- check file existence ---
        # if not os.path.exists(inifile):
        #     print(f"[Warning] Snapshot {numsp:03d}: missing {inifile}")
        #     return

        dw_particles = GalaxySimulation.get_dw_num(inifile)
        simulation = GalaxySimulation(folder_path=folder_path, snapshot_num=numsp)
        df, tsnap = simulation.load_snapshot()

        mw_mask = simulation.find_MW_particles(
            dw_particles_num=dw_particles,
            cube_file=cube_path
        )
        dw_mask = simulation.find_dwarf_particles(dw_particles_num=dw_particles)
        df = simulation.convert_coordinates(df=df)

        # --- Define Masks ---
        mw_star_mask = mw_mask & ((df.tp == 2) | (df.tp == 3))
        dw_star_mask = dw_mask & ((df.tp == 2) | (df.tp == 3))
        dw_gas_mask = dw_mask & (df.tp == 0)

        dw_hot_gas_mask = dw_gas_mask & (df.temp >= 20000)
        dw_cold_gas_mask = dw_gas_mask & (df.temp < 20000)

        new_dw_star_mask, new_mw_star_mask = simulation.classify_new_dwarf_stars(dw_star_mask)
        total_dw_star_mask = dw_star_mask | new_dw_star_mask
        total_dw_gas_mask = dw_cold_gas_mask | dw_hot_gas_mask
        total_dw_mask = total_dw_star_mask | total_dw_gas_mask

        # --- Compute Centers of Mass ---
        x_star, y_star, z_star = compute_center(df, total_dw_star_mask)
        x_gas, y_gas, z_gas = compute_center(df, total_dw_gas_mask)
        x_all, y_all, z_all = compute_center(df, total_dw_mask)

        print(f"[Snapshot {numsp:03d}] Centers (kpc):")
        print(f"  Stars:    ({x_star:.3f}, {y_star:.3f}, {z_star:.3f})")
        print(f"  Gas:      ({x_gas:.3f}, {y_gas:.3f}, {z_gas:.3f})")
        print(f"  Star+Gas: ({x_all:.3f}, {y_all:.3f}, {z_all:.3f})")

        # --- Recenter positions & velocities ---
        df["x"] -= x_star
        df["y"] -= y_star
        df["z"] -= z_star
        df["vx"] -= df.loc[total_dw_mask, "vx"].mean()
        df["vy"] -= df.loc[total_dw_mask, "vy"].mean()
        df["vz"] -= df.loc[total_dw_mask, "vz"].mean()

        # ======================================================
        # --- Compute Surface Density ---
        # ======================================================
        x, y, z, m = df[dw_cold_gas_mask][["x", "y", "z", "m"]].values.T
        r3d = np.sqrt(x**2 + y**2 + z**2)
        mask = r3d <= 10.0
        x, y, z, m = x[mask], y[mask], z[mask], m[mask]
        r = np.sqrt(x**2 + y**2)

        xs, ys, zs, ms = df[dw_star_mask][["x", "y", "z", "m"]].values.T
        r3d_star = np.sqrt(xs**2 + ys**2 + zs**2)
        mask_star = r3d_star <= 10.0
        xs, ys, zs, ms = xs[mask_star], ys[mask_star], zs[mask_star], ms[mask_star]
        r_star = np.sqrt(xs**2 + ys**2)

        # Use same bins as gas

        r_max = 10.0
        nbin = 100
        bins = np.linspace(0, r_max, nbin + 1)
        r_cent = 0.5 * (bins[:-1] + bins[1:])
        mass_in_bin, _ = np.histogram(r, bins=bins, weights=m)
        area = np.pi * (bins[1:]**2 - bins[:-1]**2)
        Sigma = mass_in_bin / area

        Sigma_H = Sigma * 1.25e-6 

        mass_in_bin_star, _ = np.histogram(r_star, bins=bins, weights=ms)
        Sigma_star = mass_in_bin_star / area
        Sigma_H_star = Sigma_star * 1.25e-6  # Same conversion factor for comparison

        # ======================================================
        # --- Plot & Save Figure ---
        # ======================================================
        fig, axes = plt.subplots(1, 2, figsize=(9, 4), dpi=150)

        # Left: spatial map
        ax = axes[0]
        ax.scatter(df[dw_star_mask].x, df[dw_star_mask].y, s=0.1, c="black", label="Stars")
        ax.scatter(df[dw_cold_gas_mask].x, df[dw_cold_gas_mask].y, s=0.2, c="cyan", alpha=0.6, zorder=0, label="Cold gas")

        r_xy = np.sqrt(df[dw_star_mask].x**2 + df[dw_star_mask].y**2)
        mask_xy = r_xy < 1.0
        sigma_z = np.std(df[dw_star_mask][mask_xy].vz)

        # XZ plane (looking along y) -> use r_xz for selecting particles
        r_xz = np.sqrt(df[dw_star_mask].x**2 + df[dw_star_mask].z**2)
        mask_xz = r_xz < 1.0
        sigma_y = np.std(df[dw_star_mask][mask_xz].vy)

        # YZ plane (looking along x) -> use r_yz for selecting particles
        r_yz = np.sqrt(df[dw_star_mask].y**2 + df[dw_star_mask].z**2)
        mask_yz = r_yz < 1.0
        sigma_x = np.std(df[dw_star_mask][mask_yz].vx)

        # Add velocity dispersion text
        disp_text = f'r < 1 kpc:\n'
        disp_text += f'σx (yz-proj) = {sigma_x:.1f} km/s\n'
        disp_text += f'σy (xz-proj) = {sigma_y:.1f} km/s\n'
        disp_text += f'σz (xy-proj) = {sigma_z:.1f} km/s'
        
        ax.text(0.05, 0.95, disp_text,
                transform=ax.transAxes,
                verticalalignment='top',
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'),
                fontsize=8)

        # Rest of the plotting code remains the same
        ax.set_xlim(-20, 20)
        ax.set_ylim(-20, 20)


        ax.set_xlabel("x [kpc]")
        ax.set_ylabel("y [kpc]")
        ax.set_title(f"T = {tsnap:.2f} Gyr")
        ax.grid(True, ls="--", alpha=0.4)
        ax.legend(markerscale=5, fontsize=7)
        ax.tick_params(axis='both', which='both', direction='in')

        # Right: Σ(R)
        ax = axes[1]
        ax.plot(r_cent, Sigma_H, lw=2., c="cyan", label=f"Gas")
        ax.plot(r_cent, Sigma_H_star, lw=2., c="black", label="Stars")
        ax.plot(wlm_r, wlm_f, label="WLM (Neel et al.)", color="orange", lw=2, ls="--")
        ax.set_yscale("log")
        ax.set_xlabel(r"$R\ [{\rm kpc}]$")
        ax.set_ylabel(r"$\Sigma(R)\ [10^{20}\ {\rm H\ atoms\ cm^{-2}}]$")
        ax.set_xlim(0, 10)
        ax.set_ylim(0.01, 55)
        ax.set_title("Radial Surface Density Profile")
        ax.grid(True, ls="--", alpha=0.4)
        ax.legend()
        ax.tick_params(axis='both', which='both', direction='in')

        fig.suptitle(f"Denstiy Profile {modelname}")

        plt.tight_layout()

        outname = output_dir / f"snapshot_{numsp:03d}_density_profile.png"
        plt.savefig(outname, dpi=200)
        plt.close(fig)

        print(f"[Snapshot {numsp:03d}] Saved → {outname}")

    except Exception as e:
        print(f"[Error] Snapshot {numsp:03d}: {e}")


# ======================================================
# --- Parallel Execution ---
# ======================================================
def get_max_snapshot(folder_path="./output"):
    """Find the highest snapshot number in output directory"""
    try:
        snapshot_files = [f for f in os.listdir(folder_path) 
                         if f.startswith('snapshot') and f.endswith('.hdf5')]
        if not snapshot_files:
            return 400  # fallback default
            
        numbers = []
        for f in snapshot_files:
            match = re.search(r'snapshot_(\d+)\.hdf5', f)
            if match:
                numbers.append(int(match.group(1)))
                
        return max(numbers) if numbers else 400
    except Exception as e:
        print(f"Warning: Error detecting max snapshot: {e}")
        return 400

def main():
    parser = argparse.ArgumentParser(description='Process Fornax simulation snapshots.')
    
    # Get max snapshot number
    max_snap = get_max_snapshot()
    default_range = f"0,{max_snap}"
    
    # Add arguments
    parser.add_argument('--range', type=str, default=default_range,
                      help='Snapshot range to process (e.g. "0,400")')
    parser.add_argument('--processes', type=int, default=8,
                      help='Number of parallel processes')
    
    args = parser.parse_args()
    
    # Parse snapshot range
    try:
        start, end = map(int, args.range.split(','))
    except ValueError:
        print(f"Invalid range format. Using default: {default_range}")
        start, end = 0, max_snap
        
    print(f"\n=== Processing snapshots {start}–{end} with {args.processes} processes ===\n")
    
    # Process snapshots in parallel
    snapshots = list(range(start, end + 1))
    with ProcessPoolExecutor(max_workers=args.processes) as executor:
        futures = {executor.submit(process_snapshot, numsp): numsp 
                  for numsp in snapshots}
        for future in as_completed(futures):
            num = futures[future]
            try:
                future.result()
            except Exception as exc:
                print(f"[Error] Snapshot {num:03d}: {exc}")

    print("\n=== All snapshots processed ===\n")

if __name__ == "__main__":
    main()

    