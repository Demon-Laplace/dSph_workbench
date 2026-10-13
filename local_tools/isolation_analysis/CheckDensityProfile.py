from basefunc import *
from variable import *
from pathlib import Path
from model_paths import model_id
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
modelname = model_id()
match = re.match(r'([A-Za-z]+)(\d+)', modelname)
if not match:
    raise ValueError(f"can't get '{modelname}' to dwarf name and number")

dwarf_name, num_str = match.groups()
num = int(num_str)
inifile = f'./{dwarf_name}.ini'

output_dir = Path("./density_output")
output_dir.mkdir(exist_ok=True)

# Full evolution CSV file (merged snapshot data)
EVO_FULL = output_dir / f"IC_{modelname}_evolution.csv"

# max_snap = 400           # process snapshots 0–400
# num_threads = 8          # number of parallel workers

sfr_eff = DataProcessor.read_sfreff_value(gizmo_config)
fd_value = DataProcessor.read_feedback_value(gizmo_config)

# Observational comparison data
try: 
    wlm_r = np.load("../radiusbin.npy")
    wlm_f = np.load("../flux.npy")
except:
    print("no wlm documents")
    wlm_r, wlm_f = 0, 0


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
    """Compute and save density profile figure for one snapshot. Returns data row dict."""
    try:
        folder_path = f"./output/"
        # inifile = folder_path + f"{dwarf_name}.ini"

        # --- check file existence ---
        # if not os.path.exists(inifile):
        #     print(f"[Warning] Snapshot {numsp:03d}: missing {inifile}")
        #     return

        dw_particles = GalaxySimulation.get_dw_num(inifile)
        simulation = GalaxySimulation(folder_path=folder_path, snapshot_num=numsp)
        df, tsnap = simulation.load_snapshot(box=False)

        if dw_particles != len(df):
            
            mw_mask = simulation.find_MW_particles(
                dw_particles_num=dw_particles,
                cube_file=cube_path
            )

        dw_mask = simulation.find_dwarf_particles(dw_particles_num=dw_particles)
        df = simulation.convert_coordinates(df=df)

        # --- Define Masks ---
        # mw_star_mask = mw_mask & ((df.tp == 2) | (df.tp == 3))
        # dw_star_mask = dw_mask & ((df.tp == 2) | (df.tp == 4))
        # dw_gas_mask = dw_mask & (df.tp == 0)

        dw_star_mask = ((df.tp == 2) | (df.tp == 4))
        dw_gas_mask = (df.tp == 0)

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
        xcg, ycg, zcg, mcg = df[dw_cold_gas_mask][["x", "y", "z", "m"]].values.T
        r3d = np.sqrt(xcg**2 + ycg**2 + zcg**2)
        mask = r3d <= 10.0
        xcg, ycg, zcg, mcg = xcg[mask], ycg[mask], zcg[mask], mcg[mask]
        r_cg = np.sqrt(xcg**2 + ycg**2)

        xhg, yhg, zhg, mhg = df[dw_hot_gas_mask][["x", "y", "z", "m"]].values.T
        r3d = np.sqrt(xhg**2 + yhg**2 + zhg**2)
        mask = r3d <= 10.0
        xhg, yhg, zhg, mhg = xhg[mask], yhg[mask], zhg[mask], mhg[mask]
        r_hg = np.sqrt(xhg**2 + yhg**2)

        xs, ys, zs, ms = df[total_dw_star_mask][["x", "y", "z", "m"]].values.T
        r3d_star = np.sqrt(xs**2 + ys**2 + zs**2)
        mask_star = r3d_star <= 10.0
        xs, ys, zs, ms = xs[mask_star], ys[mask_star], zs[mask_star], ms[mask_star]
        r_star = np.sqrt(xs**2 + ys**2)

        ellipticity, pa = Analysis.calculate_ellipticity(
            xs, zs, mass=ms, n_neighbors=30
        )

        Rotate = np.array([[np.cos(-pa), -np.sin(-pa)],
                        [np.sin(-pa),  np.cos(-pa)]])

        coords = np.vstack((xs, zs))  # shape: (2, N)
        rotated_coords = Rotate @ coords        # shape: (2, N)

        x_rot, z_rot = rotated_coords[0], rotated_coords[1]

        r_half = Analysis.half_light_radius(x_rot, z_rot, df[total_dw_star_mask].m, ep=ellipticity)

        star_points = np.column_stack([xs, ys])
        boundary_path = Analysis.estimate_galaxy_boundary(xs, ys)
        star_in = boundary_path.contains_points(star_points)

        subset_star = df.loc[total_dw_star_mask].copy()
        r3d_star = np.sqrt(subset_star["x"]**2 + subset_star["y"]**2 + subset_star["z"]**2)
        mask_star = r3d_star <= 5.0
        subset_star = subset_star.loc[mask_star]
        star_points = np.column_stack([subset_star["x"], subset_star["y"]])
        star_in = boundary_path.contains_points(star_points)
        M_star = subset_star.loc[star_in, "m"].sum()

        # For gas, only use r<=5kpc, do not use boundary
        subset_cold = df.loc[dw_cold_gas_mask].copy()
        r3d_cold = np.sqrt(subset_cold["x"]**2 + subset_cold["y"]**2 + subset_cold["z"]**2)
        mask_cold = r3d_cold <= 5.0
        subset_cold = subset_cold.loc[mask_cold]
        M_cold_gas = subset_cold["m"].sum()

        subset_hot = df.loc[dw_hot_gas_mask].copy()
        r3d_hot = np.sqrt(subset_hot["x"]**2 + subset_hot["y"]**2 + subset_hot["z"]**2)
        mask_hot = r3d_hot <= 5.0
        subset_hot = subset_hot.loc[mask_hot]
        M_hot_gas = subset_hot["m"].sum()

        # 星
        subset_star_proj = subset_star.copy()
        r_proj_star = np.sqrt(subset_star_proj["x"]**2 + subset_star_proj["z"]**2)
        mask_star_proj = r_proj_star <= r_half
        M_star_rhalf = subset_star_proj.loc[mask_star_proj, "m"].sum()

        # 冷气体
        subset_cold_proj = subset_cold.copy()
        r_proj_cold = np.sqrt(subset_cold_proj["x"]**2 + subset_cold_proj["z"]**2)
        mask_cold_proj = r_proj_cold <= r_half
        M_cold_rhalf = subset_cold_proj.loc[mask_cold_proj, "m"].sum()

        # 热气体
        subset_hot_proj = subset_hot.copy()
        r_proj_hot = np.sqrt(subset_hot_proj["x"]**2 + subset_hot_proj["z"]**2)
        mask_hot_proj = r_proj_hot <= r_half
        M_hot_rhalf = subset_hot_proj.loc[mask_hot_proj, "m"].sum()


        # Use same bins as gas

        r_max = 10.0
        nbin = 100
        bins = np.linspace(0, r_max, nbin + 1)
        r_cent = 0.5 * (bins[:-1] + bins[1:])
        mass_in_bin, _ = np.histogram(r_cg, bins=bins, weights=mcg)
        area = np.pi * (bins[1:]**2 - bins[:-1]**2)
        Sigma = mass_in_bin / area

        Sigma_H = Sigma * 1.25e-6 

        mass_in_bin_star, _ = np.histogram(r_star, bins=bins, weights=ms)
        Sigma_star = mass_in_bin_star / area
        Sigma_H_star = Sigma_star * 1.25e-6  # Same conversion factor for comparison

        # ======================================================
        # --- Plot & Save Figure ---
        # ======================================================
        fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=150)

        # Top Left: spatial map (xy plane)
        ax = axes[0,0]
        ax.scatter(df[total_dw_star_mask].x, df[total_dw_star_mask].y, s=0.1, c="black", label="Stars")
        ax.scatter(df[dw_cold_gas_mask].x, df[dw_cold_gas_mask].y, s=0.2, c="cyan", alpha=0.6, zorder=0, label="Cold gas")

        r_xy = np.sqrt(df[total_dw_star_mask].x**2 + df[total_dw_star_mask].y**2)
        mask_xy = r_xy < 1.0
        sigma_z = np.std(df[total_dw_star_mask][mask_xy].vz)

        # XZ plane (looking along y) -> use r_xz for selecting particles
        r_xz = np.sqrt(df[total_dw_star_mask].x**2 + df[total_dw_star_mask].z**2)
        mask_xz = r_xz < 1.0
        sigma_y = np.std(df[total_dw_star_mask][mask_xz].vy)

        # YZ plane (looking along x) -> use r_yz for selecting particles
        r_yz = np.sqrt(df[total_dw_star_mask].y**2 + df[total_dw_star_mask].z**2)
        mask_yz = r_yz < 1.0
        sigma_x = np.std(df[total_dw_star_mask][mask_yz].vx)

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

        ax.set_xlim(-20, 20)
        ax.set_ylim(-20, 20)
        ax.set_xlabel("x [kpc]")
        ax.set_ylabel("y [kpc]")
        ax.set_title(f"T = {tsnap:.2f} Gyr")
        ax.grid(True, ls="--", alpha=0.4)
        ax.legend(markerscale=5, fontsize=7, loc='lower right')
        ax.tick_params(axis='both', which='both', direction='in')

        # Top Right: Σ(R)
        ax = axes[0,1]
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

        # Bottom Left: xz projection
        ax = axes[1,0]
        ax.scatter(df[total_dw_star_mask].x, df[total_dw_star_mask].z, s=0.1, c="black", label="Stars")
        ax.scatter(df[dw_cold_gas_mask].x, df[dw_cold_gas_mask].z, s=0.2, c="cyan", alpha=0.6, zorder=0, label="Cold gas")
        
        # Calculate stellar mass within 5kpc
        r_5kpc = np.sqrt(df[dw_star_mask].x**2 + df[dw_star_mask].y**2 + df[dw_star_mask].z**2)
        stellar_mass_5kpc = df[dw_star_mask][r_5kpc < 5.0].m.sum()
        
        # 新增total_gas和粒子数信息
        total_gas = df.loc[df.tp == 0, "m"].sum()
        n_gas = (df.tp == 0).sum()
        n_star = (df.tp >= 2).sum()
        ax.text(0.05, 0.95,
            f'Stellar Mass (<5kpc):{stellar_mass_5kpc:.2e} M⊙\n'
            f'Stellar Half-light Radius: {r_half:.2f} kpc\n'
            f'Stellar Mass within R_half: {M_star_rhalf:.2e} M⊙\n'
            f'Cold Gas Mass within R_half: {M_cold_rhalf:.2e} M⊙\n'
            f'Total Gas Mass (tp==0): {total_gas:.2e} M⊙\n'
            f'Gas Particles (tp==0): {n_gas}\n'
            f'Star Particles (tp>=2): {n_star}\n'
            f'Gas fraction: {M_cold_gas / (M_cold_gas + M_hot_gas + stellar_mass_5kpc):.2f}\n'
            f'Gas fraction within R_half: {M_cold_rhalf / (M_cold_rhalf + M_hot_rhalf + M_star_rhalf):.2f}',
            transform=ax.transAxes,
            verticalalignment='top',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'),
            fontsize=8)
        
        ax.set_xlim(-20, 20)
        ax.set_ylim(-20, 20)
        ax.set_xlabel("x [kpc]")
        ax.set_ylabel("z [kpc]")
        ax.grid(True, ls="--", alpha=0.4)
        ax.legend(markerscale=5, fontsize=7, loc='lower right')
        ax.tick_params(axis='both', which='both', direction='in')

        # Bottom Right: Density vs Time
        ax = axes[1,1]
        
        # Compute average surface density within 1 kpc (projected)
        area_inner = np.pi * (1.0**2)  # kpc^2

        # Gas: use projected radius r
        mask_inner_gas = r_cg <= 1.0
        mass_inner_gas = mcg[mask_inner_gas].sum() if mask_inner_gas.size else 0.0
        Sigma_inner = mass_inner_gas / area_inner
        Sigma_H_inner = Sigma_inner * 1.25e-6

        # Stars: use projected radius r_star
        mask_inner_star = r_star <= 1.0
        mass_inner_star = ms[mask_inner_star].sum() if mask_inner_star.size else 0.0
        Sigma_inner_star = mass_inner_star / area_inner
        Sigma_H_inner_star = Sigma_inner_star * 1.25e-6

        # Save current inner densities and time to separate file for this snapshot
        density_file = output_dir / f"density_snapshot_{numsp:03d}.txt"
        np.savetxt(density_file,
                  [[tsnap, Sigma_H_inner, Sigma_H_inner_star]],
                  header="time gas_density stellar_density",
                  fmt='%.6f %.6e %.6e')

        # Bottom Right: Density vs Time (read all existing density files and plot)
        ax = axes[1,1]
        
        try:
            # Collect all density snapshot files created so far
            density_files = sorted(list(output_dir.glob("density_snapshot_*.txt")))
            if density_files:
                all_data = []
                for f in density_files:
                    try:
                        data = np.loadtxt(f, skiprows=1)
                        if data.size > 0:
                            all_data.append(data)
                    except:
                        continue
                
                if all_data:
                    # Combine and sort by time
                    all_data = [np.atleast_1d(d) for d in all_data]
                    all_data = np.vstack(all_data)
                    sort_idx = np.argsort(all_data[:, 0])
                    times_plot = all_data[sort_idx, 0]
                    gas_dens_plot = all_data[sort_idx, 1]
                    star_dens_plot = all_data[sort_idx, 2]
                    
                    # Handle zeros for log scale
                    pos_vals = np.concatenate((gas_dens_plot[gas_dens_plot > 0], 
                                               star_dens_plot[star_dens_plot > 0]))
                    use_log = pos_vals.size > 0
                    
                    if use_log:
                        minpos = pos_vals.min()
                        eps = minpos * 1e-3 if minpos > 0 else 1e-12
                        gas_plot = np.where(gas_dens_plot <= 0, eps, gas_dens_plot)
                        star_plot = np.where(star_dens_plot <= 0, eps, star_dens_plot)
                        ax.plot(times_plot, gas_plot, lw=2.0, color='cyan', label='Gas')
                        ax.plot(times_plot, star_plot, lw=2.0, color='black', label='Stars')
                        ax.set_yscale('log')
                    else:
                        ax.plot(times_plot, gas_dens_plot, lw=2.0, color='cyan', label='Gas')
                        ax.plot(times_plot, star_dens_plot, lw=2.0, color='black', label='Stars')
                    
                    ax.legend(fontsize=8, loc='lower right')
        except Exception as e:
            # If plotting fails, just show empty plot with message
            ax.text(0.5, 0.5, f'Density plot error:\n{str(e)[:40]}', 
                    ha='center', va='center', transform=ax.transAxes, fontsize=8)
        
        ax.set_xlabel("Time [Gyr]")
        ax.set_ylim(0.01, 55)
        ax.set_xlim(0, 8)
        ax.set_ylabel(r"$\Sigma(R)\ [10^{20}\ {\rm H\ atoms\ cm^{-2}}]$")
        ax.grid(True, ls="--", alpha=0.4)
        ax.tick_params(axis='both', which='both', direction='in')

        fig.suptitle(f"IC_{modelname}\n"
                     f"feedback = {fd_value / 8.3:.2f} median, sfr_eff = {sfr_eff:.3f}")

        plt.tight_layout()

        outname = output_dir / f"snapshot_{numsp:03d}_density_profile.png"
        plt.savefig(outname, dpi=200)
        plt.close('all')

        print(f"[Snapshot {numsp:03d}] Saved → {outname}")
        
        # Save analysis data to individual JSON file (parallel-safe, one per snapshot)
        import json
        try:
            # 1. total_gas为tp==0所有粒子的总质量
            total_gas = df.loc[df.tp == 0, "m"].sum()
            # 2. gas_fraction算法保持不变
            gas_fraction = M_cold_gas / (M_cold_gas + M_hot_gas + stellar_mass_5kpc) if (M_cold_gas + M_hot_gas + stellar_mass_5kpc) > 0 else 0
            gas_fraction_rhalf = M_cold_rhalf / (M_cold_rhalf + M_hot_rhalf + M_star_rhalf) if (M_cold_rhalf + M_hot_rhalf + M_star_rhalf) > 0 else 0

            # 3. 输出tp==0和tp>=2粒子的数量
            n_gas = (df.tp == 0).sum()
            n_star = (df.tp >= 2).sum()
            print(f"[Snapshot {numsp:03d}] Number of gas particles (tp==0): {n_gas}")
            print(f"[Snapshot {numsp:03d}] Number of star particles (tp>=2): {n_star}")

            data_row = {
                'snapshot': int(numsp),
                'time_gyr': float(tsnap),
                'sigma_x': float(sigma_x),
                'sigma_y': float(sigma_y),
                'sigma_z': float(sigma_z),
                'stellar_mass_5kpc': float(stellar_mass_5kpc),
                'half_light_radius': float(r_half),
                'stellar_mass_rhalf': float(M_star_rhalf),
                'cold_gas_mass_rhalf': float(M_cold_rhalf),
                'hot_gas_mass_rhalf': float(M_hot_rhalf),
                'cold_gas_mass': float(M_cold_gas),
                'hot_gas_mass': float(M_hot_gas),
                'total_gas_mass': float(total_gas),
                'gas_fraction': float(gas_fraction),
                'gas_fraction_rhalf': float(gas_fraction_rhalf),
                'n_gas_particles': int(n_gas),
                'n_star_particles': int(n_star),
            }

            print(f"[Snapshot {numsp:03d}] Collected data row")
            return data_row
        except Exception as e:
            print(f"[Warning] Failed to collect data for snapshot {numsp:03d}: {e}")
            return None

    except Exception as e:
        print(f"[Error] Snapshot {numsp:03d}: {e}")
        return None


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

def create_density_evolution_csv():
    """Create separate density evolution CSV from density_snapshot_*.txt files."""
    import numpy as np
    
    density_files = sorted(list(output_dir.glob("density_snapshot_*.txt")))
    if not density_files:
        print("No density data files found!")
        return
    
    try:
        all_data = []
        for f in density_files:
            try:
                data = np.loadtxt(f, skiprows=1)
                if data.size > 0:
                    all_data.append(data)
            except:
                continue
        
        if all_data:
            # Combine and sort by time
            all_data = [np.atleast_1d(d) for d in all_data]
            all_data = np.vstack(all_data)
            sort_idx = np.argsort(all_data[:, 0])
            times = all_data[sort_idx, 0]
            gas_densities = all_data[sort_idx, 1]
            stellar_densities = all_data[sort_idx, 2]
            
            # Save to CSV
            evo_out = output_dir / "density_evolution.csv"
            data_to_save = np.column_stack((times, gas_densities, stellar_densities))
            with evo_out.open('w') as f:
                f.write(f"#Timescale_fd: {fd_value}\n")
                f.write(f"#sfr_effiency: {sfr_eff}\n")
                f.write("time,gas_density,stellar_density\n")
                np.savetxt(f, data_to_save, delimiter=",", fmt="%.6f,%.6e,%.6e")
            print(f"✓ Saved density evolution CSV → {evo_out}")
    except Exception as e:
        print(f"Warning: failed to create density evolution CSV: {e}")

def main():
    import csv
    
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
    
    # Process snapshots in parallel and collect results
    snapshots = list(range(start, end + 1))
    snapshot_data = {}  # dict to store snapshot number -> data row
    
    with ProcessPoolExecutor(max_workers=args.processes) as executor:
        futures = {executor.submit(process_snapshot, numsp): numsp 
                  for numsp in snapshots}
        for future in as_completed(futures):
            num = futures[future]
            try:
                result = future.result()
                if result is not None:
                    snapshot_data[num] = result
            except Exception as exc:
                print(f"[Error] Snapshot {num:03d}: {exc}")

    print("\n=== All snapshots processed ===\n")
    
    # Write collected data to EVO_FULL CSV
    if snapshot_data:
        try:
            # Sort by snapshot number
            rows = [snapshot_data[k] for k in sorted(snapshot_data.keys())]
            
            with EVO_FULL.open('w') as f:
                f.write(f"#Timescale_fd: {fd_value}\n")
                f.write(f"#sfr_effiency: {sfr_eff}\n")
                
                writer = csv.DictWriter(f, fieldnames=rows[0].keys())
                writer.writeheader()
                writer.writerows(rows)
            
            print(f"✓ Saved {len(rows)} snapshot rows to EVO_FULL → {EVO_FULL}\n")
        except Exception as e:
            print(f"[ERROR] Failed to write EVO_FULL: {e}")
            import traceback
            traceback.print_exc()
    
    # Create density evolution CSV
    print("Creating density evolution CSV...")
    create_density_evolution_csv()
    
    # Clean up per-snapshot density files
    try:
        density_files = list(output_dir.glob("density_snapshot_*.txt"))
        for f in density_files:
            try:
                f.unlink()
            except Exception:
                pass
        if density_files:
            print(f"Removed {len(density_files)} density snapshot files.")
    except Exception as e:
        print(f"Warning: failed to remove density snapshot files: {e}")

if __name__ == "__main__":
    main()

    