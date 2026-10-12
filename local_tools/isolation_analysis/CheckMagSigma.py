#!/usr/bin/env python3
"""
CheckMagSigma.py - Analyze magnitude profile and velocity dispersion evolution
Saves combined plots (magnitude + sigma_z) at each snapshot
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import binned_statistic
from scipy.interpolate import UnivariateSpline
import argparse
import sys

from basefunc import *
from variable import *

# ======================================================
# --- Configuration ---
# ======================================================
def get_config():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(description='Analyze magnitude and velocity dispersion profiles')
    parser.add_argument('--dwarf', type=str, default='Fornax', help='Dwarf galaxy name')
    parser.add_argument('--modelnum', type=int, default=45, help='Model number')
    parser.add_argument('--snapshots', type=str, default='0,300', help='Snapshot range (e.g., "0,300")')
    parser.add_argument('--output-dir', type=str, default='./mag_sigma_output', help='Output directory')
    args = parser.parse_args()
    return args

# ======================================================
# --- Helper Functions ---
# ======================================================
def analyze_snapshot(folder_path, snapshot_num, output_dir):
    """
    Analyze a single snapshot: compute magnitude profile and velocity dispersion
    
    Returns:
        dict with keys: 'snapshot', 'time_gyr', 'mag_profile', 'sigma_profile', 'r_centers', 'r_mid'
        or None if analysis fails
    """
    try:
        # Load snapshot
        simulation = GalaxySimulation(folder_path=folder_path, snapshot_num=snapshot_num)
        df, tsnap = simulation.load_snapshot(box=False)
        
        # Define dwarf stellar mask (simplified: all particles with tp >= 2)
        dw_star_mask = df.tp >= 2
        
        if dw_star_mask.sum() < 10:
            print(f"[Snapshot {snapshot_num:03d}] Skipped: too few stars ({dw_star_mask.sum()})")
            return None
        
        # Calculate ellipticity and position angle
        ellipticity, pa = Analysis.calculate_ellipticity(
            df[dw_star_mask].x, df[dw_star_mask].y, 
            mass=df[dw_star_mask].m.values, n_neighbors=30
        )
        
        # Rotate to principal axes
        Rotate = np.array([[np.cos(-pa), -np.sin(-pa)],
                           [np.sin(-pa),  np.cos(-pa)]])
        coords = np.vstack((df[dw_star_mask].x, df[dw_star_mask].y))
        rotated_coords = Rotate @ coords
        x_rot, y_rot = rotated_coords[0], rotated_coords[1]
        
        # --- Magnitude Profile ---
        r_centers, mag_profile = Analysis.radial_magnitude_profile(
            x_rot, y_rot, mass=df[dw_star_mask].m, 
            d_kpc=130, ep=ellipticity, bins=80
        )
        r_centers = np.array(r_centers, dtype=float)
        mag_profile = np.array(mag_profile, dtype=float)
        
        # Clean magnitude data
        mask_mag = np.isfinite(r_centers) & np.isfinite(mag_profile)
        r_centers_clean = r_centers[mask_mag]
        mag_profile_clean = mag_profile[mask_mag]
        
        # Spline fit for magnitude
        if len(r_centers_clean) > 3:
            spline_mag = UnivariateSpline(r_centers_clean, mag_profile_clean, k=3, s=15)
            r_smooth_mag = np.linspace(r_centers_clean.min(), r_centers_clean.max(), 50)
            mag_smooth = spline_mag(r_smooth_mag)
        else:
            r_smooth_mag = r_centers_clean
            mag_smooth = mag_profile_clean
        
        # --- Velocity Dispersion Profile (σz vs R) ---
        r_xy = np.sqrt(df[dw_star_mask].x**2 + df[dw_star_mask].y**2)
        bin_edges = np.arange(0, r_xy.max() + 0.5, 0.5)
        
        sigma_z, _, _ = binned_statistic(
            r_xy, df[dw_star_mask].vz, statistic='std', bins=bin_edges
        )
        r_mid = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        
        # Clean sigma data
        mask_sigma = np.isfinite(sigma_z)
        r_mid_clean = r_mid[mask_sigma]
        sigma_z_clean = sigma_z[mask_sigma]
        
        # Spline fit for sigma
        if len(r_mid_clean) > 2:
            spline_sigma = UnivariateSpline(r_mid_clean, sigma_z_clean, k=2, s=5)
            r_smooth_sigma = np.linspace(r_mid_clean.min(), r_mid_clean.max(), 30)
            sigma_smooth = spline_sigma(r_smooth_sigma)
        else:
            r_smooth_sigma = r_mid_clean
            sigma_smooth = sigma_z_clean
        
        print(f"[Snapshot {snapshot_num:03d}] T={tsnap:.4f} Gyr: Analyzed successfully")
        
        return {
            'snapshot': snapshot_num,
            'time_gyr': tsnap,
            'r_smooth_mag': r_smooth_mag,
            'mag_smooth': mag_smooth,
            'r_smooth_sigma': r_smooth_sigma,
            'sigma_smooth': sigma_smooth,
            'r_mid_clean': r_mid_clean,
            'sigma_z_clean': sigma_z_clean,
        }
        
    except Exception as e:
        print(f"[Error] Snapshot {snapshot_num:03d}: {e}")
        return None

# ======================================================
# --- Main Processing ---
# ======================================================
def main():
    args = get_config()
    
    # Setup paths and directories
    folder_path = f"./{args.dwarf}{args.modelnum}/"
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Parse snapshot range
    try:
        start, end = map(int, args.snapshots.split(','))
    except ValueError:
        print(f"Invalid snapshot range format. Using default: 0,300")
        start, end = 0, 300
    
    print(f"\n=== Analyzing {args.dwarf}{args.modelnum} ===")
    print(f"Snapshots: {start} – {end}")
    print(f"Output: {output_dir}\n")
    
    # Process snapshots
    results = []
    for numsp in range(start, end + 1):
        result = analyze_snapshot(folder_path, numsp, output_dir)
        if result is not None:
            results.append(result)
    
    if not results:
        print("No snapshots analyzed successfully!")
        return
    
    # Sort by time
    results.sort(key=lambda x: x['time_gyr'])
    
    print(f"\n=== Generated plots for {len(results)} snapshots ===\n")
    
    # Generate combined plots for each snapshot
    for result in results:
        numsp = result['snapshot']
        tsnap = result['time_gyr']
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), dpi=150)
        
        # Left: Magnitude profile
        ax1.plot(result['r_smooth_mag'], result['mag_smooth'], lw=2.5, color='navy')
        ax1.invert_yaxis()
        ax1.set_xlabel("Radius (kpc)", fontsize=12)
        ax1.set_ylabel("Magnitude", fontsize=12)
        ax1.set_title(f"Magnitude Profile (T={tsnap:.4f} Gyr)", fontsize=13)
        ax1.grid(True, ls='--', alpha=0.3)
        ax1.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=11)
        ax1.minorticks_on()
        
        # Right: Velocity dispersion
        ax2.plot(result['r_smooth_sigma'], result['sigma_smooth'], lw=2.5, color='darkred', label='Spline fit')
        ax2.scatter(result['r_mid_clean'], result['sigma_z_clean'], s=50, alpha=0.6, color='darkred', label='Binned data')
        ax2.set_xlabel("Radius (kpc)", fontsize=12)
        ax2.set_ylabel("σ_z (km/s)", fontsize=12)
        ax2.set_title(f"Velocity Dispersion (T={tsnap:.4f} Gyr)", fontsize=13)
        ax2.grid(True, ls='--', alpha=0.3)
        ax2.legend(fontsize=10, loc='best')
        ax2.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=11)
        ax2.minorticks_on()
        
        plt.suptitle(f"IC_{args.dwarf}{args.modelnum} | Snapshot {numsp:03d}", 
                     fontsize=14, fontweight='bold', y=1.00)
        plt.tight_layout()
        
        # Save figure
        outname = output_dir / f"snapshot_{numsp:03d}_mag_sigma.png"
        plt.savefig(outname, dpi=200, bbox_inches='tight')
        plt.close(fig)
        
        print(f"[Snapshot {numsp:03d}] Saved → {outname}")
    
    # Create summary plot showing evolution (all snapshots on one figure)
    if len(results) > 1:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), dpi=150)
        
        times = [r['time_gyr'] for r in results]
        
        # Plot mag profiles
        for i, result in enumerate(results):
            alpha = 0.3 + 0.7 * (i / len(results))  # Gradual fade
            ax1.plot(result['r_smooth_mag'], result['mag_smooth'], 
                    lw=1.5, alpha=alpha, 
                    label=f"T={result['time_gyr']:.3f} Gyr")
        ax1.invert_yaxis()
        ax1.set_xlabel("Radius (kpc)", fontsize=12)
        ax1.set_ylabel("Magnitude", fontsize=12)
        ax1.set_title("Magnitude Profile Evolution", fontsize=13)
        ax1.grid(True, ls='--', alpha=0.3)
        ax1.legend(fontsize=9, loc='best', ncol=2)
        ax1.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=11)
        ax1.minorticks_on()
        
        # Plot sigma profiles
        for i, result in enumerate(results):
            alpha = 0.3 + 0.7 * (i / len(results))  # Gradual fade
            ax2.plot(result['r_smooth_sigma'], result['sigma_smooth'], 
                    lw=1.5, alpha=alpha, 
                    label=f"T={result['time_gyr']:.3f} Gyr")
        ax2.set_xlabel("Radius (kpc)", fontsize=12)
        ax2.set_ylabel("σ_z (km/s)", fontsize=12)
        ax2.set_title("Velocity Dispersion Evolution", fontsize=13)
        ax2.grid(True, ls='--', alpha=0.3)
        ax2.legend(fontsize=9, loc='best', ncol=2)
        ax2.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=11)
        ax2.minorticks_on()
        
        plt.suptitle(f"IC_{args.dwarf}{args.modelnum} | Evolution from T={times[0]:.3f} to {times[-1]:.3f} Gyr", 
                     fontsize=14, fontweight='bold', y=0.995)
        plt.tight_layout()
        
        evo_name = output_dir / "evolution_mag_sigma.png"
        plt.savefig(evo_name, dpi=200, bbox_inches='tight')
        plt.close(fig)
        
        print(f"\n[Evolution] Saved → {evo_name}")
    
    print(f"\n✓ All done! Results saved to {output_dir}")

if __name__ == "__main__":
    main()
