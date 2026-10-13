from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from basefunc import *
from variable import *
from scipy.interpolate import UnivariateSpline
from pathlib import Path
from model_paths import model_id

modelname = model_id()
match = re.match(r'([A-Za-z]+)(\d+)', modelname)
if not match:
    raise ValueError(f"can't get '{modelname}' to dwarf name and number")

dwarf_name, num_str = match.groups()
num = int(num_str)

numsp_list = [50, 100, 150, 200, 300, 400, 500, 600, 700, 800]

print(f"Folder: {folder_path}")
print(f"Snapshots to analyze: {numsp_list}")


# Load and analyze all snapshots
results = []

for numsp in numsp_list:
    try:
        simulation = GalaxySimulation(folder_path=folder_path, snapshot_num=numsp)
        df, tsnap = simulation.load_snapshot(box=False)
        results.append({
            'numsp': numsp,
            'df': df,
            'tsnap': tsnap
        })
        print(f"✓ Snapshot {numsp}: T={tsnap:.4f} Gyr")
    except Exception as e:
        print(f"✗ Snapshot {numsp}: {e}")

print(f"\nSuccessfully loaded {len(results)} snapshots")


# Define dwarf stellar mask for each snapshot
for result in results:
    result['dw_star_mask'] = result['df'].tp >= 2
    print(f"Snapshot {result['numsp']}: {result['dw_star_mask'].sum()} stars")


# Calculate ellipticity and position angle for each snapshot
for result in results:
    df = result['df']
    mask = result['dw_star_mask']
    
    ellipticity, pa = Analysis.calculate_ellipticity(
        df[mask].x, df[mask].y, mass=df[mask].m.values, n_neighbors=30
    )
    result['ellipticity'] = ellipticity
    result['pa'] = pa
    print(f"Snapshot {result['numsp']}: ellipticity={ellipticity:.3f}, pa={pa:.3f}")


# Rotate coordinates to principal axes
for result in results:
    df = result['df']
    mask = result['dw_star_mask']
    pa = result['pa']
    
    Rotate = np.array([[np.cos(-pa), -np.sin(-pa)],
                       [np.sin(-pa),  np.cos(-pa)]])
    
    coords = np.vstack((df[mask].x, df[mask].y))  # shape: (2, N)
    rotated_coords = Rotate @ coords        # shape: (2, N)
    
    x_rot, y_rot = rotated_coords[0], rotated_coords[1]
    result['x_rot'] = x_rot
    result['y_rot'] = y_rot


# Calculate magnitude profiles
for result in results:
    df = result['df']
    mask = result['dw_star_mask']
    ellipticity = result['ellipticity']
    x_rot = result['x_rot']
    y_rot = result['y_rot']
    
    r_centers, mag_profile = Analysis.radial_magnitude_profile(
        x_rot, y_rot, mass=df[mask].m, 
        d_kpc=130, ep=ellipticity, bins=20
    )
    
    r_centers = np.array(r_centers, dtype=float)
    mag_profile = np.array(mag_profile, dtype=float)
    
    # 去掉 NaN 和 Inf
    mask_valid = np.isfinite(r_centers) & np.isfinite(mag_profile)
    r_centers_clean = r_centers[mask_valid]
    mag_profile_clean = mag_profile[mask_valid]
    
    result['r_centers_clean'] = r_centers_clean
    result['mag_profile_clean'] = mag_profile_clean
    
    print(f"Snapshot {result['numsp']}: {len(r_centers_clean)} valid points")


# Fit splines for magnitude profiles
for result in results:
    r_centers_clean = result['r_centers_clean']
    mag_profile_clean = result['mag_profile_clean']
    
    spline = UnivariateSpline(r_centers_clean, mag_profile_clean, k=3, s=15)
    
    r_smooth = np.linspace(r_centers_clean.min(), r_centers_clean.max(), 50)
    mag_smooth = spline(r_smooth)
    
    result['r_smooth_mag'] = r_smooth
    result['mag_smooth'] = mag_smooth


# This cell is now for reference only - plotting happens in the next cell after velocity dispersion is calculated
print("Data preparation complete. Running velocity dispersion calculation next...")


# Calculate velocity dispersion profiles
from scipy.stats import binned_statistic

for result in results:
    df = result['df']
    mask = result['dw_star_mask']
    
    r_xy = np.sqrt(df[mask].x**2 + df[mask].y**2)
    
    bin_edges = np.arange(0, r_xy.max() + 0.2, 0.2)
    
    sigma_z, _, _ = binned_statistic(
        r_xy, df[mask].vz, statistic='std', bins=bin_edges
    )
    
    r_mid = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    
    mask_valid = np.isfinite(sigma_z)
    r_mid_clean = r_mid[mask_valid]
    sigma_z_clean = sigma_z[mask_valid]
    
    spline_sigma = UnivariateSpline(r_mid_clean, sigma_z_clean, k=2, s=5)
    
    r_smooth_sigma = np.linspace(r_mid_clean.min(), r_mid_clean.max(), 30)
    sigma_smooth = spline_sigma(r_smooth_sigma)
    
    result['r_mid_clean'] = r_mid_clean
    result['sigma_z_clean'] = sigma_z_clean
    result['r_smooth_sigma'] = r_smooth_sigma
    result['sigma_smooth'] = sigma_smooth
    
    print(f"Snapshot {result['numsp']}: σ_z computed")


# Optional: Save combined figure
import os

output_dir = "./"

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), dpi=150)

colors = plt.cm.viridis(np.linspace(0, 1, len(results)))

# Magnitude profiles
for i, result in enumerate(results):
    tsnap = result['tsnap']
    r_smooth = result['r_smooth_mag']
    mag_smooth = result['mag_smooth']
    
    ax1.plot(r_smooth, mag_smooth, lw=1.5, color=colors[i], 
            label=f'T={tsnap:.2f} Gyr')
    
    # 补充散点
    r_centers_clean = result['r_centers_clean']
    mag_profile_clean = result['mag_profile_clean']
    ax1.scatter(r_centers_clean, mag_profile_clean, s=40, alpha=0.5, color=colors[i])

ax1.invert_yaxis()
ax1.set_xlabel("Radius (kpc)", fontsize=12)
ax1.set_ylabel("Magnitude", fontsize=12)
ax1.set_title("Magnitude Profile Evolution", fontsize=13)
ax1.grid(True, ls='--', alpha=0.3)
ax1.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=11)
ax1.minorticks_on()
ax1.legend(fontsize=10, loc='best')

# Velocity dispersion profiles
for i, result in enumerate(results):
    tsnap = result['tsnap']
    r_smooth = result['r_smooth_sigma']
    sigma_smooth = result['sigma_smooth']
    r_mid_clean = result['r_mid_clean']
    sigma_z_clean = result['sigma_z_clean']
    
    ax2.plot(r_smooth, sigma_smooth, lw=2.5, color=colors[i], 
            label=f'T={tsnap:.2f} Gyr')
    ax2.scatter(r_mid_clean, sigma_z_clean, s=40, alpha=0.5, color=colors[i])

ax2.set_xlabel("Radius (kpc)", fontsize=12)
ax2.set_ylabel("σ_z (km/s)", fontsize=12)
ax2.set_title("Velocity Dispersion Evolution", fontsize=13)
ax2.grid(True, ls='--', alpha=0.3)
ax2.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=11)
ax2.minorticks_on()
ax2.legend(fontsize=10, loc='best')

plt.suptitle(f"IC{dwarf_name}{num_str}", 
            fontsize=14, fontweight='bold', y=1.00)
plt.tight_layout()

# Save
# snapshot_str = "_".join(str(n) for n in numsp_list)
outname = os.path.join(output_dir, f"mag_sigma.png")
plt.savefig(outname, dpi=200, bbox_inches='tight')
print(f"✓ Saved to {outname}")
plt.show()


