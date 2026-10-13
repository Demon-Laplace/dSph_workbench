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

from mpl_toolkits.axes_grid1 import make_axes_locatable

import argparse
# from multiprocessing import Pool
from multiprocessing import get_context

warnings.filterwarnings("ignore")
plt.switch_backend('Agg')

modelname = model_id()
match = re.match(r'([A-Za-z]+)(\d+)', modelname)
if not match:
    raise ValueError(f"can't get '{modelname}' to dwarf name and number")

dwarf_name, num_str = match.groups()
num = int(num_str)
inifile = f'./{dwarf_name}.ini'
dw_particles = GalaxySimulation.get_dw_num(inifile)

modelfolder = './output'
elinfo_path = f'./elinfo_{modelname}.csv'

core_radius = core_radius_dict.get(dwarf_name, None)
if core_radius is None:
    raise ValueError(f"please input '{dwarf_name}''s radius")

dw_elinfo = pd.read_csv(elinfo_path, skiprows=2)

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    print(f"Created output folder: {output_folder}")

ini_x, ini_y, ini_z, ini_vx, ini_vy, ini_vz, fd_value = DataProcessor.get_ini_data(elinfo_path=elinfo_path)
fd_factor_raw = fd_value / 8.3
if fd_factor_raw >= 0.90:
    fd_factor_str = f"{round(fd_factor_raw):.0f}"
else:
    fd_factor_str = f"{fd_factor_raw:.1f}"

dw_sfrmass = Analysis.get_sfh(dw_particles=dw_particles, folder_path=folder_path)

def safe_plot(numsp):
    try:
        plot(numsp)
    except Exception as e:
        print(f"Error in numsp={numsp}: {e}")
        raise

def plot(numsp, dw_particles=dw_particles, dw_elinfo=dw_elinfo):

    star_mass = dw_elinfo.loc[dw_elinfo['numsp'] == numsp, 'star_mass'].iloc[0]
    star_half_mass = dw_elinfo.loc[dw_elinfo['numsp'] == numsp, 'star_half_mass'].iloc[0]
    coldgas_mass = dw_elinfo.loc[dw_elinfo['numsp'] == numsp, 'coldgas_mass'].iloc[0]
    coldgas_half_mass = dw_elinfo.loc[dw_elinfo['numsp'] == numsp, 'coldgas_half_mass'].iloc[0]

    gas_frac = coldgas_mass / star_mass
    gas_half_frac = coldgas_half_mass / star_half_mass

    sm_str = DataProcessor.sci_notation(star_mass)
    shm_str = DataProcessor.sci_notation(star_half_mass)
    cg_str = DataProcessor.sci_notation(coldgas_mass)
    cghm_str = DataProcessor.sci_notation(coldgas_half_mass)

    simulation = GalaxySimulation(folder_path=folder_path, snapshot_num=numsp)
    df, tsnap = simulation.load_snapshot()
    mw_mask = simulation.find_MW_particles(dw_particles_num=dw_particles, cube_file=cube_path)
    dw_mask = simulation.find_dwarf_particles(dw_particles_num=dw_particles)
    df = simulation.center_on_MW(mw_mask, iterations=2, radius=100)
    df = simulation.convert_coordinates(df=df)

    # mw_star_mask = mw_mask & ((df.tp == 2) | (df.tp == 3))
    dw_star_mask = dw_mask & ((df.tp == 2) | (df.tp == 3))
    dw_gas_mask = dw_mask & (df.tp == 0)
    # mw_gas_mask = mw_mask & (df.tp == 0)
    # new_star_mask = df.tp == 4
    dw_hot_gas_mask = dw_gas_mask & (df.temp >= 20000)
    dw_cold_gas_mask = dw_gas_mask & (df.temp < 20000)
    new_dw_star_mask, new_mw_star_mask = simulation.classify_new_dwarf_stars(dw_star_mask)
    total_dw_star_mask = dw_star_mask | new_dw_star_mask
    # total_mw_star_mask = mw_star_mask | new_mw_star_mask 

    dw_cra, dw_cdec = Analysis.find_center(df[total_dw_star_mask].ra.to_numpy(), df[total_dw_star_mask].dec.to_numpy())
    rotra_all, rotdec_all = Analysis.rotate_to_sky(df.ra.to_numpy(), df.dec.to_numpy(), center=[dw_cra, dw_cdec])
    # mw_center_ra, mw_center_dec = CoordinateTransform.galactic_to_equatorial([0], [0])
    # mw_rotra_center, mw_rotdec_center = Analysis.rotate_to_sky(mw_center_ra, mw_center_dec, center=[dw_cra, dw_cdec])

    rotra_dw_star, rotdec_dw_star = rotra_all[total_dw_star_mask], rotdec_all[total_dw_star_mask]
    rotra_dw_hot_gas, rotdec_dw_hot_gas = rotra_all[dw_hot_gas_mask], rotdec_all[dw_hot_gas_mask]
    rotra_dw_cold_gas, rotdec_dw_cold_gas = rotra_all[dw_cold_gas_mask], rotdec_all[dw_cold_gas_mask]

    d_mean = np.mean(df[total_dw_star_mask].rh)
    d_mean_gal = np.mean(df[total_dw_star_mask].r)
    x_kpc = np.radians(np.array(rotra_dw_star)) * d_mean
    y_kpc = np.radians(np.array(rotdec_dw_star)) * d_mean
    hot_gas_x_kpc = np.radians(np.array(rotra_dw_hot_gas)) * d_mean
    hot_gas_y_kpc = np.radians(np.array(rotdec_dw_hot_gas)) * d_mean
    cold_gas_x_kpc = np.radians(np.array(rotra_dw_cold_gas)) * d_mean
    cold_gas_y_kpc = np.radians(np.array(rotdec_dw_cold_gas)) * d_mean

    img, extent = Analysis.surface_brightness_map(rotra_dw_star, rotdec_dw_star, df[total_dw_star_mask].m, d_mean,
                                                  npix=130, size_kpc=4, mag_sys='vega')
    
    df_sub = df[total_dw_star_mask]
    pos = CartesianRepresentation(df_sub['x'].values * u.kpc,
                                df_sub['y'].values * u.kpc,
                                df_sub['z'].values * u.kpc)

    vel = CartesianDifferential(df_sub['vx'].values * u.km/u.s,
                                df_sub['vy'].values * u.km/u.s,
                                df_sub['vz'].values * u.km/u.s)

    coord = SkyCoord(pos.with_differentials(vel), frame='galactocentric')
    coord_icrs = coord.transform_to('icrs')

    mean_pmra = np.mean(coord_icrs.pm_ra_cosdec).value
    mean_pmdec = np.mean(coord_icrs.pm_dec).value
    
    ellipticity, pa = Analysis.calculate_ellipticity(x_kpc, y_kpc, mass=df[total_dw_star_mask].m.to_numpy(), n_neighbors=30)

    Rotate = np.array([[np.cos(-pa), -np.sin(-pa)],
                   [np.sin(-pa),  np.cos(-pa)]])

    coords = np.vstack((x_kpc, y_kpc))  # shape: (2, N)
    rotated_coords = Rotate @ coords        # shape: (2, N)

    x_rot, y_rot = rotated_coords[0], rotated_coords[1]

    r_centers, mag_profile = Analysis.radial_magnitude_profile(x_rot, y_rot, mass=df[total_dw_star_mask].m,
                                                           d_kpc=d_mean, ep=ellipticity, bins=40)
    r_half = Analysis.half_light_radius(x_rot, y_rot, df[total_dw_star_mask].m.to_numpy(), ep=ellipticity)
    vlos_array = Analysis.cal_vlos(df[total_dw_star_mask].xh.to_numpy(), df[total_dw_star_mask].yh.to_numpy(), 
                                   df[total_dw_star_mask].zh.to_numpy(),
                                   df[total_dw_star_mask].vxh.to_numpy(), df[total_dw_star_mask].vyh.to_numpy(), 
                                   df[total_dw_star_mask].vzh.to_numpy())
    
    bin_centers, vlos_dispersion = Analysis.get_vlos_dispersion(x_kpc, y_kpc, vlos_array)

    img, extent = Analysis.surface_brightness_map(rotra_dw_star, rotdec_dw_star, df[total_dw_star_mask].m, d_mean,
                                                npix=120, size_kpc=4, mag_sys='vega')
    vf, extent_vf = Analysis.get_vlos_field(x_kpc, y_kpc, vlos_array, npix=40)

    dprof = pd.read_csv(fndens)
    factor = 1.0385305633697762e-11
    I_obs = dprof['density'] * factor
    I_obs_err = dprof['density_error'] * factor
    mag_obs = -2.5 * np.log10(I_obs)
    mag_obs_err = 1.0857 * (I_obs_err / I_obs)

    sigma_mw = np.where(dw_elinfo.sigma >= dw_elinfo.tsigma, 
                        np.sqrt(dw_elinfo.sigma**2 - dw_elinfo.tsigma**2), 
                        0)
    
    sfh_dw = pd.read_csv(f'../SFH_{dwarf_name}.csv')
    Tage = dw_elinfo.loc[
        dw_elinfo.loc[:dw_elinfo['distance_gal'].idxmin(), 'distance_gal']
        .sub(d_today)
        .abs()
        .idxmin(),
        'age'
    ]
    sfh_dw = DataProcessor.calsfh(sfh_dw, Tage)

    tidal_eff, dw_energy = Analysis.cal_tidal_eff(dw_elinfo=dw_elinfo)
    tidal_eff_numsp, dw_energy_nusmp = Analysis.cal_tidal_eff(dw_elinfo=dw_elinfo[dw_elinfo['numsp'].isin([numsp])])

    index_peri = dw_elinfo['distance_gal'].idxmin()
    d_peri = dw_elinfo.loc[index_peri, 'distance_gal']
    T_peri = dw_elinfo.loc[index_peri, 'age']

    vr = dw_elinfo[dw_elinfo['numsp'].isin([numsp])].vr.to_numpy()[0]
    vtheta = dw_elinfo[dw_elinfo['numsp'].isin([numsp])].vtheta.to_numpy()[0]
    vphi = dw_elinfo[dw_elinfo['numsp'].isin([numsp])].vphi.to_numpy()[0]

    vgsr = np.sqrt(vr**2 + vtheta**2 + vphi**2)

    ### --- ###

    fig = plt.figure(figsize=(21, 15), dpi=150)
    gs = gridspec.GridSpec(3, 4)

    # Create axes excluding (2,3) position
    axes = [[None]*4 for _ in range(3)]
    for i in range(3):
        for j in range(4):
            if not (i == 2 and j == 3):  # Skip the (2,3) position
                axes[i][j] = fig.add_subplot(gs[i, j])

    axes = np.array(axes) 

    plt.subplots_adjust(
        left=0.05,   
        right=0.95, 
        top=0.95,   
        bottom=0.08,
        wspace=0.35, 
        hspace=0.22
    )

    #--- Title ---#

    fig.suptitle(f'{modelname}', fontsize=22)

    #--- Gas sky map (kpc) ---#

    axes[1, 1].scatter(20, 20, s=18, color='aquamarine', label='HI Gas')
    axes[1, 1].scatter(20, 20, s=18, color='orange', label='HII Gas')
    axes[1, 1].scatter(cold_gas_x_kpc, cold_gas_y_kpc, s=4, color='aquamarine')
    axes[1, 1].scatter(hot_gas_x_kpc, hot_gas_y_kpc, s=5, color='orange')
    axes[1, 1].set_xlim(-5, 5)
    axes[1, 1].set_ylim(-5, 5)
    axes[1, 1].set_xlabel('RA (kpc)', fontsize=xylabel_size)
    axes[1, 1].set_ylabel('DEC (kpc)', fontsize=xylabel_size)
    axes[1, 1].minorticks_on()
    axes[1, 1].tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=15)
    axes[1, 1].legend(loc='upper left', fontsize=10)
    axes[1, 1].invert_xaxis()
    axes[1, 1].grid(True)

    #--- Stellar sky map (deg) ---#

    arrow_length = 1.5
    dx = mean_pmra
    dy = mean_pmdec
    norm = np.sqrt(dx**2 + dy**2)
    dx_scaled = dx / norm * arrow_length
    dy_scaled = dy / norm * arrow_length

    axes[0, 0].arrow(0, 0, dx_scaled, dy_scaled, 
                    head_width=0.1, head_length=0.1,
                    fc='black', ec='black', alpha=0.7, lw=2)

    im = axes[0, 0].imshow(img.T, origin='lower', extent=extent, cmap='bone_r',
                        vmin=22, vmax=34)
    axes[0, 0].set_xlim(-2, 2)
    axes[0, 0].set_ylim(-2, 2)
    axes[0, 0].set_xlabel("RA (deg)", fontsize=xylabel_size)
    axes[0, 0].set_ylabel("DEC (deg)", fontsize=xylabel_size)
    axes[0, 0].invert_xaxis()
    axes[0, 0].minorticks_on()
    axes[0, 0].tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=15)

    divider = make_axes_locatable(axes[0, 0])
    cax = divider.append_axes("right", size="4%", pad=0.05)
    cbar = fig.colorbar(im, cax=cax, label=r'$\mu_V$ (mag arcsec$^{-2}$)')
    cbar.set_label(r'$\mu_V$ (mag arcsec$^{-2}$)', fontsize=14)
    cax.tick_params(axis='both', which='both', direction='in')
    cax.minorticks_on()

    #--- vlos sky map (deg) ---#

    sc = axes[0, 1].imshow(vf.T, origin='lower', extent=extent_vf, cmap='coolwarm', vmin=-12, vmax=12)
    axes[0, 1].set_xlim(-5, 5)
    axes[0, 1].set_ylim(-5, 5)
    axes[0, 1].set_xlabel("RA (kpc)", fontsize=xylabel_size)
    axes[0, 1].set_ylabel("DEC (kpc)", fontsize=xylabel_size)
    axes[0, 1].invert_xaxis()
    axes[0, 1].minorticks_on()
    axes[0, 1].tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=15)

    divider = make_axes_locatable(axes[0, 1])
    cax = divider.append_axes("right", size="4%", pad=0.05)
    cbar = fig.colorbar(sc, cax=cax)
    cbar.set_label(r'<vlos> (km/s)', fontsize=14)  
    cax.tick_params(axis='both', which='both', direction='in')
    cax.minorticks_on()

    #--- 3D Scene in Y-Z ---#

    cy, cz = Analysis.find_center(df[total_dw_star_mask].y, df[total_dw_star_mask].z, units='kpc')
    arrow_length = 7
    total_distance = np.sqrt(cy**2 + cz**2)
    if total_distance > 0:  

        normalized_cz = (-cy/total_distance) * arrow_length
        normalized_cy = (-cz/total_distance) * arrow_length
    else:
        normalized_cz = normalized_cy = 0

    axes[0, 3].arrow(0, 0, normalized_cz, normalized_cy, 
                    head_width=1.2, head_length=1.2,
                    fc='black', ec='black', zorder=6)

    axes[0, 3].scatter(df[dw_cold_gas_mask].y.to_numpy() - cy, df[dw_cold_gas_mask].z.to_numpy() - cz, c='aquamarine', s=1.5, lw=0, alpha=0.5)
    axes[0, 3].scatter(df[dw_hot_gas_mask].y.to_numpy() - cy, df[dw_hot_gas_mask].z.to_numpy() - cz, c='orange', s=3, lw=0)
    axes[0, 3].scatter(df[total_dw_star_mask].y.to_numpy() - cy, df[total_dw_star_mask].z.to_numpy() - cz, c='gray', s=2, lw=0)
    axes[0, 3].scatter(50, 50, s=70, color='aquamarine', label='HI Gas')
    axes[0, 3].scatter(50, 50, s=70, color='orange', label='HII Gas')
    axes[0, 3].scatter(50, 50, s=70, color='gray', label='Stars')
    axes[0, 3].set_xlim(-20, 20)
    axes[0, 3].set_ylim(-20, 20)
    axes[0, 3].set_xlabel("Y (kpc)", fontsize=xylabel_size)
    axes[0, 3].set_ylabel("Z (kpc)", fontsize=xylabel_size)
    axes[0, 3].minorticks_on()
    axes[0, 3].tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=15)
    axes[0, 3].legend(loc='lower left', fontsize=14)

    #--- velocity dispersion ---#

    axes[1, 0].scatter(bin_centers, vlos_dispersion, marker='o', color='black', label='Simulation data', s=8, zorder=3)
    axes[1, 0].errorbar(r_pc/1000, sigma, yerr=err_sigma, fmt='o', color='darkorchid', ecolor='gray',
                        elinewidth=1, label='Walker 2009', capsize=4, alpha=0.6, markersize=4, zorder=2)
    axes[1, 0].hlines(dw_elinfo[dw_elinfo['numsp'].isin([numsp])].tsigma, 0, 3, lw=2, ls='--', color='gray', zorder=0)
    axes[1, 0].vlines(r_half, 0, 15, lw=2, ls='--', color='gray', zorder=1)
    axes[1, 0].scatter(r_half, dw_elinfo[dw_elinfo['numsp'].isin([numsp])].tsigma, c='red', s=80, zorder=10)

    axes[1, 0].set_xlabel("Radius (kpc)", fontsize=xylabel_size)
    axes[1, 0].set_ylabel("Velocity Dispersion (km/s)", fontsize=xylabel_size)
    # axes[1, 0].set_title("LOS Velocity Dispersion vs. Radius")
    axes[1, 0].set_xlim(0, 2.2)
    axes[1, 0].set_ylim(4, 13.5)
    axes[1, 0].grid(True)
    axes[1, 0].legend(fontsize=12, loc='upper right')
    axes[1, 0].minorticks_on()
    axes[1, 0].tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=15)

    #--- Surface luminosity profile ---#

    axes[2, 0].errorbar(np.radians(dprof['xr']) * d_today , mag_obs,
                        yerr = mag_obs_err, fmt='o', capsize=4, color='darkorchid', ecolor='gray',
                        markersize=4, elinewidth=1, alpha=0.5, zorder=1, label='Yang 2022')

    axes[2, 0].scatter(r_centers, mag_profile, label='Simulation data', c='black', s=10, zorder=3)
    axes[2, 0].set_xlabel("Radius (kpc)", fontsize=xylabel_size)
    axes[2, 0].set_ylabel("Magnitude", fontsize=xylabel_size)
    axes[2, 0].legend(fontsize=12, loc='upper right')
    axes[2, 0].set_ylim(36, 23)
    axes[2, 0].set_xlim(0, 5)
    axes[2, 0].minorticks_on()
    axes[2, 0].tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=15)
    axes[2, 0].grid(True)

    #--- Gas fraction ---#

    axes[1, 2].plot(dw_elinfo.age.to_numpy(), 
                    dw_elinfo.coldgas_half_mass.to_numpy() / dw_elinfo.star_half_mass.to_numpy(), lw=3, c='black', 
                    label=r'within $R_{\mathrm{half}}$')
    axes[1, 2].plot(dw_elinfo.age.to_numpy(), 
                    dw_elinfo.coldgas_mass.to_numpy() / dw_elinfo.star_mass.to_numpy(), lw=3, c='gray',
                    label=r'within stellar region')
    axes[1, 2].scatter(tsnap, 
                       coldgas_half_mass / star_half_mass, 
                    c='red', marker='o', s=100, zorder=5)
    axes[1, 2].scatter(tsnap, 
                       dw_elinfo[dw_elinfo['numsp'].isin([numsp])].coldgas_mass.to_numpy() / dw_elinfo[dw_elinfo['numsp'].isin([numsp])].star_mass.to_numpy(), 
                    c='salmon', marker='o', s=100, zorder=5)
    axes[1, 2].vlines(tsnap, 0, 10, color='gray', lw=2, ls='--', zorder=0)
    axes[1, 2].minorticks_on()
    axes[1, 2].tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=15)
    axes[1, 2].grid(True)
    axes[1, 2].legend(loc='lower left', fontsize=15)
    axes[1, 2].set_yscale('log')
    axes[1, 2].set_xlim(0, T_peri)
    axes[1, 2].set_xlabel("Simulation Time (Gyr)", fontsize=xylabel_size)
    axes[1, 2].set_ylabel(r'$M_{\mathrm{HI}} / M_{stars}$', fontsize=xylabel_size)
    axes[1, 2].set_ylim(0, 10.8)

    #--- Sigma history ---#

    axes[0, 2].plot(dw_elinfo.age.to_numpy(), dw_elinfo.sigma.to_numpy(), lw=2, color='black', label=r'$\sigma_{los}$')
    axes[0, 2].plot(dw_elinfo.age.to_numpy(), dw_elinfo.tsigma.to_numpy(), lw=2, color='gray', label=r'$\sigma_{stars}$')
    axes[0, 2].plot(dw_elinfo.age.to_numpy(), sigma_mw, lw=2, color='c', label=r'$\sigma_{MW}$')
    axes[0, 2].scatter(tsnap, dw_elinfo[dw_elinfo['numsp'].isin([numsp])].sigma, s=70, zorder=5, color='red')
    axes[0, 2].vlines(tsnap, 0, 20, color='gray', lw=2, ls='--', zorder=0)
    axes[0, 2].minorticks_on()
    axes[0, 2].tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=15)
    axes[0, 2].set_ylim(-0.08, 13.5)
    axes[0, 2].set_xlim(0, T_peri)
    axes[0, 2].grid(True)
    axes[0, 2].legend(loc='upper left', fontsize=15)
    axes[0, 2].set_xlabel("Simulation Time (Gyr)", fontsize=xylabel_size)
    axes[0, 2].set_ylabel("velocity dispersion (km/s)", fontsize=xylabel_size)

    #------#

    #--- Evolution Histroy ---#

    axes[2, 2].plot(dw_elinfo.age.to_numpy(), dw_elinfo.coldgas_mass.to_numpy(), lw=2.5, c='aquamarine', label='HI mass within Stellar region')
    axes[2, 2].plot(dw_elinfo.age.to_numpy(), dw_elinfo.coldgas_half_mass.to_numpy(), lw=2.5, c='teal', label='HI mass within $R_{{half}}$')
    axes[2, 2].plot(dw_elinfo.age.to_numpy(), dw_elinfo.star_mass.to_numpy(), lw=2, c='lightblue', alpha=0.9, label='Stellar mass')
    axes[2, 2].plot(dw_elinfo.age.to_numpy(), dw_elinfo.star_half_mass.to_numpy(), lw=2, c='darkblue', label='Stellar mass within $R_{{half}}$')

    axes[2, 2].step(sfh_dw['cosmic_time_shifted'].to_numpy(), 1e9*sfh_dw['cumulative_mass'].to_numpy(), 
                    lw=2.5, color='maroon', ls='--', label='SFH in Observation')
    axes[2, 2].step(sfh_dw['cosmic_time_shifted'].to_numpy(), 1e9*sfh_dw['pure_sfh'].to_numpy(), 
                    lw=2.5, color='coral', ls='--', label='Pure SFH in Observation')
    if dw_sfrmass is not None:
        axes[2, 2].step(dw_sfrmass['birth'].to_numpy(), dw_sfrmass['cumulative_mass'].to_numpy(), lw=2, color='indigo', label='SFH in Simulation')

    axes[2, 2].vlines(tsnap, 0, 1e11, color='gray', lw=2, ls='--')
    axes[2, 2].hlines(dw_cold_gas_lower_limit, -0.1, T_peri*2, lw=2, ls='--', color='paleturquoise', label='lower limit of HI mass')
    axes[2, 2].minorticks_on()
    axes[2, 2].tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=15)
    axes[2, 2].set_yscale('log')
    axes[2, 2].set_xlim(0, T_peri)
    axes[2, 2].set_ylim(1e5, 8e10)
    axes[2, 2].grid(True)
    axes[2, 2].legend(loc='upper right')
    axes[2, 2].set_ylabel(r'Mass ($M_{\odot}$)', fontsize=xylabel_size)
    axes[2, 2].set_xlabel('Simulation Time (Gyr)', fontsize=xylabel_size)

    #--- oribt ---#

    axes[1, 3].plot(dw_elinfo.age.to_numpy(), dw_elinfo.distance_gal.to_numpy(), color='black', lw=3, zorder=1)
    axes[1, 3].scatter(tsnap, dw_elinfo[dw_elinfo['numsp'].isin([numsp])].distance_gal.to_numpy(), color='red', s=80, zorder=3, label=f'{tsnap:.2f} Gyr')
    axes[1, 3].hlines(d_today, -0.1, T_peri * 3, color='gray', lw=2, ls='--', zorder=0)
    axes[1, 3].minorticks_on()
    axes[1, 3].tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=15)
    axes[1, 3].set_xlim(0, T_peri*1.2)
    axes[1, 3].set_ylim(50, 820)
    axes[1, 3].legend(loc='lower left', fontsize=16)
    axes[1, 3].set_xlabel('Simulation Time (Gyr)', fontsize=xylabel_size)
    axes[1, 3].set_ylabel('Distance (kpc)', fontsize=xylabel_size)

    #--- Tidal Effect ---#

    te = axes[2, 1].scatter(np.log10(tidal_eff), np.log10(dw_energy), c=dw_elinfo['age'].to_numpy(), s=10, vmin=0, vmax=max(dw_elinfo['age']), cmap='PuBu')
    axes[2, 1].scatter(np.log10(tidal_eff_numsp), np.log10(dw_energy_nusmp), c='orangered', s=260, marker='+', lw=3)
    axes[2, 1].minorticks_on()
    axes[2, 1].tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=10)
    axes[2, 1].grid(True)
    axes[2, 1].set_xlim(-33.9, -29.1)
    axes[2, 1].set_ylim(-33.9, -29.1)
    axes[2, 1].set_ylabel(r'$log_{10} (\frac{\sigma^2_{los}}{r^2_{half}}) ~ (s^{-2})$', fontsize=xylabel_size)
    axes[2, 1].set_xlabel(r'$log_{10} (\frac{128}{27} \frac{g^{2}_{MW}}{v^2}) ~ (s^{-2})$', fontsize=xylabel_size)

    x_tidal = np.linspace(-100, 0, 20)
    y_tidal = x_tidal

    axes[2, 1].plot(x_tidal, y_tidal, lw=2, ls='--', c='orchid')

    divider = make_axes_locatable(axes[2, 1])
    cax = divider.append_axes("right", size="4%", pad=0.05)
    cbar = fig.colorbar(te, cax=cax)
    cbar.set_label(r'Simulation Time (Gyr)', fontsize=12)
    cax.tick_params(axis='both', which='both', direction='in')
    cax.minorticks_on()

    #--- Infomation ---#

    info_orbit_text = (
        f"Orbit : {ini_x:.2f}, {ini_y:.2f}, {ini_z:.2f}\n"
        f"        {ini_vx:.2f}, {ini_vy:.2f}, {ini_vz:.2f}\n"
    )

    # info_mass_text = (
    #     f"Star Mass: {sm_str}\n"
    #     f"Star Half-Mass: {shm_str}\n"
    #     f"Cold Gas Mass: {cg_str}\n"
    #     f"Cold Gas Half-Mass: {cghm_str}"
    # )

    axes[1, 3].text(
        0.05, 0.95,
        f'D : {d_mean_gal:.2f} kpc\n'
        f'$v_{{r}}$ : {vr:.2f} km/s\n'
        f'$v_{{\\theta}}$ : {vtheta:.2f} km/s\n'
        f'$v_{{\\phi}}$ : {vphi:.2f} km/s\n'
        f'$v_{{gsr}}$ : {vgsr:.2f} km/s',
        transform=axes[1, 3].transAxes,
        fontsize=20,
        verticalalignment='top',
        horizontalalignment='left'
    )

    fig.text(
        0.935, 0.031,
        info_orbit_text, 
        fontsize=20,
        verticalalignment='bottom',
        horizontalalignment='right'
    )

    fig.text(
        0.74, 0.32,               
        # f'{dwarf_name} stellar mass:\n'
        f'$M_{{stars}}$: {sm_str} $M_{{\odot}}$ \n'  
        f'$M_{{R_{{half}}}}$: {shm_str} $M_{{\odot}}$ \n'
        f'$Gas_{{stars}}$: {cg_str} $M_{{\odot}}$ \n'   
        f'$Gas_{{R_{{half}}}}$: {cghm_str} $M_{{\odot}}$ \n'
        f'Gas Fraction in Stellar Region: {gas_frac:.3f} \n'
        f'Gas Fraction in $R_{{half}}$: {gas_half_frac:.3f} \n'
        f'\n'
        f'$R_{{half}}$: {r_half:.2f} kpc \n'  
        f'$D_{{peri}}$: {d_peri:.2f} kpc \n'
        f'Feedback: {fd_factor_str} median feedback \n',
        fontsize=20,
        ha='left', va='top'
    )

    output_file = os.path.join(output_folder, f'check_nature_{numsp:03d}.png')
    plt.savefig(output_file)
    plt.close(fig)

def main():

    snapshot_files = [f for f in os.listdir(folder_path) if f.startswith('snapshot')]
    max_snapshot_num = len(snapshot_files) - 1

    parser = argparse.ArgumentParser(description='Process some integers.')
    
    snapshot_files = [f for f in os.listdir(folder_path) if f.startswith('snapshot')]
    max_snapshot_num = len(snapshot_files) - 1
    
    default_range = f'0,{max_snapshot_num}'
    
    parser.add_argument('--range', type=str, default=default_range, help='the range of snapshots, e.g., "1,20"')
    parser.add_argument('--processes', type=int, default=1, help='number of processes to use')
    args = parser.parse_args()

    start, end = map(int, args.range.split(','))
    processes = args.processes

    numsp_list = range(start, end + 1)

    with get_context("spawn").Pool(processes=processes) as pool:
        pool.starmap(safe_plot, [(numsp,) for numsp in numsp_list])

if __name__ == "__main__":
    main() 