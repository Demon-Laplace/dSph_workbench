from basefunc import *
from variable import *
from pathlib import Path
import matplotlib.pyplot as plt
import warnings
from matplotlib.widgets import Slider
import plotly.graph_objs as go

warnings.filterwarnings("ignore")
# 用户输入
num_str = input('input num_str: ')
numsp = int(input('input numsp: '))

dwarf_name = Path().resolve().name
folder_path = f"./{dwarf_name}{num_str}/"

inifile = folder_path + f'{dwarf_name}.ini'
dw_particles = GalaxySimulation.get_dw_num(inifile)

core_radius = core_radius_dict.get(dwarf_name, None)
if core_radius is None:
    raise ValueError(f"please input '{dwarf_name}''s radius")

simulation = GalaxySimulation(folder_path=folder_path, snapshot_num=numsp)
df, tsnap = simulation.load_snapshot()
mw_mask = simulation.find_MW_particles(dw_particles_num=dw_particles, cube_file=folder_path+cube_path)
dw_mask = simulation.find_dwarf_particles(dw_particles_num=dw_particles)
df = simulation.center_on_MW(mw_mask, iterations=2, radius=100)
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
mw_center_ra, mw_center_dec = CoordinateTransform.galactic_to_equatorial([0], [0])
mw_rotra_center, mw_rotdec_center = Analysis.rotate_to_sky(mw_center_ra, mw_center_dec, center=[dw_cra, dw_cdec])

rotra_dw_star, rotdec_dw_star = rotra_all[total_dw_star_mask], rotdec_all[total_dw_star_mask]
rotra_dw_hot_gas, rotdec_dw_hot_gas = rotra_all[dw_hot_gas_mask], rotdec_all[dw_hot_gas_mask]
rotra_dw_cold_gas, rotdec_dw_cold_gas = rotra_all[dw_cold_gas_mask], rotdec_all[dw_cold_gas_mask]

d_mean = np.mean(df[total_dw_star_mask].rh)
d_mean_gal = np.mean(df[total_dw_star_mask].r)
# x_kpc = np.radians(rotra_dw_star) * d_mean
# y_kpc = np.radians(rotdec_dw_star) * d_mean
# hot_gas_x_kpc = np.radians(rotra_dw_hot_gas) * d_mean
# hot_gas_y_kpc = np.radians(rotdec_dw_hot_gas) * d_mean
# cold_gas_x_kpc = np.radians(rotra_dw_cold_gas) * d_mean
# cold_gas_y_kpc = np.radians(rotdec_dw_cold_gas) * d_mean

img, extent = Analysis.surface_brightness_map(rotra_dw_star, rotdec_dw_star, df[total_dw_star_mask].m, d_mean,
                                              npix=130, size_kpc=4, mag_sys='vega')

df_sub = df[total_dw_star_mask]
# pos = CartesianRepresentation(df_sub['x'].values * u.kpc,
#                               df_sub['y'].values * u.kpc,
#                               df_sub['z'].values * u.kpc)

# vel = CartesianDifferential(df_sub['vx'].values * u.km/u.s,
#                             df_sub['vy'].values * u.km/u.s,
#                             df_sub['vz'].values * u.km/u.s)

# coord = SkyCoord(pos.with_differentials(vel), frame='galactocentric')
# coord_icrs = coord.transform_to('icrs')

# mean_pmra = np.mean(coord_icrs.pm_ra_cosdec).value
# mean_pmdec = np.mean(coord_icrs.pm_dec).value

# ellipticity, pa = Analysis.calculate_ellipticity(
#     x_kpc, y_kpc, mass=df[total_dw_star_mask].m.values, n_neighbors=30
# )

# Rotate = np.array([[np.cos(-pa), -np.sin(-pa)],
#                    [np.sin(-pa),  np.cos(-pa)]])

# coords = np.vstack((x_kpc, y_kpc))  # shape: (2, N)
# rotated_coords = Rotate @ coords        # shape: (2, N)

# x_rot, y_rot = rotated_coords[0], rotated_coords[1]

# r_centers, mag_profile = Analysis.radial_magnitude_profile(x_rot, y_rot, mass=df[total_dw_star_mask].m,
#                                                            d_kpc=d_mean, ep=ellipticity, bins=80)
# r_half = Analysis.half_light_radius(x_rot, y_rot, df[total_dw_star_mask].m, ep=ellipticity)
# vlos_array = Analysis.cal_vlos(df[total_dw_star_mask].xh, df[total_dw_star_mask].yh, df[total_dw_star_mask].zh,
#                                    df[total_dw_star_mask].vxh, df[total_dw_star_mask].vyh, df[total_dw_star_mask].vzh)
# bin_centers, vlos_dispersion = Analysis.get_vlos_dispersion(x_kpc, y_kpc, vlos_array)

# img, extent = Analysis.surface_brightness_map(rotra_dw_star, rotdec_dw_star, df[total_dw_star_mask].m, d_mean,
#                                               npix=120, size_kpc=4, mag_sys='vega')
# vf, extent_vf = Analysis.get_vlos_field(x_kpc, y_kpc, vlos_array, npix=40)

df_sub = df[total_dw_star_mask]
xh, yh, zh = df_sub['xh'], df_sub['yh'], df_sub['zh']

vlos_array = Analysis.cal_vlos(df[total_dw_star_mask].xh, df[total_dw_star_mask].yh, df[total_dw_star_mask].zh,
                                   df[total_dw_star_mask].vxh, df[total_dw_star_mask].vyh, df[total_dw_star_mask].vzh)

mean_xh, mean_yh, mean_zh = df_sub['xh'].mean(), df_sub['yh'].mean(), df_sub['zh'].mean()
mean_vxh, mean_vyh, mean_vzh = df_sub['vxh'].mean(), df_sub['vyh'].mean(), df_sub['vzh'].mean()
# 箭头长度缩放
arrow_scale = 1e4  # 调节线段长度，单位kpc
km_s_to_kpc_yr = 1.0227e-6
arrow_dx = mean_vxh * km_s_to_kpc_yr * arrow_scale
arrow_dy = mean_vyh * km_s_to_kpc_yr * arrow_scale
arrow_dz = mean_vzh * km_s_to_kpc_yr * arrow_scale
# 线段起点为质心，终点为速度方向
arrow_x = [mean_xh, mean_xh + arrow_dx]
arrow_y = [mean_yh, mean_yh + arrow_dy]
arrow_z = [mean_zh, mean_zh + arrow_dz]
# 指向太阳方向的短线段（单位向量，长度为sun_dir_scale）
sun_dir_scale = 5  # kpc, 可调节
r = np.sqrt(mean_xh**2 + mean_yh**2 + mean_zh**2)
if r > 0:
    sun_dir_x = [mean_xh, mean_xh - mean_xh/r * sun_dir_scale]
    sun_dir_y = [mean_yh, mean_yh - mean_yh/r * sun_dir_scale]
    sun_dir_z = [mean_zh, mean_zh - mean_zh/r * sun_dir_scale]
else:
    sun_dir_x = [mean_xh, mean_xh]
    sun_dir_y = [mean_yh, mean_yh]
    sun_dir_z = [mean_zh, mean_zh]

fig = plt.figure(figsize=(8, 7))
ax = fig.add_subplot(111, projection='3d')
sc = ax.scatter(xh[::10], yh[::10], zh[::10], s=2, alpha=0.7, c=zh[::10], cmap='viridis')
# 平均速度矢量
ax.plot([mean_xh, mean_xh + arrow_dx], [mean_yh, mean_yh + arrow_dy], [mean_zh, mean_zh + arrow_dz], 
        color='r', linewidth=3, label='Mean 3D Velocity')
# 速度矢量两端圆点
ax.scatter([mean_xh, mean_xh + arrow_dx], [mean_yh, mean_yh + arrow_dy], [mean_zh, mean_zh + arrow_dz], 
           color='r', s=60, zorder=10)

# 太阳方向线
ax.plot([mean_xh, sun_dir_x[1]], [mean_yh, sun_dir_y[1]], [mean_zh, sun_dir_z[1]], 
        color='orange', linewidth=2, linestyle='dashed', label='To Sun')
# 太阳方向两端圆点
ax.scatter([mean_xh, sun_dir_x[1]], [mean_yh, sun_dir_y[1]], [mean_zh, sun_dir_z[1]], 
           color='orange', s=60, zorder=10)

ax.set_xlabel('Xh (kpc)')
ax.set_ylabel('Yh (kpc)')
ax.set_zlabel('Zh (kpc)')
ax.set_title('Dwarf Galaxy 3D (Heliocentric)')
ax.legend()

axcolor = 'lightgoldenrodyellow'
ax_elev = plt.axes([0.15, 0.02, 0.65, 0.03], facecolor=axcolor)
ax_azim = plt.axes([0.15, 0.06, 0.65, 0.03], facecolor=axcolor)
slider_elev = Slider(ax_elev, 'Elev', 0, 90, valinit=30)
slider_azim = Slider(ax_azim, 'Azim', 0, 360, valinit=45)

def update(val):
    ax.view_init(elev=slider_elev.val, azim=slider_azim.val)
    fig.canvas.draw_idle()

slider_elev.on_changed(update)
slider_azim.on_changed(update)

plt.show()
