import numpy as np
import pandas as pd
import zim
import readsnap as giz
from astropy.coordinates import Galactocentric, Galactic
from astropy import units as u
from astropy.constants import G
from sklearn.neighbors import KDTree, NearestNeighbors
# from sklearn.mixture import GaussianMixture
from scipy.stats import binned_statistic, linregress
from scipy.ndimage import gaussian_filter
# from galpy.orbit import Orbit
from scipy.stats import gaussian_kde
from matplotlib.path import Path
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
# import configparser
import re
# import csv
# from scipy.interpolate import UnivariateSpline
import gc
import os
from variable import fornax_core_radius, gizmo_config, orbit_file, folder_path, cd_threshold

core_radius = fornax_core_radius

class Particle:
    """Base class for particle data handling"""
    def __init__(self, x, y, z, vx, vy, vz, mass=None):
        self.x = x
        self.y = y 
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.mass = mass if mass is not None else np.ones_like(x)
        
    @property
    def radius(self):
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)
    
    def to_cartesian(self):
        return np.column_stack([self.x, self.y, self.z])
    
    def to_velocity(self):
        return np.column_stack([self.vx, self.vy, self.vz])
    
    def to_mass(self):
        return self.mass
    
    def select_particles(self, mask):
        """Return a new Particle instance with only the selected particles"""
        return Particle(
            x=self.x[mask],
            y=self.y[mask],
            z=self.z[mask],
            vx=self.vx[mask],
            vy=self.vy[mask],
            vz=self.vz[mask],
            mass=self.mass[mask]
        )

class GalaxySimulation:
    """Main class for galaxy simulation data analysis"""
    def __init__(self, folder_path, snapshot_num, boxsize='Auto'):
        self.folder_path = folder_path
        self.snapshot_num = snapshot_num
        self.boxsize = boxsize
        self.df = None
        self.time = None    

    def load_snapshot(self, dftype=3, box=True):
        """Load snapshot data into DataFrame"""
        try:
            # Read header
            hdr = giz.readsnap(self.folder_path, self.snapshot_num, 0, header_only=1)
            if hdr is None:
                raise ValueError(f"Failed to read snapshot {self.snapshot_num} from {self.folder_path}")
                
            self.time = hdr['time']
            npartTotal = hdr['npartTotal']
            Ntot = np.sum(npartTotal[[0,1,2,3,4]])

            # Define columns based on dftype
            base_columns = ['id', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'm', 'tp']
            gas_columns = ['temp', 'rho', 'nh', 'sfr', 'ugas']
            star_columns = ['birth']
            
            columns = base_columns
            if dftype >= 2:
                columns.extend(gas_columns)
            if dftype >= 3:
                columns.extend(star_columns)

            # Initialize DataFrame with zeros instead of NaN
            df_data = {
                'id': np.zeros(Ntot, dtype=np.int32),
                'x': np.zeros(Ntot, dtype=np.float32),
                'y': np.zeros(Ntot, dtype=np.float32),
                'z': np.zeros(Ntot, dtype=np.float32),
                'vx': np.zeros(Ntot, dtype=np.float32),
                'vy': np.zeros(Ntot, dtype=np.float32),
                'vz': np.zeros(Ntot, dtype=np.float32),
                'm': np.zeros(Ntot, dtype=np.float32),
                'tp': np.zeros(Ntot, dtype=np.int32)
            }
            
            # Add additional columns based on dftype
            if dftype >= 2:
                for col in gas_columns:
                    df_data[col] = np.zeros(Ntot, dtype=np.float32)
            # 仅在dftype==4时添加ags列
            if dftype == 4:
                df_data['ags'] = np.zeros(Ntot, dtype=np.float32)
            if dftype >= 3:
                for col in star_columns:
                    df_data[col] = np.zeros(Ntot, dtype=np.float32)
                    
            df = pd.DataFrame(df_data)

            kids = 0
            for tp in [0,1,2,3,4]:
                # 仅在dftype==4时读取ags
                if tp == 0 and dftype == 4:
                    data = giz.readsnap(self.folder_path, self.snapshot_num, tp, readags=True)
                else:
                    data = giz.readsnap(self.folder_path, self.snapshot_num, tp)
                npart = npartTotal[tp]
                
                if npart == 0:
                    continue
                
                if data is None or len(data['p'][:,0]) != npart:
                    raise RuntimeError(f"Particle count mismatch for type {tp}")
                
                slice_idx = slice(kids, kids + npart)
                
                # Basic particle properties
                df.iloc[slice_idx, df.columns.get_loc('id')] = data['id'].astype(np.int32)
                df.iloc[slice_idx, [df.columns.get_loc(c) for c in ['x','y','z']]] = data['p'].astype(np.float32)
                df.iloc[slice_idx, [df.columns.get_loc(c) for c in ['vx','vy','vz']]] = data['v'].astype(np.float32)
                df.iloc[slice_idx, df.columns.get_loc('m')] = (data['m'] * 1e10).astype(np.float32)  # Convert to solar masses
                df.iloc[slice_idx, df.columns.get_loc('tp')] = tp
                
                # Gas properties (type 0)
                if tp == 0 and dftype >= 2:
                    temp_values = zim.unitconversion(data['u'], data['ne'])
                    df.iloc[slice_idx, df.columns.get_loc('ugas')] = data['u'].astype(np.float32)
                    df.iloc[slice_idx, df.columns.get_loc('temp')] = temp_values.astype(np.float32)
                    df.iloc[slice_idx, df.columns.get_loc('rho')] = (data['rho'] * 404.7).astype(np.float32)
                    df.iloc[slice_idx, df.columns.get_loc('nh')] = data['nh'].astype(np.float32)
                    df.iloc[slice_idx, df.columns.get_loc('sfr')] = data['sfr'].astype(np.float32)
                    # 仅在dftype==4时赋值ags
                    if dftype == 4 and 'ags' in data:
                        df.iloc[slice_idx, df.columns.get_loc('ags')] = data['ags'].astype(np.float32)
                
                # Star properties (type 4)
                if tp == 4 and dftype >= 3:
                    df.iloc[slice_idx, df.columns.get_loc('birth')] = data['age'].astype(np.float32)
                
                kids += npart

            # Auto-calculate boxsize if needed
            if self.boxsize == 'Auto':
                self.boxsize = round(np.ptp(df['x']), -1)
                print(f'Using auto-calculated boxsize: {self.boxsize}')
                
            # Center coordinates and calculate radii
            if box == True:
                df[['x','y','z']] -= self.boxsize / 2
            else:
                print("Box centering skipped. Ensure coordinates are already centered if needed.")
                pass

            # df[['x','y','z']] -= self.boxsize / 2
            
            self.df = df

            return df, self.time
            
        except Exception as e:
            raise RuntimeError(f"Error loading snapshot: {str(e)}")
        
    def find_MW_particles(self, dw_particles_num, cube_file):
        """Find Milky Way particles"""
        if cube_file is None:
            cube_count = 0
        else:
            with open(cube_file, 'r') as f:
                cube_count = sum(1 for _ in f) - 1
            
        df_sorted = self.df.sort_values(by='id')
        front_ids = df_sorted.head(cube_count)['id']
        back_ids = df_sorted.tail(dw_particles_num)['id']
        return ~self.df['id'].isin(front_ids) & ~self.df['id'].isin(back_ids)
    
    def find_dwarf_particles(self, dw_particles_num):
        """Find dwarf galaxy particles"""
        if dw_particles_num <= 0:
            return pd.Series([False] * len(self.df), index=self.df.index)
        dw_all = self.df.nlargest(dw_particles_num, 'id')
        return self.df['id'].isin(dw_all['id'])
 
    def classify_mw_dwarf(self,
                          r_exclude=5,
                          dwarf_radius=5,
                          k_density=32,
                          df=None):
      
        if df is None:
            df = self.df
        star_sel = df.tp >= 2
        df_star = df[star_sel]

        x = df_star["x"].to_numpy()
        y = df_star["y"].to_numpy()
        z = df_star["z"].to_numpy()

        r = np.sqrt(x*x + y*y + z*z)

        density = Analysis.estimate_local_density(x, y, z, k=k_density)

        outer = r > r_exclude

        idx_peak = np.argmax(density * outer)

        xd = x[idx_peak]
        yd = y[idx_peak]
        zd = z[idx_peak]

        dx = x - xd
        dy = y - yd
        dz = z - zd

        r_d = np.sqrt(dx*dx + dy*dy + dz*dz)

        dw_star_mask = r_d < dwarf_radius

        mw_star_mask = ~dw_star_mask

        dw_mask = np.zeros(len(df), dtype=bool)
        mw_mask = np.zeros(len(df), dtype=bool)

        dw_mask[star_sel] = dw_star_mask
        mw_mask[star_sel] = mw_star_mask

        return mw_mask, dw_mask

    def center_on_MW(self, r_init=300.0, r_min=5.0, shrink=0.8):

        df = self.df

        x = df["x"].to_numpy()
        y = df["y"].to_numpy()
        z = df["z"].to_numpy()
        m = df["m"].to_numpy()
        tp = df["tp"].to_numpy()

        # use stars
        sel = (tp == 2) | (tp == 4)

        x = x[sel]
        y = y[sel]
        z = z[sel]
        m = m[sel]

        # initial center
        xc = np.average(x, weights=m)
        yc = np.average(y, weights=m)
        zc = np.average(z, weights=m)

        r_current = r_init

        # shrinking sphere loop
        while r_current > r_min:

            dx = x - xc
            dy = y - yc
            dz = z - zc

            r2 = dx*dx + dy*dy + dz*dz
            inside = r2 < r_current*r_current

            if inside.sum() < 100:
                break

            xc = np.average(x[inside], weights=m[inside])
            yc = np.average(y[inside], weights=m[inside])
            zc = np.average(z[inside], weights=m[inside])

            r_current *= shrink

        # apply final center shift
    
        df['x'] -= xc
        df['y'] -= yc
        df['z'] -= zc

        # compute radius AFTER centering
        df['r'] = np.sqrt(
            df['x']**2 +
            df['y']**2 +
            df['z']**2
        )

        return df

    
    def convert_coordinates(self, df=None):

        df = (self.df if df is None else df).copy()

        # Extract arrays once
        x  = df['x'].to_numpy()
        y  = df['y'].to_numpy()
        z  = df['z'].to_numpy()
        vx = df['vx'].to_numpy()
        vy = df['vy'].to_numpy()
        vz = df['vz'].to_numpy()

        # Heliocentric
        xh, yh, zh, vxh, vyh, vzh = CoordinateTransform.to_heliocentric(
            x, y, z, vx, vy, vz
        )
        df['xh'], df['yh'], df['zh'] = xh, yh, zh
        df['vxh'], df['vyh'], df['vzh'] = vxh, vyh, vzh

        # Distance
        df['rh'] = CoordinateTransform.calculate_distances(xh, yh, zh)

        # Galactic (use heliocentric if physically valid)
        l, b = CoordinateTransform.to_galactic(xh, yh, zh)
        df['l'], df['b'] = l, b

        # Equatorial
        ra, dec = CoordinateTransform.galactic_to_equatorial(l, b)
        df['ra'], df['dec'] = ra, dec

        return df

    def classify_new_dwarf_stars(self, dw_star_mask, r_cut=core_radius):
        try:
            # Get positions of new stars (type 4)
            new_stars_mask = self.df.tp == 4

            if not new_stars_mask.any() or not dw_star_mask.any():
                empty_mask = pd.Series(False, index=self.df.index, dtype=bool)
                return empty_mask, empty_mask

            new_pos = self.df.loc[new_stars_mask, ['x', 'y', 'z']].values
            dw_pos = self.df.loc[dw_star_mask, ['x', 'y', 'z']].values

            dw_tree = KDTree(dw_pos)
            dw_dist, _ = dw_tree.query(new_pos, k=1)
            # dw_dist = np.atleast_1d(dw_dist)
            threshold = r_cut
            is_dw = (dw_dist.squeeze() < threshold)

            new_dw_stars = pd.Series(False, index=self.df.index, dtype=bool)
            new_mw_stars = pd.Series(False, index=self.df.index, dtype=bool)

            new_dw_stars.loc[new_stars_mask] = is_dw.astype(bool)
            new_mw_stars.loc[new_stars_mask] = (~is_dw).astype(bool)

            return new_dw_stars, new_mw_stars

        except Exception as e:
            raise RuntimeError(f"Error classifying new stars: {str(e)}")

    def find_dwarf_gas(self,
                   dw_star_mask,
                   radius=30.0):

        df = self.df

        x = df["x"].to_numpy()
        y = df["y"].to_numpy()
        z = df["z"].to_numpy()
        m = df["m"].to_numpy()

        if dw_star_mask.sum() == 0:
            return np.zeros(len(df), dtype=bool)

        xc = np.average(x[dw_star_mask], weights=m[dw_star_mask])
        yc = np.average(y[dw_star_mask], weights=m[dw_star_mask])
        zc = np.average(z[dw_star_mask], weights=m[dw_star_mask])

        gas_sel = df["tp"].to_numpy() == 0

        dx = x - xc
        dy = y - yc
        dz = z - zc

        r2 = dx*dx + dy*dy + dz*dz

        dw_gas_mask = gas_sel & (r2 < radius**2)

        return dw_gas_mask    

    def find_MW_gas(self,
                    radius=200.0):

        df = self.df

        x = df["x"].to_numpy()
        y = df["y"].to_numpy()
        z = df["z"].to_numpy()
        tp = df["tp"].to_numpy()

        gas_sel = tp == 0
        r2 = x*x + y*y + z*z

        mw_gas_mask = gas_sel & (r2 < radius**2)

        return mw_gas_mask

    @staticmethod
    def get_dw_num(file_path):

        try:
            with open(file_path, 'r') as file:
                # Read first line and sum all numeric values
                first_line = file.readline().strip()
                # Handle both space and comma separated values
                numbers = [float(num) for num in first_line.replace(',', ' ').split() 
                        if num.replace('.', '', 1).isdigit()]
                total_sum = int(sum(numbers))
                
                if total_sum <= 0:
                    raise ValueError("No valid particle count found")
                    
                return total_sum
                
        except FileNotFoundError:
            print(f"Error: File not found - {file_path}")
            return 0
        except (ValueError, TypeError) as e:
            print(f"Error parsing dwarf particle count: {e}")
            return 0
        except Exception as e:
            print(f"Unexpected error reading dwarf config: {e}")
            return 0

class CoordinateTransform:
    """Handle coordinate transformations"""
    @staticmethod
    def to_heliocentric(x, y, z, vx=None, vy=None, vz=None):

        galactocentric_frame = Galactocentric()
        X_sun = galactocentric_frame.galcen_distance.to(u.kpc).value
        Z_sun = galactocentric_frame.z_sun.to(u.kpc).value
        
        xh = x - X_sun
        yh = y
        zh = z - Z_sun
        
        if all(v is not None for v in [vx, vy, vz]):
            vsun = galactocentric_frame.galcen_v_sun.d_xyz.to(u.km/u.s).value
            vxh = vx - vsun[0]
            vyh = vy - vsun[1]
            vzh = vz - vsun[2]
            return xh, yh, zh, vxh, vyh, vzh
            
        return xh, yh, zh
    
    @staticmethod
    def to_galactic(x, y, z, chunk_size=100000):

        l_list = []
        b_list = []
        
        for i in range(0, len(x), chunk_size):
            chunk_x = x[i:i + chunk_size]
            chunk_y = y[i:i + chunk_size]
            chunk_z = z[i:i + chunk_size]
            
            galactic_coords = Galactocentric(
                x=chunk_x * u.kpc,
                y=chunk_y * u.kpc,
                z=chunk_z * u.kpc
            ).transform_to(Galactic())
            
            l_list.append(galactic_coords.l.deg)
            b_list.append(galactic_coords.b.deg)
        
        return np.concatenate(l_list), np.concatenate(b_list)
    
    @staticmethod
    def calculate_distances(x, y, z):

        return np.sqrt(
            x.astype(np.float32)**2 + 
            y.astype(np.float32)**2 + 
            z.astype(np.float32)**2
        )

    @staticmethod  
    def galactic_to_equatorial(l, b, chunk_size=100000):
 
        from astropy.coordinates import SkyCoord
        
        ra_list = []
        dec_list = []
        
        for i in range(0, len(l), chunk_size):
            chunk_l = l[i:i + chunk_size]
            chunk_b = b[i:i + chunk_size]
            
            # Convert to SkyCoord object
            coords = SkyCoord(l=chunk_l*u.degree, 
                            b=chunk_b*u.degree, 
                            frame='galactic')
            
            # Transform to ICRS (equatorial) frame
            eq_coords = coords.transform_to('icrs')
            
            ra_list.append(eq_coords.ra.deg)
            dec_list.append(eq_coords.dec.deg)
        
        return np.concatenate(ra_list), np.concatenate(dec_list)

class Analysis:
    """Analysis tools for galaxy data"""
    @staticmethod
    def cal_vlos(x, y, z, vx, vy, vz):

        r_vec = np.stack((x, y, z), axis=1)
        v_vec = np.stack((vx, vy, vz), axis=1)

        r_norm = np.linalg.norm(r_vec, axis=1)
        epsilon = 1e-8 
        r_norm_safe = np.where(r_norm < epsilon, epsilon, r_norm)
        r_hat = r_vec / r_norm_safe[:, np.newaxis]
        vlos = np.sum(v_vec * r_hat, axis=1)

        return vlos

    @staticmethod
    def mass_weighted_mean_and_dispersion(values, mass):
        """
        Compute the mass-weighted mean velocity and dispersion.

        This is the global second-moment definition used for both the full
        stellar sample and the projected sample inside R_e. We do not average
        dispersions from radial bins, because that does not preserve the total
        second moment.
        """
        values = np.asarray(values, dtype=float)
        mass = np.asarray(mass, dtype=float)

        if values.size == 0 or mass.size == 0:
            return np.nan, np.nan

        total_mass = np.sum(mass)
        if total_mass <= 0:
            return np.nan, np.nan

        v_mean = np.sum(mass * values) / total_mass
        sigma2 = np.sum(mass * (values - v_mean)**2) / total_mass

        return v_mean, np.sqrt(np.maximum(sigma2, 0.0))

    @staticmethod
    def calculate_half_mass_radius(x, y, z, m, center_x=None, center_y=None, center_z=None):
        """
        Calculate half-mass radius from particle positions and masses.
        
        Parameters:
        -----------
        x, y, z : array-like
            3D positions of particles
        m : array-like
            Masses of particles
        center_x, center_y, center_z : float, optional
            Center coordinates for the calculation. If provided, measurements are made
            relative to this center. If not provided, assumes particles are already 
            centered at the origin.
            
        Returns:
        --------
        float
            Half-mass radius in the same units as input coordinates
        """
        x = np.asarray(x)
        y = np.asarray(y)
        z = np.asarray(z)
        m = np.asarray(m)
        
        # Apply center offset if provided
        if center_x is not None:
            x = x - center_x
        if center_y is not None:
            y = y - center_y
        if center_z is not None:
            z = z - center_z
        
        r = np.sqrt(x**2 + y**2 + z**2)
        idx = np.argsort(r)
        r_sorted = r[idx]
        m_sorted = m[idx]
        m_cum = np.cumsum(m_sorted)
        half_mass = m_cum[-1] / 2
        idx_half = np.searchsorted(m_cum, half_mass)
        return r_sorted[idx_half] if idx_half < len(r_sorted) else r_sorted[-1]
    
    @staticmethod
    def find_center(x, y, units='degree', shrink_rate=0.95, threshold=0.1):

        x_arr = np.array(x)
        y_arr = np.array(y)

        if units == 'degree':
            hist, xedges, yedges = np.histogram2d(x_arr, y_arr, bins=[360, 180], range=[[0, 360], [-90, 90]])
        elif units == 'kpc':
            hist, xedges, yedges = np.histogram2d(x_arr, y_arr, bins=[1000, 1000], range=[[-1000, 1000], [-1000, 1000]])
        else:
            raise ValueError("units must be 'degree' or 'kpc'")

        i, j = np.unravel_index(np.argmax(hist), hist.shape)
        center_x = np.mean(xedges[i:i+2])
        center_y = np.mean(yedges[j:j+2])

        range_x, range_y = core_radius / 3, core_radius / 3

        while range_x > threshold and range_y > threshold:
            mask = ((x_arr >= center_x - range_x) & (x_arr <= center_x + range_x) &
                    (y_arr >= center_y - range_y) & (y_arr <= center_y + range_y))
            filtered_l_arr = x_arr[mask]
            filtered_b_arr = y_arr[mask]

            if filtered_l_arr.size == 0 or filtered_b_arr.size == 0:
                break

            center_x = np.median(filtered_l_arr)
            center_y = np.median(filtered_b_arr)

            range_x *= shrink_rate
            range_y *= shrink_rate

        return center_x, center_y
    
    @staticmethod
    def compute_sigma_mw(m, r, rhalf, f=0.25):

        # 引力常数 G (kpc * (km/s)^2 / Msun)
        G = 4.302e-6
        # g_MW = GM / r^2
        g_mw = G * m / (r**2)
        # Δσ² = 2√2 g rhalf f
        delta_sigma2 = 2 * np.sqrt(2) * g_mw * rhalf * f
        sigma_mw = np.sqrt(delta_sigma2)

        return sigma_mw
    
    @staticmethod
    def estimate_local_density(x, y, z, k=32):

        pos = np.vstack((x, y, z)).T
        tree = cKDTree(pos)

        dist, _ = tree.query(pos, k=k)

        # 第 k 近邻距离
        rk = dist[:, -1]

        volume = (4/3) * np.pi * rk**3

        density = k / volume

        return density
    
    @staticmethod
    def rotate_to_sky(l, b, center=None):
        l_arr = np.array(l)
        b_arr = np.array(b)

        if center == None:
            center_l, center_b = Analysis.find_center(l, b)
        elif len(center) == 2:
            center_l, center_b = center[0], center[1]
        else:
            print('Please input correct center coordinate')

        l_rad_arr, b_rad_arr = np.radians(l_arr), np.radians(b_arr)
        center_rad_l, center_rad_b = np.radians(center_l), np.radians(center_b)

        deltal = np.cos(b_rad_arr) * np.sin(l_rad_arr - center_rad_l)
        deltab = np.sin(b_rad_arr) * np.cos(center_rad_b) - np.cos(b_rad_arr) * np.sin(center_rad_b) * np.cos(l_rad_arr-center_rad_l)

        dl_deg = np.rad2deg(deltal)
        db_deg = np.rad2deg(deltab)

        return dl_deg, db_deg
    
    @staticmethod
    def L_to_ABmag(mass, d_kpc, mag_sys='vega', m_l_relation=2.5):

        L_d = mass/m_l_relation
        m_AB = -2.5 * np.log10(L_d) + 5 * np.log10(d_kpc*1000) + 4.83 - 5
        m_vega = m_AB - 0.02 # V band

        if mag_sys == 'AB':
            return m_AB
        elif mag_sys == 'vega':
            return m_vega

    @staticmethod
    def surface_brightness_map(x_kpc, y_kpc, mass, d_mean,
                            npix=100, size_kpc=10.8, mag_sys='vega', m_l_relation=2.5):

        H, xedges, yedges = np.histogram2d(
            x_kpc, y_kpc,
            bins=npix,
            range=[[-size_kpc/2, size_kpc/2], [-size_kpc/2, size_kpc/2]],
            weights=mass
        )

        img_mag = np.full_like(H, np.inf, dtype=float)

        pix_size_kpc = size_kpc / npix
        pix_size_rad = (pix_size_kpc / d_mean)  
        pix_area_arcsec2 = (pix_size_rad * 206265)**2

        for i in range(npix):
            for j in range(npix):
                m_pix = H[i, j]
                if m_pix > 0:
                    m_tot = Analysis.L_to_ABmag(m_pix, d_mean, mag_sys=mag_sys, m_l_relation=m_l_relation)
                    mu = m_tot + 2.5*np.log10(pix_area_arcsec2)
                    img_mag[i, j] = mu

        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        return img_mag, extent
    
    @staticmethod
    def calculate_gas_loss_rate(age, mass, window=None):

        age = np.array(age)
        mass = np.array(mass)
        df = pd.DataFrame({'age': age, 'mass': mass})
        if window is None:
            window = max(3, len(df) // 6)
        
        def fit_slope(x):
            if x.isnull().any():
                return np.nan
            return np.polyfit(range(len(x)), x, 1)[0]
        gas_loss_rate_smooth = (
            df['mass']
            .rolling(window, center=True)
            .apply(fit_slope, raw=False)
        )
       
        dt_mean = df['age'].diff().mean()
        gas_loss_rate_smooth = -gas_loss_rate_smooth / dt_mean
        return -gas_loss_rate_smooth.values

    @staticmethod
    def calculate_ellipticity(x, y, mass, n_neighbors=20):
       
        coords = np.vstack([x, y]).T

        nn = NearestNeighbors(n_neighbors=n_neighbors).fit(coords)
        distances, _ = nn.kneighbors(coords)
        local_density = 1 / (np.mean(distances, axis=1) + 1e-5) 

        weights = mass * local_density
        weights /= np.sum(weights)  
        x_c = np.sum(x * weights)
        y_c = np.sum(y * weights)

        x_rel = x - x_c
        y_rel = y - y_c

        Ixx = np.sum(x_rel**2 * weights)
        Iyy = np.sum(y_rel**2 * weights)
        Ixy = np.sum(x_rel * y_rel * weights)

        tmp = np.sqrt((Ixx - Iyy)**2 + 4 * Ixy**2)
        a2 = (Ixx + Iyy + tmp) / 2
        b2 = (Ixx + Iyy - tmp) / 2

        eps = 1 - np.sqrt(b2 / a2)
        pa = 0.5 * np.arctan2(2 * Ixy, Ixx - Iyy)

        return eps, pa
        
    @staticmethod
    def radial_magnitude_profile(x_rot, y_rot, mass, d_kpc, ep, r_cut=core_radius, bins=20, mag_sys='vega', m_l_relation=2.5):

        q = 1-ep
        r_elliptical = np.sqrt(x_rot**2 + (y_rot / q)**2)
        factor = (206265 / d_kpc)**2

        r_edges = np.linspace(0, r_cut, bins + 1)
        r_centers = 0.5 * (r_edges[1:] + r_edges[:-1])

        mag_profile = []

        for r1, r2 in zip(r_edges[:-1], r_edges[1:]):
            mask = (r_elliptical >= r1) & (r_elliptical < r2)
            if np.any(mask):
                m_bin = np.sum(mass[mask])
                mag_r = Analysis.L_to_ABmag(mass=m_bin, d_kpc=d_kpc, mag_sys=mag_sys, m_l_relation=m_l_relation)
    
                area_arcsec2 = np.pi * (r2**2 - r1**2) * q * factor
                mag = mag_r + 2.5 * np.log10(area_arcsec2)
                mag_profile.append(mag) 
            else:
                mag = np.nan
                mag_profile.append(mag)

        return r_centers, mag_profile

    @staticmethod
    def half_light_radius(x, y, mass, ep):

        mass = np.asarray(mass)
        q = 1 - ep  # b/a
        r = np.sqrt(x**2 + (y / q)**2)

        idx = np.argsort(r)
        r_sorted = r[idx]
        m_sorted = mass[idx]

        cum_mass = np.cumsum(m_sorted)
        half_mass = 0.5 * cum_mass[-1]

        Re = np.interp(half_mass, cum_mass, r_sorted)

        return Re
    
    @staticmethod
    def get_vlos_dispersion(x, y, vlos, r_range=core_radius, bins=50):

        r_proj = np.sqrt(x ** 2 + y ** 2)
        mask_r = r_proj < r_range

        r_bins = np.linspace(0, r_proj[mask_r].max(), bins)
        bin_centers = 0.5 * (r_bins[:-1] + r_bins[1:])

        vlos_dispersion, _, _ = binned_statistic(r_proj[mask_r], vlos[mask_r], statistic='std', bins=r_bins)

        return bin_centers, vlos_dispersion

    @staticmethod
    def est_gas_contour(x, y, m, 
                        box=core_radius * 0.8, 
                        nbins=50, 
                        cd_threshold=cd_threshold,
                        smooth_sigma=2.0,
                        unit="kpc",
                        d_mean=None):

        x = np.asarray(x)
        y = np.asarray(y)
        m = np.asarray(m)

        mask = (np.abs(x) < box) & (np.abs(y) < box)
        x_cut = x[mask]
        y_cut = y[mask]
        m_cut = m[mask]

        # 2D histogram
        H, xedges, yedges = np.histogram2d(x_cut, y_cut, bins=nbins, weights=m_cut)

        xc = 0.5 * (xedges[:-1] + xedges[1:])
        yc = 0.5 * (yedges[:-1] + yedges[1:])
        X, Y = np.meshgrid(xc, yc)

        dx = xedges[1] - xedges[0]
        dy = yedges[1] - yedges[0]
        area = dx * dy  # kpc^2

        cd = H / area * 1.25e14  # M_sun/kpc^2 → atoms/cm^2
        cd_smooth = gaussian_filter(cd, sigma=smooth_sigma)

        mass_high = H[cd_smooth > cd_threshold].sum()

        # --- New: coordinate unit conversion using d_mean ---
        if unit == "deg":
            if d_mean is None:
                raise ValueError("unit='deg' requires d_mean (distance in kpc).")
            kpc_to_deg = (180.0 / np.pi) / d_mean
            X = X * kpc_to_deg
            Y = Y * kpc_to_deg

        return mass_high, X, Y, cd_smooth
    
    @staticmethod
    def get_vlos_field(x_kpc, y_kpc, vlos_array, size_kpc=10.8, npix=100):

        vlos_residual = vlos_array - np.mean(vlos_array)
        
        H, xedges, yedges = np.histogram2d(
            x_kpc, y_kpc,
            bins=npix,
            range=[[-size_kpc/2, size_kpc/2], [-size_kpc/2, size_kpc/2]],
            weights=vlos_residual
        )

        N, _, _ = np.histogram2d(
            x_kpc, y_kpc,
            bins=npix,
            range=[[-size_kpc/2, size_kpc/2], [-size_kpc/2, size_kpc/2]]
        )

        with np.errstate(invalid='ignore', divide='ignore'):
            vlos_mean = np.where(N > 0, H / N, np.nan)

        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

        return vlos_mean, extent

    @staticmethod
    def get_radial_mask(x, y, radius):
        r = np.sqrt(x**2 + y**2)
        return r < radius
    
    @staticmethod
    def get_3d_radial_mask(x, y, z, radius):
        r = np.sqrt(x**2 + y**2 + z**2)
        return r < radius

    @staticmethod
    def estimate_galaxy_boundary(x, y, grid_resolution=15, density_percentile=0.998):

        kde = gaussian_kde(np.vstack([x, y]))
        xmin, xmax = x.min(), x.max()
        ymin, ymax = y.min(), y.max()

        xx, yy = np.mgrid[xmin:xmax:complex(grid_resolution), ymin:ymax:complex(grid_resolution)]
        zz = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)

        zz_flat = zz.flatten()
        sorted_z = np.sort(zz_flat)
        cdf = np.cumsum(sorted_z) / np.sum(sorted_z)
        level = sorted_z[np.searchsorted(cdf, 1 - density_percentile)]

        fig = plt.figure()
        ax = fig.add_subplot()
        cs = ax.contour(xx, yy, zz, levels=[level])
        plt.close(fig) 

        valid_contours = [c for c in cs.allsegs[0] if c.ndim == 2 and c.shape[1] == 2 and len(c) >= 3]

        if not valid_contours:
            raise ValueError("No valid contours found for galaxy boundary.")

        largest_contour = max(
            valid_contours, 
            key=lambda pts: 0.5 * abs(
                np.dot(pts[:, 0], np.roll(pts[:, 1], 1)) - np.dot(pts[:, 1], np.roll(pts[:, 0], 1))
                )
        )

        boundary_path = Path(largest_contour)

        return boundary_path

    @staticmethod
    # @jit(nopython=True)
    def cartesian_to_spherical(x, y, z, vx, vy, vz):

        r = np.sqrt(x**2 + y**2 + z**2)
        theta = np.arccos(z / r)
        phi = np.arctan2(y, x)
        
        sinθ = np.sin(theta)
        cosθ = np.cos(theta)
        sinφ = np.sin(phi)
        cosφ = np.cos(phi)

        e_r = np.stack([sinθ * cosφ,
                        sinθ * sinφ,
                        cosθ], axis=1)

        e_theta = np.stack([cosθ * cosφ,
                            cosθ * sinφ,
                            -sinθ], axis=1)

        e_phi = np.stack([-sinφ,
                        cosφ,
                        np.zeros_like(x)], axis=1)

        v_vec = np.stack([vx, vy, vz], axis=1)

        v_r = np.sum(v_vec * e_r, axis=1)
        v_theta = np.sum(v_vec * e_theta, axis=1)
        v_phi = np.sum(v_vec * e_phi, axis=1)

        return v_r, v_theta, v_phi

    @staticmethod
    def fit_velocity_gradient_binned(x, y, n_bins=30, method='uniform'):
   
        df_bin = pd.DataFrame({'x': x, 'y': y})

        if method == 'quantile':
            df_bin['bin'] = pd.qcut(df_bin['x'], q=n_bins)
        else:
            df_bin['bin'] = pd.cut(df_bin['x'], bins=n_bins)
        
        bin_stats = df_bin.groupby('bin').agg({'x': 'mean', 'y': 'mean'}).dropna()
        slope, intercept, r_value, p_value, std_err = linregress(bin_stats['x'], bin_stats['y'])
        
        return slope, intercept, bin_stats['x'].values, bin_stats['y'].values

    @staticmethod
    def cal_tidal_eff(dw_elinfo):

        d_array = dw_elinfo['distance_gal'].values * u.kpc
        v_array = u.Quantity(np.sqrt(dw_elinfo['vr']**2 + dw_elinfo['vtheta']**2 + dw_elinfo['vphi']**2), u.km / u.s)
        Mw_array = dw_elinfo['mw_mass_r'].values * u.Msun
        r_half_km2 = (dw_elinfo['rhalf'].values * u.kpc).to(u.km).value ** 2
        sigma_array = dw_elinfo['sigma'].values

        expr = (G * Mw_array) / (d_array**2 * v_array)
        tidal_eff = u.Quantity((128.0 / 27.0) * expr**2, 1 / u.s**2).value

        dw_energy = sigma_array**2 / r_half_km2

        g_mw = (G * Mw_array / d_array**2).to(u.km / u.s**2)
        sigma_mw = np.sqrt((np.sqrt(2) * 1 * g_mw * dw_elinfo['rhalf'].values * u.km)).value

        return tidal_eff, dw_energy, sigma_mw

    @staticmethod
    def GetSFH_sim(folder_path=folder_path, Tage=None):

        if Tage is None:
            print('Please input the age of the universe at the snapshot (Tage) in Gyr')
        
        numsp = DataProcessor.find_snapshot_range(folder_path)
        print(f'Choose snapshot{numsp} for SFH calculation')
        simulation = GalaxySimulation(folder_path=folder_path, snapshot_num=numsp)
        df, tsnap = simulation.load_snapshot()
        df = simulation.center_on_MW()
        _, dw_star_mask = simulation.classify_mw_dwarf(df=df, r_exclude=5, dwarf_radius=3*core_radius, k_density=16)

        df_dw = df.loc[dw_star_mask, ['birth', 'm']].copy()
        df_dw = df_dw[df_dw['birth'] != 0]
        df_dw['age'] = Tage - df_dw['birth']
        dw_sfrmass = df_dw[['age', 'm']].astype(float)

        return dw_sfrmass

    @staticmethod
    def get_sfh(folder_path=folder_path, range=core_radius):

        numsp = DataProcessor.find_snapshot_range(folder_path)
        print(f'Choose snapshot{numsp}')
        simulation = GalaxySimulation(folder_path=folder_path, snapshot_num=numsp)
        df, tsnap = simulation.load_snapshot()

        _, dw_star_mask = simulation.classify_mw_dwarf(r_exclude=5, dwarf_radius=core_radius, k_density=32)

        x_arr = df.loc[dw_star_mask, 'x'].to_numpy()
        y_arr = df.loc[dw_star_mask, 'y'].to_numpy()
        xc, yc = Analysis.find_center(x_arr, y_arr, units='kpc')

        dx = df.x.to_numpy() - xc
        dy = df.y.to_numpy() - yc
        r = np.sqrt(dx**2 + dy**2)

        radial_mask = r < range
        tp4_mask = df.tp.to_numpy() == 4
        final_mask = radial_mask & tp4_mask

        sfrmass = df[final_mask]

        if sfrmass.empty:
            print(f"snapshot{numsp} no new formation stars")
            return None

        dw_sfrmass = (
            sfrmass[['birth', 'm']]
            .sort_values(by='birth')
            .assign(cumulative_mass=lambda df: df['m'].cumsum())
        )

        last_row = dw_sfrmass.iloc[-1].copy()
        last_row['birth'] = tsnap
        dw_sfrmass = pd.concat([dw_sfrmass, pd.DataFrame([last_row])], ignore_index=True)

        return dw_sfrmass

    @staticmethod
    def GetGasdensity(x, y, m, r_half):

        R = np.sqrt(x**2 + y**2)          # kpc
        mask = R <= r_half

        M_tot = m[mask].sum()             # Msun
        area_kpc2 = np.pi * r_half**2     # kpc^2

        Sigma_kpc2 = M_tot / area_kpc2    # Msun/kpc^2

        # 转成 Msun/pc^2
        Sigma_pc2 = Sigma_kpc2 / 1e6

        # 再转成 column density
        N_H = Sigma_pc2 * 1.25e20         # atoms/cm^2

        return N_H

    @staticmethod
    def particles_to_image(x, y, m, npix=200, size_kpc=core_radius):
        """
        Convert particle positions to a 2D image for photutils.
        """
        xmin, xmax = -size_kpc, size_kpc
        ymin, ymax = -size_kpc, size_kpc

        image, _, _ = np.histogram2d(
            y, x,  
            bins=npix,
            range=[[ymin, ymax], [xmin, xmax]],
            weights=m
        )
        return image

class DataProcessor:

    @staticmethod
    def read_feedback_value(file_path=gizmo_config):
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith('%'):
                    continue
                if line.startswith("timescale_fd"):
                    parts = line.split()
                    if len(parts) > 1:
                        return float(parts[1].split(";")[0])
        return None
    
    @staticmethod
    def read_sfreff_value(file_path=gizmo_config):
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith('%'):
                    continue
                if line.startswith("sfr_effiency"):
                    parts = line.split()
                    if len(parts) > 1:
                        return float(parts[1].split(";")[0])
        return None

    @staticmethod
    def get_ini_data(elinfo_path):
        with open(elinfo_path, 'r') as f:

            line1 = f.readline().strip().removeprefix("# Orbit Data:").strip()
            x, y, z, vx, vy, vz = map(float, line1.split(','))

            line2 = f.readline().strip().removeprefix("# Timescale_fd:").strip()
            fd_value = float(line2)

        return x, y, z, vx, vy, vz, fd_value

    @staticmethod
    def read_orbit_parameters(orbit_file=orbit_file):
        if not os.path.exists(orbit_file):
            print(f"Orbit file not found: {orbit_file}")
            return None

        with open(orbit_file, 'r') as f:
            f.readline()
            second_line = f.readline().strip()
            columns = second_line.split()
            if len(columns) < 6:
                print("Orbit file does not have enough columns.")
                return None
            return [float(val) for val in columns[:6]]
            
    @staticmethod
    def calsfh(sfr_csv, Tage):

        sfr_idw = sfr_csv.copy()
        sfr_idw['sfr'] = sfr_idw['sfr'] * 1e-4
        sfr_idw['sfr0'] = sfr_idw['sfr'] / sfr_idw['sfr'].sum()

        sfr_idw['dage'] = sfr_idw['age'].shift(-1) - sfr_idw['age']
        sfr_idw.iloc[-1, sfr_idw.columns.get_loc('dage')] = sfr_idw.iloc[-2]['dage']
        sfr_idw['mt'] = sfr_idw['sfr'] * sfr_idw['dage']

        max_age = sfr_idw['age'].max()
        sfr_idw['cosmic_time'] = max_age - sfr_idw['age']
        sfr_idw['cosmic_time_shifted'] = sfr_idw['cosmic_time'] - sfr_idw['cosmic_time'].max() + Tage

        sfr_idw.sort_values('cosmic_time_shifted', inplace=True)
        sfr_idw['cumulative_mass'] = sfr_idw['mt'].cumsum()

        mask_positive = sfr_idw['cosmic_time_shifted'] > 0
        if mask_positive.any():
            min_positive_time = sfr_idw.loc[mask_positive, 'cosmic_time_shifted'].min()
            mt_at_min = sfr_idw.loc[sfr_idw['cosmic_time_shifted'] == min_positive_time, 'cumulative_mass'].values[0]
        else:
            mt_at_min = 0.0  # fallback if no positive time exists

        sfr_idw['pure_sfh'] = sfr_idw['cumulative_mass'] - mt_at_min

        return sfr_idw

    @staticmethod
    def GetSFR(sfr_csv, bins=None):

        sfr = sfr_csv.copy()

        sfr['dage'] = sfr['age'].shift(-1) - sfr['age']
        sfr.iloc[-1, sfr.columns.get_loc('dage')] = sfr.iloc[-2]['dage']

        lookback = sfr['age'].to_numpy()
        sfr_values = sfr['sfr'].to_numpy()

        if bins is None:
            bin_width = 0.5
            bins = np.arange(0 - bin_width, np.max(lookback), bin_width)

        sfr_in_bin, edges = np.histogram(
            lookback,
            bins=bins,
            weights=sfr_values
        )

        count_in_bin, _ = np.histogram(lookback, bins=bins)
        avg_sfr = np.divide(
            sfr_in_bin,
            count_in_bin,
            out=np.zeros_like(sfr_in_bin, dtype=float),
            where=count_in_bin > 0
        )

        return edges[:-1], avg_sfr

    @staticmethod
    def find_snapshot_range(folder) :
        files = os.listdir(folder)
        numbers = []
        for f in files:
            match = re.match(r"snapshot_(\d+)\.hdf5", f)
            if match:
                numbers.append(int(match.group(1)))
        if not numbers:
            return None, None
        return max(numbers)

    @staticmethod
    def sci_notation(value, precision=2):
        if value == 0:
            return f"{0:.{precision}f} × 10⁰"
        exponent = int(np.floor(np.log10(abs(value))))
        base = value / (10 ** exponent)

        superscript_map = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
        exponent_str = str(exponent).translate(superscript_map)
        
        return f"{base:.{precision}f} × 10{exponent_str}"
