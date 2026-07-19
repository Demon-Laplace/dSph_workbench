import numpy as np
import pandas as pd
from astropy import units as u
from astropy.constants import G
from matplotlib.path import Path
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree
from scipy.stats import binned_statistic, gaussian_kde, linregress
from sklearn.neighbors import NearestNeighbors

from variable import cd_threshold, fornax_core_radius, sim_plot_mass_to_light_ratio


core_radius = fornax_core_radius


class Analysis:
    """Analysis tools for galaxy data."""

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
        values = np.asarray(values, dtype=float)
        mass = np.asarray(mass, dtype=float)

        if values.size == 0 or mass.size == 0:
            return np.nan, np.nan

        total_mass = np.sum(mass)
        if total_mass <= 0:
            return np.nan, np.nan

        v_mean = np.sum(mass * values) / total_mass
        sigma2 = np.sum(mass * (values - v_mean) ** 2) / total_mass

        return v_mean, np.sqrt(np.maximum(sigma2, 0.0))

    @staticmethod
    def calculate_half_mass_radius(*args, **kwargs):
        if len(args) == 1 and hasattr(args[0], "radius") and hasattr(args[0], "mass"):
            particle = args[0]
            r = np.asarray(particle.radius, dtype=float)
            m_sorted = np.asarray(particle.mass, dtype=float)
        elif len(args) >= 4:
            x = np.asarray(args[0], dtype=float)
            y = np.asarray(args[1], dtype=float)
            z = np.asarray(args[2], dtype=float)
            m_sorted = np.asarray(args[3], dtype=float)

            center_x = kwargs.get("center_x", 0.0)
            center_y = kwargs.get("center_y", 0.0)
            center_z = kwargs.get("center_z", 0.0)

            r = np.sqrt((x - center_x) ** 2 + (y - center_y) ** 2 + (z - center_z) ** 2)
        else:
            raise TypeError(
                "calculate_half_mass_radius expects either a Particle-like object "
                "or x, y, z, mass arrays."
            )

        if r.size == 0 or m_sorted.size == 0:
            return np.nan

        idx = np.argsort(r)
        r_sorted = r[idx]
        m_sorted = m_sorted[idx]
        m_cum = np.cumsum(m_sorted)
        if m_cum.size == 0 or m_cum[-1] <= 0:
            return np.nan
        half_mass = m_cum[-1] / 2
        idx_half = np.searchsorted(m_cum, half_mass)
        return r_sorted[idx_half] if idx_half < len(r_sorted) else r_sorted[-1]

    @staticmethod
    def find_center_2d(x, y, units='degree', shrink_rate=0.95, threshold=0.1):
        x_arr = np.array(x)
        y_arr = np.array(y)

        if x_arr.size == 0 or y_arr.size == 0:
            return np.nan, np.nan

        if units == 'degree':
            # Absolute sky coordinates keep the legacy full-sky histogram.
            hist_range = [[0, 360], [-90, 90]]
            bins = [360, 180]
        elif units == 'offset':
            # Local tangent-plane offsets can be negative and should use data-driven ranges.
            x_pad = max(np.ptp(x_arr) * 0.05, threshold)
            y_pad = max(np.ptp(y_arr) * 0.05, threshold)
            hist_range = [
                [x_arr.min() - x_pad, x_arr.max() + x_pad],
                [y_arr.min() - y_pad, y_arr.max() + y_pad],
            ]
            bins = [200, 200]
        elif units == 'kpc':
            x_pad = max(np.ptp(x_arr) * 0.05, threshold)
            y_pad = max(np.ptp(y_arr) * 0.05, threshold)
            hist_range = [
                [x_arr.min() - x_pad, x_arr.max() + x_pad],
                [y_arr.min() - y_pad, y_arr.max() + y_pad],
            ]
            bins = [200, 200]
        else:
            raise ValueError("units must be 'degree', 'offset', or 'kpc'")

        hist, xedges, yedges = np.histogram2d(
            x_arr,
            y_arr,
            bins=bins,
            range=hist_range,
        )

        i, j = np.unravel_index(np.argmax(hist), hist.shape)
        center_x = np.mean(xedges[i:i + 2])
        center_y = np.mean(yedges[j:j + 2])

        if units == 'degree':
            range_x, range_y = core_radius / 3, core_radius / 3
        elif units == 'kpc':
            range_x = max(np.ptp(x_arr) / 2, core_radius / 3, threshold)
            range_y = max(np.ptp(y_arr) / 2, core_radius / 3, threshold)
        else:
            range_x = max(np.ptp(x_arr) / 2, threshold)
            range_y = max(np.ptp(y_arr) / 2, threshold)

        while range_x > threshold and range_y > threshold:
            mask = (
                (x_arr >= center_x - range_x) & (x_arr <= center_x + range_x) &
                (y_arr >= center_y - range_y) & (y_arr <= center_y + range_y)
            )
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
    def find_center_3d(x, y, z, mass=None, shrink_rate=0.95, threshold=0.1):
        x_arr = np.asarray(x, dtype=float)
        y_arr = np.asarray(y, dtype=float)
        z_arr = np.asarray(z, dtype=float)

        if x_arr.size == 0 or y_arr.size == 0 or z_arr.size == 0:
            return np.nan, np.nan, np.nan

        if mass is not None:
            mass_arr = np.asarray(mass, dtype=float)
            if mass_arr.size != x_arr.size:
                raise ValueError("mass must have the same length as x, y, z")
        else:
            mass_arr = None

        if mass_arr is not None and np.isfinite(mass_arr).any() and np.sum(mass_arr) > 0:
            center_x = np.average(x_arr, weights=mass_arr)
            center_y = np.average(y_arr, weights=mass_arr)
            center_z = np.average(z_arr, weights=mass_arr)
        else:
            center_x = np.median(x_arr)
            center_y = np.median(y_arr)
            center_z = np.median(z_arr)

        radius = np.max(
            np.sqrt((x_arr - center_x) ** 2 + (y_arr - center_y) ** 2 + (z_arr - center_z) ** 2)
        )
        if not np.isfinite(radius):
            return np.nan, np.nan, np.nan

        while radius > threshold:
            r = np.sqrt((x_arr - center_x) ** 2 + (y_arr - center_y) ** 2 + (z_arr - center_z) ** 2)
            mask = r <= radius
            if not np.any(mask):
                break

            if mass_arr is not None:
                mass_masked = mass_arr[mask]
                if np.isfinite(mass_masked).any() and np.sum(mass_masked) > 0:
                    center_x = np.average(x_arr[mask], weights=mass_masked)
                    center_y = np.average(y_arr[mask], weights=mass_masked)
                    center_z = np.average(z_arr[mask], weights=mass_masked)
                else:
                    center_x = np.median(x_arr[mask])
                    center_y = np.median(y_arr[mask])
                    center_z = np.median(z_arr[mask])
            else:
                center_x = np.median(x_arr[mask])
                center_y = np.median(y_arr[mask])
                center_z = np.median(z_arr[mask])

            radius *= shrink_rate

        return center_x, center_y, center_z

    @staticmethod
    def find_center(*args, **kwargs):
        return Analysis.find_center_2d(*args, **kwargs)

    @staticmethod
    def compute_sigma_mw(m, r, rhalf, f=0.25):
        m = np.asarray(m, dtype=float)
        r = np.asarray(r, dtype=float)
        rhalf = np.asarray(rhalf, dtype=float)
        with np.errstate(divide='ignore', invalid='ignore'):
            g_mw = np.where(r > 0, 4.302e-6 * m / (r ** 2), np.nan)
        delta_sigma2 = 2 * np.sqrt(2) * g_mw * rhalf * f
        sigma_mw = np.sqrt(np.maximum(delta_sigma2, 0.0))

        return sigma_mw

    @staticmethod
    def estimate_local_density(x, y, z, k=32):
        pos = np.vstack((x, y, z)).T
        tree = cKDTree(pos)

        dist, _ = tree.query(pos, k=k)
        rk = dist[:, -1]
        volume = (4 / 3) * np.pi * rk ** 3

        density = k / volume
        return density

    @staticmethod
    def rotate_to_sky(l, b, center=None):
        l_arr = np.array(l)
        b_arr = np.array(b)

        if l_arr.size == 0 or b_arr.size == 0:
            return np.array([], dtype=float), np.array([], dtype=float)

        if center is None:
            center_l, center_b = Analysis.find_center_2d(l, b, units='degree')
        elif len(center) == 2:
            center_l, center_b = center[0], center[1]
        else:
            raise ValueError('Please input correct center coordinate')

        l_rad_arr, b_rad_arr = np.radians(l_arr), np.radians(b_arr)
        center_rad_l, center_rad_b = np.radians(center_l), np.radians(center_b)

        deltal = np.cos(b_rad_arr) * np.sin(l_rad_arr - center_rad_l)
        deltab = (
            np.sin(b_rad_arr) * np.cos(center_rad_b)
            - np.cos(b_rad_arr) * np.sin(center_rad_b) * np.cos(l_rad_arr - center_rad_l)
        )

        dl_deg = np.rad2deg(deltal)
        db_deg = np.rad2deg(deltab)

        return dl_deg, db_deg

    @staticmethod
    def L_to_ABmag(mass, d_kpc, mag_sys='vega', m_l_relation=sim_plot_mass_to_light_ratio):
        l_d = mass / m_l_relation
        m_ab = -2.5 * np.log10(l_d) + 5 * np.log10(d_kpc * 1000) + 4.83 - 5
        m_vega = m_ab - 0.02

        if mag_sys == 'AB':
            return m_ab
        if mag_sys == 'vega':
            return m_vega
        raise ValueError(f"Unsupported magnitude system: {mag_sys}")

    @staticmethod
    def surface_brightness_map(
        x_kpc,
        y_kpc,
        mass,
        d_mean,
        npix=100,
        size_kpc=10.8,
        mag_sys='vega',
        m_l_relation=sim_plot_mass_to_light_ratio,
        coordinate_unit='kpc',
        coordinate_size=None,
    ):
        coordinate_size = float(size_kpc if coordinate_size is None else coordinate_size)
        h, xedges, yedges = np.histogram2d(
            x_kpc,
            y_kpc,
            bins=npix,
            range=[
                [-coordinate_size / 2, coordinate_size / 2],
                [-coordinate_size / 2, coordinate_size / 2],
            ],
            weights=mass,
        )

        img_mag = np.full_like(h, np.inf, dtype=float)

        pix_size = coordinate_size / npix
        if coordinate_unit == 'kpc':
            pix_size_rad = pix_size / d_mean
        elif coordinate_unit in ('deg', 'degree'):
            pix_size_rad = np.radians(pix_size)
        elif coordinate_unit in ('rad', 'radian'):
            pix_size_rad = pix_size
        else:
            raise ValueError(f"Unsupported coordinate_unit: {coordinate_unit}")
        pix_area_arcsec2 = (pix_size_rad * 206265) ** 2

        positive = h > 0
        if np.any(positive):
            m_tot = Analysis.L_to_ABmag(h[positive], d_mean, mag_sys=mag_sys, m_l_relation=m_l_relation)
            img_mag[positive] = m_tot + 2.5 * np.log10(pix_area_arcsec2)

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

        gas_loss_rate_smooth = df['mass'].rolling(window, center=True).apply(fit_slope, raw=False)

        dt_mean = df['age'].diff().mean()
        gas_loss_rate_smooth = -gas_loss_rate_smooth / dt_mean
        return -gas_loss_rate_smooth.values

    @staticmethod
    def calculate_ellipticity(x, y, mass, n_neighbors=20):
        shape = Analysis.calculate_projected_shape(
            x,
            y,
            mass,
            n_neighbors=n_neighbors,
        )
        return shape['eps'], shape['pa']

    @staticmethod
    def calculate_projected_shape(x, y, mass, n_neighbors=20):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        mass = np.asarray(mass, dtype=float)
        invalid_shape = {
            'center_x': np.nan,
            'center_y': np.nan,
            'eps': np.nan,
            'pa': np.nan,
            'q': np.nan,
        }
        if x.size == 0 or y.size == 0 or mass.size == 0:
            return invalid_shape

        finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(mass) & (mass > 0)
        if not np.any(finite):
            return invalid_shape
        x = x[finite]
        y = y[finite]
        mass = mass[finite]

        coords = np.vstack([x, y]).T
        if coords.shape[0] < 2:
            return {
                'center_x': float(x[0]),
                'center_y': float(y[0]),
                'eps': 0.0,
                'pa': 0.0,
                'q': 1.0,
            }

        effective_neighbors = min(max(int(n_neighbors), 1), coords.shape[0])
        nn = NearestNeighbors(n_neighbors=effective_neighbors).fit(coords)
        distances, _ = nn.kneighbors(coords)
        local_density = 1 / (np.mean(distances, axis=1) + 1e-5)

        weights = mass * local_density
        weight_sum = np.sum(weights)
        if weight_sum <= 0:
            return invalid_shape
        weights /= weight_sum
        x_c = np.sum(x * weights)
        y_c = np.sum(y * weights)

        x_rel = x - x_c
        y_rel = y - y_c

        ixx = np.sum(x_rel ** 2 * weights)
        iyy = np.sum(y_rel ** 2 * weights)
        ixy = np.sum(x_rel * y_rel * weights)

        tmp = np.sqrt((ixx - iyy) ** 2 + 4 * ixy ** 2)
        a2 = (ixx + iyy + tmp) / 2
        b2 = (ixx + iyy - tmp) / 2
        if a2 <= 0:
            return invalid_shape

        eps = 1 - np.sqrt(np.clip(b2 / a2, 0.0, None))
        pa = 0.5 * np.arctan2(2 * ixy, ixx - iyy)
        q = 1.0 - eps

        return {
            'center_x': float(x_c),
            'center_y': float(y_c),
            'eps': float(eps),
            'pa': float(pa),
            'q': float(q),
        }

    @staticmethod
    def circularized_half_light_radius(r_half_major, ep):
        r_half_major = np.asarray(r_half_major, dtype=float)
        ep = np.asarray(ep, dtype=float)
        q = np.clip(1.0 - ep, np.finfo(float).eps, None)
        with np.errstate(invalid='ignore'):
            r_half_circularized = r_half_major * np.sqrt(q)
        return r_half_circularized

    @staticmethod
    def projected_elliptical_radius(x, y, ep, pa=0.0, center_x=0.0, center_y=0.0):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        x_centered = x - center_x
        y_centered = y - center_y

        cos_pa = np.cos(-pa)
        sin_pa = np.sin(-pa)
        x_rot = cos_pa * x_centered - sin_pa * y_centered
        y_rot = sin_pa * x_centered + cos_pa * y_centered

        q = 1 - ep
        if q <= 0:
            q = np.finfo(float).eps
        return np.sqrt(x_rot ** 2 + (y_rot / q) ** 2)

    @staticmethod
    def get_elliptical_radial_mask(x, y, radius, ep, pa=0.0, center_x=0.0, center_y=0.0):
        if not np.isfinite(radius) or radius <= 0:
            return np.zeros_like(np.asarray(x, dtype=float), dtype=bool)
        r = Analysis.projected_elliptical_radius(
            x,
            y,
            ep,
            pa=pa,
            center_x=center_x,
            center_y=center_y,
        )
        return r < radius

    @staticmethod
    def radial_magnitude_profile(
        x_rot,
        y_rot,
        mass,
        d_kpc,
        ep,
        r_cut=core_radius,
        bins=20,
        mag_sys='vega',
        m_l_relation=sim_plot_mass_to_light_ratio,
        pa=0.0,
        center_x=0.0,
        center_y=0.0,
        r_edges=None,
        r_intervals=None,
    ):
        x_rot = np.asarray(x_rot, dtype=float)
        y_rot = np.asarray(y_rot, dtype=float)
        mass = np.asarray(mass, dtype=float)

        if r_intervals is not None:
            intervals = np.asarray(r_intervals, dtype=float)
            if intervals.ndim != 2 or intervals.shape[1] != 2:
                raise ValueError("r_intervals must have shape (n, 2)")
            r_inner = intervals[:, 0]
            r_outer = intervals[:, 1]
        elif r_edges is not None:
            edges = np.asarray(r_edges, dtype=float)
            if edges.ndim != 1 or edges.size < 2:
                raise ValueError("r_edges must be a one-dimensional array with at least two entries")
            r_inner = edges[:-1]
            r_outer = edges[1:]
        else:
            edges = np.linspace(0, r_cut, bins + 1)
            r_inner = edges[:-1]
            r_outer = edges[1:]

        r_centers = 0.5 * (r_inner + r_outer)
        if x_rot.size == 0 or y_rot.size == 0 or mass.size == 0 or d_kpc <= 0:
            return r_centers, np.full_like(r_centers, np.nan, dtype=float)

        q = 1 - ep
        if q <= 0:
            q = np.finfo(float).eps
        r_elliptical = Analysis.projected_elliptical_radius(
            x_rot,
            y_rot,
            ep,
            pa=pa,
            center_x=center_x,
            center_y=center_y,
        )
        factor = (206265 / d_kpc) ** 2

        if r_intervals is not None:
            mass_per_bin = np.array(
                [
                    np.sum(mass[(r_elliptical >= inner) & (r_elliptical < outer)])
                    if np.isfinite(inner) and np.isfinite(outer) and outer > inner
                    else np.nan
                    for inner, outer in zip(r_inner, r_outer)
                ],
                dtype=float,
            )
        else:
            mass_per_bin, _ = np.histogram(r_elliptical, bins=np.r_[r_inner, r_outer[-1]], weights=mass)
        area_arcsec2 = np.pi * (r_outer ** 2 - r_inner ** 2) * q * factor

        mag_profile = np.full_like(r_centers, np.nan, dtype=float)
        positive = (mass_per_bin > 0) & np.isfinite(area_arcsec2) & (area_arcsec2 > 0)
        if np.any(positive):
            mag_r = Analysis.L_to_ABmag(
                mass=mass_per_bin[positive],
                d_kpc=d_kpc,
                mag_sys=mag_sys,
                m_l_relation=m_l_relation,
            )
            mag_profile[positive] = mag_r + 2.5 * np.log10(area_arcsec2[positive])

        return r_centers, mag_profile

    @staticmethod
    def half_light_radius(x, y, mass, ep, pa=0.0, center_x=0.0, center_y=0.0):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        mass = np.asarray(mass, dtype=float)
        if x.size == 0 or y.size == 0 or mass.size == 0:
            return np.nan
        r = Analysis.projected_elliptical_radius(
            x,
            y,
            ep,
            pa=pa,
            center_x=center_x,
            center_y=center_y,
        )

        idx = np.argsort(r)
        r_sorted = r[idx]
        m_sorted = mass[idx]

        cum_mass = np.cumsum(m_sorted)
        if cum_mass.size == 0 or cum_mass[-1] <= 0:
            return np.nan
        half_mass = 0.5 * cum_mass[-1]

        re = np.interp(half_mass, cum_mass, r_sorted)
        return re

    @staticmethod
    def get_vlos_dispersion(x, y, vlos, r_range=core_radius, bins=50):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        vlos = np.asarray(vlos, dtype=float)
        r_proj = np.sqrt(x ** 2 + y ** 2)
        mask_r = r_proj < r_range
        if not np.any(mask_r):
            return np.array([], dtype=float), np.array([], dtype=float)

        r_max = r_proj[mask_r].max()
        if not np.isfinite(r_max) or r_max <= 0:
            return np.array([], dtype=float), np.array([], dtype=float)

        r_bins = np.linspace(0, r_max, bins)
        bin_centers = 0.5 * (r_bins[:-1] + r_bins[1:])

        vlos_dispersion, _, _ = binned_statistic(r_proj[mask_r], vlos[mask_r], statistic='std', bins=r_bins)

        return bin_centers, vlos_dispersion

    @staticmethod
    def est_gas_contour(
        x,
        y,
        m,
        box=core_radius * 0.8,
        nbins=50,
        cd_threshold=cd_threshold,
        smooth_sigma=2.0,
        unit="kpc",
        d_mean=None,
        histogram_range=None,
        smooth_mode="reflect",
    ):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        m = np.asarray(m, dtype=float)

        mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(m)
        if box is not None:
            mask &= (np.abs(x) < box) & (np.abs(y) < box)

        x_cut = x[mask]
        y_cut = y[mask]
        m_cut = m[mask]

        if histogram_range is None and x_cut.size == 0:
            empty = np.empty((0, 0), dtype=float)
            return 0.0, empty, empty, empty

        h, xedges, yedges = np.histogram2d(
            x_cut,
            y_cut,
            bins=nbins,
            range=histogram_range,
            weights=m_cut,
        )

        xc = 0.5 * (xedges[:-1] + xedges[1:])
        yc = 0.5 * (yedges[:-1] + yedges[1:])
        x_grid, y_grid = np.meshgrid(xc, yc)

        dx = xedges[1] - xedges[0]
        dy = yedges[1] - yedges[0]
        area = dx * dy

        cd = h / area * 1.25e14
        cd_smooth = gaussian_filter(cd, sigma=smooth_sigma, mode=smooth_mode, cval=0.0)

        mass_high = h[cd_smooth > cd_threshold].sum()

        if unit == "deg":
            if d_mean is None:
                raise ValueError("unit='deg' requires d_mean (distance in kpc).")
            kpc_to_deg = (180.0 / np.pi) / d_mean
            x_grid = x_grid * kpc_to_deg
            y_grid = y_grid * kpc_to_deg

        return mass_high, x_grid, y_grid, cd_smooth

    @staticmethod
    def get_vlos_field(x_kpc, y_kpc, vlos_array, size_kpc=10.8, npix=100):
        vlos_residual = vlos_array - np.mean(vlos_array)

        h, xedges, yedges = np.histogram2d(
            x_kpc,
            y_kpc,
            bins=npix,
            range=[[-size_kpc / 2, size_kpc / 2], [-size_kpc / 2, size_kpc / 2]],
            weights=vlos_residual,
        )

        n, _, _ = np.histogram2d(
            x_kpc,
            y_kpc,
            bins=npix,
            range=[[-size_kpc / 2, size_kpc / 2], [-size_kpc / 2, size_kpc / 2]],
        )

        with np.errstate(invalid='ignore', divide='ignore'):
            vlos_mean = np.where(n > 0, h / n, np.nan)

        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        return vlos_mean, extent

    @staticmethod
    def get_radial_mask(x, y, radius):
        r = np.sqrt(x ** 2 + y ** 2)
        return r < radius

    @staticmethod
    def get_3d_radial_mask(x, y, z, radius):
        r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
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
            ),
        )

        boundary_path = Path(largest_contour)
        return boundary_path

    @staticmethod
    def cartesian_to_spherical(x, y, z, vx, vy, vz):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)
        vx = np.asarray(vx, dtype=float)
        vy = np.asarray(vy, dtype=float)
        vz = np.asarray(vz, dtype=float)
        if x.size == 0:
            empty = np.array([], dtype=float)
            return empty, empty, empty

        r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        r_safe = np.where(r > 0, r, np.finfo(float).eps)
        theta = np.arccos(np.clip(z / r_safe, -1.0, 1.0))
        phi = np.arctan2(y, x)

        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        sin_phi = np.sin(phi)
        cos_phi = np.cos(phi)

        e_r = np.stack([sin_theta * cos_phi, sin_theta * sin_phi, cos_theta], axis=1)
        e_theta = np.stack([cos_theta * cos_phi, cos_theta * sin_phi, -sin_theta], axis=1)
        e_phi = np.stack([-sin_phi, cos_phi, np.zeros_like(x)], axis=1)

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
        v_array = u.Quantity(
            np.sqrt(dw_elinfo['vr'] ** 2 + dw_elinfo['vtheta'] ** 2 + dw_elinfo['vphi'] ** 2),
            u.km / u.s,
        )
        mw_array = dw_elinfo['mw_mass_r'].values * u.Msun
        rhalf_values = dw_elinfo['rhalf_circularized'].to_numpy(dtype=float)
        r_half_km2 = (rhalf_values * u.kpc).to(u.km).value ** 2
        sigma_array = dw_elinfo['sigma'].values

        expr = (G * mw_array) / (d_array ** 2 * v_array)
        tidal_eff = u.Quantity((128.0 / 27.0) * expr ** 2, 1 / u.s ** 2).value

        dw_energy = sigma_array ** 2 / r_half_km2

        g_mw = (G * mw_array / d_array ** 2).to(u.km / u.s ** 2)
        sigma_mw = np.sqrt((np.sqrt(2) * 1 * g_mw * rhalf_values * u.km)).value

        return tidal_eff, dw_energy, sigma_mw

    @staticmethod
    def GetGasdensity(x, y, m, r_half, ep=0.0, pa=0.0, center_x=0.0, center_y=0.0):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        m = np.asarray(m, dtype=float)
        if x.size == 0 or y.size == 0 or m.size == 0 or not np.isfinite(r_half) or r_half <= 0:
            return 0.0
        r = Analysis.projected_elliptical_radius(
            x,
            y,
            ep,
            pa=pa,
            center_x=center_x,
            center_y=center_y,
        )
        mask = r <= r_half

        m_tot = m[mask].sum()
        q = 1 - ep
        if q <= 0:
            q = np.finfo(float).eps
        area_kpc2 = np.pi * r_half ** 2 * q
        sigma_kpc2 = m_tot / area_kpc2
        sigma_pc2 = sigma_kpc2 / 1e6
        n_h = sigma_pc2 * 1.25e20

        return n_h

    @staticmethod
    def particles_to_image(x, y, m, npix=200, size_kpc=core_radius):
        xmin, xmax = -size_kpc, size_kpc
        ymin, ymax = -size_kpc, size_kpc

        image, _, _ = np.histogram2d(
            y,
            x,
            bins=npix,
            range=[[ymin, ymax], [xmin, xmax]],
            weights=m,
        )
        return image
