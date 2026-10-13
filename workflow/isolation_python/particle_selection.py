import numpy as np
import pandas as pd
from sklearn.neighbors import KDTree

from analysis_core import Analysis as AnalysisCore
from coordinate_transform import CoordinateTransform
from variable import fornax_core_radius


core_radius = fornax_core_radius


class ParticleSelectionMixin:
    def find_MW_particles(self, dw_particles_num, cube_file):
        """Find Milky Way particles."""
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
        """Find dwarf galaxy particles."""
        if dw_particles_num <= 0:
            return pd.Series([False] * len(self.df), index=self.df.index)
        dw_all = self.df.nlargest(dw_particles_num, 'id')
        return self.df['id'].isin(dw_all['id'])

    def classify_mw_dwarf(self, r_exclude=5, dwarf_radius=5, k_density=32, df=None):
        df = self.df if df is None else df
        star_sel = df.tp >= 2
        df_star = df[star_sel]

        x = df_star["x"].to_numpy()
        y = df_star["y"].to_numpy()
        z = df_star["z"].to_numpy()

        r = np.sqrt(x * x + y * y + z * z)
        density = AnalysisCore.estimate_local_density(x, y, z, k=k_density)
        outer = r > r_exclude
        idx_peak = np.argmax(density * outer)

        xd = x[idx_peak]
        yd = y[idx_peak]
        zd = z[idx_peak]

        dx = x - xd
        dy = y - yd
        dz = z - zd

        r_d = np.sqrt(dx * dx + dy * dy + dz * dz)

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

        sel = (tp == 2) | (tp == 4)
        x = x[sel]
        y = y[sel]
        z = z[sel]
        m = m[sel]

        xc = np.average(x, weights=m)
        yc = np.average(y, weights=m)
        zc = np.average(z, weights=m)

        r_current = r_init
        while r_current > r_min:
            dx = x - xc
            dy = y - yc
            dz = z - zc

            r2 = dx * dx + dy * dy + dz * dz
            inside = r2 < r_current * r_current

            if inside.sum() < 100:
                break

            xc = np.average(x[inside], weights=m[inside])
            yc = np.average(y[inside], weights=m[inside])
            zc = np.average(z[inside], weights=m[inside])

            r_current *= shrink

        df['x'] -= xc
        df['y'] -= yc
        df['z'] -= zc
        df['r'] = np.sqrt(df['x'] ** 2 + df['y'] ** 2 + df['z'] ** 2)

        return df

    def convert_coordinates(self, df=None, copy_df=True, include_galactic=True):
        df = self.df if df is None else df
        if copy_df:
            df = df.copy()

        x = df['x'].to_numpy()
        y = df['y'].to_numpy()
        z = df['z'].to_numpy()
        vx = df['vx'].to_numpy()
        vy = df['vy'].to_numpy()
        vz = df['vz'].to_numpy()

        xh, yh, zh, vxh, vyh, vzh = CoordinateTransform.to_heliocentric(x, y, z, vx, vy, vz)
        df['xh'], df['yh'], df['zh'] = xh, yh, zh
        df['vxh'], df['vyh'], df['vzh'] = vxh, vyh, vzh

        df['rh'] = CoordinateTransform.calculate_distances(xh, yh, zh)

        l, b = CoordinateTransform.to_galactic(xh, yh, zh)
        if include_galactic:
            df['l'], df['b'] = l, b

        ra, dec = CoordinateTransform.galactic_to_equatorial(l, b)
        df['ra'], df['dec'] = ra, dec

        return df

    def convert_coordinates_for_mask(self, mask, df=None, with_velocity=True, include_galactic=False):
        df = self.df if df is None else df
        mask = np.asarray(mask, dtype=bool)

        x = df.loc[mask, 'x'].to_numpy()
        y = df.loc[mask, 'y'].to_numpy()
        z = df.loc[mask, 'z'].to_numpy()

        result = {}
        if with_velocity:
            vx = df.loc[mask, 'vx'].to_numpy()
            vy = df.loc[mask, 'vy'].to_numpy()
            vz = df.loc[mask, 'vz'].to_numpy()
            xh, yh, zh, vxh, vyh, vzh = CoordinateTransform.to_heliocentric(x, y, z, vx, vy, vz)
            result['vxh'] = vxh
            result['vyh'] = vyh
            result['vzh'] = vzh
        else:
            xh, yh, zh = CoordinateTransform.to_heliocentric(x, y, z)

        result['x'] = x
        result['y'] = y
        result['z'] = z
        result['xh'] = xh
        result['yh'] = yh
        result['zh'] = zh
        result['rh'] = CoordinateTransform.calculate_distances(xh, yh, zh)

        l, b = CoordinateTransform.to_galactic(xh, yh, zh)
        if include_galactic:
            result['l'] = l
            result['b'] = b

        ra, dec = CoordinateTransform.galactic_to_equatorial(l, b)
        result['ra'] = ra
        result['dec'] = dec

        return result

    def classify_new_dwarf_stars(self, dw_star_mask, r_cut=core_radius):
        try:
            new_stars_mask = self.df.tp == 4

            if not new_stars_mask.any() or not dw_star_mask.any():
                empty_mask = pd.Series(False, index=self.df.index, dtype=bool)
                return empty_mask, empty_mask

            new_pos = self.df.loc[new_stars_mask, ['x', 'y', 'z']].values
            dw_pos = self.df.loc[dw_star_mask, ['x', 'y', 'z']].values

            dw_tree = KDTree(dw_pos)
            dw_dist, _ = dw_tree.query(new_pos, k=1)
            is_dw = dw_dist.squeeze() < r_cut

            new_dw_stars = pd.Series(False, index=self.df.index, dtype=bool)
            new_mw_stars = pd.Series(False, index=self.df.index, dtype=bool)

            new_dw_stars.loc[new_stars_mask] = is_dw.astype(bool)
            new_mw_stars.loc[new_stars_mask] = (~is_dw).astype(bool)

            return new_dw_stars, new_mw_stars

        except Exception as e:
            raise RuntimeError(f"Error classifying new stars: {str(e)}")

    def find_dwarf_gas(self, dw_star_mask, radius=25.0):
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
        r2 = dx * dx + dy * dy + dz * dz

        dw_gas_mask = gas_sel & (r2 < radius ** 2)
        return dw_gas_mask

    def find_MW_gas(self, radius=200.0):
        df = self.df

        x = df["x"].to_numpy()
        y = df["y"].to_numpy()
        z = df["z"].to_numpy()
        tp = df["tp"].to_numpy()

        gas_sel = tp == 0
        r2 = x * x + y * y + z * z

        mw_gas_mask = gas_sel & (r2 < radius ** 2)
        return mw_gas_mask

    @staticmethod
    def get_dw_num(file_path):
        try:
            with open(file_path, 'r') as file:
                first_line = file.readline().strip()
                numbers = [
                    float(num)
                    for num in first_line.replace(',', ' ').split()
                    if num.replace('.', '', 1).isdigit()
                ]
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
