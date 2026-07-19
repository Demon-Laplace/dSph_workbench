import numpy as np
import pandas as pd
from analysis_core import Analysis as AnalysisCore
from coordinate_transform import CoordinateTransform
from data_processing import DataProcessor
from particle_selection import ParticleSelectionMixin
from snapshot_io import SnapshotLoaderMixin
from variable import fornax_core_radius, folder_path

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

class GalaxySimulation(SnapshotLoaderMixin, ParticleSelectionMixin):
    """Main class for galaxy simulation data analysis"""
    def __init__(self, folder_path, snapshot_num, boxsize='Auto'):
        self.folder_path = folder_path
        self.snapshot_num = snapshot_num
        self.boxsize = boxsize
        self.df = None
        self.time = None    

class Analysis(AnalysisCore):
    """Simulation-aware analysis helpers."""

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
        xc, yc = Analysis.find_center_2d(x_arr, y_arr, units='kpc')

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
