import io
import os
import numpy as np
import pandas as pd
from contextlib import redirect_stdout
import h5py

import readsnap as giz
import zim

AUTO_BOXSIZE_PRINTED_LOCAL = False


class SnapshotLoaderMixin:
    GAS_FIELD_ORDER = ('temp', 'rho', 'nh', 'sfr', 'ugas')
    STAR_FIELD_ORDER = ('birth',)
    DEFAULT_PARTICLE_TYPES = (0, 1, 2, 3, 4)

    @staticmethod
    def _readsnap_quiet(*args, **kwargs):
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            return giz.readsnap(*args, **kwargs)

    @classmethod
    def _resolve_requested_fields(cls, dftype, gas_fields=None, star_fields=None):
        if dftype >= 2:
            if gas_fields is None:
                gas_columns = list(cls.GAS_FIELD_ORDER)
            else:
                gas_columns = [field for field in cls.GAS_FIELD_ORDER if field in set(gas_fields)]
        else:
            gas_columns = []

        if dftype >= 3:
            if star_fields is None:
                star_columns = list(cls.STAR_FIELD_ORDER)
            else:
                star_columns = [field for field in cls.STAR_FIELD_ORDER if field in set(star_fields)]
        else:
            star_columns = []

        return gas_columns, star_columns

    def _resolve_snapshot_hdf5(self):
        fname, fname_base, fname_ext = giz.check_if_filename_exists(
            self.folder_path,
            self.snapshot_num,
            snapshot_name='snapshot',
            extension='.hdf5',
            four_char=0,
        )
        if fname == 'NULL' or fname_ext != '.hdf5':
            return None
        return fname, fname_base, fname_ext

    @classmethod
    def _resolve_particle_types(cls, particle_types=None):
        selected_types = cls.DEFAULT_PARTICLE_TYPES if particle_types is None else tuple(particle_types)
        if not selected_types:
            raise ValueError("particle_types must contain at least one particle type")
        invalid_types = [tp for tp in selected_types if tp not in cls.DEFAULT_PARTICLE_TYPES]
        if invalid_types:
            raise ValueError(f"Unsupported particle types requested: {invalid_types}")
        return selected_types

    @staticmethod
    def _iter_hdf5_parts(fname_base, fname_ext, numfiles):
        if numfiles <= 1:
            yield fname_base + fname_ext
            return
        for i_file in range(numfiles):
            yield f"{fname_base}.{i_file}{fname_ext}"

    @staticmethod
    def _fill_constant(target, slice_idx, value):
        target[slice_idx] = np.float32(value)

    @staticmethod
    def _announce_auto_boxsize_once(boxsize):
        global AUTO_BOXSIZE_PRINTED_LOCAL

        run_id = os.environ.get("DSPH_RUN_ID")
        if run_id:
            marker_path = f"/tmp/dsph_boxsize_{run_id}.flag"
            try:
                fd = os.open(marker_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.close(fd)
                print(f'Using auto-calculated boxsize: {boxsize}')
            except FileExistsError:
                pass
            return

        if not AUTO_BOXSIZE_PRINTED_LOCAL:
            AUTO_BOXSIZE_PRINTED_LOCAL = True
            print(f'Using auto-calculated boxsize: {boxsize}')

    def _load_snapshot_direct_hdf5(self, dftype=3, box=True, gas_fields=None, star_fields=None, particle_types=None):
        resolved = self._resolve_snapshot_hdf5()
        if resolved is None:
            return None
        selected_types = self._resolve_particle_types(particle_types)

        fname, fname_base, fname_ext = resolved
        with h5py.File(fname, 'r') as first_file:
            header = first_file['Header'].attrs
            self.time = float(header['Time'])
            npart_total = np.asarray(header['NumPart_Total'], dtype=np.int64)
            mass_table = np.asarray(header['MassTable'], dtype=np.float64)
            numfiles = int(header['NumFilesPerSnapshot'])
            flag_cooling = int(header['Flag_Cooling'])

        total_particles = int(sum(int(npart_total[tp]) for tp in selected_types))
        gas_columns, star_columns = self._resolve_requested_fields(
            dftype,
            gas_fields=gas_fields,
            star_fields=star_fields,
        )

        df_data = {
            'id': np.zeros(total_particles, dtype=np.int32),
            'x': np.zeros(total_particles, dtype=np.float32),
            'y': np.zeros(total_particles, dtype=np.float32),
            'z': np.zeros(total_particles, dtype=np.float32),
            'vx': np.zeros(total_particles, dtype=np.float32),
            'vy': np.zeros(total_particles, dtype=np.float32),
            'vz': np.zeros(total_particles, dtype=np.float32),
            'm': np.zeros(total_particles, dtype=np.float32),
            'tp': np.zeros(total_particles, dtype=np.int32),
        }

        if dftype >= 2:
            for col in gas_columns:
                df_data[col] = np.zeros(total_particles, dtype=np.float32)
        if dftype == 4:
            df_data['ags'] = np.zeros(total_particles, dtype=np.float32)
        if dftype >= 3:
            for col in star_columns:
                df_data[col] = np.zeros(total_particles, dtype=np.float32)

        offsets = {tp: 0 for tp in selected_types}
        running_offset = 0
        for tp in selected_types:
            offsets[tp] = running_offset
            running_offset += int(npart_total[tp])
        for snapshot_file in self._iter_hdf5_parts(fname_base, fname_ext, numfiles):
            with h5py.File(snapshot_file, 'r') as file:
                npart_file = np.asarray(file['Header'].attrs['NumPart_ThisFile'], dtype=np.int64)

                for tp in selected_types:
                    npart = int(npart_file[tp])
                    if npart <= 0:
                        continue

                    group = file[f'PartType{tp}']
                    start = int(offsets[tp])
                    stop = start + npart
                    slice_idx = slice(start, stop)

                    coords = group['Coordinates']
                    vels = group['Velocities']
                    df_data['x'][slice_idx] = coords[:, 0]
                    df_data['y'][slice_idx] = coords[:, 1]
                    df_data['z'][slice_idx] = coords[:, 2]
                    df_data['vx'][slice_idx] = vels[:, 0]
                    df_data['vy'][slice_idx] = vels[:, 1]
                    df_data['vz'][slice_idx] = vels[:, 2]
                    df_data['id'][slice_idx] = group['ParticleIDs'][:]
                    df_data['tp'][slice_idx] = tp

                    if mass_table[tp] > 0.0:
                        self._fill_constant(df_data['m'], slice_idx, mass_table[tp] * 1e10)
                    else:
                        df_data['m'][slice_idx] = group['Masses'][:] * 1e10

                    if tp == 0 and dftype >= 2:
                        need_ugas = ('ugas' in gas_columns) or ('temp' in gas_columns)
                        ugas = group['InternalEnergy'][:] if need_ugas else None

                        need_ne = ('temp' in gas_columns)
                        need_nh = ('nh' in gas_columns)
                        if (need_ne or need_nh) and flag_cooling > 0 and 'ElectronAbundance' in group:
                            ne = group['ElectronAbundance'][:] if need_ne else None
                            nh = group['NeutralHydrogenAbundance'][:] if need_nh else None
                        else:
                            ne = np.zeros(npart, dtype=np.float32) if need_ne else None
                            nh = np.zeros(npart, dtype=np.float32) if need_nh else None

                        if 'ugas' in gas_columns:
                            df_data['ugas'][slice_idx] = ugas
                        if 'temp' in gas_columns:
                            df_data['temp'][slice_idx] = zim.unitconversion(ugas, ne)
                        if 'rho' in gas_columns:
                            df_data['rho'][slice_idx] = group['Density'][:] * 404.7
                        if 'nh' in gas_columns:
                            df_data['nh'][slice_idx] = nh
                        if 'sfr' in gas_columns and 'StarFormationRate' in group:
                            df_data['sfr'][slice_idx] = group['StarFormationRate'][:]
                        if dftype == 4 and 'AGS-Softening' in group:
                            df_data['ags'][slice_idx] = group['AGS-Softening'][:]

                    if 'birth' in star_columns and tp == 4 and dftype >= 3 and 'StellarFormationTime' in group:
                        df_data['birth'][slice_idx] = group['StellarFormationTime'][:]

                    offsets[tp] = stop

        if self.boxsize == 'Auto':
            self.boxsize = round(np.ptp(df_data['x']), -1)
            self._announce_auto_boxsize_once(self.boxsize)

        if box:
            box_center = self.boxsize / 2
            df_data['x'] -= box_center
            df_data['y'] -= box_center
            df_data['z'] -= box_center
        else:
            print("Box centering skipped. Ensure coordinates are already centered if needed.")

        df = pd.DataFrame(df_data, copy=False)
        self.df = df
        return df, self.time

    def load_snapshot(self, dftype=3, box=True, gas_fields=None, star_fields=None, particle_types=None):
        """Load snapshot data into a DataFrame."""
        try:
            selected_types = self._resolve_particle_types(particle_types)
            try:
                direct_result = self._load_snapshot_direct_hdf5(
                    dftype=dftype,
                    box=box,
                    gas_fields=gas_fields,
                    star_fields=star_fields,
                    particle_types=selected_types,
                )
            except Exception:
                direct_result = None
            if direct_result is not None:
                return direct_result

            header_ptype = selected_types[0]
            hdr = self._readsnap_quiet(self.folder_path, self.snapshot_num, header_ptype, header_only=1)
            if hdr is None:
                raise ValueError(f"Failed to read snapshot {self.snapshot_num} from {self.folder_path}")

            self.time = hdr['time']
            npart_total = hdr['npartTotal']
            total_particles = int(sum(int(npart_total[tp]) for tp in selected_types))

            gas_columns, star_columns = self._resolve_requested_fields(
                dftype,
                gas_fields=gas_fields,
                star_fields=star_fields,
            )

            df_data = {
                'id': np.zeros(total_particles, dtype=np.int32),
                'x': np.zeros(total_particles, dtype=np.float32),
                'y': np.zeros(total_particles, dtype=np.float32),
                'z': np.zeros(total_particles, dtype=np.float32),
                'vx': np.zeros(total_particles, dtype=np.float32),
                'vy': np.zeros(total_particles, dtype=np.float32),
                'vz': np.zeros(total_particles, dtype=np.float32),
                'm': np.zeros(total_particles, dtype=np.float32),
                'tp': np.zeros(total_particles, dtype=np.int32),
            }

            if dftype >= 2:
                for col in gas_columns:
                    df_data[col] = np.zeros(total_particles, dtype=np.float32)
            if dftype == 4:
                df_data['ags'] = np.zeros(total_particles, dtype=np.float32)
            if dftype >= 3:
                for col in star_columns:
                    df_data[col] = np.zeros(total_particles, dtype=np.float32)

            offset = 0
            for tp in selected_types:
                data = self._readsnap_quiet(
                    self.folder_path,
                    self.snapshot_num,
                    tp,
                    readags=(tp == 0 and dftype == 4),
                )
                npart = npart_total[tp]

                if npart == 0:
                    continue
                if data is None or len(data['p'][:, 0]) != npart:
                    raise RuntimeError(f"Particle count mismatch for type {tp}")

                slice_idx = slice(offset, offset + npart)
                pos = data['p'].astype(np.float32, copy=False)
                vel = data['v'].astype(np.float32, copy=False)

                df_data['id'][slice_idx] = data['id'].astype(np.int32, copy=False)
                df_data['x'][slice_idx] = pos[:, 0]
                df_data['y'][slice_idx] = pos[:, 1]
                df_data['z'][slice_idx] = pos[:, 2]
                df_data['vx'][slice_idx] = vel[:, 0]
                df_data['vy'][slice_idx] = vel[:, 1]
                df_data['vz'][slice_idx] = vel[:, 2]
                df_data['m'][slice_idx] = (data['m'] * 1e10).astype(np.float32, copy=False)
                df_data['tp'][slice_idx] = tp

                if tp == 0 and dftype >= 2:
                    if 'ugas' in gas_columns:
                        df_data['ugas'][slice_idx] = data['u'].astype(np.float32, copy=False)
                    if 'temp' in gas_columns:
                        temp_values = zim.unitconversion(data['u'], data['ne'])
                        df_data['temp'][slice_idx] = temp_values.astype(np.float32, copy=False)
                    if 'rho' in gas_columns:
                        df_data['rho'][slice_idx] = (data['rho'] * 404.7).astype(np.float32, copy=False)
                    if 'nh' in gas_columns:
                        df_data['nh'][slice_idx] = data['nh'].astype(np.float32, copy=False)
                    if 'sfr' in gas_columns:
                        df_data['sfr'][slice_idx] = data['sfr'].astype(np.float32, copy=False)
                    if dftype == 4 and 'ags' in data:
                        df_data['ags'][slice_idx] = data['ags'].astype(np.float32, copy=False)

                if 'birth' in star_columns and tp == 4 and dftype >= 3:
                    df_data['birth'][slice_idx] = data['age'].astype(np.float32, copy=False)

                offset += npart

            if self.boxsize == 'Auto':
                self.boxsize = round(np.ptp(df_data['x']), -1)
                self._announce_auto_boxsize_once(self.boxsize)

            if box:
                box_center = self.boxsize / 2
                df_data['x'] -= box_center
                df_data['y'] -= box_center
                df_data['z'] -= box_center
            else:
                print("Box centering skipped. Ensure coordinates are already centered if needed.")

            df = pd.DataFrame(df_data, copy=False)
            self.df = df
            return df, self.time

        except Exception as e:
            raise RuntimeError(f"Error loading snapshot: {str(e)}")
