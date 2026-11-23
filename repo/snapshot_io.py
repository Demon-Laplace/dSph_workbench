import io
import os
import numpy as np
import pandas as pd
from contextlib import redirect_stdout
import h5py

import readsnap as giz
import zim
from canonical_observer import LENGTH_CODE_TO_KPC, VELOCITY_CODE_TO_KMS

AUTO_BOXSIZE_PRINTED_LOCAL = False


class SnapshotLoaderMixin:
    GAS_FIELD_ORDER = ('temp', 'rho', 'nh', 'sfr', 'ugas')
    STAR_FIELD_ORDER = ('birth', 'stellar_age_gyr')
    DEFAULT_PARTICLE_TYPES = (0, 1, 2, 3, 4)
    # Noncosmological production units, matching the explicit position/velocity
    # conversion. Dataset time-unit attributes take precedence for stored fields.
    stellar_time_code_to_gyr = LENGTH_CODE_TO_KPC * 3.0856775814913673e16 / VELOCITY_CODE_TO_KMS / (365.25 * 86400 * 1.e9)
    stellar_age_default_unit = 'code_time'

    @staticmethod
    def _source_particle_ids(source, key, count):
        """Optional source labels: a metadata read cannot invalidate physical rows."""
        if key not in source:
            return None, False, 'field_absent'
        try:
            values = np.asarray(source[key][:])
            if values.ndim != 1 or len(values) != count:
                return None, False, 'source_label_shape_does_not_align_with_particle_rows'
            if values.dtype.kind == 'O':
                if not all(isinstance(value, (str, bytes)) for value in values):
                    return None, False, 'unsupported_object_source_label_dtype'
                values = np.asarray([value.decode('utf-8') if isinstance(value, bytes) else value
                                     for value in values], dtype=str)
            return values, True, 'recorded_source_labels'
        except (OSError, TypeError, ValueError, KeyError, RuntimeError) as error:
            return None, False, 'source_label_read_error: ' + str(error)

    def _attach_source_ids(self, frame):
        """Keep optional source arrays in their original non-object dtype."""
        blocks = self.source_particle_id_blocks
        complete = sum(len(block['values']) for block in blocks) == len(frame)
        same_dtype = bool(blocks) and all(block['values'].dtype == blocks[0]['values'].dtype for block in blocks)
        self.source_particle_ids = None
        if complete and same_dtype:
            labels = np.empty(len(frame), dtype=blocks[0]['values'].dtype)
            for block in blocks:
                labels[block['start']:block['stop']] = block['values']
            self.source_particle_ids = labels
            frame['source_particle_id'] = labels
            # Historical alias remains actual source metadata, never row labels.
            frame['id'] = labels

    def _configure_stellar_time_units(self, header):
        """Resolve physical time units from explicit input or the actual run parameters."""
        from pathlib import Path
        seconds_per_gyr = 365.25 * 86400 * 1.e9
        wind_tunnel = getattr(self, 'wind_tunnel', None)
        if wind_tunnel is not None:
            self.stellar_time_code_to_gyr = wind_tunnel.time_to_myr / 1000.
            self.stellar_time_unit_source = str(wind_tunnel.path)
        elif 'DSPH_STELLAR_TIME_CODE_TO_GYR' in os.environ:
            self.stellar_time_code_to_gyr = float(os.environ['DSPH_STELLAR_TIME_CODE_TO_GYR'])
            self.stellar_time_unit_source = 'explicit DSPH_STELLAR_TIME_CODE_TO_GYR'
        elif 'stellar_time_code_to_gyr' in self.__dict__:
            self.stellar_time_unit_source = getattr(self, 'stellar_time_unit_source', 'explicit loader time-unit parameter')
        elif 'UnitTime_in_s' in header:
            self.stellar_time_code_to_gyr = float(header['UnitTime_in_s']) / seconds_per_gyr
            self.stellar_time_unit_source = 'snapshot Header/UnitTime_in_s'
        else:
            parameter_units = {}
            folder = getattr(self, 'folder_path', None)
            if folder is not None:
                run = Path(folder)
                candidates = [directory / name for name in ('GZWJL.PARAM-usedvalues', 'GZWJL.PARAM')
                              for directory in (run, run.parent)]
                for path in candidates:
                    if not path.is_file():
                        continue
                    parameter_units = {}
                    for line in path.read_text(encoding='utf-8').splitlines():
                        parts = line.split()
                        if len(parts) >= 2 and parts[0] in ('UnitLength_in_cm', 'UnitVelocity_in_cm_per_s'):
                            parameter_units[parts[0]] = float(parts[1])
                    if len(parameter_units) == 2:
                        self.stellar_time_code_to_gyr = (parameter_units['UnitLength_in_cm'] /
                            parameter_units['UnitVelocity_in_cm_per_s'] / seconds_per_gyr)
                        self.stellar_time_unit_source = str(path)
                        break
            if not hasattr(self, 'stellar_time_unit_source'):
                self.stellar_time_unit_source = 'explicit canonical production length/velocity units'
        if not np.isfinite(self.stellar_time_code_to_gyr) or self.stellar_time_code_to_gyr <= 0:
            raise ValueError('Stellar time unit must be a finite positive number of Gyr per code time')

    def _stellar_time_factor(self, dataset, default_unit):
        annotated = next((dataset.attrs[key] for key in ('Units', 'Unit', 'units', 'unit')
                          if key in dataset.attrs), default_unit)
        if isinstance(annotated, bytes):
            annotated = annotated.decode('utf-8')
        unit = str(annotated).lower().replace(' ', '').replace('_', '')
        factors = {'gyr': 1., 'myr': 1.e-3, 'kyr': 1.e-6, 'yr': 1.e-9,
                   'year': 1.e-9, 'years': 1.e-9, 's': 1. / (365.25 * 86400 * 1.e9),
                   'sec': 1. / (365.25 * 86400 * 1.e9),
                   'codetime': self.stellar_time_code_to_gyr,
                   'codeunits': self.stellar_time_code_to_gyr}
        return factors.get(unit), str(annotated)

    def _stellar_age_fields(self, group, count):
        """Read recorded ages/formation times independently of type, ID or origin."""
        birth = np.full(count, np.nan, dtype=np.float32)
        age = np.full(count, np.nan, dtype=np.float32)
        formation = next((name for name in ('StellarFormationTime', 'BirthTime', 'FormationTime')
                          if name in group), None)
        explicit_age = next((name for name in ('StellarAge', 'StellarAges', 'Age', 'Ages')
                             if name in group), None)
        info = dict(formation_field=formation, age_field=explicit_age,
                    snapshot_time_code=float(self.time),
                    time_code_to_Gyr=float(self.stellar_time_code_to_gyr),
                    time_unit_source=getattr(self, 'stellar_time_unit_source', 'explicit loader time-unit parameter'))
        if formation is not None:
            factor, unit = self._stellar_time_factor(group[formation], 'code_time')
            info['formation_unit'] = unit
            if factor is not None:
                birth[:] = np.asarray(group[formation][:], dtype=np.float32) * (factor / self.stellar_time_code_to_gyr)
        if explicit_age is not None:
            factor, unit = self._stellar_time_factor(group[explicit_age], self.stellar_age_default_unit)
            info.update(age_unit=unit, method='recorded_age' if factor is not None else 'unsupported_recorded_age_unit')
            if factor is not None:
                age[:] = np.asarray(group[explicit_age][:], dtype=np.float32) * factor
        elif formation is not None:
            # Zero and negative formation times are actual recorded values;
            # negative times can describe formation before simulation start.
            age[:] = (float(self.time) - birth) * self.stellar_time_code_to_gyr
            info['method'] = 'snapshot_time_minus_recorded_formation_time'
        else:
            info['method'] = 'recorded_age_and_formation_time_absent'
        info['available_N'] = int(np.count_nonzero(np.isfinite(age) & (age >= 0)))
        return birth, age, info

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

        if star_columns:
            # Keep both field definitions available whenever stellar times are requested.
            star_columns = list(cls.STAR_FIELD_ORDER)

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
        target[slice_idx] = value

    @staticmethod
    def _announce_auto_boxsize_once(boxsize):
        global AUTO_BOXSIZE_PRINTED_LOCAL

        run_id = os.environ.get("DSPH_RUN_ID")
        if run_id:
            marker_dir = os.environ.get('TMPDIR') or os.path.join(os.getcwd(), '.dsph_tmp')
            os.makedirs(marker_dir, exist_ok=True)
            marker_path = os.path.join(marker_dir, f"dsph_boxsize_{run_id}.flag")
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
        wind_tunnel = getattr(self, 'wind_tunnel', None)

        fname, fname_base, fname_ext = resolved
        with h5py.File(fname, 'r') as first_file:
            header = first_file['Header'].attrs
            self.time = float(header['Time'])
            self._configure_stellar_time_units(header)
            npart_total = np.asarray(header['NumPart_Total'], dtype=np.int64)
            mass_table = np.asarray(header['MassTable'], dtype=np.float64)
            numfiles = int(header['NumFilesPerSnapshot'])
            flag_cooling = int(header['Flag_Cooling'])
            if wind_tunnel is not None:
                self.boxsize = float(header['BoxSize'])
                if not np.isfinite(self.boxsize) or self.boxsize <= 0:
                    raise ValueError('A wind tunnel snapshot needs a positive Header/BoxSize')

        mass_to_msun = wind_tunnel.mass_to_msun if wind_tunnel else 1e10
        energy_scale = wind_tunnel.velocity_to_kms ** 2 if wind_tunnel else 1.
        density_scale = 404.7
        if wind_tunnel:
            density_scale *= (mass_to_msun / 1e10) * (LENGTH_CODE_TO_KPC / wind_tunnel.length_to_kpc) ** 3

        total_particles = int(sum(int(npart_total[tp]) for tp in selected_types))
        gas_columns, star_columns = self._resolve_requested_fields(
            dftype,
            gas_fields=gas_fields,
            star_fields=star_fields,
        )

        df_data = {
            'analysis_row': np.arange(total_particles, dtype=np.uint64),
            'source_id_available': np.zeros(total_particles, dtype=bool),
            'x': np.zeros(total_particles, dtype=np.float32),
            'y': np.zeros(total_particles, dtype=np.float32),
            'z': np.zeros(total_particles, dtype=np.float32),
            'vx': np.zeros(total_particles, dtype=np.float32),
            'vy': np.zeros(total_particles, dtype=np.float32),
            'vz': np.zeros(total_particles, dtype=np.float32),
            'm': np.zeros(total_particles, dtype=np.float64),
            'tp': np.zeros(total_particles, dtype=np.int32),
        }

        if dftype >= 2:
            for col in gas_columns:
                df_data[col] = np.zeros(total_particles, dtype=np.float32)
        if dftype == 4:
            df_data['ags'] = np.zeros(total_particles, dtype=np.float32)
        if dftype >= 3:
            for col in star_columns:
                df_data[col] = np.full(total_particles, np.nan, dtype=np.float32)

        self.stellar_age_metadata = []
        self.source_id_metadata = []
        self.source_particle_id_blocks = []
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

                    # Read each vector dataset once. Column-wise HDF5 reads
                    # otherwise scan the same interleaved rows three times.
                    for dataset, columns in (
                        ('Coordinates', ('x', 'y', 'z')),
                        ('Velocities', ('vx', 'vy', 'vz')),
                    ):
                        vectors = group[dataset][:]
                        for axis, column in enumerate(columns):
                            df_data[column][slice_idx] = vectors[:, axis]
                        del vectors
                    labels, available, id_status = self._source_particle_ids(group, 'ParticleIDs', npart)
                    df_data['source_id_available'][slice_idx] = available
                    self.source_id_metadata.append(dict(particle_type=tp, status=id_status))
                    if labels is not None:
                        self.source_particle_id_blocks.append(dict(start=start, stop=stop, values=labels))
                    df_data['tp'][slice_idx] = tp

                    if mass_table[tp] > 0.0:
                        self._fill_constant(df_data['m'], slice_idx, mass_table[tp] * mass_to_msun)
                    else:
                        df_data['m'][slice_idx] = np.asarray(group['Masses'][:],dtype=np.float64) * mass_to_msun

                    if tp == 0 and dftype >= 2:
                        need_ugas = ('ugas' in gas_columns) or ('temp' in gas_columns)
                        ugas = group['InternalEnergy'][:] if need_ugas else None

                        need_ne = ('temp' in gas_columns)
                        need_nh = ('nh' in gas_columns)
                        ne = (group['ElectronAbundance'][:] if flag_cooling > 0 and 'ElectronAbundance' in group
                              else np.zeros(npart, dtype=np.float32)) if need_ne else None
                        if need_nh:
                            if flag_cooling > 0 and 'NeutralHydrogenAbundance' not in group:
                                raise ValueError("Snapshot lacks NeutralHydrogenAbundance required for H I statistics")
                            nh = (group['NeutralHydrogenAbundance'][:] if 'NeutralHydrogenAbundance' in group
                                  else np.zeros(npart, dtype=np.float32))
                        else:
                            nh = None

                        if 'ugas' in gas_columns:
                            df_data['ugas'][slice_idx] = ugas
                        if 'temp' in gas_columns:
                            df_data['temp'][slice_idx] = zim.unitconversion(
                                ugas if energy_scale == 1. else ugas * energy_scale, ne)
                        if 'rho' in gas_columns:
                            df_data['rho'][slice_idx] = group['Density'][:] * density_scale
                        if 'nh' in gas_columns:
                            df_data['nh'][slice_idx] = nh
                        if 'sfr' in gas_columns and 'StarFormationRate' in group:
                            df_data['sfr'][slice_idx] = group['StarFormationRate'][:]
                        if dftype == 4 and 'AGS-Softening' in group:
                            df_data['ags'][slice_idx] = group['AGS-Softening'][:]

                    if star_columns and tp in (2, 3, 4):
                        birth, age, info = self._stellar_age_fields(group, npart)
                        df_data['birth'][slice_idx] = birth
                        df_data['stellar_age_gyr'][slice_idx] = age
                        self.stellar_age_metadata.append(dict(particle_type=tp, **info))

                    offsets[tp] = stop

        if self.boxsize == 'Auto':
            self.boxsize = round(np.ptp(df_data['x']), -1)
            self._announce_auto_boxsize_once(self.boxsize)

        if wind_tunnel:
            # Keep raw periodic coordinates for the moving-frame transform.
            pass
        elif box:
            box_center = self.boxsize / 2
            df_data['x'] -= box_center
            df_data['y'] -= box_center
            df_data['z'] -= box_center
        else:
            print("Box centering skipped. Ensure coordinates are already centered if needed.")

        df = pd.DataFrame(df_data, copy=False)
        self._attach_source_ids(df)
        self.df = df
        return df, self.time

    def load_snapshot(self, dftype=3, box=True, gas_fields=None, star_fields=None, particle_types=None):
        """Load snapshot data into a DataFrame."""
        try:
            selected_types = self._resolve_particle_types(particle_types)
            direct_result = self._load_snapshot_direct_hdf5(
                dftype=dftype,
                box=box,
                gas_fields=gas_fields,
                star_fields=star_fields,
                particle_types=selected_types,
            )
            if direct_result is not None:
                return direct_result
            if getattr(self, 'wind_tunnel', None) is not None:
                raise FileNotFoundError(f'No wind tunnel HDF5 snapshot {self.snapshot_num} under {self.folder_path}')

            header_ptype = selected_types[0]
            hdr = self._readsnap_quiet(self.folder_path, self.snapshot_num, header_ptype, header_only=1)
            if hdr is None:
                raise ValueError(f"Failed to read snapshot {self.snapshot_num} from {self.folder_path}")

            self.time = hdr['time']
            self._configure_stellar_time_units(hdr)
            npart_total = hdr['npartTotal']
            total_particles = int(sum(int(npart_total[tp]) for tp in selected_types))

            gas_columns, star_columns = self._resolve_requested_fields(
                dftype,
                gas_fields=gas_fields,
                star_fields=star_fields,
            )

            df_data = {
                'analysis_row': np.arange(total_particles, dtype=np.uint64),
                'source_id_available': np.zeros(total_particles, dtype=bool),
                'x': np.zeros(total_particles, dtype=np.float32),
                'y': np.zeros(total_particles, dtype=np.float32),
                'z': np.zeros(total_particles, dtype=np.float32),
                'vx': np.zeros(total_particles, dtype=np.float32),
                'vy': np.zeros(total_particles, dtype=np.float32),
                'vz': np.zeros(total_particles, dtype=np.float32),
                'm': np.zeros(total_particles, dtype=np.float64),
                'tp': np.zeros(total_particles, dtype=np.int32),
            }

            if dftype >= 2:
                for col in gas_columns:
                    df_data[col] = np.zeros(total_particles, dtype=np.float32)
            if dftype == 4:
                df_data['ags'] = np.zeros(total_particles, dtype=np.float32)
            if dftype >= 3:
                for col in star_columns:
                    df_data[col] = np.full(total_particles, np.nan, dtype=np.float32)

            self.stellar_age_metadata = []
            self.source_id_metadata = []
            self.source_particle_id_blocks = []
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

                labels, available, id_status = self._source_particle_ids(data, 'id', npart)
                df_data['source_id_available'][slice_idx] = available
                self.source_id_metadata.append(dict(particle_type=tp, status=id_status))
                if labels is not None:
                    self.source_particle_id_blocks.append(dict(start=offset, stop=offset + npart, values=labels))
                df_data['x'][slice_idx] = pos[:, 0]
                df_data['y'][slice_idx] = pos[:, 1]
                df_data['z'][slice_idx] = pos[:, 2]
                df_data['vx'][slice_idx] = vel[:, 0]
                df_data['vy'][slice_idx] = vel[:, 1]
                df_data['vz'][slice_idx] = vel[:, 2]
                df_data['m'][slice_idx] = np.asarray(data['m'],dtype=np.float64) * 1e10
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

                if star_columns and tp in (2, 3, 4):
                    # readsnap's ``age`` key denotes StellarFormationTime, not age.
                    if 'age' in data:
                        formed = np.asarray(data['age'], dtype=np.float32)
                        df_data['birth'][slice_idx] = formed
                        df_data['stellar_age_gyr'][slice_idx] = (float(self.time) - formed) * self.stellar_time_code_to_gyr
                    self.stellar_age_metadata.append(dict(particle_type=tp,
                        method='readsnap_recorded_formation_time' if 'age' in data else 'recorded_time_absent',
                        time_code_to_Gyr=float(self.stellar_time_code_to_gyr)))

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
            self._attach_source_ids(df)
            self.df = df
            return df, self.time

        except Exception as e:
            raise RuntimeError(f"Error loading snapshot: {str(e)}")
