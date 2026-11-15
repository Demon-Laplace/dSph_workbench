"""Stream standalone MW snapshots into versioned radial gas/mass archives.

No snapshot, IC, launcher or existing analysis is modified. All PartType0 gas
is retained in the radial statistics, including ambient cube gas. ``CGM`` is
therefore a radial diagnostic, not a guessed initial-particle membership.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor
import contextlib
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import time
import h5py
import numpy as np

SCHEMA = 'mw_cgm_shell_history_v1'
PHASES = ('all_gas', 'T_lt_20000K', 'T_ge_20000K')
GAS_FIELDS = ('mass', 'm_temperature', 'm_temperature2', 'm_density_cgs',
              'm_density_cgs2', 'm_pressure_cgs', 'volume_kpc3',
              'm_vr', 'm_vr2', 'm_vt2', 'm_vx', 'm_vy', 'm_vz',
              'inflow_m_vr', 'outflow_m_vr', 'HI_mass', 'm_internal_energy',
              'm_electron_abundance', 'm_metallicity', 'density_valid_mass')
MOMENTS = ('xx', 'xy', 'xz', 'yy', 'yz', 'zz')
RADIAL_EDGES = np.unique(np.r_[0., np.geomspace(0.1, 2000., 256),
                              1., 5., 10., 20., 30., 50., 100., 150., 200., 250., 300., 500., 1000., 1500.])
TEMP_EDGES = np.unique(np.r_[0., np.geomspace(1., 1.e10, 51), 2.e4, 1.e5, 1.e6, np.inf])

def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8*1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def fingerprint(path):
    s = Path(path).stat()
    return dict(bytes=s.st_size, mtime_ns=s.st_mtime_ns, inode=s.st_ino)


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def read_parameters(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        fields = line.split('%', 1)[0].split('#', 1)[0].split()
        if len(fields) >= 2:
            result[fields[0]] = fields[1]
    return result


def units(parameters):
    """Match current production analysis's documented legacy unit convention."""
    from canonical_observer import LENGTH_CODE_TO_KPC, VELOCITY_CODE_TO_KMS
    expected = dict(UnitLength_in_cm=3.085678e21, UnitMass_in_g=1.989e43,
                    UnitVelocity_in_cm_per_s=1.e5)
    recorded = {}
    for name, value in expected.items():
        recorded[name] = float(parameters[name])
        # GIZMO's parameters-usedvalues prints %g: the standard length becomes
        # 3.08568e21. Accept that exact serialization, not arbitrary near units.
        if recorded[name] not in (value, float(format(value, '.6g'))):
            raise ValueError('Unsupported units; an explicit reviewed conversion is required: '+name)
    if int(parameters['ComovingIntegrationOn']) != 0:
        raise ValueError('Archive currently supports non-comoving MW runs only')
    return dict(length_to_kpc=LENGTH_CODE_TO_KPC, velocity_to_kms=VELOCITY_CODE_TO_KMS,
                mass_to_Msun=1.e10, density_to_cgs=expected['UnitMass_in_g']/expected['UnitLength_in_cm']**3,
                time_to_Gyr=expected['UnitLength_in_cm']/expected['UnitVelocity_in_cm_per_s']/31557600.e9,
                velocity_kms_to_kpc_yr=1.022712165045695e-9,
                X_H=0.76, gamma=5./3., boltzmann_cgs=1.3806e-16, proton_mass_g=1.6726e-24,
                native_units=expected, recorded_unit_values=recorded,
                unit_serialization='exact standard units or their six-significant-digit usedvalues serialization',
                convention='current_production_legacy_units_v1')


def unwrap_stellar_positions(positions, weights, selected, box):
    """Choose a periodic reference without averaging stars across opposite faces."""
    positions=np.asarray(positions,dtype=float);weights=np.asarray(weights,dtype=float)
    phase=np.exp((2j*np.pi/box)*positions[selected])
    circular=np.average(phase,axis=0,weights=weights[selected])
    concentration=np.abs(circular)
    if np.any(concentration<1.e-6):
        raise ValueError('Stellar periodic reference is unconstrained')
    reference=(np.angle(circular)%(2*np.pi))*box/(2*np.pi)
    delta=positions-reference
    delta-=box*np.rint(delta/box)
    return reference+delta,reference,concentration


def host_frame(path):
    """Use the existing TP2/4 shrinking centre and TP2/3/4 R5 bulk exactly."""
    from basefunc import GalaxySimulation
    path = Path(path)
    match = re.fullmatch(r'snapshot_(\d+)\.hdf5', path.name)
    if not match:
        raise ValueError('Expected a single-file HDF5 snapshot: '+str(path))
    simulation = GalaxySimulation(str(path.parent), int(match[1]))
    with contextlib.redirect_stdout(io.StringIO()):
        simulation.load_snapshot(dftype=1, box=False, gas_fields=(), star_fields=(), particle_types=(2,3,4))
        with h5py.File(path,'r') as snapshot:box=float(snapshot['Header'].attrs['BoxSize'])
        df=simulation.df
        unwrapped,reference,concentration=unwrap_stellar_positions(
            df[['x','y','z']].to_numpy(),df['m'].to_numpy(),
            np.isin(df['tp'].to_numpy(),(2,4)),box)
        for axis,name in enumerate(('x','y','z')):df[name]=unwrapped[:,axis]
        simulation.center_on_MW()
    frame=simulation.observer_host_state
    frame['unwrapped_host_center_input_code']=frame['host_center_input_code']
    frame['host_center_input_code']=(np.asarray(frame['host_center_input_code'])%box).tolist()
    frame.update(periodic_reference_input_code=reference.tolist(),periodic_concentration=concentration.tolist(),
                 centering_method='periodic_circular_unwrap_then_existing_TP24_shrink_R5_bulk')
    return frame


def temperature(internal_energy, electrons):
    # Same function/constants as the current loader; not an ad-hoc temperature cut.
    from zim import unitconversion
    return unitconversion(internal_energy, electrons)


def add_bins(target, bins, weights, nbin):
    target += np.bincount(bins, weights=weights, minlength=nbin)


def reduce_snapshot(path, output, parameter_values, chunk_size=262144, frame=None):
    path, output = Path(path), Path(output)
    before = fingerprint(path)
    scale = units(parameter_values)
    frame = host_frame(path) if frame is None else frame
    center = np.asarray(frame['host_center_input_code'], dtype=float)
    bulk = np.asarray(frame['host_bulk_input_code'], dtype=float)
    nr, nt = len(RADIAL_EDGES)-1, len(TEMP_EDGES)-1
    data = dict(radial_edges_kpc=RADIAL_EDGES, temperature_edges_K=TEMP_EDGES,
                mass_by_type_Msun=np.zeros((6,nr)), count_by_type=np.zeros((6,nr), dtype=np.int64),
                mass_second_moments_Msun_kpc2=np.zeros((6,nr,6)),
                outer_quadrupole_Msun_per_kpc3=np.zeros((6,nr,6)),
                angular_momentum_Msun_kpc_kms=np.zeros((6,nr,3)),
                gas_sums=np.zeros((3,nr,len(GAS_FIELDS))), gas_count=np.zeros((3,nr), dtype=np.int64),
                gas_temperature_histogram_Msun=np.zeros((nr,nt)))
    native_mass = np.zeros(6)
    field_hashes = {}
    def read(group, name, start, stop):
        values = group[name][start:stop]
        key = group.name+'/'+name
        field_hashes.setdefault(key, hashlib.sha256()).update(np.ascontiguousarray(values).tobytes())
        return values
    with h5py.File(path, 'r') as snapshot:
        header = snapshot['Header'].attrs
        if int(header['NumFilesPerSnapshot']) != 1:
            raise ValueError('Multipart snapshot must be handled explicitly before retirement')
        counts = np.asarray(header['NumPart_Total'], dtype=np.int64)
        mass_table = np.asarray(header['MassTable'], dtype=float)
        box = float(header['BoxSize'])
        if not np.isfinite(box) or box <= 0:
            raise ValueError('A positive periodic box size is required')
        time_code = float(header['Time'])
        for kind in range(6):
            if not counts[kind]:
                continue
            group = snapshot['PartType%d'%kind]
            if len(group['Coordinates']) != counts[kind]:
                raise ValueError('Particle count disagrees with Header')
            for start in range(0, int(counts[kind]), chunk_size):
                stop = min(start+chunk_size, int(counts[kind]))
                pos = read(group,'Coordinates',start,stop)
                vel = read(group,'Velocities',start,stop)
                read(group,'ParticleIDs',start,stop)
                if mass_table[kind] > 0:
                    m = np.full(stop-start, mass_table[kind], dtype=float)
                else:
                    m = np.asarray(read(group,'Masses',start,stop),dtype=float)
                native_mass[kind] += m.sum()
                m *= scale['mass_to_Msun']
                # Bounded chunks only; float64 products/accumulators retain small velocity moments.
                delta = np.asarray(pos,dtype=float)-center
                delta -= box*np.rint(delta/box)
                delta *= scale['length_to_kpc']
                velocity = (np.asarray(vel,dtype=float)-bulk)*scale['velocity_to_kms']
                radius2 = np.einsum('ij,ij->i',delta,delta)
                radius = np.sqrt(radius2)
                if (not np.isfinite(delta).all() or not np.isfinite(velocity).all()
                        or not np.isfinite(m).all() or np.any(m<=0)):
                    raise ValueError('Invalid position/velocity/particle mass')
                bins = np.searchsorted(RADIAL_EDGES, radius, side='right')-1
                if np.any(bins<0) or np.any(bins>=nr):
                    raise ValueError('Radial archive grid does not cover every particle')
                add_bins(data['mass_by_type_Msun'][kind], bins, m, nr)
                data['count_by_type'][kind] += np.bincount(bins,minlength=nr)
                for axis,(a,b) in enumerate(((0,0),(0,1),(0,2),(1,1),(1,2),(2,2))):
                    product = delta[:,a]*delta[:,b]
                    add_bins(data['mass_second_moments_Msun_kpc2'][kind,:,axis],bins,m*product,nr)
                    quadrupole = m*(3*product-(radius2 if a==b else 0.))
                    outer = np.divide(quadrupole,radius**5,out=np.zeros_like(m),where=radius>0)
                    add_bins(data['outer_quadrupole_Msun_per_kpc3'][kind,:,axis],bins,outer,nr)
                angular = np.cross(delta, velocity)
                for axis in range(3):
                    add_bins(data['angular_momentum_Msun_kpc_kms'][kind,:,axis],bins,m*angular[:,axis],nr)
                if kind != 0:
                    continue
                rho = np.asarray(read(group,'Density',start,stop),dtype=float)*scale['density_to_cgs']
                energy = np.asarray(read(group,'InternalEnergy',start,stop),dtype=float)
                if 'ElectronAbundance' in group:
                    electrons = np.asarray(read(group,'ElectronAbundance',start,stop),dtype=float)
                elif int(header['Flag_Cooling']) == 0:
                    electrons = np.zeros_like(energy)
                else:
                    raise ValueError('Cooling snapshot lacks electron abundance')
                temp = temperature(energy, electrons)
                if (not np.isfinite(temp).all() or np.any(temp<0) or not np.isfinite(rho).all() or np.any(rho<0)):
                    raise ValueError('Invalid gas thermodynamic values')
                vr = np.divide(np.einsum('ij,ij->i',delta,velocity),radius,
                               out=np.zeros_like(radius),where=radius>0)
                vt2 = np.maximum(0.,np.einsum('ij,ij->i',velocity,velocity)-vr*vr)
                density_valid = rho>0
                rho_msun_kpc3 = rho/scale['density_to_cgs']*scale['mass_to_Msun']/scale['length_to_kpc']**3
                volume = np.divide(m,rho_msun_kpc3,out=np.zeros_like(m),where=density_valid)
                neutral = np.asarray(read(group,'NeutralHydrogenAbundance',start,stop),dtype=float) if 'NeutralHydrogenAbundance' in group else np.zeros_like(m)
                metallicity = read(group,'Metallicity',start,stop) if 'Metallicity' in group else np.zeros_like(m)
                if metallicity.ndim==2:metallicity=metallicity[:,0]
                values = (m,m*temp,m*temp*temp,m*rho,m*rho*rho,
                          m*((scale['gamma']-1)*rho*energy*1.e10),volume,
                          m*vr,m*vr*vr,m*vt2,m*velocity[:,0],m*velocity[:,1],m*velocity[:,2],
                          m*np.minimum(vr,0.),m*np.maximum(vr,0.),m*scale['X_H']*neutral,
                          m*energy,m*electrons,m*metallicity,m*density_valid)
                for phase,selection in enumerate((np.ones(len(m),dtype=bool),temp<2.e4,temp>=2.e4)):
                    selected_bins=bins[selection]
                    data['gas_count'][phase] += np.bincount(selected_bins,minlength=nr)
                    for field,weights in enumerate(values):
                        add_bins(data['gas_sums'][phase,:,field],selected_bins,weights[selection],nr)
                tbins=np.searchsorted(TEMP_EDGES,temp,side='right')-1
                histogram=np.bincount(bins*nt+tbins,weights=m,minlength=nr*nt).reshape(nr,nt)
                data['gas_temperature_histogram_Msun'] += histogram
        header_saved = {k:np.asarray(v).tolist() for k,v in header.items()}
    if fingerprint(path)!=before:
        raise ValueError('Snapshot changed while being reduced')
    total = data['mass_by_type_Msun'].sum(axis=1)
    np.testing.assert_allclose(total,native_mass*scale['mass_to_Msun'],rtol=2.e-12,atol=1.e-6)
    np.testing.assert_array_equal(data['count_by_type'].sum(axis=1),counts)
    np.testing.assert_allclose(data['gas_sums'][0,:,0],data['gas_sums'][1:,:,0].sum(axis=0),rtol=2.e-12,atol=1.e-6)
    np.testing.assert_allclose(data['gas_temperature_histogram_Msun'].sum(axis=1),data['gas_sums'][0,:,0],rtol=2.e-12,atol=1.e-6)
    for name,value in data.items():
        if name!='temperature_edges_K' and not np.isfinite(value).all():
            raise ValueError('Non-finite additive archive field: '+name)
    metadata = dict(schema=SCHEMA,source=str(path),source_fingerprint=before,header=header_saved,
                    source_module_sha256=sha(Path(__file__)),
                    scale=scale,host_frame=frame,phase_labels=PHASES,gas_sum_fields=GAS_FIELDS,
                    tensor_components=MOMENTS,gas_selection='all PartType0; no guessed source membership',
                    thermodynamics='zim.unitconversion; X_H=0.76; gamma=5/3; hot split 20000 K',
                    boundary='nearest periodic image relative to per-snapshot stellar centre',
                    field_sha256={k:v.hexdigest() for k,v in field_hashes.items()},mass_balance_passed=True)
    data.update(time_code=np.asarray(time_code),time_Gyr=np.asarray(time_code*scale['time_to_Gyr']),
                host_center_raw=center,host_bulk_raw=bulk,native_total_mass_code=native_mass,
                metadata_json=np.asarray(json.dumps(metadata,allow_nan=False)))
    output.parent.mkdir(parents=True,exist_ok=True)
    pending=output.with_name(output.name+'.tmp.npz')
    np.savez_compressed(pending,**data);pending.replace(output)
    return dict(path=str(output),sha256=sha(output),source_fingerprint=before,time_code=time_code,
                time_Gyr=time_code*scale['time_to_Gyr'],mass_balance_passed=True)


def worker(task):
    run,source,destination=task
    os.chdir(run['model_path'])
    parameter_values=read_parameters(run['parameter_file'])
    destination=Path(destination)
    if destination.exists():
        with np.load(destination,allow_pickle=False) as saved:
            metadata=json.loads(str(saved['metadata_json']))
            if (metadata['schema']==SCHEMA and metadata['source_fingerprint']==fingerprint(source)
                    and metadata['mass_balance_passed']):
                return dict(path=str(destination),sha256=sha(destination),time_code=float(saved['time_code']),
                            time_Gyr=float(saved['time_Gyr']),mass_balance_passed=True)
        raise ValueError('An existing archive has a different source identity')
    try:
        return reduce_snapshot(source,destination,parameter_values)
    except Exception as error:
        raise RuntimeError(str(source)+': '+str(error)) from error


def check_reusable_frame(metadata):
    source=metadata['source'];os.chdir(Path(source).parent.parent)
    if fingerprint(source)!=metadata['source_fingerprint']:
        raise ValueError('Reusable archive source changed: '+source)
    frame=host_frame(source);old=metadata['host_frame'];box=float(metadata['header']['BoxSize'])
    delta=np.asarray(frame['host_center_input_code'])-old['host_center_input_code']
    delta-=box*np.rint(delta/box)
    dr=float(np.max(np.abs(delta)))
    dv=float(np.max(np.abs(np.asarray(frame['host_bulk_input_code'])-old['host_bulk_input_code'])))
    return dict(source=source,center_difference_code=dr,bulk_difference_code=dv,equivalent=dr<1.e-8 and dv<1.e-8)


def combine_run(run,outputs,destination):
    """One compact time series; keep per-frame metadata/source hashes in a manifest."""
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    arrays={};metadata=[];fixed=None
    for record in outputs:
        with np.load(record['path'],allow_pickle=False) as frame:
            current={k:frame[k].copy() for k in ('radial_edges_kpc','temperature_edges_K')}
            if fixed is None:fixed=current
            else:
                for key in fixed:np.testing.assert_array_equal(current[key],fixed[key])
            metadata.append(json.loads(str(frame['metadata_json'])))
            for key in frame.files:
                if key not in (*current,'metadata_json'):arrays.setdefault(key,[]).append(frame[key].copy())
    combined={key:np.stack(values) for key,values in arrays.items()}
    combined.update(fixed)
    combined['schema_id']=np.asarray(SCHEMA)
    combined['gas_sum_fields']=np.asarray(GAS_FIELDS)
    combined['phase_labels']=np.asarray(PHASES)
    path=destination/(run['key']+'_cgm_shells.npz')
    np.savez_compressed(path,**combined)
    summary=destination/(run['key']+'_stability.csv')
    with summary.open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(('snapshot','time_Gyr','gas_mass_Msun','cold_gas_mass_Msun','hot_gas_mass_Msun','HI_mass_Msun','mass_r_ge_10kpc_Msun','host_x_code','host_y_code','host_z_code'))
        nr=len(RADIAL_EDGES)-1
        for i,record in enumerate(outputs):
            gas=combined['gas_sums'][i]
            # Complete shells outside 10 kpc; exact shell edges are retained in NPZ.
            outside=RADIAL_EDGES[:-1]>=10.
            writer.writerow((Path(metadata[i]['source']).name,record['time_Gyr'],gas[0,:,0].sum(),gas[1,:,0].sum(),gas[2,:,0].sum(),gas[0,:,15].sum(),gas[0,outside,0].sum(),*combined['host_center_raw'][i]))
    # A plain schema + full source identities survive even when snapshots are retired.
    figure=destination/(run['key']+'_cgm_stability.png')
    plot_history(combined,figure,run['key'])
    manifest=destination/(run['key']+'_cgm_manifest.json')
    document=dict(schema=SCHEMA,run=run,snapshots=len(outputs),snapshot_metadata=metadata,
                  output_files=[dict(name=p.name,bytes=p.stat().st_size,sha256=sha(p)) for p in (path,summary,figure)],
                  limitations='Binned moments cannot restore individual particles or refine radial/temperature bins after raw snapshot removal.')
    write_json(manifest,document)
    return dict(run=run['key'],files=[str(path),str(summary),str(figure),str(manifest)],snapshots=len(outputs),passed=True)


def plot_history(data,path,title):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    times=data['time_Gyr'];edges=data['radial_edges_kpc'];radius=(edges[1:]+edges[:-1])/2
    volume=4*np.pi/3*np.diff(edges**3)
    gas=data['gas_sums'][:,0]
    mass=gas[:,:,0]
    temp=np.divide(gas[:,:,1],mass,out=np.full_like(mass,np.nan),where=mass>0)
    vr=np.divide(gas[:,:,7],mass,out=np.full_like(mass,np.nan),where=mass>0)
    fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    view=(radius>=1)&(radius<=500)
    for axis,values,label,norm,cmap in (
        (axes[0,0],mass/volume,'Mean gas density [Msun/kpc^3]',LogNorm(),'viridis'),
        (axes[0,1],temp,'Mass-weighted temperature [K]',LogNorm(),'inferno'),
        (axes[1,0],vr,'Mass-weighted radial velocity [km/s]',None,'coolwarm')):
        kwargs=dict(shading='nearest',cmap=cmap)
        if norm is not None:kwargs['norm']=norm
        mesh=axis.pcolormesh(times,radius[view],values[:,view].T,**kwargs)
        axis.set_yscale('log');axis.set(xlabel='Time [Gyr]',ylabel='Radius [kpc]')
        fig.colorbar(mesh,ax=axis,label=label)
    for cutoff in (10.,50.,150.,300.):
        inside=edges[1:]<=cutoff
        axes[1,1].plot(times,mass[:,inside].sum(axis=1),label='r < %g kpc'%cutoff)
    axes[1,1].set(xlabel='Time [Gyr]',ylabel='Enclosed gas mass [Msun]');axes[1,1].legend()
    fig.suptitle(title+'; all gas including ambient cube')
    fig.savefig(path,dpi=140);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=16)
    parser.add_argument('--validate-first',action='store_true',help='Compare two chunk sizes on the first real snapshot before the full batch')
    args=parser.parse_args();plan=json.loads(args.plan.read_text());root=Path(plan['archive_root'])
    if args.validate_first and not (root/'representative-validation.json').exists():
        run=plan['runs'][0];source=run['snapshots'][0]['path'];os.chdir(run['model_path'])
        values=read_parameters(run['parameter_file']);paths=[]
        for size in (131071,262144):
            path=root/('validation-%d.npz'%size)
            reduce_snapshot(source,path,values,chunk_size=size);paths.append(path)
        with np.load(paths[0],allow_pickle=False) as a,np.load(paths[1],allow_pickle=False) as b:
            for key in a.files:
                if key=='metadata_json':continue
                left,right=a[key],b[key]
                if key=='temperature_edges_K':np.testing.assert_array_equal(left,right);continue
                scale=max(float(np.max(np.abs(left))),1.e-100)
                np.testing.assert_allclose(left,right,rtol=2.e-10,atol=scale*2.e-12)
            assert json.loads(str(a['metadata_json']))['field_sha256']==json.loads(str(b['metadata_json']))['field_sha256']
        write_json(root/'representative-validation.json',dict(passed=True,source=source,chunk_sizes=[131071,262144],
                   additive_arrays_match=True,source_field_hashes_match=True))
        print('REPRESENTATIVE_VALIDATION_PASSED',source,flush=True)
    completed=[];pending=list(plan['runs'])
    for reuse in plan.get('reuse_completed_runs',[]):
        manifest=json.loads(Path(reuse['manifest']).read_text())
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            checks=list(pool.map(check_reusable_frame,manifest['snapshot_metadata']))
        passed=all(row['equivalent'] for row in checks)
        write_json(root/(reuse['result']['run']+'_reuse-validation.json'),dict(passed=passed,checks=checks))
        if passed:
            completed.append(dict(reuse['result'],reused_after_periodic_center_check=True))
            print('REUSED_COMPLETE_RUN',reuse['result']['run'],len(checks),flush=True)
        else:
            pending.insert(0,reuse['original_run'])
            print('CHANGED_CENTER_REQUIRES_REDUCTION',reuse['result']['run'],flush=True)
    for run in pending:
        started=time.monotonic();tasks=[]
        for source in run['snapshots']:
            target=root/'frames'/run['key']/(Path(source['path']).stem+'.npz')
            if fingerprint(source['path'])!={k:source[k] for k in ('bytes','mtime_ns','inode')}:raise ValueError('Frozen source changed')
            tasks.append((run,source['path'],str(target)))
        print('RUN_START',run['key'],len(tasks),flush=True)
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            outputs=[]
            for record in pool.map(worker,tasks):
                outputs.append(record)
                if len(outputs)%10==0:print('FRAMES',run['key'],len(outputs),len(tasks),flush=True)
        result=combine_run(run,outputs,root/'results'/run['model'])
        result['seconds']=time.monotonic()-started;completed.append(result)
        write_json(root/'progress.json',dict(runs=completed))
        print('RUN_COMPLETE',run['key'],len(outputs),result['seconds'],flush=True)
    write_json(root/'archive-completed.json',dict(schema=SCHEMA,passed=True,runs=completed,
               snapshots=sum(r['snapshots'] for r in completed),deletion_performed=False))
    print('ARCHIVE_COMPLETE',flush=True)


if __name__=='__main__':main()
