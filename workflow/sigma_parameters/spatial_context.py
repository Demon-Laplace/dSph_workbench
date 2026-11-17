"""Property-only snapshot context and spatial centre seeds; no particle-ID I/O."""
from pathlib import Path
import hashlib
import json
import shlex
import numpy as np
import h5py


def parameters(path):
    result={}
    for line in Path(path).read_text(encoding='utf-8',errors='replace').splitlines():
        fields=line.split('#',1)[0].split('%',1)[0].split()
        if len(fields)>1:result[fields[0]]=fields[1]
    return result


def merger_parameters(path):
    """Parse two source paths followed by twelve transform numbers."""
    rows=[]
    for line in Path(path).read_text(encoding='utf-8',errors='strict').splitlines():
        fields=shlex.split(line.split('#',1)[0],comments=False)
        if fields:rows.append(fields)
        if len(rows)==2:break
    if len(rows)!=2 or len(rows[0])!=2 or len(rows[1])!=12:
        raise ValueError('Unsupported merger parameter layout: '+str(path))
    v=np.asarray([float(x) for x in rows[1]],dtype=float)
    if not np.isfinite(v).all():raise ValueError('Nonfinite merger transform')
    return rows[0],v


def merger_rotation(v):
    """Legacy transform of the second input: Rx then Ry then Rz.

    ic_merge.c uses columns m[i,j] with m[i,j]=+sin for cyclic i,j,
    m[j,i]=-sin and legacy radians=degrees*0.017453.
    """
    rotation=np.eye(3)
    for axis in range(3):
        theta=float(v[9+axis])*.017453
        sine,cosine=np.sin(theta),np.cos(theta)
        i,j=(axis+1)%3,(axis+2)%3
        step=np.eye(3);step[i,i]=step[j,j]=cosine
        step[i,j]=sine;step[j,i]=-sine
        rotation=step@rotation
    return rotation


def initial_seed_from_mergers(run):
    run=Path(run)
    _,galaxy=merger_parameters(run/'mergegalaxy.par')
    _,cube=merger_parameters(run/'mergecube.par')
    rotation=merger_rotation(cube)
    return rotation@galaxy[:3]+cube[:3],dict(
        method='nominal second-galaxy origin transformed by the outer cube merger',
        galaxy_translation_raw=galaxy[:3].tolist(),cube_translation_raw=cube[:3].tolist(),
        cube_second_input_rotation=rotation.tolist(),
        exact_rule_source='/travail/xhuang/workspace/code/github-sync/workflow/ic_merge.c',
        caveat='This is a search seed; original internal IC centre need not equal its origin.')


def prepare_context(entry,output_dir):
    e=dict(entry);run=Path(e.get('run_dir',e.get('path'))).resolve()
    pp=run/'GZWJL.PARAM';used=run/'GZWJL.PARAM-usedvalues'
    effective=parameters(pp)
    if used.is_file():effective.update(parameters(used))
    if int(effective.get('ComovingIntegrationOn','0')):
        raise ValueError('Property-only context supports noncosmological snapshots')
    ul,um,uv=[float(effective[k]) for k in
      ('UnitLength_in_cm','UnitMass_in_g','UnitVelocity_in_cm_per_s')]
    if not all(np.isfinite(v) and v>0 for v in (ul,um,uv)):
        raise ValueError('Explicit valid physical units required')
    kpc,msun,gyr=3.0856775814913673e21,1.989e33,365.25*86400*1e9
    gi=float(effective.get('GravityConstantInternal','0'))
    units=dict(length_to_kpc=ul/kpc,velocity_to_kms=uv/1e5,mass_to_msun=um/msun,
      time_to_gyr=ul/uv/gyr,UnitLength_in_cm=ul,UnitMass_in_g=um,
      UnitVelocity_in_cm_per_s=uv,IAU_kpc_cm=kpc,adopted_Msun_g=msun,
      seconds_per_Gyr=gyr,GravityConstantInternal=gi,
      G_KPC_KMS2_MSUN=(gi*(ul/kpc)*(uv/1e5)**2/(um/msun) if gi else 6.672e-8*msun/(kpc*1e10)),
      unit_source=str(used if used.is_file() else pp))
    seed,seed_info=initial_seed_from_mergers(run)
    if e.get('existing_context'):
        old=json.loads(Path(e['existing_context']).read_text(encoding='utf-8'))
        candidate=np.asarray(old.get('initial_center_raw',[]),dtype=float)
        if candidate.shape==(3,) and np.isfinite(candidate).all() and Path(old['run_dir']).resolve()==run:
            seed=candidate
            seed_info=dict(method='existing geometric initial centre used only as a spatial search seed',
              path=e['existing_context'],source_fields_reused=['initial_center_raw','run_dir'],
              member_selection_reused=False)
    paths={str(s['snapshot']):s['paths'] for s in e['snapshots']}
    fingerprints={str(s['snapshot']):s['files'] for s in e['snapshots']}
    source_sha={}
    for p in (pp,used,run/'mergegalaxy.par',run/'mergecube.par'):
        if p.is_file():source_sha[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    return dict(contract_id='sigma_parameters_property_selection_v1',
      canonical_observer_id='astropy_native_direct_v1',model=e.get('model',run.name),
      run_id=e.get('run_id',run.name),run_dir=str(run),output_dir=str(Path(output_dir).resolve()),
      orbit_label=e.get('orbit_label',''),units=units,params=effective,config=e,
      snapshot_paths=paths,input_snapshot_fingerprints=fingerprints,
      initial_seed_raw=seed.tolist(),initial_seed_method=seed_info,
      center_method='local mass-weighted stellar shrinking sphere; no ID membership',
      initial_search_radius_kpc=30.,subsequent_search_radius_kpc=5.,
      source_sha256=source_sha,methods_path=e.get('methods_path'),jeans_path=e.get('jeans_path'),
      particle_IDs_read=False,origin_status='PROPERTY_SELECTION',
      population_definition='pre-existing birth<=0 versus formed birth>0; absent birth falls back to simulation stellar type')


def shrinking_center(p,m,seed,length_to_kpc,initial_radius_kpc=30.,
                     min_radius_kpc=.15,min_count=64,shrink=.8):
    """Return position in raw units, anchored to a local predicted seed.

    Never use a global mean that could jump to the Milky Way. Shrinking stops
    before losing the required sample. The returned point is the mass centre
    of the last supported local shrinking sphere.
    """
    p=np.asarray(p);m=np.asarray(m,dtype=float);center=np.asarray(seed,dtype=float).copy()
    if p.ndim!=2 or p.shape[1]!=3 or m.shape!=(len(p),) or center.shape!=(3,):
        raise ValueError('Invalid spatial-centre array dimensions')
    valid=np.isfinite(p).all(axis=1)&np.isfinite(m)&(m>0)
    radius=float(initial_radius_kpc)
    q=valid&(np.linalg.norm((p-center)*length_to_kpc,axis=1)<radius)
    if int(q.sum())<min_count:
        raise ValueError('Too few local stars around predicted dwarf centre; no global fallback')
    for _ in range(60):
        candidate=np.average(p[q],weights=m[q],axis=0)
        proposed=max(min_radius_kpc,radius*shrink)
        trial=valid&(np.linalg.norm((p-candidate)*length_to_kpc,axis=1)<proposed)
        center=candidate
        if proposed<=min_radius_kpc or int(trial.sum())<min_count:break
        radius=proposed;q=trial
    return center


def read_stellar_properties(ctx,n,types=(2,3,4),birth=True):
    """Load only stellar physical properties; do not access ParticleIDs."""
    pieces={k:[] for k in ('p','v','m','tp','birth','age_Gyr')}
    stats=[];digest=hashlib.sha256();stamp=None;age_sources=[]
    for name in ctx['snapshot_paths'][str(int(n))]:
        path=Path(name);st=path.stat()
        expected=ctx.get('input_snapshot_fingerprints',{}).get(str(int(n)),[])
        record=next((x for x in expected if Path(x['path']).resolve()==path.resolve()),None)
        if record and (st.st_size,st.st_mtime_ns)!=(int(record['size_bytes']),int(record['mtime_ns'])):
            raise RuntimeError('Source snapshot changed: '+str(path))
        stats.append((str(path),st.st_size,st.st_mtime_ns,str(st.st_dev),str(st.st_ino)))
        with h5py.File(str(path),'r') as f:
            header=f['Header'].attrs;time=float(header['Time'])*ctx['units']['time_to_gyr']
            if stamp is not None and time!=stamp:raise ValueError('Snapshot parts have differing times')
            stamp=time;table=header['MassTable']
            for tp in types:
                group='PartType'+str(tp)
                if group not in f:continue
                g=f[group];p=g['Coordinates'][:];v=g['Velocities'][:];count=len(p)
                mass=(g['Masses'][:] if 'Masses' in g else np.full(count,table[tp]))*ctx['units']['mass_to_msun']
                if p.shape!=(count,3) or v.shape!=(count,3):raise ValueError('Invalid phase arrays')
                if not all(np.isfinite(a).all() for a in (p,v,mass)):raise ValueError('Nonfinite phase arrays')
                for key,array in (('p',p),('v',v),('m',mass)):
                    pieces[key].append(array)
                    digest.update((str(path.resolve())+':'+group+':'+key+str(array.dtype)+str(array.shape)).encode())
                    digest.update(np.ascontiguousarray(array).tobytes())
                pieces['tp'].append(np.full(count,tp,dtype=np.int8))
                age,formed,age_info=stellar_age_properties(g,time,ctx['units'],count,
                    ctx.get('explicit_age_default_unit','code_time'))
                pieces['age_Gyr'].append(age)
                age_sources.append(dict(path=str(path),particle_type=int(tp),**age_info))
                digest.update(json.dumps(age_info,sort_keys=True).encode())
                digest.update(np.ascontiguousarray(age).tobytes())
                if birth:
                    pieces['birth'].append(formed);digest.update(np.ascontiguousarray(formed).tobytes())
        after=path.stat()
        if (after.st_size,after.st_mtime_ns)!=(st.st_size,st.st_mtime_ns):
            raise RuntimeError('Source changed while reading: '+str(path))
    if not pieces['p']:raise ValueError('No requested stellar properties')
    result={k:np.concatenate(v) for k,v in pieces.items() if v}
    result.update(time=stamp,source_stats=stats,source_fields_sha256=digest.hexdigest(),age_source=age_sources)
    return result


def _time_field_factor(dataset,units,default_unit):
    annotated=None
    for key in ('units','unit','Units','Unit'):
        if key in dataset.attrs:
            annotated=dataset.attrs[key]
            if isinstance(annotated,bytes):annotated=annotated.decode('utf-8')
            annotated=str(annotated).strip()
            break
    label=str(annotated if annotated is not None else default_unit).lower().replace(' ','').replace('_','')
    factors={'gyr':1.,'gigayear':1.,'gigayears':1.,'myr':1.e-3,'megayear':1.e-3,
             'yr':1.e-9,'year':1.e-9,'years':1.e-9,'kyr':1.e-6,
             's':1./units.get('seconds_per_Gyr',365.25*86400*1.e9),
             'sec':1./units.get('seconds_per_Gyr',365.25*86400*1.e9),
             'codetime':units['time_to_gyr'],'codeunits':units['time_to_gyr']}
    if label not in factors:
        return None,dict(unit=annotated,unit_source='unsupported_dataset_time_unit')
    return factors[label],dict(unit=annotated if annotated is not None else default_unit,
        unit_source='dataset_attribute' if annotated is not None else 'explicit_loader_default',
        factor_to_Gyr=factors[label])


def stellar_age_properties(group,time_Gyr,units,count,default_explicit_unit='code_time'):
    """Prefer explicit ages; else derive from the recorded formation-time field.

    A recorded formation time of zero is a real value in this noncosmological
    contract and yields age=current time. Missing fields, rather than old-star
    membership or any ID, determine age unavailability.
    """
    formed=np.full(count,np.nan,dtype=float)
    formation_name=next((key for key in ('StellarFormationTime','BirthTime','FormationTime') if key in group),None)
    formation_info={}
    if formation_name is not None:
        factor,formation_info=_time_field_factor(group[formation_name],units,'code_time')
        if factor is not None:formed=np.asarray(group[formation_name][:],dtype=float)*factor
    explicit=next((key for key in ('StellarAge','Age','StellarAges','Ages') if key in group),None)
    if explicit is not None:
        factor,info=_time_field_factor(group[explicit],units,default_explicit_unit)
        if factor is None:
            return np.full(count,np.nan),formed,dict(field=explicit,method='unavailable_explicit_age_unit',**info)
        age=np.asarray(group[explicit][:],dtype=float)*factor
        method='recorded_age_field'
    elif formation_name is not None and formation_info.get('factor_to_Gyr') is not None:
        age=float(time_Gyr)-formed
        info=formation_info;explicit=formation_name;method='snapshot_time_minus_recorded_formation_time_including_zero'
    else:
        return np.full(count,np.nan),formed,dict(field=formation_name,method='unavailable_missing_or_unsupported_age_field',
            reason='No explicit age field and no supported recorded formation-time field; no inference from particle IDs or population')
    if age.shape!=(count,):raise ValueError('Age field has an unsupported shape')
    invalid=~np.isfinite(age)|(age<0)
    age=age.copy();age[invalid]=np.nan
    return age,formed,dict(field=explicit,method=method,invalid_count=int(invalid.sum()),**info)


# Kept explicit for a caller that wants a short property-loader name.
load_stellar_properties=read_stellar_properties
