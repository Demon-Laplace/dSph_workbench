"""One allocation, sequential models and one 16-worker pool per model."""
import csv,gc,hashlib,json,os,re,sys,traceback
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import physical_helpers as C
import stellar_parameters as S
from gas_parameters import analyze_gas
from spatial_context import prepare_context

TASK=Path(__file__).resolve().parent
PLAN=json.loads((TASK/'plan.json').read_text(encoding='utf8'))
CTX=None;REFERENCE={};FIXED={};CURRENT_ENTRY=None
def now():return datetime.now(timezone.utc).isoformat()
def clean(obj):
 if isinstance(obj,dict):return {str(k):clean(v) for k,v in obj.items()}
 if isinstance(obj,(list,tuple)):return [clean(v) for v in obj]
 if isinstance(obj,np.ndarray):return clean(obj.tolist())
 if isinstance(obj,np.generic):return clean(obj.item())
 if isinstance(obj,float) and not np.isfinite(obj):return None
 return obj
def save(p,d):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name('.'+p.name+'.partial')
 tmp.write_text(json.dumps(clean(d),ensure_ascii=False,allow_nan=False,indent=2)+'\n',encoding='utf8');os.replace(tmp,p)
def append(p,d):
 with Path(p).open('a',encoding='utf8') as f:f.write(json.dumps(clean(d),ensure_ascii=False,allow_nan=False)+'\n');f.flush()
def fingerprint_matches(record,ctx,n):
 current=ctx['input_snapshot_fingerprints'][str(n)];old=record.get('input_snapshot_fingerprints',[])
 return bool(old) and {(str(Path(s['path']).resolve()),s['size_bytes'],s['mtime_ns']) for s in old}=={(str(Path(s['path']).resolve()),s['size_bytes'],s['mtime_ns']) for s in current}

def setup_context(e):
 import spatial_context as P
 out=Path(e['output_dir']);out.mkdir(parents=True,exist_ok=True)
 ctx=P.prepare_context(e,out);ctx['task_source_sha256']=PLAN['source_sha256'];ctx['task_contract_id']='sigma_spatial_parameters_v1'
 references={}
 if e.get('existing_epochs'):
  with Path(e['existing_epochs']).open(encoding='utf8') as f:
   for line in f:
    d=json.loads(line);references[int(d['snapshot'])]=d.get('row',d)
 seeds={};previous=np.asarray(ctx['initial_seed_raw']);velocity=None;lasttime=None
 for item in e['snapshots']:
  n=item['snapshot'];ref=references.get(n)
  if ref and fingerprint_matches(ref,ctx,n) and ref.get('center_raw'):
   center=np.asarray(ref['center_raw']);time=ref['time_Gyr']
   velocity=np.asarray(ref['center_velocity_native_kms']) if ref.get('center_velocity_native_kms') else None
  else:
   a=P.read_stellar_properties(ctx,n);time=a['time'];seed=previous
   if velocity is not None and lasttime is not None:
    seed=previous+velocity*(time-lasttime)*1.022712165/ctx['units']['length_to_kpc']
   center=P.shrinking_center(a['p'],a['m'],seed,ctx['units']['length_to_kpc'],initial_radius_kpc=30. if n==0 else 5.,min_radius_kpc=.3,min_count=128)
   r=np.linalg.norm((a['p']-center)*ctx['units']['length_to_kpc'],axis=1);q=r<.5
   if q.sum()<8:q=r<2
   velocity=np.average(a['v'][q],weights=a['m'][q],axis=0)*ctx['units']['velocity_to_kms']
   del a
  seeds[str(n)]=center.tolist();previous=center;lasttime=time
 ctx['center_seeds']=seeds
 save(out/(e['run_id']+'_sigma_parameters_context.json'),ctx)
 return ctx
def init_worker(ctx,reference):
 global CTX,REFERENCE
 CTX=ctx;REFERENCE=reference;C.CACHE.clear()

def worker(n):
 out=Path(CTX['output_dir']);p=out/(CTX['run_id']+'_sigma_parameters_s%04d.json'%n)
 try:
  if p.is_file():
   old=json.loads(p.read_text())
   if old.get('task_source_sha256')==PLAN['source_sha256'] and old.get('source_snapshot_fingerprints')==CTX['input_snapshot_fingerprints'][str(n)] and old.get('status') in ('COMPLETE','PARTIAL'):return dict(snapshot=n,status=old['status'],file=str(p),reused=True)
  ref=REFERENCE.get(n)
  if ref:
   ref=dict(ref,member_reference_matches=fingerprint_matches(ref,CTX,n))
  row,center,bulk=S.snapshot(CTX,n,ref)
  roster={}
  gasctx=dict(CTX,snapshot_host_center_raw=row['host_center_raw'],snapshot_Rgc_kpc=row['Rgc_kpc'])
  try:gas=analyze_gas(gasctx,n,center,bulk,roster,CTX['units'])
  except Exception as exc:gas=dict(status='FAILED',reason=repr(exc),traceback=traceback.format_exc())
  row['gas']=gas
  # Keep reusable old measurements with their own unchanged definitions.
  row['existing_reference']=ref if ref and ref['member_reference_matches'] else None
  G=CTX['units']['G_KPC_KMS2_MSUN'];gravity={}
  for rad in S.RADII:
   k=str(rad);gm=gas.get('gas',{}).get(k,{})
   massgas=gm.get('all',{}).get('Mgas')
   stellar=row['stellar_mass_spheres'][k]['all'];dm=row['local_DM_spheres_Msun'][k]
   total=stellar+dm+massgas if all(v is not None for v in (stellar,dm,massgas)) else None
   gravity[k]=dict(Mtotal_local_Msun=total,GM_over_r_kms2=G*total/rad if total is not None else None)
  row['gravity_spheres']=gravity
  row['HI_fractions']={}
  for label,rad in [('strict_r5',5.),('loose_r15',15.)]:
   hi=gas.get('gas',{}).get(str(rad),{}).get('all',{}).get('MHI')
   stars=row['stellar_mass_spheres'][str(rad)]['all']
   row['HI_fractions'][label]=dict(radius_3d_kpc=rad,MHI_Msun=hi,Mstar_Msun=stars,
     fraction=hi/(hi+stars) if hi is not None and hi+stars>0 else None,
     definition='MHI/(Mstar+MHI); both masses in the same3Dspatial aperture; no30kpc preselection mass included')
  gcpos=(center-np.asarray(row['host_center_raw']))*CTX['units']['length_to_kpc'];Rgc=np.linalg.norm(gcpos)
  # Enclosed-host mass and spherical tide proxy, supplied by certified old results only when present.
  host_components=[row.get('host_enclosed_stellar_mass_Msun'),row.get('host_enclosed_DM_mass_Msun'),gas.get('host_enclosed_gas_mass_Msun')]
  totalhost=sum(host_components) if all(v is not None for v in host_components) else None
  row['host_tide_proxy']=dict(status='supported' if totalhost is not None else 'partial_missing_host_component',
     host_enclosed_mass_Msun=totalhost,GM_over_Rgc3_kms2_kpc2=G*totalhost/Rgc**3 if totalhost is not None else None,
     definition='spherical enclosed-host GM/Rgc^3 proxy; not full tidal tensor',components_Msun=host_components)
  row.update(status='PARTIAL' if gas.get('status')!='COMPLETE' or CTX.get('origin_status')=='OLD_STARS_ONLY' else 'COMPLETE',
    task_source_sha256=PLAN['source_sha256'],source_snapshot_fingerprints=CTX['input_snapshot_fingerprints'][str(n)],finished_utc=now())
  save(p,row)
  return dict(snapshot=n,status=row['status'],file=str(p))
 except Exception as exc:
  result=dict(snapshot=n,status='FAILED',error=repr(exc),traceback=traceback.format_exc());save(p,result);return result

def scalar_flat(d,prefix='',result=None):
 if result is None:result={}
 for k,v in d.items():
  key=prefix+k
  if k in ('existing_reference','task_source_sha256','source_snapshot_fingerprints','source_stellar_fields_sha256','fixed_initial_selections'):continue
  if isinstance(v,dict):scalar_flat(v,key+'.',result)
  elif not isinstance(v,(list,tuple)):result[key]=v
 return result

def finalize(e,ctx,records):
 out=Path(e['output_dir']);stem=e['run_id']+'_sigma_parameters'
 rows=[]
 for s in e['snapshots']:
  p=out/(e['run_id']+'_sigma_parameters_s%04d.json'%s['snapshot'])
  if p.is_file():rows.append(json.loads(p.read_text(encoding='utf8')))
 full=out/(stem+'_epochs.jsonl')
 with full.open('w',encoding='utf8') as f:
  for r in rows:f.write(json.dumps(r,separators=(',',':'),ensure_ascii=False)+'\n')
 scalars=[scalar_flat(r) for r in rows];keys=sorted(set().union(*(r.keys() for r in scalars)))
 table=out/(stem+'_timeseries.csv')
 with table.open('w',newline='',encoding='utf8') as f:
  w=csv.DictWriter(f,keys);w.writeheader();w.writerows(scalars)
 valid=[r for r in rows if 'profiles' in r];arrays={}
 for pop in ('old','new','all'):
  for geom in ('sphere','projected'):
   for field in ('N','mass_Msun','sigma_los_raw_kms','sigma_los_Walker_kms','sigma_3d_kms','sigma_1d_equiv_kms','beta_sigma','velocity_covariance_kms2'):
    shape=(len(valid),len(S.EDGES)-1,3,3) if field=='velocity_covariance_kms2' else (len(valid),len(S.EDGES)-1)
    arr=np.full(shape,np.nan)
    for i,r in enumerate(valid):
     for j,item in enumerate(r['profiles'].get(pop+'_'+geom,[])):
      if item.get(field) is not None:arr[i,j]=item[field]
    arrays[pop+'_'+geom+'_'+field]=arr
 arrays.update(snapshot=np.array([r['snapshot'] for r in valid]),time_Gyr=np.array([r['time_Gyr'] for r in valid]),radial_edges_kpc=S.EDGES)
 profile=out/(stem+'_profiles.npz');np.savez_compressed(profile,**arrays)
 files=[full,table,profile,out/(stem+'_context.json')]
 meta=out/(stem+'_metadata.json')
 counts={s:sum(r.get('status')==s for r in rows) for s in ('COMPLETE','PARTIAL','FAILED')}
 metadata=dict(contract_id='sigma_spatial_parameters_v1',observer_id=C.OBSERVER,run_id=e['run_id'],model_directory=e['model_directory'],
   source=e['run_dir'],snapshot_count=len(e['snapshots']),snapshot_results=counts,particle_IDs_read=False,
   membership_status=ctx.get('origin_status'),source_parameters=e['selected_parameters'],input_documents=e.get('input_documents',{}),
   inputs=ctx.get('source_sha256'),source_code_sha256=PLAN['source_sha256'],units=ctx['units'],
   definitions=dict(LOS='per-star canonical heliocentric ray; equal stars ddof0',
     Walker='current spatial old stars within3D15kpc and forward R<5kpc;61directions;one common fit',
     sigma_3d='sqrt(trace covariance) after aperture mean subtraction; spatial ordered motions remain',
     sigma_1d_equiv='sqrt(trace covariance/3)',R='circular projected D0*tan(rho), kpc',r='3D dwarf-centric kpc',
     centre='R<0.5kpc',outskirt='1<R<2kpc',Re='direct circular half-mass radius of spatial pre-existing stars at each epoch',
     beta_sigma='spherical component variances about component means in nativeZ basis; beta_second saved separately',
     mass='all local particles by physical position; gas strict r<5 and loose r<15; no ancestry exclusions',
     existing_reference='historical results retained with original definitions; new primary fields use physical spatial membership',
     theoretical_sigma='existing Jeans predictions retained if input fingerprints match; not interpreted as causal heating',
     unavailable='null plus field availability; repeated ID never makes a physical field unavailable',
     particle_duplicates='same ID is not evidence of duplicate physical particles; no deduplication by ID',
     gas_membership='strict3D r<5kpc; loose3D r<15kpc; annulus5<=r<15kpc, regardless of source label'),
   delivery_status='calculated_with_explicit_gaps' if counts['PARTIAL'] or counts['FAILED'] else 'complete',
   deferred=['correlation/lag/event inference and PDFs are local follow-up','precise tides/binding/potential and local ordered-motion decomposition require deeper analysis'],created_utc=now())
 save(meta,metadata);files.append(meta)
 delivery=dict(run_id=e['run_id'],model_directory=e['model_directory'],status=metadata['delivery_status'],
   files=[dict(path=str(p),name=p.name,bytes=p.stat().st_size,sha256=C.sha(p)) for p in files if p.is_file()])
 # Prefix every exported filename with this run's unique identity.
 delivery['files']=[f for f in delivery['files'] if f['name'].startswith(e['run_id'])]
 return delivery,metadata

def main():
 deliveries=[];results=[]
 for i,e in enumerate(PLAN['entries']):
  print(now(),i+1,len(PLAN['entries']),e['run_id'],flush=True)
  try:
   ctx=setup_context(e);ref={}
   if e.get('existing_epochs'):
    with Path(e['existing_epochs']).open(encoding='utf8') as f:
     for line in f:
      d=json.loads(line);ref[int(d['snapshot'])]=d.get('row',d)
   records=[]
   with ProcessPoolExecutor(max_workers=PLAN['workers'],initializer=init_worker,initargs=(ctx,ref)) as pool:
    pending={pool.submit(worker,s['snapshot']):s['snapshot'] for s in e['snapshots']}
    for future in as_completed(pending):
     result=future.result();records.append(result)
     append(TASK/'epoch-results.jsonl',dict(run_id=e['run_id'],**result))
   delivery,meta=finalize(e,ctx,records);deliveries.append(delivery)
   results.append(dict(run_id=e['run_id'],status=meta['delivery_status'],counts=meta['snapshot_results']))
  except Exception as exc:
   error=dict(run_id=e['run_id'],status='FAILED',error=repr(exc),traceback=traceback.format_exc());results.append(error)
   out=Path(e['output_dir']);out.mkdir(parents=True,exist_ok=True);failure=out/(e['run_id']+'_sigma_parameters_failure.json');save(failure,error)
   deliveries.append(dict(run_id=e['run_id'],model_directory=e['model_directory'],files=[dict(path=str(failure),name=failure.name,bytes=failure.stat().st_size,sha256=C.sha(failure))]))
   print(error['traceback'],flush=True)
  save(TASK/'ready-deliveries.json',dict(deliveries=deliveries,updated_utc=now()))
  save(TASK/'status.json',dict(status='running',completed=i+1,total=len(PLAN['entries']),results=results))
  C.CACHE.clear();gc.collect()
 save(TASK/'status.json',dict(status='complete_with_gaps' if any(r['status']!='complete' for r in results) else 'complete',results=results,excluded=PLAN['excluded'],finished_utc=now()))
 save(TASK/'transfer-final.json',dict(terminal=True,final=True,expected_batches=1,completed_batches=1,finished_utc=now()))
if __name__=='__main__':main()
