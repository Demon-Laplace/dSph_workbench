"""Freeze current input identities, existing-cache routes and single-batch plan."""
import hashlib,json,re,shutil,sys
from pathlib import Path
ROOT=Path('/travail/xhuang');TASK=ROOT/'workspace/jobs/sigma-parameters'
BASE=ROOT/'workspace/analysis/canonical_model_comparison'
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def lines(p):
 d={}
 if p.is_file():
  for s in p.read_text(errors='replace').splitlines():
   t=s.split('%')[0].split('#')[0].split()
   if len(t)>1:d[t[0]]=t[1]
 return d
def main():
 TASK.mkdir(exist_ok=True)
 entries={}
 keep=('model','run_id','orbit_label','run_dir','snapshot_dir','dwarf_ini','mw_ini','cube_ini')
 for name in ('run_manifest.json','run_manifest_pre2047.json'):
  manifest=read(BASE/name)
  for e in manifest['runs']:entries[str(Path(e['run_dir']).resolve())]={k:e[k] for k in keep if k in e}
  del manifest
 for e in read(ROOT/'workspace/jobs/video-maintenance/plan.json')['entries']:
  source=str(Path(e['source']).resolve())
  if source not in entries:entries[source]=dict(model=e['model'],run_id=e['model'] if Path(source).parent==ROOT else e['model']+'__'+Path(source).parent.name,run_dir=source,snapshot_dir=e['snapshot_dir'],parameter_source=str(Path(source)/'GZWJL.PARAM-usedvalues'))
 active_snapshot_dirs=[]
 work=[];excluded=[]
 for source,e0 in entries.items():
  print('Register',e0.get('run_id',source),file=sys.stderr,flush=True)
  e=dict(e0);run=Path(source);out=Path(e.get('snapshot_dir',run/'output'))
  grouped={}
  for p in out.glob('snapshot_*.hdf5'):
   match=re.fullmatch(r'snapshot_(\d+)(?:\.(\d+))?\.hdf5',p.name)
   if match and p.is_file():grouped.setdefault(int(match[1]),[]).append(p.resolve())
  if not grouped or 0 not in grouped:
   excluded.append(dict(run_id=e.get('run_id',run.name),source=source,reason='no_current_snapshot000_or_no_snapshots'));continue
  snapshots=[]
  for n,paths in sorted(grouped.items()):
   single=[p for p in paths if re.fullmatch(r'snapshot_\d+\.hdf5',p.name)]
   paths=single[:1] if single else sorted(paths)
   files=[dict(path=str(p),size_bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in paths]
   snapshots.append(dict(snapshot=n,paths=[str(p) for p in paths],files=files))
  run_id=e.get('run_id',run.name)
  # The actual server top-level model directory determines the local archive.
  model_directory=run.relative_to(ROOT).parts[0]
  if not re.fullmatch(r'Fornax\d+[A-Za-z0-9_.-]*',model_directory):raise ValueError(model_directory)
  parameter=run/'GZWJL.PARAM-usedvalues'
  if not parameter.is_file():parameter=run/'GZWJL.PARAM'
  output=out/'analysis/sigma_spatial_parameters'
  e.update(run_id=run_id,model_directory=model_directory,run_dir=str(run),snapshot_dir=str(out),snapshots=snapshots,
   parameter_source=str(parameter),parameter_file=str(run/'GZWJL.PARAM'),output_dir=str(output),
   methods_path=str(TASK/'observation.py'),jeans_path=str(TASK/'jeans_reference.py'),selected_parameters=lines(parameter))
  if not e.get('dwarf_ini'):
   import shlex
   s=next(s for s in (run/'mergegalaxy.par').read_text().splitlines() if s.strip() and not s.startswith('#'))
   g=shlex.split(s);e['dwarf_ini']=str((run/g[1]).resolve());e['mw_ini']=str((run/g[0]).resolve())
   c=shlex.split(next(s for s in (run/'mergecube.par').read_text().splitlines() if s.strip() and not s.startswith('#')))
   e['cube_ini']=str((run/c[0]).resolve())
  cp=BASE/'results'/run_id/(run_id+'_context.json')
  ep=BASE/'results'/run_id/(run_id+'_epochs.jsonl')
  e['existing_context']=str(cp) if cp.is_file() else None
  e['existing_epochs']=str(ep) if ep.is_file() else None
  documents={}
  for parent in dict.fromkeys([run,run.parent]):
   for name in ('case.json','simulation-details.json','IC_PARAMETERS.json','MODEL_CHANGES.txt'):
    p=parent/name
    if p.is_file() and p.stat().st_size<2000000:
     documents[str(p)]=dict(sha256=sha(p),text=p.read_text(encoding='utf8',errors='replace'))
  for name in ('mergegalaxy.par','mergecube.par'):
   p=run/name
   if p.is_file():documents[str(p)]=dict(sha256=sha(p),text=p.read_text(encoding='utf8',errors='replace'))
  e['input_documents']=documents
  if run_id in ('Fornax5028','Fornax4001','Fornax5028__anchor139p4_orbitonly'):
   e.update(effective_sf_density_threshold_code=.00171,effective_sf_threshold_evidence=str(BASE/'environment_timing/Fornax_environment_method_audit.json'))
  work.append(e)
 work.sort(key=lambda e:(e['model_directory']!='Fornax3002',e['model_directory'],e['run_id']))
 plan=dict(contract_id='sigma_spatial_parameters_v1',scope='all available production runs; one batch; physical spatial selection; repeated/missing particle IDs do not fail computation',
   observer_id='astropy_native_direct_v1',entries=work,excluded=excluded,workers=16,
   local_root='F:/Fornax_data/Fornax处理',video_dependency=None,
   video_dependency_reason='All core fields are measured from stable raw snapshots or existing completed structured caches; media and incomplete elinfo are unnecessary.',
   source_sha256={p.name:sha(p) for p in TASK.glob('*.py')},
   particle_IDs_read=False,gas_membership=dict(strict_radius_3d_kpc=5,loose_radius_3d_kpc=15),
   particle_age='read recorded physical age/formation time, independent of IDs')
 (TASK/'plan.json').write_text(json.dumps(plan,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
 (TASK/'ready-deliveries.json').write_text('{"deliveries":[]}\n')
 print(json.dumps(dict(runs=len(work),snapshots=sum(len(e['snapshots']) for e in work),reuse_contexts=sum(bool(e['existing_context']) for e in work),excluded=excluded),indent=2))
if __name__=='__main__':main()
