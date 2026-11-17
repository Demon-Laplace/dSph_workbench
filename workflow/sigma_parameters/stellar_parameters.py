"""Same-observer, explicitly named population/aperture moments; no rendering."""
import gc,hashlib,importlib.util,json,os
from pathlib import Path
import numpy as np
import physical_helpers as C

CONTRACT='sigma_spatial_parameters_v1'
EDGES=np.array([0.,.25,.5,1.,1.5,2.,3.,5.])
RADII=np.array([.25,.3,.5,1.,1.5,2.,3.,5.,15.])

def moments(row_labels,x,v,m,los,resid,q,axis):
 n=int(q.sum());out=dict(N=n,mass_Msun=float(m[q].sum()))
 if n<2:return dict(out,supported=False,reason='fewer_than_two_stars')
 z=v[q].astype('f8');xx=x[q].astype('f8');mass=m[q].astype('f8');a=los[q];b=resid[q]
 avg=z.mean(axis=0);dv=z-avg;cov=dv.T@dv/n
 eig,vec=np.linalg.eigh(cov)
 sph_r=np.linalg.norm(xx,axis=1);er=np.divide(xx,sph_r[:,None],out=np.zeros_like(xx),where=sph_r[:,None]>0)
 ep=np.c_[-er[:,1],er[:,0],np.zeros(n)];norm=np.linalg.norm(ep,axis=1)
 ep=np.divide(ep,norm[:,None],out=np.zeros_like(ep),where=norm[:,None]>0);et=np.cross(er,ep)
 sph=np.c_[np.sum(z*er,axis=1),np.sum(z*et,axis=1),np.sum(z*ep,axis=1)]
 vari=np.var(sph,axis=0);sec=np.mean(sph*sph,axis=0)
 w=mass/mass.sum();massmean=np.sum(z*w[:,None],axis=0);dwm=z-massmean;covm=(dwm*w[:,None]).T@dwm
 shape=xx-xx.mean(axis=0);shape_cov=shape.T@shape/n;se=np.linalg.eigvalsh(shape_cov)
 out.update(supported=True,sigma_los_raw_kms=float(np.std(a)),sigma_los_Walker_kms=float(np.std(b)),
  mean_Walker_residual_kms=float(np.mean(b)),rms_Walker_residual_kms=float(np.sqrt(np.mean(b*b))),
  mean_los_kms=float(np.mean(a)),sigma_3d_kms=float(np.sqrt(np.trace(cov))),
  sigma_1d_equiv_kms=float(np.sqrt(np.trace(cov)/3)),velocity_covariance_kms2=cov.tolist(),
  velocity_mean_relative_core_kms=avg.tolist(),velocity_eigenvalues_kms2=eig.tolist(),
  velocity_principal_axes_native=vec.tolist(),sigma_common_LOS_tensor_kms=float(np.sqrt(max(0,axis@cov@axis))),
  sigma_xyz_kms=np.sqrt(np.diag(cov)).tolist(),sigma_spherical_kms=np.sqrt(vari).tolist(),
  spherical_mean_kms=sph.mean(axis=0).tolist(),beta_sigma=float(1-(vari[1]+vari[2])/(2*vari[0])) if vari[0]>0 else None,
  beta_second=float(1-(sec[1]+sec[2])/(2*sec[0])) if sec[0]>0 else None,
  angular_momentum_number_kpc_kms=np.mean(np.cross(xx,z),axis=0).tolist(),
  velocity_covariance_mass_kms2=covm.tolist(),sigma_3d_mass_kms=float(np.sqrt(np.trace(covm))),
  sigma_los_raw_mass_kms=float(np.sqrt(np.sum(w*(a-np.sum(w*a))**2))),
  N_eff_mass=float(1/np.sum(w*w)),spatial_covariance_kpc2=shape_cov.tolist(),
  rms_axis_c_over_a=float(np.sqrt(max(0,se[0])/se[-1])) if se[-1]>0 else None,
  los_depth_std_kpc=float(np.std(xx@axis)),mean_r_kpc=float(np.mean(sph_r)))
 return out

def projected_shape(sky,m,q):
 if q.sum()<32:return dict(supported=False,N=int(q.sum()))
 x=np.c_[sky['east_kpc'][q],sky['north_kpc'][q]];mass=m[q];com=np.average(x,weights=mass,axis=0)
 d=x-com;cov=(d*mass[:,None]).T@d/mass.sum();eig,vec=np.linalg.eigh(cov)
 return dict(supported=True,N=int(q.sum()),center_EN_kpc=com.tolist(),covariance_kpc2=cov.tolist(),
             q=float(np.sqrt(max(0,eig[0])/eig[1])) if eig[1]>0 else None,
             PA_EofN_deg=float(np.degrees(np.arctan2(vec[0,1],vec[1,1]))%180))

def snapshot(ctx,n,reference=None):
 import spatial_context as P
 rt=ctx['_runtime'] if '_runtime' in ctx else None
 if rt is None:
  O=C.module('spatial_observation',ctx['methods_path'])
  probe=O.astropy_native_direct_observables(np.array([[0.,0.,0.],[1.,2.,3.]]),np.zeros((2,3)),np.arange(2,dtype='u8'),
      host_center_raw=[0,0,0],host_bulk_raw=[0,0,0],length_to_kpc=1.,velocity_to_kms=1.,**C.FRAME)
  rt=dict(O=O,M=probe['physical_to_icrs']);ctx['_runtime']=rt
 O=rt['O'];L=ctx['units']['length_to_kpc'];V=ctx['units']['velocity_to_kms']
 a=P.read_stellar_properties(ctx,n)
 seed=np.asarray(ctx['center_seeds'][str(n)])
 center=P.shrinking_center(a['p'],a['m'],seed,L,initial_radius_kpc=5,min_radius_kpc=.3,min_count=128)
 r0=np.linalg.norm((a['p']-center)*L,axis=1);core=r0<.5
 if core.sum()<8:core=r0<2
 bulk=np.average(a['v'][core],weights=a['m'][core],axis=0)*V
 H,B,nMW=C._mw_center_bulk(a,L)
 # Membership is spatial. Stored ParticleIDs are never read or used.
 keep=r0<15.;b={k:a[k][keep] for k in ('p','v','m','tp','birth')}
 born=np.isfinite(b['birth']);old=b['tp']!=4
 x=((b['p']-center)*L).astype('f4');v=(b['v']*V-bulk).astype('f4');m=b['m']
 rows=np.arange(len(m),dtype='u8')  # temporary catalogue row labels, not physical particle IDs
 native=O.native_observables(b['p'],b['v']*V,rows,host_center_raw=H,host_bulk_raw=B*V,
       length_to_kpc=L,velocity_to_kms=1.,sun_position_kpc=C.SUN,sun_velocity_kms=C.VS)
 h=(center-H)*L-C.SUN;axis=h/np.linalg.norm(h)
 sky=O.project_sky(native['h_native_kpc'],rows,center_h_native_kpc=h,physical_to_icrs=rt['M'])
 fit=old&sky['fit_R5_mask']
 if fit.sum()<3:raise ValueError('Too few spatial old stars for Walker fit')
 walker=O.walker61_zero_error(sky,native['vlos_native_kms'],fit_mask=fit)
 theta=np.radians(walker['theta0_best_deg'])
 resid=native['vlos_native_kms']-walker['v0_kms']-walker['k_signed']*(sky['east_deg']*np.sin(theta)+sky['north_deg']*np.cos(theta))
 r=np.linalg.norm(x,axis=1);R=sky['projected_radius_kpc'];front=sky['front_mask']
 Re=C._quantile(R[front&old],m[front&old]);rh=C._quantile(r[old],m[old])
 populations={'old':old,'new':~old,'all':np.ones(len(m),bool)}
 stats={};profiles={}
 for pop,popmask in populations.items():
  regions={}
  for rad in RADII:regions['Rlt'+str(rad)]=front&(R<rad);regions['rlt'+str(rad)]=r<rad
  regions.update(R1_2=front&(R>=1)&(R<2),r1_2=(r>=1)&(r<2),RltRe=front&(R<Re),rltrhalf=r<rh)
  stats[pop]={key:moments(rows,x,v,m,native['vlos_native_kms'],resid,popmask&q,axis) for key,q in regions.items()}
  for geom,distance,mask in [('projected',R,front),('sphere',r,np.ones(len(r),bool))]:
   result=[]
   for lo,hi in zip(EDGES[:-1],EDGES[1:]):
    q=popmask&mask&(distance>=lo)&(distance<hi)
    moment=moments(rows,x,v,m,native['vlos_native_kms'],resid,q,axis)
    moment.update(lo_kpc=float(lo),hi_kpc=float(hi),geometry=geom)
    area=np.pi*(hi*hi-lo*lo);volume=4*np.pi*(hi**3-lo**3)/3
    moment['surface_density_Msun_kpc2' if geom=='projected' else 'density_Msun_kpc3']=float(m[q].sum()/(area if geom=='projected' else volume))
    result.append(moment)
   profiles[pop+'_'+geom]=result
 star_mass={str(rad):dict(old=float(m[old&(r<rad)].sum()),new=float(m[(~old)&(r<rad)].sum()),all=float(m[r<rad].sum())) for rad in RADII}
 gcpos=(center-H)*L;Rgc=np.linalg.norm(gcpos);vGC=bulk-B*V
 row=dict(snapshot=int(n),time_Gyr=float(a['time']),canonical_observer_id=C.OBSERVER,measurement_contract_id=CONTRACT,
   Dsun_kpc=float(np.linalg.norm(h)),Rgc_kpc=float(Rgc),Vr_gc_kms=float(vGC@gcpos/Rgc),Vgc_kms=float(np.linalg.norm(vGC)),
   Vtan_gc_kms=float(np.linalg.norm(np.cross(vGC,gcpos/Rgc))),core_center_raw=center.tolist(),core_bulk_native_kms=bulk.tolist(),
   host_center_raw=H.tolist(),host_bulk_raw=B.tolist(),N_MW_bulk=int(nMW),
   Re_circular_old_kpc=Re,rhalf_mass_old_kpc=rh,shape_old_Rlt2=projected_shape(sky,m,front&old&(R<2)),
   stats=stats,profiles=profiles,stellar_mass_spheres=star_mass,
   Walker=dict(N=walker['N_fit'],theta_deg=walker['theta0_best_deg'],k_kms_per_deg=walker['k_abs'],v0_kms=walker['v0_kms'],
       fit_population='spatial old stars, current forward R<5kpc within3D15kpc'),
   source_stellar_fields_sha256=a['source_fields_sha256'],particle_IDs_read=False,
   membership=dict(stars='3D r<15kpc, then physical projected/3D apertures; no origin-ID exclusion',
        old_new='simulation stellar type2/3 pre-existing and type4 formed; recorded age/formation time handled independently'),
   new_formation_ages_Gyr=dict(count_known=int(np.count_nonzero((~old)&born)),count_unknown=int(np.count_nonzero((~old)&~born))))
 age=a['age_Gyr'][keep] if 'age_Gyr' in a else a['time']-b['birth'];row['age_groups']={}
 row['age_source']=a.get('age_source','snapshot_time_minus_recorded_formation_time')
 row['age_availability']={}
 for pop,q in populations.items():
  valid=q&np.isfinite(age)&(age>=0)&front&(R<5)
  row['age_availability'][pop]=dict(N_available=int(valid.sum()),N_missing=int((q&~np.isfinite(age)&front&(R<5)).sum()),
    age_mean_Gyr=float(np.mean(age[valid])) if valid.any() else None,
    age_median_Gyr=float(np.median(age[valid])) if valid.any() else None,
    source='recorded particle properties, independent of particle IDs')
 for label,lo,hi in [('0_0p1',0,.1),('0p1_0p5',.1,.5),('0p5_1',.5,1),('gte1',1,np.inf)]:
  q=np.isfinite(age)&(age>=lo)&(age<hi)&front&(R<5)
  row['age_groups'][label]=moments(rows,x,v,m,native['vlos_native_kms'],resid,q,axis)
 hoststar=(np.linalg.norm((a['p']-H)*L,axis=1)<Rgc)&(r0>=15)
 row['host_enclosed_stellar_mass_Msun']=float(a['m'][hoststar].sum())
 dm={str(rad):0. for rad in RADII};host_dm=0.
 for pp in ctx['snapshot_paths'][str(n)]:
  with C.h5py.File(pp,'r') as stream:
   if 'PartType1' not in stream:continue
   group=stream['PartType1'];count=len(group['Coordinates']);table=stream['Header'].attrs['MassTable']
   for start in range(0,count,262144):
    sl=slice(start,min(start+262144,count));position=group['Coordinates'][sl].astype('f8')
    mass=(group['Masses'][sl].astype('f8') if 'Masses' in group else np.full(len(position),table[1]))*ctx['units']['mass_to_msun']
    rr=np.linalg.norm((position-center)*L,axis=1);rhost=np.linalg.norm((position-H)*L,axis=1)
    host_dm+=float(mass[(rr>=15)&(rhost<Rgc)].sum())
    for rad in RADII:dm[str(rad)]+=float(mass[rr<rad].sum())
 row['local_DM_spheres_Msun']=dm;row['DM_status']='all_local_DM_by_position_includes_host_background'
 row['host_enclosed_DM_mass_Msun']=host_dm
 C._unchanged(a['source_stats']);del a,b,native,sky;gc.collect()
 return row,center,bulk
