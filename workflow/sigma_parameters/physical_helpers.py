"""Shared physical units and host geometry. No particle-ID interface."""
import hashlib,importlib.util,json
from pathlib import Path
import numpy as np
import h5py
OBSERVER='astropy_native_direct_v1'
FRAME=dict(galcen_distance_kpc=8.122,z_sun_kpc=.0208,galcen_v_sun_kms=[12.9,245.6,7.78],
           galcen_icrs_deg=[266.4051,-28.936175],roll_deg=0.,native_to_GC=np.eye(3).tolist())
SUN=np.array([-np.sqrt(8.122**2-.0208**2),0.,.0208]);VS=np.array([12.9,245.6,7.78]);CACHE={}
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def module(name,path):
 spec=importlib.util.spec_from_file_location(name,path);obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj);return obj
def _quantile(radius,mass,f=.5):
 if not len(radius):return float('nan')
 order=np.argsort(radius);r=np.asarray(radius)[order];w=np.asarray(mass)[order]
 if w.sum()<=0:return float('nan')
 return float(np.interp(f*w.sum(),np.cumsum(w),r))
def _mw_center_bulk(a,L):
 p,m=a['p'],a['m'];center=np.average(p,weights=m,axis=0);radius=300.
 while radius>5.:
  q=np.linalg.norm((p-center)*L,axis=1)<radius
  if q.sum()<100:raise ValueError('Insufficient MW-centre stars')
  center=np.average(p[q],weights=m[q],axis=0);radius*=.8
 q=np.linalg.norm((p-center)*L,axis=1)<5
 if q.sum()<100:raise ValueError('Insufficient stellar MW bulk sample')
 return center,np.average(a['v'][q],weights=m[q],axis=0),int(q.sum())
def _unchanged(stats):
 for p,size,mtime,*_ in stats:
  st=Path(p).stat()
  if (st.st_size,st.st_mtime_ns)!=(size,mtime):raise RuntimeError('Physical source changed during calculation: '+p)
