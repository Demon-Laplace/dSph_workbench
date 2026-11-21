#!/usr/bin/env python3
"""Build an offline 3-D viewer from local GIZMO HDF5 snapshots or particle INIs."""
import argparse
import json
from pathlib import Path
import re
import sys
import webbrowser

import numpy as np


LABELS = ('Gas', 'Type 1', 'Old stars', 'Type 3', 'New stars', 'Type 5')
COLORS = ('#3498db', '#9877b0', '#f39c12', '#807060', '#e65c72', '#56a772')
SNAPSHOT = re.compile(r'((?:snapshot|snap)_(\d+))(?:\.(\d+))?\.h(?:df5|5)$', re.I)
CHUNK = 262144


def discover(paths, frame_limit, skip_latest=False):
    """Keep split files in one frame; order directories by snapshot number."""
    frames = {}
    for path in paths:
        path = path.expanduser().resolve(strict=True)
        if path.is_dir():
            entries = list(path.iterdir())
            for sub in (p for p in tuple(entries) if re.fullmatch(r'snapdir_\d+', p.name) and p.is_dir()):
                entries.extend(sub.iterdir())
            groups = {}
            for candidate in entries:
                if m := SNAPSHOT.fullmatch(candidate.name):
                    if candidate.is_file():
                        groups.setdefault((int(m[2]), candidate.parent, m[1]), []).append(candidate)
            order = sorted(groups)
            if skip_latest and order:
                latest = order[-1][0]
                order = [k for k in order if k[0] < latest]
            for number, parent, stem in order:
                frames[parent/stem] = part_order(groups[(number, parent, stem)])
        elif m := SNAPSHOT.fullmatch(path.name):
            if m[3] is None:
                frames[path.parent/m[1]] = (path,)
            else:
                files = [p for p in path.parent.iterdir() if (v := SNAPSHOT.fullmatch(p.name))
                         and v[1] == m[1] and v[3] is not None and p.is_file()]
                frames[path.parent/m[1]] = part_order(files)
        else:
            frames[path] = (path,)
    items = list(frames.items())
    if not items:
        raise ValueError('No snapshots found. Pass individual snapshot files or a snapshot directory.')
    if len(items) > frame_limit:
        indices = (len(items)-1,) if frame_limit == 1 else np.linspace(0, len(items)-1, frame_limit, dtype=int)
        items = [items[i] for i in indices]
    return items


def part_order(paths):
    suffixes = [SNAPSHOT.fullmatch(p.name)[3] for p in paths]
    if len(paths) > 1 and None in suffixes:
        raise ValueError(f'Ambiguous single and split files for the same snapshot: {paths[0].parent}')
    if None not in suffixes and sorted(map(int, suffixes)) != list(range(len(paths))):
        raise ValueError(f'Missing or duplicate snapshot part numbers: {paths[0].parent}')
    return tuple(sorted(paths, key=lambda p: int(SNAPSHOT.fullmatch(p.name)[3] or 0)))


def empty():
    return (np.empty((0, 3), dtype=np.float32), np.empty((0, 3), dtype=np.float32),
            np.empty(0, dtype=np.float32))


def join(parts):
    return tuple(np.concatenate(a) for a in zip(*parts)) if parts else empty()


def hdf_header(paths):
    import h5py
    time = None
    for path in paths:
        with h5py.File(path, 'r') as f:
            attrs = f['Header'].attrs
            if int(attrs.get('NumFilesPerSnapshot', 1)) != len(paths):
                raise ValueError(f'Incomplete HDF5 snapshot: {path}')
            current = float(attrs['Time']) if 'Time' in attrs else None
            if current is not None and not np.isfinite(current):
                raise ValueError(f'Non-finite snapshot time: {path}')
            if path != paths[0] and current != time:
                raise ValueError(f'Snapshot part has a different time: {path}')
            time = current
    return time


def get_hdf5(paths, kind, length, velocity, center=None, radius=None):
    import h5py
    parts = []
    for path in paths:
        with h5py.File(path, 'r') as f:
            group = f.get(f'PartType{kind}')
            if group is None:
                continue
            total = len(group['Coordinates'])
            if not total:
                continue
            masses = group.get('Masses')
            table = f['Header'].attrs.get('MassTable', np.zeros(6))
            if masses is None and table[kind] <= 0:
                raise ValueError(f'{path}: PartType{kind} has neither Masses nor a positive MassTable entry')
            for start in range(0, total, CHUNK):
                rows = slice(start, start+CHUNK)
                x = group['Coordinates'][rows].astype(np.float32)*length
                mask = np.isfinite(x).all(axis=1)
                if center is not None:
                    dx = x-center
                    mask &= np.einsum('ij,ij->i', dx, dx) < radius**2
                if not np.any(mask):
                    continue
                v = group['Velocities'][rows].astype(np.float32)[mask]*velocity
                m = (masses[rows].astype(np.float32)[mask] if masses is not None
                     else np.full(mask.sum(), table[kind], dtype=np.float32))
                valid = np.isfinite(v).all(axis=1) & np.isfinite(m) & (m > 0)
                if np.any(valid):
                    parts.append((x[mask][valid], v[valid], m[valid]))
    return join(parts)


def read_ini(path, length, velocity):
    """Particle export layout is id, xyz, vxyz, mass, extra, type; skip ID."""
    with path.open() as stream:
        head = stream.readline().split()
        if len(head) != 4 or not all(s.isdigit() for s in head):
            raise ValueError(f'{path} is not a particle INI (parameter INIs are not supported)')
        data = np.loadtxt(stream, usecols=(*range(1, 8), 9), dtype=np.float32, ndmin=2)
    data = data[np.isfinite(data).all(axis=1) & (data[:, 6] > 0)]
    groups = {}
    for k in range(len(LABELS)):
        d = data[data[:, 7] == k]
        groups[k] = (d[:, :3]*length, d[:, 3:6]*velocity, d[:, 6])
    return groups


def weighted_mean(x, w):
    return np.einsum('i,ij->j', w, x, dtype=np.float64) / w.sum(dtype=np.float64)


def within(x, center, radius):
    dx = x-center
    return np.einsum('ij,ij->i', dx, dx) < radius**2


def define_frame(stars, radius, axis_radius, axis_mode, center=None, bulk=None):
    x, v, w = stars
    if center is None:
        if len(x) < 3:
            raise ValueError('Auto centering needs stars; choose --center-type or supply --center and --bulk')
        center = np.median(x, axis=0)
        for width in (radius*2, radius, axis_radius):
            member = within(x, center, width)
            if member.sum() >= 3:
                center = weighted_mean(x[member], w[member])
    center = np.asarray(center)
    core = within(x, center, axis_radius)
    if bulk is None:
        use = core if core.sum() >= 3 else within(x, center, radius)
        if use.sum() < 3:
            raise ValueError('Too few anchor particles for bulk velocity; use --bulk or adjust center/aperture')
        bulk = weighted_mean(v[use], w[use])
    bulk = np.asarray(bulk)
    z = np.array([0., 0., 1.])
    shape = None
    if core.sum() >= 3:
        xx, mm = x[core]-center, w[core]
        tensor = (xx.T*mm) @ xx / mm.sum(dtype=np.float64)
        eigen, vectors = np.linalg.eigh(tensor)
        if eigen[-1] > 0:
            shape = np.sqrt(np.maximum(eigen[[1, 0]], 0)/eigen[-1]).tolist()
        if axis_mode == 'short':
            z = vectors[:, 0]
            if z[2] < 0:
                z = -z
        elif axis_mode == 'angular-momentum':
            z = weighted_mean(np.cross(xx, v[core]-bulk), mm)
            norm = np.linalg.norm(z)
            if norm < 1e-12:
                raise ValueError('Anchor angular momentum is negligible; use --axis native-z or short')
            z = z/norm
    elif axis_mode != 'native-z':
        raise ValueError('Too few anchor particles inside --axis-radius to measure the requested axis')
    # Rotate coordinates and velocities by the same proper rotation.
    ref = np.array([1., 0., 0.]) if abs(z[0]) < .95 else np.array([0., 1., 0.])
    y = np.cross(z, ref)
    y /= np.linalg.norm(y)
    return center, bulk, np.column_stack((np.cross(y, z), y, z)), shape


def group_display(x, v, m, radius, factor, point_limit, seed, profile_limit=0):
    # Statistics use the entire aperture. Only the 3-D marker payload is sampled.
    R = np.hypot(x[:, 0], x[:, 1])
    jz = x[:, 0]*v[:, 1]-x[:, 1]*v[:, 0]
    ok = R > 0
    vp = np.divide(jz, R, out=np.zeros_like(jz), where=ok)
    vr = np.divide((x[:, :2]*v[:, :2]).sum(axis=1), R, out=np.zeros_like(jz), where=ok)
    edges = np.linspace(0, radius, 17)
    bins = np.minimum(np.digitize(R[ok], edges)-1, len(edges)-2)
    # Bincount accumulates weighted moments in double precision.
    weight = np.bincount(bins, weights=m[ok], minlength=16)
    take = weight > 0

    def radial(vv):
        a = np.bincount(bins, weights=m[ok]*vv[ok], minlength=16)
        b = np.bincount(bins, weights=m[ok]*vv[ok]**2, minlength=16)
        mu = a[take]/weight[take]
        return mu.tolist(), np.sqrt(np.maximum(0, b[take]/weight[take]-mu**2)).tolist()

    mean, phi_sigma = radial(vp)
    _, radial_sigma = radial(vr)
    total = m.sum(dtype=np.float64)
    occupied = m[ok].sum(dtype=np.float64)
    dims = {'xy':(0,1), 'xz':(0,2), 'yz':(1,2)}
    grid = {}
    for key, (i, j) in dims.items():
        hist, _, _ = np.histogram2d(x[:, i], x[:, j], bins=56,
                                     range=((-radius,radius),)*2, weights=m)
        grid[key] = np.round(hist.T*factor/(2*radius/56)**2, 6).tolist()
    if len(x) > point_limit:
        sample = np.sort(np.random.default_rng(seed).choice(len(x), point_limit, replace=False))
    else:
        sample = slice(None)

    def marker(vals):
        return np.round(np.asarray(vals[sample], dtype=np.float64), 5).tolist()

    rr = np.linalg.norm(x, axis=1)
    order = np.argsort(rr)
    r50 = np.interp(.5*total, np.cumsum(m[order],dtype=np.float64)-m[order]*.5, rr[order])
    light = None
    if profile_limit:
        pick = np.sort(np.random.default_rng(seed+17).choice(len(x), profile_limit, replace=False)) if len(x) > profile_limit else slice(None)
        selected_mass = m[pick].astype(float)
        scaled = selected_mass * (factor*total/selected_mass.sum())
        light = {'xyz':np.round(x[pick], 5).tolist(), 'mass':scaled.tolist()}
    return {'n':len(x), 'mass':float(total*factor), 'mean':float(np.dot(m[ok].astype(float), vp[ok])/occupied) if occupied else None,
            'prograde':float(m[jz > 0].sum(dtype=np.float64)/occupied) if occupied else None,
            'j':float(np.dot(m.astype(float), jz)/total) if total else None,
            'r50':float(r50), 'zrms':float(np.sqrt(np.dot(m.astype(float),x[:,2]**2)/total)),
            'light':light, 'maps':grid,
            'profile':{'r':((edges[:-1]+edges[1:])/2)[take].tolist(), 'mean':mean,
                       'radial_sigma':radial_sigma, 'phi_sigma':phi_sigma},
            'points':{**{k:marker(x[:,i]) for i,k in enumerate(('x','y','z'))},
                      **{k:marker(v[:,i]) for i,k in enumerate(('vx','vy','vz'))},
                      'vphi':marker(vp), 'speed':marker(np.linalg.norm(v, axis=1)),
                      'height':marker(np.abs(x[:, 2]))}}


def build_frame(key, paths, args, length, velocity, mass):
    hdf = paths[0].suffix.lower() in ('.hdf5', '.h5')
    if hdf:
        time = hdf_header(paths)
        # Skip reading an entire anchor component when a fixed frame is supplied.
        direct = args.center is not None and args.bulk is not None and args.axis == 'native-z'
        stars = empty() if direct else get_hdf5(paths, args.center_type, length, velocity, args.center,
                                               args.radius if args.center is not None else None)
    else:
        if paths[0].suffix.lower() != '.ini':
            raise ValueError(f'Unsupported input: {paths[0]} (use a particle INI or HDF5 snapshot)')
        groups = read_ini(paths[0], length, velocity)
        stars, time = groups[args.center_type], 0.
    center, bulk, basis, shape = define_frame(stars, args.radius, args.axis_radius, args.axis, args.center, args.bulk)
    result = {}
    for k in args.types:
        if hdf:
            part = stars if k == args.center_type and not direct else get_hdf5(paths, k, length, velocity, center, args.radius)
        else:
            part = groups[k]
        x, v, m = part
        member = within(x, center, args.radius)
        if member.any():
            result[str(k)] = group_display((x[member]-center)@basis, (v[member]-bulk)@basis,
                                          m[member], args.radius, mass, args.points, args.seed+k,
                                          args.profile_points if k in args.light_types else 0)
    if not result:
        raise ValueError(f'No particles from the selected types inside the display aperture: {key}')
    return {'name':key.name, 'directory':str(key.parent), 'time':time, 'center':center.tolist(),
            'bulk':bulk.tolist(), 'axis':basis[:, 2].tolist(), 'shape':shape, 'groups':result}


HTML = r"""<!doctype html>
<html lang="en"><meta charset="utf-8"><title>Fornax 3D viewer</title>
<style>
body{font:15px system-ui;margin:18px;background:#f6f8fb;color:#192435}
main{max-width:1600px;margin:auto} aside{margin:9px 0;display:flex;flex-wrap:wrap;gap:12px;align-items:center}
button,select,input{font:inherit} #scene{height:735px} .pair{display:grid;grid-template-columns:1fr 1fr}
.overview{display:grid;grid-template-columns:minmax(0,1.5fr) minmax(0,1fr);gap:8px}
#lightmap,#lightcurve{height:337px} #projection,#profiles{height:430px}
table{border-collapse:collapse;margin:10px 0}td,th{padding:5px 12px;text-align:right;border:1px solid #ccd3dc}
.panel{background:white} small{overflow-wrap:anywhere}
@media(max-width:900px){.pair,.overview{grid-template-columns:1fr}#scene{height:515px}}
</style><main><h1>Fornax · simulation explorer</h1>
<aside><label>Frame <input id="frame" type="range" min="0" value="0"></label><button id="play">Play</button>
<label>Color <select id="color"><option value="type">Component</option><option value="vphi">Signed azimuthal velocity</option>
<option value="speed">Speed</option><option value="vx">Vx</option><option value="vy">Vy</option><option value="vz">Vz</option>
<option value="height">Absolute height</option></select></label>
<label>Projection <select id="plane"><option>xy</option><option>xz</option><option>yz</option></select></label>
<button id="face">Look down Z</button><button id="edge">Look along Y</button><button id="reset">Reset view</button></aside>
<aside id="types"></aside><small id="meta"></small>
<div class="overview"><div id="scene" class="panel"></div>
<div><h2 style="margin:5px 0">Starlight along the camera</h2>
<small id="lightmeta"></small>
<div id="lightmap" class="panel"></div><div id="lightcurve" class="panel"></div></div></div>
<div class="pair"><div id="projection" class="panel"></div><div id="profiles" class="panel"></div></div>
<table id="stats"></table>
<p>Drag to rotate, scroll to zoom. Camera changes only the view. Coordinates and velocities use the declared dwarf frame;
 Vx/Vy/Vz are not heliocentric observables. The spherical aperture defines plotted membership, not a binding or HI measurement.
 The 3-D markers and live stellar-light preview use separate samples. The fixed-axis mass maps and kinematics use the full aperture.</p>
</main><script>__PLOTLY__</script>
<script id="payload" type="application/json">__DATA__</script>
<script>
'use strict';
const D=JSON.parse(document.getElementById('payload').textContent);
const $=id=>document.getElementById(id), cfg={responsive:true,scrollZoom:true};
const selection={}; let timer=null, queuedLight=null, cameraHooked=false;
const oblique={eye:{x:1.45,y:1.45,z:1.18},center:{x:0,y:0,z:0},up:{x:0,y:0,z:1}};
let camera=oblique;
for(const k of D.types){
  selection[k]=true;
  const input=document.createElement('input'), label=document.createElement('label');
  input.type='checkbox'; input.checked=true;
  input.onchange=()=>{selection[k]=input.checked; draw();};
  label.append(input,document.createTextNode(D.labels[k]+' · TP'+k)); $('types').append(label);
}
$('frame').max=D.frames.length-1;
function fmt(x){return x===null?'—':Number(x).toPrecision(4);}
const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
function unit(a){const n=Math.hypot(...a); return n ? a.map(v=>v/n) : null;}
function vector(a,fallback){return a?[a.x,a.y,a.z]:fallback;}
function cameraLight(){
 const frame=D.frames[Number($('frame').value)], radius=D.radius, nbin=32, pixels=64;
 const eye=vector(camera.eye,[1.45,1.45,1.18]), center=vector(camera.center,[0,0,0]);
 const n=unit(eye.map((v,i)=>v-center[i]))||[0,0,1];
 const east=unit(cross(vector(camera.up,[0,0,1]),n))||unit(cross([0,1,0],n))||[1,0,0];
 const north=cross(n,east), inclination=Math.acos(Math.min(1,Math.abs(n[2])))*180/Math.PI;
 const map=Array.from({length:pixels},()=>Array(pixels).fill(0)), radial=Array(nbin).fill(0);
 let sample=0;
 for(const k of D.types){
  const p=selection[k] && frame.groups[k] && frame.groups[k].light;
  if(!p)continue;
  sample+=p.xyz.length;
  for(let i=0;i<p.xyz.length;i++){
   const x=dot(p.xyz[i],east),y=dot(p.xyz[i],north);
   const a=Math.floor((x/radius+1)*pixels/2),b=Math.floor((y/radius+1)*pixels/2);
   const q=Math.floor(Math.hypot(x,y)*nbin/radius);
   if(a>=0&&a<pixels&&b>=0&&b<pixels)map[b][a]+=p.mass[i];
   if(q<nbin)radial[q]+=p.mass[i];
  }
 }
 const ph=D.mass==='Msun'&&D.length==='kpc';
 const toValue=(mass,area)=>mass>0?(ph?4.83+21.572-2.5*Math.log10(mass/D.ml/area/1e6):Math.log10(mass/area)):null;
 const x=Array.from({length:pixels},(_,i)=>radius*(2*(i+.5)/pixels-1));
 const r=Array.from({length:nbin},(_,i)=>radius*(i+.5)/nbin);
 const mu=radial.map((m,i)=>toValue(m,Math.PI*((i+1)**2-i**2)*(radius/nbin)**2));
 const ylabel=ph?'μV (mag arcsec⁻²)':'log₁₀ stellar surface mass density';
 $('lightmeta').textContent='Inclination to Z: '+inclination.toFixed(1)+'° · '+sample+
  ' sampled stars · '+(ph?'V-band preview, fixed M/L='+D.ml+' solar units':'Native mass projection (no physical brightness conversion)');
 Plotly.react('lightmap',[{type:'heatmap',x,y:x,z:map.map(row=>row.map(m=>toValue(m,(2*radius/pixels)**2))),
    colorscale:'Viridis',reversescale:ph,colorbar:{title:{text:ylabel}}}],
    {title:{text:'Camera-aligned stellar projection'},margin:{t:50},
     xaxis:{title:'Sky-plane x ('+D.length+')'},yaxis:{title:'Sky-plane y ('+D.length+')',scaleanchor:'x'}},cfg);
 Plotly.react('lightcurve',[{type:'scatter',mode:'lines+markers',x:r,y:mu,connectgaps:false,
    name:ph?'V band':'Stellar mass',line:{color:'#973d2a'}}],
    {title:{text:'Circular-annulus stellar profile'},margin:{t:50},xaxis:{title:'Projected R ('+D.length+')',range:[0,radius]},
     yaxis:{title:ylabel,autorange:ph?'reversed':true}},cfg);
}
function trackCamera(event){
 const newCamera=event&&event['scene.camera'];
 if(newCamera) camera=newCamera;
 else if(event&&Object.keys(event).some(k=>k.startsWith('scene.camera')))camera=$('scene').layout.scene.camera||camera;
 else return;
 if(queuedLight===null)queuedLight=setTimeout(()=>{queuedLight=null;cameraLight();},90);
}
function setCamera(next){camera={center:{x:0,y:0,z:0},...next};Plotly.relayout('scene',{'scene.camera':camera});cameraLight();}
function draw(){
 const frame=D.frames[Number($('frame').value)], color=$('color').value, plane=$('plane').value;
 $('meta').textContent=frame.name+' | t (code): '+fmt(frame.time)+' | '+frame.directory+
  ' | center: '+frame.center.map(fmt).join(', ')+' | axis: '+frame.axis.map(fmt).join(', ')+
  ' | b/a, c/a: '+(frame.shape?frame.shape.map(fmt).join(', '):'—');
 const scene=[], prof=[], view=Array.from({length:56},()=>Array(56).fill(0));
 const table=$('stats'); table.replaceChildren();
 function row(labels,head){const tr=document.createElement('tr'); for(const v of labels){
   const td=document.createElement(head?'th':'td'); td.textContent=v; tr.append(td);
 } table.append(tr);}
 row(['Component','Particles in sphere','Mass ('+D.mass+')','Vφ ('+D.velocity+')','Prograde %','jZ ('+D.length+' × '+D.velocity+')',
      'r₅₀ ('+D.length+')','z rms ('+D.length+')'],true);
 const signed=['vphi','vx','vy','vz'].includes(color), radius=D.radius;
 let showScale=true;
 for (const k of D.types) {
  const g=frame.groups[k]; if(!g) continue;
  row([D.labels[k],g.n,fmt(g.mass),fmt(g.mean),g.prograde===null?'—':fmt(g.prograde*100),fmt(g.j),fmt(g.r50),fmt(g.zrms)],false);
  if(!selection[k]) continue;
  const p=g.points, style={size:3,opacity:.66,color:D.colors[k]};
  if(color!=='type'){
   style.color=p[color]; style.colorscale=signed?'RdBu':'Viridis'; style.reversescale=signed;
   style.cmin=signed?-D.limits[color]:0; style.cmax=D.limits[color];
   style.showscale=showScale; style.colorbar={title:{text:color+' ('+(color==='height'?D.length:D.velocity)+')'}};
   showScale=false;
  }
  scene.push({type:'scatter3d',mode:'markers',name:D.labels[k], x:p.x,y:p.y,z:p.z,marker:style,customdata:p.vphi,
     hovertemplate:'x=%{x:.3f}<br>y=%{y:.3f}<br>z=%{z:.3f}<br>Vφ=%{customdata:.3f}<extra>%{fullData.name}</extra>'});
  for(let i=0;i<56;i++)for(let j=0;j<56;j++) view[i][j]+=g.maps[plane][i][j];
  for(const [metric,tag,dash] of [['mean','Vφ','solid'],['radial_sigma','σR','dash'],['phi_sigma','σφ','dot']])
    prof.push({type:'scatter',mode:'lines+markers',name:D.labels[k]+' '+tag,x:g.profile.r,y:g.profile[metric],
       marker:{size:4,color:D.colors[k]},line:{color:D.colors[k],dash}});
 }
 const axis=name=>({title:name+' ('+D.length+')',range:[-radius,radius]});
 const grid=Array.from({length:56},(_,i)=>radius*((i+.5)*2/56-1));
 const logmap=view.map(row=>row.map(v=>v>0?Math.log10(v):null));
 Plotly.react('scene',scene,{title:{text:'Particle structure · '+frame.name},
    uirevision:'keep-camera',margin:{t:52,b:15,l:0,r:0},showlegend:false,
    scene:{xaxis:axis('x'),yaxis:axis('y'),zaxis:axis('z'),aspectmode:'cube',camera}},cfg).then(()=>{
   if (!cameraHooked){$('scene').on('plotly_relayout', trackCamera);$('scene').on('plotly_relayouting',trackCamera);cameraHooked=true;}
 });
 Plotly.react('projection',[{type:'heatmap',x:grid,y:grid,z:logmap,colorscale:'Viridis',
    colorbar:{title:{text:'log₁₀ Σ ('+D.mass+'/'+D.length+'²)'}}}],
    {title:{text:'Projected mass density · '+plane.toUpperCase()},margin:{t:52},
    xaxis:axis(plane[0]),yaxis:{...axis(plane[1]),scaleanchor:'x',constrain:'domain'}},cfg);
 Plotly.react('profiles',prof,{title:{text:'Rotation and dispersion profiles'},margin:{t:52},
    xaxis:{title:'Cylindrical R ('+D.length+')',range:[0,radius]},yaxis:{title:D.velocity},
    legend:{orientation:'h'}},cfg);
 cameraLight();
}
$('frame').oninput=draw; $('color').onchange=draw; $('plane').onchange=draw;
$('face').onclick=()=>setCamera({eye:{x:0,y:0,z:2.4},up:{x:0,y:1,z:0}});
$('edge').onclick=()=>setCamera({eye:{x:0,y:-2.4,z:0},up:{x:0,y:0,z:1}});
$('reset').onclick=()=>setCamera(oblique);
$('play').disabled=D.frames.length<2;
$('play').onclick=()=>{
 if(timer!==null){clearInterval(timer);timer=null;$('play').textContent='Play';return;}
 $('play').textContent='Pause';
 timer=setInterval(()=>{$('frame').value=(Number($('frame').value)+1)%D.frames.length;draw();},1200);
};
draw();
</script></html>
"""


def payload(frames, args, units):
    marks = ('vphi', 'vx', 'vy', 'vz', 'speed', 'height')
    limits = {}
    for key in marks:
        options = [p['points'][key] for frame in frames for p in frame['groups'].values()]
        limits[key] = max(float(np.quantile(np.abs(np.concatenate(options)), .98)), 1e-12)
    return {'frames':frames, 'types':args.types, 'labels':LABELS, 'colors':COLORS, 'radius':args.radius,
            'length':units[0], 'velocity':units[1], 'mass':units[2], 'limits':limits, 'ml':args.mass_to_light}


def types_arg(value):
    try:
        types = list(dict.fromkeys(map(int, value.split(','))))
    except ValueError as error:
        raise argparse.ArgumentTypeError('Use comma-separated Gadget types, e.g. 0,2,4') from error
    if not types or min(types) < 0 or max(types) >= len(LABELS):
        raise argparse.ArgumentTypeError('Particle types must be 0–5')
    return types


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    arguments = list(sys.argv[1:] if argv is None else argv)
    if not arguments or arguments == ['--gui']:
        try:
            from fornax_viewer_gui import launch
            return launch()
        except (ImportError, RuntimeError) as error:
            parser.exit(1, f'{error}\nUse --help for offline HTML export without a desktop.\n')
    parser.add_argument('inputs', type=Path, nargs='+', help='snapshot file(s), directory of snapshots, or particle INI(s)')
    parser.add_argument('--output', type=Path, required=True, help='new HTML file (existing output is never overwritten)')
    parser.add_argument('--open', action='store_true', help='open the generated page in the default browser')
    parser.add_argument('--types', type=types_arg, default=[0,2,4])
    parser.add_argument('--center-type', choices=range(6), type=int, default=2, help='anchor particle type; default TP2')
    parser.add_argument('--center', type=float, nargs=3, help='fixed center in displayed coordinates')
    parser.add_argument('--bulk', type=float, nargs=3, help='fixed velocity in displayed velocity units')
    parser.add_argument('--axis', choices=('native-z','short','angular-momentum'), default='short')
    parser.add_argument('--radius', type=float, default=5., help='spherical aperture in displayed length units')
    parser.add_argument('--axis-radius', type=float, default=2., help='inner sphere for shape and axis')
    parser.add_argument('--units', choices=('native','fornax'), default='native',
                        help='use raw numbers, or explicitly use kpc, km/s, 1e10 Msun mass unit; no automatic h factor')
    parser.add_argument('--length-factor', type=float, default=1.)
    parser.add_argument('--velocity-factor', type=float, default=1.)
    parser.add_argument('--mass-factor', type=float, default=1., help='multiplier in addition to selected unit convention')
    parser.add_argument('--points', type=int, default=2500, help='maximum preview points per type per frame')
    parser.add_argument('--light-types', type=types_arg, default=[2,4], help='particle types contributing to stellar light')
    parser.add_argument('--profile-points', type=int, default=20000, help='stellar light samples per type per frame')
    parser.add_argument('--mass-to-light', type=float, default=2.6, help='constant V-band M/L, in solar units')
    parser.add_argument('--frames', type=int, default=8, help='evenly space at most this many frames; 1 selects the last')
    parser.add_argument('--skip-latest', action='store_true', help='omit highest snapshot number when reading a directory')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args(arguments)
    for attr in ('radius','axis_radius','length_factor','velocity_factor','mass_factor','mass_to_light'):
        val = getattr(args, attr)
        if not np.isfinite(val) or val <= 0:
            parser.error(f'--{attr.replace("_","-")} must be positive and finite')
    if args.axis_radius > args.radius:
        parser.error('--axis-radius must not exceed --radius')
    for vector in (args.center, args.bulk):
        if vector is not None and not np.isfinite(vector).all():
            parser.error('--center and --bulk must contain finite values')
    if min(args.points, args.frames, args.profile_points) < 1 or args.seed < 0:
        parser.error('--points, --profile-points and --frames must be positive, --seed must be nonnegative')
    output = args.output.expanduser()
    if output.exists() or output.is_symlink() or output.suffix.lower() != '.html':
        parser.error('--output must be a new .html file')
    units = ('kpc','km/s','Msun') if args.units == 'fornax' else ('native length','native velocity','native mass')
    mass = args.mass_factor * (1e10 if args.units == 'fornax' else 1.)
    try:
        from plotly.offline import get_plotlyjs
        frames = [build_frame(key, paths, args, args.length_factor, args.velocity_factor, mass)
                  for key, paths in discover(args.inputs, args.frames, args.skip_latest)]
        data = json.dumps(payload(frames,args,units), separators=(',',':'), allow_nan=False).replace('<', '\\u003c')
        rendered = HTML.replace('__PLOTLY__', get_plotlyjs(), 1).replace('__DATA__', data, 1)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('x') as stream:
            stream.write(rendered)
    except (ImportError, OSError, ValueError, KeyError) as error:
        parser.exit(1,f'Visualization failed: {error}\n')
    print(f'{len(frames)} frame(s) → {output.resolve()}')
    if args.open:
        webbrowser.open(output.resolve().as_uri())
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
