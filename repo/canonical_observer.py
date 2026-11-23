"""Production adapters for the single central observing implementation."""
from pathlib import Path as _Path
import importlib.util as _import_util
import sys as _sys


def _load_observing_core():
    code = next((p for p in _Path(__file__).resolve().parents if p.name == "Code"), None)
    if code is not None:
        # Normal module identity is importable by multiprocessing/spawn.
        if str(code) not in _sys.path:
            _sys.path.insert(0, str(code))
        from importlib import import_module
        return import_module("simulation.dSph_workbench.repo.observing_core")
    # A standalone released/job package carries a generated pinned side file.
    source = _Path(__file__).resolve().with_name("observing_core.py")
    if not source.is_file():
        raise FileNotFoundError("Use the central Code tree or a complete package including observing_core.py")
    key = "_fornax_observing_core_" + str(source.resolve()).replace("\\", "_").replace("/", "_").replace(":", "_")
    value = _sys.modules.get(key)
    if value is None:
        spec = _import_util.spec_from_file_location(key, source)
        value = _import_util.module_from_spec(spec)
        _sys.modules[key] = value
        spec.loader.exec_module(value)
    return value


_observing_core = _load_observing_core()
for _symbol, _value in vars(_observing_core).items():
    if not _symbol.startswith("__"):
        globals()[_symbol] = _value

def attach_frozen_vlos(native_result, frozen_particle_id, frozen_vlos_kms, *,
                       atol=1.0e-10, row_aligned=False):
    """Preserve production diagnostic aliases around the shared calculation."""
    result = _observing_core.attach_frozen_vlos(
        native_result, frozen_particle_id, frozen_vlos_kms,
        atol=atol, row_aligned=row_aligned)
    status = result.get("frozen_vlos_reference_status")
    if status == "attached":
        result["frozen_vlos_reference_status"] = "used_aligned_reference"
        result["vlos_source"] = "frozen_id_matched"
        source = np.asarray(frozen_vlos_kms)
        if source.dtype.kind == "f":
            result["vlos_native_kms"] = result["vlos_native_kms"].astype(source.dtype, copy=False)
    elif status == "not_used_los_mismatch":
        result["frozen_vlos_reference_status"] = "not_used_different_reference_velocities"
    elif status == "not_used_invalid_reference":
        result["frozen_vlos_reference_status"] = "not_used_invalid_reference_values"
    return result


# Production adapters reuse the central observing_core implementation.
CANONICAL_OBSERVER_ID = "astropy_native_direct_v1"
OBSERVER_FRAME = dict(galcen_distance_kpc=8.122, z_sun_kpc=.0208,
    galcen_v_sun_kms=[12.9,245.6,7.78],galcen_icrs_deg=[266.4051,-28.936175],
    roll_deg=0.,native_to_GC=[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
LENGTH_CODE_TO_KPC = 1.00000064815577
VELOCITY_CODE_TO_KMS = 1.0
SKY_FIELD_HALF_WIDTH_DEG = 2.1
KINEMATIC_TRACER_DEFINITION = "current_dwarf_members_old_and_new"
MEMBER_SELECTION_ID = "current_spatial_properties_seeded_v2"


def production_output_root(folder_path):
    """Separate current spatial-member products from historical frozen-ID results."""
    import os
    from pathlib import Path
    default = Path(folder_path).resolve() / 'analysis' / CANONICAL_OBSERVER_ID / MEMBER_SELECTION_ID
    return Path(os.environ.get('DSPH_CANONICAL_OUTPUT_ROOT', str(default)))

def observer_contract_sha256():
    import hashlib, json
    data=dict(observer_id=CANONICAL_OBSERVER_ID,frame=OBSERVER_FRAME,
        length_code_to_kpc=LENGTH_CODE_TO_KPC,velocity_code_to_kms=VELOCITY_CODE_TO_KMS,
        host_bulk="TP2/3/4 massmean within true 3D r<5kpc of host-seeded local shrinking MW centre",
        default_motion="Heliocentric proper motion",display_half_width_deg=2.1,
        display_projection="raw DeltaRA wrapped and raw DeltaDec",walker="61direction spherical Rproj<5kpc",
        kinematic_tracer=KINEMATIC_TRACER_DEFINITION)
    return hashlib.sha256(json.dumps(data,sort_keys=True,separators=(",",":")).encode()).hexdigest()

def require_canonical_csv(path):
    """Refuse mixed/legacy CSV instead of relabelling or silently continuing it."""
    from pathlib import Path
    p=Path(path)
    if not p.exists(): return []
    comments=[]
    with p.open("r",encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"): comments.append(line.strip())
            elif line.strip(): break
    if canonical_csv_comment().strip() not in comments:
        raise ValueError("Legacy or different-observer CSV: use a NEW analysis output path and rerun GetInfo; existing CSV is preserved: "+str(p))
    return comments

def canonical_csv_comment():
    return ("# canonical_observer_id=" + CANONICAL_OBSERVER_ID +
            " member_selection_id=" + MEMBER_SELECTION_ID + "\n")

def production_observables(position_grf_kpc, velocity_grf_kms, particle_id):
    """Inputs already have H/B and units applied exactly once by center_on_MW."""
    core=astropy_native_direct_observables(position_grf_kpc,velocity_grf_kms,particle_id,
        host_center_raw=[0,0,0],host_bulk_raw=[0,0,0],length_to_kpc=1.,velocity_to_kms=1.,**OBSERVER_FRAME)
    h,u=core["h_native_kpc"],core["u_native_kms"]
    p,w=core["host_relative_position_kpc"],core["grf_velocity_native_kms"]
    result=dict(x=p[:,0],y=p[:,1],z=p[:,2],vx=w[:,0],vy=w[:,1],vz=w[:,2],
        xh=h[:,0],yh=h[:,1],zh=h[:,2],vxh=u[:,0],vyh=u[:,1],vzh=u[:,2],
        id=core["particle_id"],ra=core["RA_deg"],dec=core["Dec_deg"],rh=core["distance_kpc"],
        vlos=core["vlos_native_kms"],pmra=core["mu_alpha_cosdec_masyr"],pmdec=core["mu_delta_masyr"],
        canonical_observer_id=CANONICAL_OBSERVER_ID,member_selection_id=MEMBER_SELECTION_ID,
        observer_config_sha256=observer_contract_sha256(),_native=core)
    return result

def production_sky(coords, center_h_native_kpc):
    native=coords["_native"]
    return project_sky(native["h_native_kpc"],native["particle_id"],
        center_h_native_kpc=center_h_native_kpc,physical_to_icrs=native["physical_to_icrs"],aperture_kpc=5.)

def production_motion(coords, sky, mass_msun, selected):
    return aperture_motion_summaries(coords["_native"],sky,mass_msun,fit_mask=selected)["common_display_center_bulk_motion"]

def angle_velocity_grid(delta_ra_deg,delta_dec_deg,vlos_kms,v0_kms,npix=40):
    """Raw DeltaRA/DeltaDec display bins; footprint never defines Walker sample."""
    x,y,v=map(np.asarray,(delta_ra_deg,delta_dec_deg,vlos_kms))
    edges=np.linspace(-SKY_FIELD_HALF_WIDTH_DEG,SKY_FIELD_HALF_WIDTH_DEG,int(npix)+1)
    count,_,_=np.histogram2d(x,y,bins=(edges,edges))
    total,_,_=np.histogram2d(x,y,bins=(edges,edges),weights=v-float(v0_kms))
    field=np.divide(total,count,out=np.full_like(total,np.nan),where=count>0)
    return field,(-SKY_FIELD_HALF_WIDTH_DEG,SKY_FIELD_HALF_WIDTH_DEG,-SKY_FIELD_HALF_WIDTH_DEG,SKY_FIELD_HALF_WIDTH_DEG)

def tangent_kpc_to_raw_offsets(projected,x_kpc,y_kpc):
    """Exact inverse gnomonic mapping for contour/shape overlays on raw sky axes."""
    a,d=np.deg2rad([projected['center_RA_deg'],projected['center_Dec_deg']])
    n0=np.array([np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)])
    e=np.array([-np.sin(a),np.cos(a),0.]);north=np.array([-np.sin(d)*np.cos(a),-np.sin(d)*np.sin(a),np.cos(d)])
    x,y=np.broadcast_arrays(np.asarray(x_kpc),np.asarray(y_kpc))
    q=projected['center_distance_kpc']*n0+x[...,None]*e+y[...,None]*north
    ra=np.degrees(np.arctan2(q[...,1],q[...,0]));dec=np.degrees(np.arctan2(q[...,2],np.hypot(q[...,0],q[...,1])))
    return (ra-projected['center_RA_deg']+180)%360-180,dec-projected['center_Dec_deg']

def raw_sky_pixel_area_arcsec2(center_dec_deg,npix,half_width_deg=2.1):
    """Exact spherical area for each row of a constant raw-RA/Dec grid."""
    edges=np.radians(center_dec_deg+np.linspace(-half_width_deg,half_width_deg,npix+1))
    return np.radians(2*half_width_deg/npix)*np.diff(np.sin(edges))*(180./np.pi*3600.)**2

def frozen_cohort_mask(snapshot_particle_id,folder_path):
    """Explicit historical ID selection; never used by current production membership."""
    import os
    from pathlib import Path
    default=Path(folder_path).resolve()/'analysis'/CANONICAL_OBSERVER_ID/'cohort_ids.npz'
    path=Path(os.environ.get('DSPH_PREEXISTING_ID_PATH',str(default))).resolve()
    with np.load(str(path),allow_pickle=False) as data:
        cohort=_ids(data['particle_id'])
    current=_ids(snapshot_particle_id)
    mask=np.isin(current,cohort)
    return mask,dict(cohort_definition='explicit_historical_ID_reference',cohort_path=str(path),cohort_N=len(cohort),
        not_claimed_to_be_all_initial_IC_stars=True)

def write_snapshot_observer_metadata(folder_path,snapshot_num,host_state,cohort,star_coords,walker,motion,
                                     *, wind_tunnel=None):
    """Small atomic sidecar; source snapshot/particle arrays are never modified."""
    import os,json,tempfile
    from pathlib import Path
    root=production_output_root(folder_path)
    root.mkdir(parents=True,exist_ok=True)
    def serial(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [serial(v) for v in x]
        return x
    result=dict(canonical_observer_id=CANONICAL_OBSERVER_ID,observer_config_sha256=observer_contract_sha256(),
        observer_frame=OBSERVER_FRAME,member_selection_id=MEMBER_SELECTION_ID,
        host=host_state,membership=cohort,snapshot_num=int(snapshot_num),snapshot_source_directory=str(Path(folder_path)),
        kinematic_tracer=KINEMATIC_TRACER_DEFINITION,
        walker_R5={k:v for k,v in walker.items() if k not in ('fit_particle_id',)},motion_R5=motion)
    if wind_tunnel is not None:
        result['wind_tunnel'] = wind_tunnel
    target=root/('observer_snapshot_'+str(int(snapshot_num)).zfill(3)+'.json')
    with tempfile.NamedTemporaryFile('w',encoding='utf-8',dir=str(root),prefix='.observer-',suffix='.json',delete=False) as f:
        temp=Path(f.name);json.dump(serial(result),f,ensure_ascii=False,indent=2)
    os.replace(str(temp),str(target))
    return str(target)
