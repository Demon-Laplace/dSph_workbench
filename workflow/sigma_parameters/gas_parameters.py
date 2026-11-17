"""Compact gas/environment parameters from existing snapshots only.

Mass and phase statistics follow the earlier environment_timing measurement,
but optional fields and effective star-formation thresholds are independent.
No model-specific WJL density threshold or efficiency is assumed here.
"""
from pathlib import Path
import hashlib
import json
import math

import h5py
import numpy as np

RADII = (0.25, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 15.0, 20.0)
FLUX_RADII = (2.0, 5.0, 15.0, 20.0)
CHUNK = 262144
KPC_CM = 3.0856775814913673e21
MSUN_G = 1.989e33
YEAR_S = 31557600.0
MP_G = 1.6726e-24
KB_CGS = 1.38066e-16
XH = 0.76
GAMMA = 5.0 / 3.0
COLD_K = 2.0e4
DENSITY_PROXIES = ((0.1, "Mdensity_proxy_nHgt0p1_Msun"),
                   (0.5, "Mdensity_proxy_nHgt0p5_Msun"),
                   (1.0, "Mdensity_proxy_nHgt1_Msun"),
                   (10.0, "Mdensity_proxy_nHgt10_Msun"))


def _fields_hash_add(digest, name, value):
    array = np.ascontiguousarray(value)
    digest.update(name.encode("utf-8") + b"\0")
    digest.update(str(array.dtype).encode("ascii") + b"\0")
    digest.update(json.dumps(array.shape).encode("ascii") + b"\0")
    digest.update(memoryview(array).cast("B"))


def _params(ctx, units):
    result = {}
    path = units.get("unit_source") or ctx.get("config", {}).get("parameter_source")
    if path and Path(path).is_file():
        for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
            line = line.split("%", 1)[0].split("#", 1)[0].strip()
            fields = line.split()
            if len(fields) > 1:
                result[fields[0]] = fields[1]
    return result, path


def _effective_threshold(ctx, units, params):
    configuration = ctx.get("config", {})
    specifications = [configuration, ctx]
    for specification in specifications:
        for key in ("effective_sf_density_threshold_code", "effective_star_formation_density_threshold_code"):
            value = specification.get(key)
            if value is not None and np.isfinite(float(value)) and float(value) > 0:
                return float(value), dict(status="explicit_effective", source=key,
                    evidence=specification.get("effective_sf_threshold_evidence"))
        value = specification.get("effective_sf_nH_threshold_cm3")
        if value is not None and np.isfinite(float(value)) and float(value) > 0:
            density_factor = units["UnitMass_in_g"] / units["UnitLength_in_cm"] ** 3
            return float(value) * MP_G / (XH * density_factor), dict(
                status="explicit_effective", source="effective_sf_nH_threshold_cm3",
                evidence=specification.get("effective_sf_threshold_evidence"))
    raw = {name: params[name] for name in ("rho_threshold_sfr", "CritPhysDensity", "sfr_effiency", "SfEffPerFreeFall") if name in params}
    return None, dict(status="unavailable", reason="No authenticated effective compiled-branch density threshold; raw parameter values are not substituted", raw_parameters=raw)


def _total():
    names = ("MHI", "Mcold", "Mdense", "Mcold_dense", "SFR", "MSFpositive", "mU", "mFd", "m_nH", "m_T") + tuple(name for _, name in DENSITY_PROXIES)
    return dict(N=0, Mgas=0.0, values={name: 0.0 for name in names}, covered={name: 0 for name in names})


def _add(total, selected, mass, arrays, dense_threshold):
    count = int(np.count_nonzero(selected))
    if not count:
        return
    m = mass[selected]
    total["N"] += count
    total["Mgas"] += float(np.sum(m, dtype=np.float64))
    temperature = arrays["T"][selected]
    density_code = arrays["density_code"][selected]
    nh = arrays["neutral"][selected]
    sfr = arrays["sfr"][selected]
    u = arrays["u"][selected]
    fd = arrays["fd"][selected]
    nH = arrays["nH"][selected]
    valid_t = np.isfinite(temperature)
    valid_rho = np.isfinite(density_code)
    dense = density_code >= dense_threshold if dense_threshold is not None else np.zeros(count, bool)
    fields = {
        "MHI": (np.isfinite(nh), XH * m * nh),
        "Mcold": (valid_t, m * (temperature < COLD_K)),
        "Mdense": (valid_rho & (dense_threshold is not None), m * dense),
        "Mcold_dense": (valid_t & valid_rho & (dense_threshold is not None), m * dense * (temperature < COLD_K)),
        "SFR": (np.isfinite(sfr), sfr),
        "MSFpositive": (np.isfinite(sfr), m * (sfr > 0)),
        "mU": (np.isfinite(u), m * u),
        "mFd": (np.isfinite(fd), m * fd),
        "m_nH": (np.isfinite(nH), m * nH),
        "m_T": (valid_t, m * temperature),
    }
    fields.update({name: (np.isfinite(nH), m * (nH > threshold)) for threshold, name in DENSITY_PROXIES})
    for name, (good, values) in fields.items():
        total["covered"][name] += int(np.count_nonzero(good))
        total["values"][name] += float(np.sum(values[good], dtype=np.float64))


def _finish(total, membership_supported=True):
    result = dict(N=total["N"], Mgas_Msun=total["Mgas"] if membership_supported else None,
                  coverage_N=dict(total["covered"]), membership_supported=bool(membership_supported))
    for name, value in total["values"].items():
        covered = total["covered"][name] == total["N"]
        result[name] = value if covered and membership_supported else None
    mass = result["Mgas_Msun"]
    for name, numerator in (("nH_massmean_cm3", "m_nH"), ("T_massmean_K", "m_T")):
        result[name] = result[numerator] / mass if mass and result[numerator] is not None else None
    result["cold_depletion_time_Gyr"] = result["Mcold"] / result["SFR"] / 1e9 if result["SFR"] and result["Mcold"] is not None else None
    result["dense_depletion_time_Gyr"] = result["Mdense"] / result["SFR"] / 1e9 if result["SFR"] and result["Mdense"] is not None else None
    result["feedback_over_thermal_stored"] = result["mFd"] / result["mU"] if result["mU"] and result["mFd"] is not None else None
    return result


def _ambient(a, lo, hi, axis=None, upstream=True):
    radius = a["r"]
    selected = (radius >= lo) & (radius < hi)
    if axis is not None:
        cosine = np.divide(a["x"] @ axis, radius, out=np.zeros(len(radius)), where=radius > 0)
        selected &= (cosine > 0.5) if upstream else (cosine < -0.5)
    omega = np.pi if axis is not None else 4 * np.pi
    volume = omega * (hi ** 3 - lo ** 3) / 3
    number = int(np.count_nonzero(selected))
    result = dict(N=number, supported=number >= 64, rlo_kpc=float(lo), rhi_kpc=float(hi), volume_kpc3=float(volume))
    if number == 0:
        result["reason"] = "No gas in the spatial ambient region"
        return result
    mass = a["m"][selected]
    velocity = a["v"][selected]
    good_v = np.isfinite(velocity).all(axis=1)
    total_mass = float(np.sum(mass, dtype=np.float64))
    rho_geom = total_mass * MSUN_G / (volume * KPC_CM ** 3)
    result.update(Mgas_Msun=total_mass, nH_geom_cm3=float(rho_geom * XH / MP_G))
    if np.all(good_v) and number >= 2:
        mean = np.sum(velocity.astype(np.float64) * mass[:, None], axis=0, dtype=np.float64) / total_mass
        relative = np.asarray(a["bulk"], dtype=float) - mean
        speed = float(np.linalg.norm(relative))
        result.update(velocity_native_kms=mean.tolist(), relative_velocity_native_kms=relative.tolist(),
                      vrel_kms=speed, P_ram_bulk_dyn_cm2=float(rho_geom * (speed * 1e5) ** 2))
    else:
        result.update(supported=False, reason="Ambient velocity field missing/invalid or fewer than two particles",
                      velocity_native_kms=None, relative_velocity_native_kms=None,
                      vrel_kms=None, P_ram_bulk_dyn_cm2=None)
    for field, name in (("u", "P_thermal_dyn_cm2"), ("fd", "P_feedback_stored_proxy_dyn_cm2")):
        values = a[field][selected]
        result[name] = float((GAMMA - 1) * np.sum(mass * values, dtype=np.float64) * MSUN_G * a["UV2"] / (volume * KPC_CM ** 3)) if np.isfinite(values).all() else None
    temperature = a["T"][selected]
    if np.isfinite(temperature).all():
        hot = temperature >= COLD_K
        result["T_massmean_K"] = float(np.sum(mass * temperature, dtype=np.float64) / total_mass)
        result["hot_mass_fraction"] = float(np.sum(mass[hot], dtype=np.float64) / total_mass)
    else:
        result["T_massmean_K"] = None
        result["hot_mass_fraction"] = None
    return result


def _shape(a, selected):
    selected &= np.isfinite(a["v"]).all(axis=1)
    number = int(np.count_nonzero(selected))
    if number < 32:
        return dict(N=number, supported=False, reason="Fewer than 32 gas particles with phase fields")
    x, velocity, mass = a["x"][selected], a["v"][selected], a["m"][selected]
    total = np.sum(mass, dtype=np.float64)
    center = np.sum(x.astype(np.float64) * mass[:, None], axis=0, dtype=np.float64) / total
    mean = np.sum(velocity.astype(np.float64) * mass[:, None], axis=0, dtype=np.float64) / total
    delta = x - center
    covariance = (delta * mass[:, None]).T @ delta / total
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    angular = np.sum(np.cross(delta, velocity - mean) * mass[:, None], axis=0, dtype=np.float64)
    norm = np.linalg.norm(angular)
    return dict(N=number, supported=True, covariance_kpc2=covariance.tolist(),
                eigenvalues_kpc2=eigenvalues.tolist(), minor_axis_native=eigenvectors[:, 0].tolist(),
                q_3d=float(np.sqrt(max(0, eigenvalues[0]) / eigenvalues[2])) if eigenvalues[2] > 0 else None,
                angular_momentum_Msun_kpc_kms=angular.tolist(),
                spin_axis_native=(angular / norm).tolist() if norm > 0 else None,
                COM_offset_kpc=center.tolist())


def analyze_gas(ctx, n, center_raw, bulk_native_kms, roster, units):
    """Return gas statistics; optional missing fields do not fail stellar work.

    ``center_raw`` is in snapshot coordinate units. ``bulk_native_kms`` is
    the dwarf velocity in the simulation's native axes (no solar subtraction).
    ``roster`` is ignored. Membership is determined entirely by spatial
    apertures. Particle IDs are never read, checked, or hashed by this module.
    """
    length = float(units["length_to_kpc"])
    velocity_factor = float(units["velocity_to_kms"])
    mass_factor = float(units["mass_to_msun"])
    length_cm = float(units.get("UnitLength_in_cm", length * KPC_CM))
    mass_g = float(units.get("UnitMass_in_g", mass_factor * MSUN_G))
    velocity_cms = float(units.get("UnitVelocity_in_cm_per_s", velocity_factor * 1e5))
    resolved_units = dict(units, UnitLength_in_cm=length_cm, UnitMass_in_g=mass_g)
    params, parameter_source = _params(ctx, resolved_units)
    dense_threshold, threshold_record = _effective_threshold(ctx, resolved_units, params)
    threshold_nh = dense_threshold * mass_g / length_cm ** 3 * XH / MP_G if dense_threshold is not None else None
    inventories = {str(radius): _total() for radius in RADII}
    annulus5_15 = _total()
    flux = {str(radius): dict(N=0, outward_Msun_yr=0.0, inward_Msun_yr=0.0, net_outward_Msun_yr=0.0, velocity_missing_N=0) for radius in FLUX_RADII}
    names = ("x", "r", "v", "m", "u", "fd", "T")
    ambient_lists = {name: [] for name in names}
    own_lists = {name: [] for name in names}
    availability = {}
    errors = []
    stats = []
    field_hash = hashlib.sha256()
    counts = dict(gas_particles_read=0, invalid_position_or_mass=0)
    total_time = None
    source_changed = False
    host_center = np.asarray(ctx.get("snapshot_host_center_raw", []), dtype=float)
    host_radius = ctx.get("snapshot_Rgc_kpc")
    host_geometry_supported = bool(host_center.shape == (3,) and np.isfinite(host_center).all()
                                   and host_radius is not None and np.isfinite(float(host_radius)) and float(host_radius) > 0)
    host_gas_mass = 0.0
    host_gas_number = 0
    paths = ctx["snapshot_paths"][str(int(n))]
    if isinstance(paths, (str, Path)):
        paths = [paths]
    for pp in paths:
        path = Path(pp)
        file_stat = dict(path=str(path), stat_available=False)
        stats.append(file_stat)
        try:
            before = path.stat()
            file_stat.update(stat_available=True, size_bytes=before.st_size, mtime_ns=before.st_mtime_ns,
                             device=str(before.st_dev), inode=str(before.st_ino))
            expected = next((record for record in ctx.get("input_snapshot_fingerprints", {}).get(str(int(n)), []) if record.get("path") == str(path)), None)
            if expected and (before.st_size, before.st_mtime_ns) != (expected["size_bytes"], expected["mtime_ns"]):
                source_changed = True
                file_stat.update(changed_before_read=True, expected_size_bytes=expected["size_bytes"], expected_mtime_ns=expected["mtime_ns"])
                errors.append(dict(path=str(path), reason="Frozen snapshot changed before gas read", expected=expected, observed=dict(file_stat)))
                continue
            with h5py.File(str(path), "r") as snapshot:
                total_time = float(snapshot["Header"].attrs["Time"]) * float(units["time_to_gyr"])
                if "PartType0" not in snapshot:
                    continue
                gas = snapshot["PartType0"]
                if "Coordinates" not in gas:
                    errors.append(dict(path=str(path), reason="Gas coordinates absent; gas statistics unavailable for this part"))
                    continue
                number = len(gas["Coordinates"])
                table_mass = float(snapshot["Header"].attrs["MassTable"][0])
                field_hash.update((str(path) + ":PartType0:").encode("utf-8"))
                for name in ("Coordinates", "Velocities", "Masses", "Density", "InternalEnergy", "ElectronAbundance", "NeutralHydrogenAbundance", "StarFormationRate", "default"):
                    count = availability.setdefault(name, dict(parts_present=0, parts_absent=0))
                    count["parts_present" if name in gas else "parts_absent"] += 1
                    field_hash.update((name + ":present=" + str(name in gas)).encode("ascii"))
                for first in range(0, number, CHUNK):
                    sl = slice(first, min(number, first + CHUNK))
                    size = sl.stop - sl.start
                    raw = {}
                    for name in ("Coordinates", "Velocities", "Masses", "Density", "InternalEnergy", "ElectronAbundance", "NeutralHydrogenAbundance", "StarFormationRate", "default"):
                        if name in gas:
                            raw[name] = gas[name][sl]
                            _fields_hash_add(field_hash, name, raw[name])
                    x = (np.asarray(raw["Coordinates"], dtype=np.float64) - np.asarray(center_raw, dtype=np.float64)) * length
                    x = x.astype(np.float32)
                    mass = (np.asarray(raw["Masses"], dtype=np.float64) if "Masses" in raw else np.full(size, table_mass)) * mass_factor
                    good = np.isfinite(x).all(axis=1) & np.isfinite(mass) & (mass > 0)
                    counts["gas_particles_read"] += size
                    counts["invalid_position_or_mass"] += int(np.count_nonzero(~good))
                    if not np.any(good):
                        continue
                    x, mass = x[good], mass[good]
                    radius = np.linalg.norm(x, axis=1)
                    if host_geometry_supported:
                        host_distance = np.linalg.norm((np.asarray(raw["Coordinates"], dtype=np.float64)[good] - host_center) * length, axis=1)
                        host_selected = (radius >= 15.0) & (host_distance < float(host_radius))
                        host_gas_mass += float(np.sum(mass[host_selected], dtype=np.float64))
                        host_gas_number += int(np.count_nonzero(host_selected))
                    vv = np.asarray(raw["Velocities"], dtype=np.float32)[good] * velocity_factor if "Velocities" in raw else np.full_like(x, np.nan)
                    def optional(name, minimum=None, maximum=None):
                        value = np.asarray(raw[name], dtype=np.float32)[good] if name in raw else np.full(len(mass), np.nan, np.float32)
                        valid = np.isfinite(value)
                        if minimum is not None:
                            valid &= value >= minimum
                        if maximum is not None:
                            valid &= value <= maximum
                        return np.where(valid, value, np.nan).astype(np.float32)
                    density = optional("Density", minimum=0.0)
                    thermal = optional("InternalEnergy", minimum=0.0)
                    electrons = optional("ElectronAbundance", minimum=0.0)
                    neutral = optional("NeutralHydrogenAbundance", minimum=0.0, maximum=1.0)
                    sfr = optional("StarFormationRate", minimum=0.0)
                    fd = optional("default", minimum=0.0)
                    mu = 4.0 / (1.0 + 3.0 * XH + 4.0 * XH * electrons)
                    temperature = ((GAMMA - 1) * mu * MP_G / KB_CGS * thermal * velocity_cms ** 2).astype(np.float32)
                    nH = (density * (mass_g / length_cm ** 3 * XH / MP_G)).astype(np.float32)
                    arrays = dict(T=temperature, density_code=density, neutral=neutral, sfr=sfr, u=thermal, fd=fd, nH=nH)
                    for rad in RADII:
                        within = radius < rad
                        _add(inventories[str(rad)], within, mass, arrays, dense_threshold)
                    _add(annulus5_15, (radius >= 5.0) & (radius < 15.0), mass, arrays, dense_threshold)
                    for rad in FLUX_RADII:
                        shell = (radius >= 0.9 * rad) & (radius < 1.1 * rad)
                        velocity_ok = np.isfinite(vv).all(axis=1)
                        flux[str(rad)]["velocity_missing_N"] += int(np.count_nonzero(shell & ~velocity_ok))
                        shell &= velocity_ok
                        vr = np.sum(x[shell].astype(float) * (vv[shell] - np.asarray(bulk_native_kms)), axis=1) / radius[shell]
                        terms = mass[shell] * vr / (0.2 * rad) * (1e5 / KPC_CM * YEAR_S)
                        flux[str(rad)]["N"] += int(np.count_nonzero(shell))
                        flux[str(rad)]["outward_Msun_yr"] += float(np.sum(terms[terms > 0], dtype=np.float64))
                        flux[str(rad)]["inward_Msun_yr"] += float(-np.sum(terms[terms < 0], dtype=np.float64))
                        flux[str(rad)]["net_outward_Msun_yr"] += float(np.sum(terms, dtype=np.float64))
                    ambient_mask = (radius >= 15) & (radius < 60)
                    own_mask = radius < 15
                    phase = dict(x=x, r=radius, v=vv, m=mass, u=thermal, fd=fd, T=temperature)
                    for name in names:
                        ambient_lists[name].append(phase[name][ambient_mask].copy())
                        own_lists[name].append(phase[name][own_mask].copy())
            after = path.stat()
            if (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns) != (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns):
                source_changed = True
                file_stat["changed_during_read"] = True
                errors.append(dict(path=str(path), reason="Snapshot changed during gas read", before=dict(file_stat),
                                   after=dict(size_bytes=after.st_size, mtime_ns=after.st_mtime_ns, device=str(after.st_dev), inode=str(after.st_ino))))
        except (OSError, KeyError, ValueError) as exc:
            file_stat["error"] = type(exc).__name__ + ": " + str(exc)
            errors.append(dict(path=str(path), reason=type(exc).__name__ + ": " + str(exc)))
    inventory_result = {}
    for radius, total in inventories.items():
        component = _finish(total)
        if dense_threshold is None:
            component["Mdense"] = None
            component["Mcold_dense"] = None
        inventory_result[radius] = dict(all=component, dwarf_spatial=component)
    annulus_result = _finish(annulus5_15)
    if dense_threshold is None:
        annulus_result["Mdense"] = None
        annulus_result["Mcold_dense"] = None
    for record in flux.values():
        record["supported"] = bool(record["N"] >= 20 and not record["velocity_missing_N"])
    def joined(parts):
        return {name: np.concatenate(values) if values else np.empty((0, 3) if name in ("x", "v") else (0,), dtype=np.float32) for name, values in parts.items()}
    ambient = joined(ambient_lists)
    ambient.update(bulk=np.asarray(bulk_native_kms), UV2=velocity_cms ** 2)
    seed = _ambient(ambient, 30, 60)
    axis = None
    if seed.get("relative_velocity_native_kms"):
        vector = np.asarray(seed["relative_velocity_native_kms"])
        norm = np.linalg.norm(vector)
        if norm > 0:
            axis = vector / norm
    environment = dict(supported=bool(axis is not None), seed_shell30_60=seed,
                       membership="ambient_spatial; all gas in spatial shells/cones outside the loose r<15 kpc aperture",
                       wind_direction_native=axis.tolist() if axis is not None else None)
    if axis is not None:
        environment.update(upstream20_40=_ambient(ambient, 20, 40, axis),
                           upstream15_25=_ambient(ambient, 15, 25, axis),
                           upstream30_60=_ambient(ambient, 30, 60, axis),
                           downstream20_40=_ambient(ambient, 20, 40, axis, False))
    own = joined(own_lists)
    shape = _shape(own, own["r"] < 5)
    coldshape = _shape(own, (own["r"] < 5) & np.isfinite(own["T"]) & (own["T"] < COLD_K))
    loose_shape = _shape(own, own["r"] < 15)
    loose_coldshape = _shape(own, (own["r"] < 15) & np.isfinite(own["T"]) & (own["T"] < COLD_K))
    for record in (shape, coldshape, loose_shape, loose_coldshape):
        if axis is not None and record.get("supported"):
            record["cos_incidence_upstream"] = float(abs(np.asarray(record["minor_axis_native"]) @ axis))
    host_supported = bool(host_geometry_supported and not errors and not counts["invalid_position_or_mass"])
    host_status = ("COMPLETE" if host_supported else "UNAVAILABLE_HOST_GEOMETRY" if not host_geometry_supported
                   else "PARTIAL_SOURCE")
    return dict(contract_id="sigma_spatial_gas_parameters_v1", snapshot=int(n), time_Gyr=total_time,
                status="INVALID_SOURCE" if source_changed else ("PARTIAL" if errors or counts["invalid_position_or_mass"] else "COMPLETE"),
                gas=inventory_result, dwarf_spatial_gas_r20=inventory_result["20.0"]["dwarf_spatial"],
                gas_flux=flux, environment=environment, gas_shape_r5=shape, cold_gas_shape_r5=coldshape,
                membership_strict_r5=inventory_result["5.0"]["dwarf_spatial"],
                membership_loose_r15=inventory_result["15.0"]["dwarf_spatial"],
                annulus5_15=annulus_result, gas_shape_r15=loose_shape, cold_gas_shape_r15=loose_coldshape,
                host_enclosed_gas_mass_Msun=host_gas_mass if host_supported else None,
                host_enclosed_gas_N=host_gas_number if host_supported else None,
                host_enclosed_gas_mass_status=host_status,
                effective_sf_density_threshold_code=dense_threshold, effective_sf_nH_threshold_cm3=threshold_nh,
                effective_sf_threshold=threshold_record, field_availability=availability, counts=counts,
                membership_status=dict(supported=True, definition="dwarf-centred spatial apertures; all gas irrespective of formation/source", particle_IDs_read=False),
                input_snapshot_stats=stats, source_gas_fields_sha256=field_hash.hexdigest(), errors=errors,
                source_unchanged=not source_changed and not any("error" in item for item in stats),
                source_module_path=str(Path(__file__).resolve()), source_module_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                definitions=dict(radii="3D dwarf-centred kpc; R is reserved for projected stellar radius",
                    gas_membership="dwarf_spatial includes all gas inside each dwarf-centred physical radius; all is an exact synonymous compatibility entry; strict r<5 kpc, loose r<15 kpc, annulus 5<=r<15 kpc; ambient_spatial is spatial gas in 15-60 kpc shells/cones",
                    HI="0.76 * gas mass * NeutralHydrogenAbundance, both temperature phases",
                    cold_temperature_K=COLD_K, dense="Only a separately authenticated effective density threshold is applied",
                    fixed_density_proxies="Mass with nH strictly >0.1,0.5,1,10 cm^-3; fixed physical diagnostics, not actual star-formation eligibility or threshold",
                    SFR="Snapshot StarFormationRate summed in Msun/yr; global-box sfr.txt is not used",
                    feedback="mFd and pressure proxy use the stored default field; interpretation as IO_FD requires run/source confirmation",
                    energy_unit="mU/mFd are Msun times stored specific-energy code units; pressure uses UnitVelocity^2",
                    ambient="20-40 kpc upstream 60-degree cone; wind from 30-60 kpc mass bulk; no solar subtraction",
                    host_enclosed_gas="All gas with |position-host_center| < current Rgc and dwarf-centred r>=15 kpc; exclusion is the loose spatial aperture, not provenance; spherical host-gravity proxy",
                    flux="Instantaneous all-local gas shell sum m*v_r/(0.2*r); 0.9r<=radius<1.1r, not permanent escape",
                    source_hash="Ordered raw physical gas fields actually read with field name/dtype/shape/bytes; ParticleIDs excluded; not whole-file HDF5 SHA"),
                parameter_source=parameter_source)
