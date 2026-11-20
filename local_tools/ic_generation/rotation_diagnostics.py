#!/usr/bin/env python3
"""Measure intrinsic rotation in a new IC or selected isolation snapshots.

This is a scientific measurement, not a transfer/integrity gate. Membership is
spatial and mass weighted; particle IDs are neither read nor required.
"""
import argparse
import csv
import json
from pathlib import Path
import re
import subprocess

import numpy as np


def average(a, m):
    return np.average(a, weights=m, axis=0)


def frame(x, v, m):
    center = average(x, m)
    for radius in (5., 3., 2., 1.5):
        use = np.sum((x-center)**2, axis=1) < radius**2
        if use.sum() >= 30:
            center = average(x[use], m[use])
    dx = x-center
    core = np.sum(dx*dx, axis=1) < 25.
    bulk = average(v[core], m[core])
    shape = np.sum(dx*dx, axis=1) < 4.
    xx = dx[shape]
    tensor = (xx.T*m[shape]) @ xx / m[shape].sum()
    values, axes = np.linalg.eigh(tensor)
    axis = axes[:, 0]
    if axis[2] < 0:
        axis = -axis
    return center, bulk, axis, np.sqrt(np.maximum(values, 0)/values[-1])


def quantile(r, m, fraction):
    order = np.argsort(r)
    cumulative = (np.cumsum(m[order])-.5*m[order])/m.sum()
    return float(np.interp(fraction, cumulative, r[order]))


def moments(x, v, m, axis, mask):
    x, v, m = x[mask], v[mask], m[mask]
    if not len(m):
        return None
    r = np.linalg.norm(x, axis=1)
    z = x @ axis
    planar = x-z[:, None]*axis
    R = np.linalg.norm(planar, axis=1)
    er = planar/np.maximum(R[:, None], 1e-15)
    ephi = np.cross(axis, er)
    vR = np.sum(v*er, axis=1)
    vp = np.sum(v*ephi, axis=1)
    vz = v @ axis
    cylindrical = np.column_stack((vR, vp, vz))
    mean = average(cylindrical, m)
    sigma = np.sqrt(average((cylindrical-mean)**2, m))
    vr = np.sum(v*x, axis=1)/np.maximum(r, 1e-15)
    vt = (z*vR-R*vz)/np.maximum(r, 1e-15)
    spherical = np.column_stack((vr, vt, vp))
    spherical_mean = average(spherical, m)
    second = average(spherical*spherical, m)
    dispersion2 = average((spherical-spherical_mean)**2, m)
    return {
        'N': len(m), 'mass_Msun': float(m.sum()*1e10),
        'mean_vphi_kms': float(mean[1]),
        'mean_abs_vphi_kms': float(average(np.abs(vp), m)),
        'mean_vR_kms': float(mean[0]), 'mean_vz_kms': float(mean[2]),
        'sigma_R_phi_z_kms': sigma.tolist(),
        'V_over_sigma_phi': float(abs(mean[1])/sigma[1]) if sigma[1] > 0 else None,
        'prograde_mass_fraction': float(m[vp > 0].sum()/m.sum()),
        'specific_L_axis_kpc_kms': float(average(R*vp, m)),
        'specific_L_xyz_kpc_kms': average(np.cross(x, v), m).tolist(),
        'signed_L_over_absL': float(np.sum(m*R*vp)/np.sum(m*R*np.abs(vp)))
            if np.sum(m*R*np.abs(vp)) > 0 else None,
        'r10_r50_r90_kpc': [quantile(r, m, f) for f in (.1, .5, .9)],
        'beta_second_spherical': float(1-(second[1]+second[2])/(2*second[0]))
            if second[0] > 0 else None,
        'beta_sigma_spherical': float(1-(dispersion2[1]+dispersion2[2])/(2*dispersion2[0]))
            if dispersion2[0] > 0 else None,
    }


def diagnose(populations, label, time_code):
    center, bulk, axis, ratios = frame(*populations['old'])
    result = {'label': label, 'time_code': time_code,
              'center_kpc': center.tolist(), 'bulk_kms': bulk.tolist(),
              'density_short_axis': axis.tolist(),
              'axis_tilt_from_Z_deg': float(np.degrees(np.arccos(np.clip(axis[2], -1, 1)))),
              'density_covariance_2kpc_axis_ratios': ratios.tolist(), 'populations': {}}
    for name, (x, v, m) in populations.items():
        x, v = x-center, v-bulk
        radius = np.linalg.norm(x, axis=1)
        R = np.linalg.norm(x-(x@axis)[:, None]*axis, axis=1)
        core = radius < 5.
        pop = {'3D5': moments(x, v, m, axis, core),
               'annuli_within_3D5': []}
        for lo, hi in zip((0., .25, .5, 1., 2.), (.25, .5, 1., 2., 5.)):
            values = moments(x, v, m, axis, core & (R >= lo) & (R < hi))
            if values is not None:
                pop['annuli_within_3D5'].append({'Rmin_kpc': lo, 'Rmax_kpc': hi, **values})
        if name == 'gas':
            pop['3D15'] = moments(x, v, m, axis, radius < 15.)
            pop['shell_5_15'] = moments(x, v, m, axis, (radius >= 5.) & (radius < 15.))
        result['populations'][name] = pop
    return result


def read_build(directory):
    build = json.loads((directory/'build.json').read_text())
    populations, libraries = {}, {}
    for name in ('stellar', 'gas'):
        data = np.loadtxt(directory/('model_component_'+name+'.nbody'), ndmin=2)
        populations['old' if name == 'stellar' else 'gas'] = (data[:, :3], data[:, 3:6], data[:, 6])
        with np.load(directory/('model_component_'+name+'.data.npz')) as archive:
            ic, weight = archive['ic'], archive['weights']
            Lz = ic[:, 0]*ic[:, 4]-ic[:, 1]*ic[:, 3]
            vp = Lz/np.maximum(np.hypot(ic[:, 0], ic[:, 1]), 1e-15)
            libraries[name] = {
                'N_orbits': len(ic),
                'equal_seed_prograde_fraction': float(np.mean(Lz > 0)),
                'optimized_prograde_weight_fraction': float(weight[Lz > 0].sum()/weight.sum()),
                'equal_seed_mean_vphi_kms': float(np.mean(vp)),
                'optimized_seed_mean_vphi_kms': float(average(vp, weight)),
                'note': 'Seed vphi is not an orbit-time-averaged rotation curve.'}
    return {'model': Path(build['source_ini']).stem.removeprefix('IC_'),
            'build_id': build.get('build_id'), 'generator_revision': build.get('generator_revision'),
            'seed': build.get('seed'), 'sampling_recipe': build.get('sampling_recipe'),
            'input_build': str(directory), 'orbit_libraries': libraries,
            'states': [diagnose(populations, 'IC', 0.)]}


def read_snapshots(directory, selected):
    import h5py
    files = {int(p.stem.split('_')[-1]): p for p in directory.glob('snapshot_*.hdf5')
             if re.fullmatch(r'snapshot_\d+', p.stem)}
    if not files:
        raise ValueError('No single-file HDF5 snapshots in '+str(directory))
    states = []
    for number in sorted(set(selected+[min(files), max(files)]).intersection(files)):
        with h5py.File(files[number], 'r') as stream:
            populations = {}
            for name, kind in [('old', 2), ('new', 4), ('gas', 0)]:
                group = stream.get('PartType'+str(kind))
                if group is None or not len(group['Coordinates']):
                    continue
                masses = group['Masses'][:] if 'Masses' in group else np.full(
                    len(group['Coordinates']), stream['Header'].attrs['MassTable'][kind])
                populations[name] = (group['Coordinates'][:], group['Velocities'][:], masses)
            states.append(diagnose(populations, 'snapshot_%03d' % number,
                                   float(stream['Header'].attrs['Time'])))
    report = {'snapshot_directory': str(directory), 'states': states}
    request = directory.parent/'request.json'
    if request.is_file():
        identity = json.loads(request.read_text())
        report['run_identity'] = {key: identity.get(key) for key in (
            'model', 'name', 'build_id', 'ic_version', 'parameter_version', 'profile')}
        report['model'] = identity.get('model')
    return report


def write_report(report, prefix, target_time=None):
    report.update(schema='intrinsic_rotation_v1',
                  units={'length': 'kpc', 'velocity': 'km/s', 'mass': 'Msun', 'time': 'solver_code_time'},
                  membership='Current physical particles, spatial apertures, mass weighted; no IDs.',
                  axis_definition='Old-star density covariance short axis within 2 kpc; sign toward +Z.',
                  frame_definition='Shrinking old-star center (5/3/2/1.5 kpc), old-star bulk within 5 kpc.',
                  validation_limit='Rotation measurement only; not proof of long-term or gas equilibrium.')
    states = report['states']
    if len(states) >= 2:
        first = states[0]['populations']['old']['3D5']
        last = states[-1]['populations']['old']['3D5']
        report['evolution_summary'] = {
            'first_time_code': states[0]['time_code'], 'last_time_code': states[-1]['time_code'],
            'old_V_first_last_kms': [first['mean_vphi_kms'], last['mean_vphi_kms']],
            'old_V_over_sigma_phi_first_last': [first['V_over_sigma_phi'], last['V_over_sigma_phi']],
            'old_r50_fractional_change': last['r10_r50_r90_kpc'][1]/first['r10_r50_r90_kpc'][1]-1,
            'old_V_fractional_change': last['mean_vphi_kms']/first['mean_vphi_kms']-1
                if abs(first['mean_vphi_kms']) > 1e-3 else None}
    if target_time is not None:
        report['requested_end_time_code'] = target_time
        report['snapshot_end_time_reached'] = states[-1]['time_code'] >= target_time*(1-1e-6)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    prefix.with_suffix('.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    fields = ['label', 'time_code', 'population', 'aperture', 'N', 'mass_Msun', 'mean_vphi_kms',
              'V_over_sigma_phi', 'prograde_mass_fraction', 'specific_L_axis_kpc_kms',
              'sigma_R_kms', 'sigma_phi_kms', 'sigma_z_kms', 'r50_kpc']
    with prefix.with_suffix('.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for state in states:
            for name, population in state['populations'].items():
                for aperture in ('3D5', '3D15', 'shell_5_15'):
                    values = population.get(aperture)
                    if values is None:
                        continue
                    row = {key: values[key] for key in fields if key in values}
                    row.update(label=state['label'], time_code=state['time_code'],
                               population=name, aperture=aperture,
                               r50_kpc=values['r10_r50_r90_kpc'][1])
                    row.update(dict(zip(('sigma_R_kms', 'sigma_phi_kms', 'sigma_z_kms'),
                                        values['sigma_R_phi_z_kms'])))
                    writer.writerow(row)
    for state in states:
        old = state['populations']['old']['3D5']
        print('%s t=%.8g: old Vphi=%.5g km/s, sigma_phi=%.5g, V/sigma_phi=%.5g, prograde=%.4f'
              % (state['label'], state['time_code'], old['mean_vphi_kms'],
                 old['sigma_R_phi_z_kms'][1], old['V_over_sigma_phi'],
                 old['prograde_mass_fraction']), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--build', type=Path)
    source.add_argument('--snapshots', type=Path)
    source.add_argument('--run-solver', type=Path, metavar='PARAMETERS',
                        help='Within the existing allocation: run MPI solver, then measure its new snapshots.')
    parser.add_argument('--solver', type=Path)
    parser.add_argument('--selected', default='0,10,25,50')
    parser.add_argument('--output-prefix', required=True, type=Path)
    args = parser.parse_args()
    target_time = None
    if args.run_solver:
        if args.solver is None:
            parser.error('--run-solver needs --solver')
        code = subprocess.run(['mpirun', str(args.solver), str(args.run_solver)]).returncode
        if code:
            return code
        args.snapshots = args.run_solver.resolve().parent/'output'
        for line in args.run_solver.read_text().splitlines():
            fields = line.split('%', 1)[0].split()
            if len(fields) >= 2 and fields[0] == 'TimeMax':
                target_time = float(fields[1])
    report = read_build(args.build) if args.build else read_snapshots(
        args.snapshots, [int(s) for s in args.selected.split(',') if s])
    write_report(report, args.output_prefix, target_time)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
