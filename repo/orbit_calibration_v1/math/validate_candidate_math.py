"""Independent numeric checks and actual current C-merger probe comparison."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import shutil
import subprocess

import numpy as np
from scipy.spatial.transform import Rotation
import candidate_math as m

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--reference-source', type=Path, required=True,
                    help='Actual current ic_merge.c source; no historical implicit source path')
parser.add_argument('--work-name', default='validation_v2')
parser.add_argument('--qa-name', default='MATH_QA.json')
args = parser.parse_args()
assert Path(args.work_name).name == args.work_name
assert Path(args.qa_name).name == args.qa_name
WORK = HERE / args.work_name
WORK.mkdir(exist_ok=False)
rng = np.random.default_rng(70047005)
checks = {}


def close(a, b, tolerance):
    error = float(np.max(abs(np.asarray(a) - np.asarray(b))))
    assert error <= tolerance, (error, tolerance)
    return error


# Algebraic calibration recovery, using independently constructed final state.
r0 = np.array([-192.3313, 164.1997, -183.4200])
v0 = np.array([34.7655, -86.0286, 0.1016])
seed = m.decompose_state(r0, v0)
ce, kl, radius = 12000.0, 0.97, 139.4
vt = kl * np.linalg.norm(np.cross(r0, v0)) / radius
speed = np.sqrt(v0 @ v0 + ce)
vr = -np.sqrt(speed ** 2 - vt ** 2)
rotation = Rotation.from_rotvec([0.2, -0.3, 0.1]).as_matrix()
rf = radius * rotation @ seed['rhat']
vf = rotation @ (vr * seed['rhat'] + vt * seed['that'])
coeff = m.empirical_coefficients(r0, v0, rf, vf)
checks['coefficient_CE_error'] = close(coeff['C_E'], ce, 2e-11)
checks['coefficient_kL_error'] = close(coeff['k_L'], kl, 1e-14)
solved = m.solve_initial_velocity(r0, v0, rf, vf, coeff['C_E'], coeff['k_L'], radial_sign=-1)
checks['recovered_initial_velocity_error'] = close(solved['velocity'], v0, 1e-12)
checks['reconstructed_target_L_error'] = close(solved['reconstructed_target_angular_momentum'], np.linalg.norm(np.cross(rf, vf)), 1e-11)
try:
    m.solve_initial_velocity(r0, v0, rf, vf, 1e9, kl)
except ValueError:
    checks['infeasible_energy_rejected'] = True
else:
    raise AssertionError('infeasible initial energy was accepted')

# Frame map must be proper and preserve all three directed basis vectors.
p = m.orbit_plane_rotation(radius * seed['rhat'], vr * seed['rhat'] + vt * seed['that'], rf, vf)
checks['orbit_frame_map_error'] = close(p, rotation, 1e-14)
checks['orbit_frame_det_error'] = close(np.linalg.det(p), 1, 1e-14)

# North->east PA sign, for a non-axis-aligned physical sky basis.
basis = Rotation.from_rotvec([0.8, -0.2, 0.3]).as_matrix()
east, north, los = basis.T
pa, target_pa = 29.46794, 43.0
major = np.sin(np.deg2rad(pa)) * east + np.cos(np.deg2rad(pa)) * north
minor = np.cos(np.deg2rad(pa)) * east - np.sin(np.deg2rad(pa)) * north
delta = m.pa_difference_deg(target_pa, pa)
a = m.pa_roll_rotation(los, delta)
checks['PA_roll_target_error_deg'] = abs(m.pa_difference_deg(target_pa, m.axis_pa_deg(a @ major, east, north)))
assert checks['PA_roll_target_error_deg'] < 1e-12
checks['PA_roll_LOS_invariant_error'] = close(a @ los, los, 1e-14)
phi = np.linspace(0, 2 * np.pi, 360, endpoint=False)
points = 3 * np.cos(phi)[:, None] * major + np.sin(phi)[:, None] * minor
before = np.cov(np.column_stack([points @ east, points @ north]), rowvar=False, ddof=0)
rotated = points @ a.T
after = np.cov(np.column_stack([rotated @ east, rotated @ north]), rowvar=False, ddof=0)
eigen_before = np.linalg.eigvalsh(before)
eigen_after, axes_after = np.linalg.eigh(after)
checks['PA_roll_shape_eigenvalues_error'] = close(eigen_after, eigen_before, 1e-13)
shape_pa = np.rad2deg(np.arctan2(axes_after[0, -1], axes_after[1, -1])) % 180
checks['PA_roll_covariance_target_error_deg'] = abs(m.pa_difference_deg(target_pa, shape_pa))
assert checks['PA_roll_covariance_target_error_deg'] < 1e-12

# Nonzero cube transform and input centering/bulk offsets cannot be omitted.
old_angles, cube_angles = [103.42, 35.84, -1.8], [21, -12, 37]
correction = m.corrected_internal_angles(old_angles, cube_angles, a, p)
checks['outer_cube_conjugation_matrix_error'] = correction['native_matrix_error']
assert correction['native_matrix_error'] < 4e-7
cm, bm = np.array([0.012, -0.017, 0.009]), np.array([0.003, -0.007, 0.004])
cd, bd = np.array([-0.004, 0.006, 0.008]), np.array([0.012, -0.014, 0.002])
target_r, target_v = p @ r0, p @ solved['velocity']
placement = m.placement_for_target(target_r, target_v, correction['angles_deg'], cube_angles, cm, bm, cd, bd)
gp = np.r_[placement, [90, 0, 0], correction['angles_deg']]
cp = np.r_[[500, 500, 500, -20, 11, 13], [41, 72, -15], cube_angles]
actual_r, actual_v = m.initial_state_from_merge(gp, cp, cm, bm, cd, bd)
checks['placement_with_offsets_r_error'] = close(actual_r, target_r, 1e-12)
checks['placement_with_offsets_v_error'] = close(actual_v, target_v, 1e-12)

inverse_errors = []
for i in range(40):
    angles = rng.uniform(-160, 160, 3)
    if i == 0:
        angles = np.array([130.2, 89.98, -157.1])
    target = m.legacy_matrix(angles, float32_trig=False)
    result = m.legacy_angles_from_matrix(target, reference_angles=angles)
    inverse_errors.append(result['max_abs_matrix_error'])
checks['inverse_angles_max_matrix_error'] = max(inverse_errors)
assert max(inverse_errors) < 4e-7
for matrix in [np.diag([-1, 1, 1]), np.diag([1, 1, 2])]:
    try:
        m.legacy_angles_from_matrix(matrix)
    except ValueError:
        pass
    else:
        raise AssertionError('reflection or non-rotation was accepted as a physical rotation')
checks['reflection_and_nonrotation_rejected'] = True

# Compile the actual current implementation; compare all roles, fields and
# both ASCII rounding steps, exactly after the C output has been re-read.
source = args.reference_source
reference = WORK / 'ic_merge_reference.c'
reference.write_bytes(source.read_bytes())
compiler = shutil.which('gcc')
if compiler is None:
    raise RuntimeError('GCC is required for the actual C-merger validation')
compile_cmd = [compiler] + (['-D__USE_MINGW_ANSI_STDIO=1'] if platform.system() == 'Windows' else [])
compile_cmd += ['-O2', '-std=c99', '-ffp-contract=off', '-fno-fast-math',
               '-Wall', '-Wextra', reference.name, '-lm', '-o', 'ic_merge_reference.exe']
subprocess.run(compile_cmd, cwd=WORK, check=True, capture_output=True)
roles = []
for role in range(3):
    count = 64
    values = rng.normal(size=(count, 8)).astype(np.float32)
    values[:, :3] *= 7
    values[:, 3:6] *= 50
    values[:, 6] = rng.uniform(1, 5e4, count)
    values[:, 7] = rng.uniform(0, 3e4, count)
    types = np.arange(count) % 4
    name = ['cube.dat', 'mw.dat', 'dwarf.dat'][role]
    with (WORK / name).open('w', encoding='ascii') as stream:
        stream.write(' '.join(str(int(np.sum(types == i))) for i in range(4)) + '\n')
        for particle_id, (row, particle_type) in enumerate(zip(values, types)):
            stream.write(str(particle_id) + ' ' + ' '.join(format(float(x), '.9g') for x in row) + ' ' + str(particle_type) + '\n')
    roles.append((values, types))

trial_records = []
for trial in range(3):
    if trial == 0:
        galaxy = np.array([-192.3313, 164.1997, -183.4200, 34.7655, -86.0286, .1016,
                           90, 0, 0, 103.42, 35.84, -1.8], dtype=np.float32)
        cube = np.array([500, 500, 500, 0, 0, 0, 90, 0, 0, 0, 0, 0], dtype=np.float32)
    elif trial == 1:
        galaxy, cube = gp.astype(np.float32), cp.astype(np.float32)
    else:
        galaxy = rng.uniform(-179, 179, 12).astype(np.float32)
        cube = rng.uniform(-179, 179, 12).astype(np.float32)
    with (WORK / 'galaxy.par').open('w', encoding='ascii') as stream:
        stream.write('mw.dat dwarf.dat\n' + ' '.join(format(float(x), '.9g') for x in galaxy) + '\n')
    with (WORK / 'cube.par').open('w', encoding='ascii') as stream:
        stream.write('cube.dat mergegalaxy.dat\n' + ' '.join(format(float(x), '.9g') for x in cube) + '\n')
    filename = 'merged_%d.dat' % trial
    process = subprocess.run([str(WORK / 'ic_merge_reference.exe'), 'galaxy.par', 'cube.par', filename],
                             cwd=WORK, check=True, capture_output=True)
    table = np.loadtxt(WORK / filename, skiprows=1)
    predicted = [m.emulate_one_pass_particle(values, role, galaxy, cube) for role, (values, _) in enumerate(roles)]
    for row in table:
        particle_id = int(row[0])
        role, local_id = divmod(particle_id, 64)
        returned = row[1:9].astype(np.float32)
        np.testing.assert_array_equal(returned.view(np.uint32), predicted[role][local_id].view(np.uint32))
    trial_records.append(dict(trial=trial, particles=len(table), float_fields_checked=len(table) * 8,
                              bitwise_equal=True, galaxy_parameters=galaxy.tolist(), cube_parameters=cube.tolist()))
checks['actual_C_merger_bitwise_comparison'] = trial_records

runtime_note = ('C probe executed on Windows with MinGW ANSI stdio selected for nearest-even decimal rounding; not a Linux deployment test.'
                if platform.system() == 'Windows' else
                'C probe executed on this platform using its default C stdio; only this reported platform is certified.')
result = dict(state='PASS', platform=platform.platform(), random_seed=70047005,
              reference_source=str(source), reference_source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
              compiler=compiler, compile_command=compile_cmd, checks=checks,
              limits=[runtime_note,
                      'No real model particles, snapshots, production ICs or evolved candidates were modified or validated.',
                      'Empirical calibration and rigid rotations are candidate construction, not a hydrodynamic re-evolution.'])
with (HERE / args.qa_name).open('x', encoding='utf-8') as stream:
    json.dump(result, stream, indent=2)
    stream.write('\n')
print(json.dumps(result, indent=2))
