# Orbit-candidate numerical helpers

`candidate_math.py` constructs explicitly physical candidates while keeping the
observer unchanged. It never reads or edits a model, snapshot, IC or production
parameter file. All calls return numeric values for the caller to inspect.

## Public API

- `empirical_coefficients(r0,v0,rf,vf)` returns `C_E=|vf|²−|v0|²` and
  `k_L=|rf×vf|/|r0×v0|`. Use the same tracer/centre definition at both epochs.
- `solve_initial_velocity(r0,v0,target_r,target_v,C_E,k_L,radial_sign=-1)` solves
  the inward initial radial and tangential velocities, preserving the seed
  radial/tangential directions. Radius must be computed from the actual target
  Galactocentric position, not copied from an old heliocentric crossing.
- `orbit_plane_rotation(current_r,current_v,target_r,target_v)` maps the three
  directed orbital basis vectors by a proper SO(3). It does not repair unequal
  radial speeds, total speeds or radii, and is not a live-host invariance claim.
- `pa_roll_rotation(los,delta_pa_deg)` increases standard north→east PA by delta.
  With `E×N=n`, its right-handed Rodrigues angle is **minus delta**.
- `corrected_internal_angles(old_angles,cube_angles,A,P=None)` solves the
  internal merger angles from `G_new=Q⁻¹ A P Q G_old` (P defaults to identity).
  An internal roll A is excluded from the target COM rotation; only P rotates
  the orbit/COM. The result reports finite-float matrix errors.
- `placement_for_target(r0_target,v0_target,new_angles,cube_angles,cM,bM,cD,bD)`
  returns the first six merger values with all input-centre/bulk offsets:
  `r_g=Q⁻¹r0_target+cM−G_new cD`,
  `v_g=Q⁻¹v0_target+bM−G_new bD`. These inputs use native **code units**; convert
  back from kpc/km/s before calling. Common cube translation/boost cancels in
  MW-relative vectors.
- `emulate_one_pass_particle(values8,role,galaxy12,cube12)` replays C float32
  operation order, `.017453` radians per parameter degree, and the intermediate
  and final `%.7g` conversions. Roles are 0 background, 1 MW, 2 dwarf. This
  excludes ID/type ordering, which the reference C merger handles separately.

The linear merger matrix is `Rz(-az*c) Ry(-ay*c) Rx(-ax*c)`, `c=.017453`.
Inverse angles use that actual scale and sign, rather than exact π/180.
Conjugation uses `Q⁻¹`, since a float-rounded Q is not exactly orthogonal.

## Verification and bounds

Latest `MATH_QA_v4.json` records analytic recovery, proper-frame mapping, a non-axis-aligned
PA/elliptical covariance check, nonzero outer cube rotation, all four input
offsets, 40 angle inversions, and bitwise comparison against the actual current
`workflow/ic_merge.c`: 3 trials × 192 particles × 8 float fields = 4608 fields.
The local checks pass. The largest inverse-angle matrix error was 1.19e-7;
the tested outer-conjugation error was 3.02e-8. Reflections and non-rotation
matrices are rejected rather than silently projected into an apparently physical
candidate. Earlier successful QA records remain for provenance.

Local execution used Windows GCC with `__USE_MINGW_ANSI_STDIO=1`, explicitly
selecting nearest-even formatting. The default Windows CRT rounds some exact
decimal halfway values differently; the observed discrepancy is retained in
`WINDOWS_CRT_ROUNDING.json` and the first probe directory. No production code
was changed. Do not call this a Linux deployment test.

The parent can run the same inexpensive probe inside the already allocated
Linux Slurm job (no separate job needed):

```bash
python3 validate_candidate_math.py \
  --reference-source /travail/xhuang/workspace/code/github-sync/workflow/ic_merge.c \
  --work-name validation_linux --qa-name MATH_QA_LINUX.json
```

The chosen work and QA names must not exist. This verifies that platform's C
runtime and source, not a new model evolution. Empirical energy/angular-momentum
calibration, rigid rotations and translated initial parameters remain candidate
construction; no pericentre or hydrodynamic outcome is certified by these helpers.
