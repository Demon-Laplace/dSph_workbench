# Production analysis agent instructions

This directory is the editable production analysis source. Preserve the physical
definitions, units, output schema, and working model paths unless the task
explicitly changes them. Inspect nearby code before editing. Changes on the
server, including temporary files, must stay inside the simulation data root.
Never alter model snapshots, particle initial conditions, or model-owned images
and videos as part of a code change.

## Required GitHub synchronization after a code change (main only)

The workspace owner explicitly authorizes automatic Git commits and pushes of
validated code changes to the configured main branch.
After changing source code, documentation, or tests here:

1. Run the smallest relevant verification; for scientific changes also compare
   the affected outputs on a representative case. Record what was not tested.
2. From the data/workspace root, run:

   ```bash
   python3 workflow.py sync-github --message "Describe the concrete code change and validation"
   ```

   When working through a local SSH connection, prefer this command in the
   local workspace instead (it uses local GitHub authentication without
   copying credentials to the server):

   ```bash
   python3 workflow.py sync-server-github --message "Describe the concrete code change and validation"
   ```

   It exports only allowed server source files, verifies their transfer
   checksums, and performs the same Git commit and push locally. Use the
   configured Python and Git modules for a direct server-side push. This command snapshots the
   allowed source files into `repo/` of the configured Git
   checkout, commits them, and pushes `main` by default.
3. Report the resulting commit and push result. Authentication failures,
   conflicting remote edits, or unavailable networking mean synchronization is
   incomplete. Preserve the local commit and report the error; never claim a
   successful push without a successful result. Do not ask for permission again
   for the already authorized scoped push.

The sync helper includes top-level Python modules and shell/Slurm launchers, this file,
`SERVER_WORKFLOW.md`, and source files in the explicit package allowlist
(`tests`, `build_production_ic`, `check_density_profile`, `evolution_pipeline`,
`paper_figures`). If introducing another
necessary source directory or file type, explicitly extend and test the
allowlist before relying on automatic synchronization. Never stage the whole
simulation workspace. Never push raw data, particle INIs, notebooks with embedded
outputs, credentials, machine configuration, Slurm logs, images, or videos.
Do not force-push or reset away user changes. Push only fast-forward updates to main.
Resolve divergent changes explicitly.

## Code publication and reproducibility

Editing source and publishing a new analysis release are separate actions.
`release` packages a clean committed revision; `deploy` uploads and verifies
checksums before activating it. Already prepared analyses retain their pinned
code. Do not edit immutable directories under `workspace/releases/`; use the
editable source directory configured by `github.source` instead. New analyses
use `analysis_code` from `workflow.json`, so a source edit requires a new release
to affect subsequent runs after release-based deployment is enabled.

## Analysis precision

Use float32 for bulk particle analysis and diagnostic mass arrays by default.
Avoid promoting full arrays to float64 solely for output comparisons. Preserve
the simulation solver precision settings. Do not use float16 for unscaled
physical mass values; its finite range is insufficient. Verify explicit
reduced-precision choices with relevant output checks.


## Canonical observing contract adopted 2026-10-18

New PlotFig/GetInfo observations use astropy_native_direct_v1. Native MW-centred
axes map directly to explicit Astropy Galactocentric axes (identity): R0=8.122
kpc, z_sun=0.0208 kpc, solar velocity=(12.9,245.6,7.78) km/s, GC ICRS=
(266.4051,-28.936175) deg, roll=0. Sun is on the negative X side. Apply the
same-snapshot MW H, TP2/3/4 massmean bulk B within true 3D r<5 kpc, and code
units once, then one complete 6D transform to ICRS. Stars and both temperature
gas phases share this contract. Camera changes never redefine bound per-ID LOS.

Default arrows are Heliocentric proper motion, using the same R5 massmean HRF
velocity at the common image-centre E/N and distance. GRF motion is separate
and must say Motion relative to the Milky Way. Display axes are raw wrapped
DeltaRA and raw DeltaDec, both +/-2.1 degrees, north up/east left; arrows start
at (0,0) with correct spherical/raw-RA endpoint geometry. The 40x40 velocity
display has colour wherever stars exist and +/-12 km/s limits, subtracting a
declared display constant only. Display window never replaces the measurement
aperture. No gradient PA or angular-momentum-axis annotations.

Walker ideal-simulation measurement uses true spherical rho in degrees, current
old+new dwarf member stars, forward Rproj=D0*tan(rho)<5 kpc, equal
stars with Pi=1 and zero measurement error, all 61 directions 0,3,...,180 deg.
Choose max abs(k) in km/s/deg. Half-light planar detrending remains a separately
named legacy diagnostic, never a Walker61/R5 replacement. Raw and residual MLE
dispersion are separately named. Do not reuse legacy observer search/probability
results as canonical without remeasurement.

Canonical analysis requires cohort_ids.npz (particle_id unique uint64) at
output/analysis/astropy_native_direct_v1 or explicit DSPH_PREEXISTING_ID_PATH.
Record selected_then_frozen, not all initial IC stars. Missing or changed
cohorts fail closed. CSV observer/config/cohort SHA guards prevent mixed-frame
or mixed-member continuation. New outputs live beside the snapshot in
output/analysis/astropy_native_direct_v1; DSPH_CANONICAL_OUTPUT_ROOT explicitly
selects a fresh validation/output directory. Legacy CSV/images/videos, raw
snapshots/ICs, GIZMO, user shell launchers, and the paper are preserved.
The frozen cohort is the old-star photometric reference. Velocity dispersions,
Walker fits, R5 motion and velocity maps/profiles include newly formed current
dwarf members; do not restrict these statistics to that initial cohort.
