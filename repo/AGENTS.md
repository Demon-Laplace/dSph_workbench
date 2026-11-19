# Production analysis agent instructions

This directory is the editable production analysis source. Preserve the physical
definitions, units, output schema, and working model paths unless the task
explicitly changes them. Inspect nearby code before editing. Changes on the
server, including temporary files, must stay inside the simulation data root.
Never alter model snapshots, particle initial conditions, or model-owned images
and videos as part of a code change.

Before implementing or running an internal-orientation search, read the
[3D orientation search and validation principles](orbit_calibration_v1/README.md#内部姿态搜索与检查原则).
Use broad SO(3) coverage with local refinement, retain per-star 3D depth in
projection, compare morphology beyond PA, and recheck projected kinematics.
Distinguish snapshot previews, initial-condition candidates and re-evolved
validation; choose reference data explicitly for each task.

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
`paper_figures`, `orbit_calibration_v1`). Orbit-processing source, examples,
tests and documentation are versioned on main; the earlier local-only orbit
source restriction was revoked by the user on 2026-10-20. Snapshot-specific
results and large inputs remain outside Git. If introducing another
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
and must say Motion relative to the Milky Way. Display uses WCS-TAN xi/eta,
both +/-2.1 plane degrees, equal aspect and square sky panels, north up/east
left. Arrows start at (0,0) with endpoint (180/pi)*tan(0.65 deg)*unit_EN;
never divide their east component by cos(Dec) again. The 40x40 velocity
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

Current production membership uses particle properties and the current spatial
density-centre aperture, under member_selection_id=current_spatial_properties_seeded_v2.
Source particle IDs and initial origin never gate the primary sample. Repeated
IDs retain every physical row. Explicit analysis_row labels align working arrays;
optional source IDs remain separately named, in their original non-object dtype.
No cohort_ids.npz or ID-origin certification is required. TP2/3 remain the
separately named light tracer; all current member stars enter the kinematics.
Formation times or recorded ages are read for every stellar type independently
of the light-tracer label. Missing fields are NaN; recorded zero or negative
formation time is not replaced by an unknown-age assumption. Units come from
the actual run parameters, an explicit override or the declared canonical units.
Gas inventories include all TP0 rows within strict 3D r<5 and r<15 kpc, plus
5<=r<15 kpc. HI fractions use stellar and HI mass in the same sphere. A larger
gas radius is preselection only; the legacy projected and half-mass estimators
keep their existing names and apertures. Observer and measurement definitions
remain unchanged. New outputs add current_spatial_properties_seeded_v2 beneath the
canonical observer directory; CSV membership identity prevents appending to
old-member results. DSPH_CANONICAL_OUTPUT_ROOT remains an explicit output option.
Legacy results, raw snapshots/ICs, frozen releases, user launchers and the paper
are preserved. Historical ID-reference helpers are outside the primary path.

Production TAN media use display_contract_id=heliocentric_tan_xi_eta_2p1_v1
and live under the canonical analysis root's display/<display_contract_id>/.
Keep physical observer/CSV identity unchanged for display-only migrations.
Rebin the frozen catalogue in TAN; use integrated 2D TAN pixel solid angles for
surface brightness and map kpc gas/shape grids linearly by (180/pi)/D0. Preserve
legacy raw CSV centre fields, transforming their same sky point for display.
PNG sidecars and display_manifest.json bind media to display/source and physical
identities. An unlabelled or changed-contract complete PNG cannot be resumed.
Old raw media and pinned runs remain untouched. Workflow release publication,
editable-source installation and Linux client synchronization are separate;
updating one does not update the others. Existing user .sh files stay unchanged.

The seeded_v2 centre follows an existing geometric trajectory only as a local search seed and remeasures physical stars per snapshot using the shared spatial helpers. DSPH_SPATIAL_CONTEXT and DSPH_SPATIAL_HELPERS_DIRECTORY provide the job-package inputs. Source IDs never gate the current sample.
