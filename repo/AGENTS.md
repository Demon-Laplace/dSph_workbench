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
