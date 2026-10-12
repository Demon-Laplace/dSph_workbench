# Production analysis agent instructions

This directory is the editable production analysis source. Preserve the physical
definitions, units, output schema, and working model paths unless the task
explicitly changes them. Inspect nearby code before editing. Changes on the
server, including temporary files, must stay inside the simulation data root.
Never alter model snapshots, particle initial conditions, or model-owned images
and videos as part of a code change.

## Required GitHub synchronization after a code change

The workspace owner explicitly authorizes automatic Git commits and pushes of
validated code changes to the configured dedicated synchronization branch.
After changing source code, documentation, or tests here:

1. Run the smallest relevant verification; for scientific changes also compare
   the affected outputs on a representative case. Record what was not tested.
2. From the data/workspace root, run:

   ```bash
   python3 workflow.py sync-github --message "Describe the concrete code change and validation"
   ```

   Use the configured Python module on the server. This command snapshots the
   allowed source files into `production_analysis/` of the configured Git
   checkout, commits them, and pushes `workflow/server-sync` by default.
3. Report the resulting commit and push result. Authentication failures,
   conflicting remote edits, or unavailable networking mean synchronization is
   incomplete. Preserve the local commit and report the error; never claim a
   successful push without a successful result. Do not ask for permission again
   for the already authorized scoped push.

The sync helper includes top-level Python modules, this file,
`SERVER_WORKFLOW.md`, and Python tests under `tests/`. If introducing another
necessary source directory or file type, explicitly extend and test the
allowlist before relying on automatic synchronization. Never stage the whole
simulation workspace. Never push raw data, particle INIs, notebooks with embedded
outputs, credentials, machine configuration, Slurm logs, images, or videos.
Do not force-push, reset away user changes, or overwrite the GitHub main branch.
Resolve divergent changes explicitly.

## Code publication and reproducibility

Editing source and publishing a new analysis release are separate actions.
`release` packages a clean committed revision; `deploy` uploads and verifies
checksums before activating it. Already prepared analyses retain their pinned
code. Do not edit immutable directories under `workspace/releases/`; use the
editable source directory configured by `github.source` instead. New analyses
use `analysis_code` from `workflow.json`, so a source edit requires a new release
to affect subsequent runs after release-based deployment is enabled.
