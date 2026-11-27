# dSph workbench agent entry point

## IC GitHub synchronization (owner correction, 2026-11-27)

Ordinary model builds retain automatic GitHub catalogue, recipe, model-note and
SIMULATION_HISTORY.md synchronization. Default local-only publication is limited
to inputs/outputs under `analysis/X_series_IC_experiments/` or tasks explicitly
marked not to publish (`--skip-catalog`, alias `--no-catalog`). Do not extend the
experiment exception to ordinary Fornax models. `--register-catalog` explicitly
requests publication. Local result publication is not GitHub synchronization.
Preserve upstream edits; never upload particle ICs, snapshots, machine settings,
credentials or the excluded X-series experiments. This authorizes metadata
synchronization, not simulation submission or evolution.

## Supabase handoff log (owner update, 2026-11-21)

Use the configured local `x10-worklog recent --project Fornax` before substantive
work. Append concise decisions/corrections (`decision`), completed milestones
(`completed`), or submitted-run/blocker/unfinished handoffs (`handoff`), normally
one milestone per task. Existing blocked events remain historical; do not rewrite
them. Exclude routine reads, polling, each tool error, full chats and raw data.
Identify production versus isolation, run and JobID; recorded status is not live.
Use disk-relative stable artifact paths. Science details remain in model notes
and SIMULATION_HISTORY.md, not duplicated into the log. A next_step is never
authorization; check current job state before any later action. Credentials stay
with the approved local/portable CLI and must not be copied to servers or Git.
The portable disk's WorkLog/README.md defines the event policy; servers without
that CLI leave a concise handoff for the authorized local client to record.

## Local IC/result layout (owner clarification, 2026-11-20)

Put IC files, initial figures and isolation results together directly in
`<model>/isolation/`. Do not create `<model>/ic/` or build subdirectories there.
Keep existing filenames; use a build-name filename prefix only when versions
would otherwise collide. Fixed-name tool views belong under
`workspace/ic-builds/<model>/<build>/` and contain links, not duplicate particles.
This supersedes the older local IC layout below. Server runs are not relocated.
See [the local layout](workflow/LOCAL_LAYOUT.md).

## Isolation run versions (owner update, 2026-11-20)

Use the automatic structured run bindings described in
[workflow/RUN_VERSIONS.md](workflow/RUN_VERSIONS.md). Compare declared IC build
versions and the small effective profile/parameter configuration, never particle
IDs, particle contents, hashes or post-copy rereads. Equal versions reuse the
original job; changed versions allocate a new run and preserve old files.
Only `--submit` authorizes job submission. Historical unversioned runs retain
their identity, not invented version labels. JSON is authoritative and TXT is a
generated view. Do not edit frozen runs or migrate historical data for this.
This supersedes conflicting older isolation transfer/preflight/hash requirements.

Read [workflow/AGENT_HANDOFF.md](workflow/AGENT_HANDOFF.md) before operating the
local/IC/server workflow. It records the directory conventions, server write
boundary, startup-retry cleanup policy, density analysis and publication steps.
Use workflow/WORKFLOW.md for further command details. Inspect installed copies
and Git status before changing code. Do not edit frozen simulation runs.

For each IC, orbit, internal-pose, physics, or consequential numerical-run change,
append a concise entry to SIMULATION_HISTORY.md: reason, actual change, explicit
IC/orbit/predecessor-run provenance, and validation status with an evidence pointer.
Keep candidates distinct from applied changes and startup checks distinct from
complete scientific validation. Append corrections; preserve earlier history.

legacy_code_archive/ is reference-only. Never install or synchronize it into
active source/model directories or releases. IC/server source checkouts exclude
it with sparse checkout. Production model compatibility links use the existing
workspace/code/production entry, not independent module copies. Preserve user
shells, frozen runs and immutable releases when removing obsolete code.

Local model results use the shallow layout requested on 2026-11-15: production
inspection images/videos/result tables directly in the model directory;
isolation results directly in its separate isolation/ directory; supporting
production analysis files in a flat analysis/ directory with descriptive
filenames. Preserve frozen IC/run data and legacy compatibility paths. This
local layout does not relocate server simulations. See workflow/WORKFLOW.md.

Submit checks for multiple simulations as one Slurm job on one node, even when
their isolation run names differ. Use repeated `--run MODEL=RUN` arguments to
identify those runs. Process models sequentially with the configured worker count
(normally 16), including video-only checks. Do not submit one job per model or
run multiple worker groups concurrently. Preserve already submitted tasks.

Production orbital simulations use astro_fat, one separate job per model,
through the existing subastrofat.sh with its current resource settings. For
Fornax7004/7005, use snap2csv exports of their completed isolation snapshot 50
and retain timescale_fd=16.6. Validate exports and merged ICs before submission;
preserve original snapshots and existing user shell scripts.
Fornax7004/7005 orbital and orientation transforms must match Fornax5028.
Use automatic startup monitoring: after 140 running seconds and nonzero
integration progress, exit the monitor and leave the simulation running.
Avoid repeated manual checks or confirmation requests.
After successful startup detection, return simulation details (IC/model and
Milky Way/background inputs, particle counts/mass resolution, orbit/orientation,
physics parameters and job resources) and end monitoring.
Integrate monitoring into the existing batch allocation as a sleeping child,
without an extra Slurm job or reserved CPU. Use a generated wrapper, preserving
user launchers. Queue time uses no monitor node resources; local disconnection
must not stop server detection or simulation.
 
Each actual model change must also append a short UTF-8 MODEL_CHANGES.txt in
that model's server directory and synchronize the complete note to
model_changes/<model>.txt (isolation uses model_changes/isolation/<model>.txt).
Append its record ID and mirror link to SIMULATION_HISTORY.md in the same task.
This includes analysis/plot changes and run lifecycle/cleanup; group related
changes from one task, and omit read-only queries. Preserve history, provenance
and validation limits. Do not claim synchronization until both copies are read
back consistently. See workflow/AGENTS.server.md and workflow/AGENTS.Fornax.md.
