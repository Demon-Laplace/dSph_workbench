# Server data workspace

All server mutations must remain inside this data directory. Existing solver
installations outside data may be read and executed, but never edited.

Keep server production models in model/ and isolation models separately in
isoFornax/model/. Do not consolidate server isolation into production folders.
New automated isolation runs use isoFornax/model/runs/run/. Only the LOCAL
workspace consolidates ICs and isolation under each production model directory.
Never move active Slurm work or overwrite snapshots, particle ICs, model images,
videos, or frozen run inputs. Existing frozen validation cases keep their paths.
Equal model numbers alone do not prove matching IC provenance.

GitHub main is the source of truth for code. Root workflow*.py files correspond
to workflow/ in the Git checkout; dSph_workbench/ corresponds to repo/. After a
validated code change, synchronize that exact change to GitHub main. Use local
authenticated synchronization when server credentials are unavailable; never
copy credentials or force-push. Stop on divergent history and report it.
Machine-specific workflow.json, particle data, logs, and private configuration
must not be committed. IC_MODELS.csv and MODEL_SOURCES.csv are small summaries.

Respect the more specific dSph_workbench/AGENTS.md. Immutable workspace/releases/
are never edited in place; publishing analysis requires a new validated release.
Use float32 for inspection/analysis where verified, and preserve simulation
precision. Run focused tests and distinguish short integration checks from
scientific equilibrium or production-scale validation.

Existing user .sh files are protected: do not rewrite, rename, move, replace,
centralize, or delete them unless the user explicitly requests that operation.
Preserve their invocation paths when organizing Python dependencies. Obsolete
standalone temporary-check shell scripts may be archived only within an
explicitly authorized cleanup; model simulation launchers are not such scripts.
The legacy isoFornax/python/ files are byte-identical shared copies tracked in
workflow/isolation_python/. Model .py compatibility links preserve existing .sh
commands. Preserve isoFornax/python/ at this depth: the density script derives
its workspace paths from __file__. Do not mix these with a different analysis
version or edit scientific behavior as part of file organization.

The user explicitly authorized adapting the unnumbered isolation template on
2026-10-13. Its five versioned sub*.sh launchers now call shared Python and
prepare-legacy. This exception does not authorize changing existing numbered
models' shells or frozen run.sh files. Future template changes still require
explicit user authorization. Template Python duplicates are archived; do not
reintroduce per-model copies. TREECOOL is materialized from the configured
template, and new manual inputs must be named after their model directory.

Local model results use a single <model>/isolation/ directory. Do not recreate
isolation/legacy/<model>/ or isolation/<run>/client/ locally. Local transport
receipts belong in workspace/receipts/isolation/<model>/<run>/. Server isolation
paths and frozen cases remain separate and are not flattened by this convention.
