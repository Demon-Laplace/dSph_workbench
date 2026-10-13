# Server data workspace

All server mutations must remain inside this data directory. Existing solver
installations outside data may be read and executed, but never edited.

Use model/ic/ for new IC builds and model/isolation/run/ for isolation runs.
Historical directories reside in model/isolation/legacy/model/ with old-path
aliases. isolation/source.json records the production IC source; equal numbers
alone do not prove a match. Never move active Slurm work or overwrite snapshots,
particle ICs, model images, videos, or frozen run inputs.

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
