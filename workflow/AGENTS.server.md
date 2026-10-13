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
