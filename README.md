# dSph workbench

`IniFile/` contains small model parameter definitions. `repo/` is the unified
analysis and IC-tool source, reconciled from the prior feature branches and
server production fixes. `workflow/` coordinates local preparation, server
Slurm jobs, result retrieval, source synchronization, and Git-based releases.

Read [the workflow guide](workflow/WORKFLOW.md) for setup and commands.
Machine-specific configuration, particle ICs, snapshots, logs, images, and
videos belong in the simulation workspace and must not enter this repository.
The remote keeps only `main`, which receives validated, scoped automatic agent
updates. Divergent history stops synchronization; force pushes are not used.

Install the `workflow/workflow*.py` modules together in the simulation data
root, copy and adapt `workflow/workflow.example.json` as `workflow.json`, and
run `python3 workflow.py doctor`. Keep the existing production source in the
configured `github.source` directory; no data migration is required.

Run the focused checks without scientific dependencies:

```bash
cd workflow
python3 -m unittest discover -s tests -v
```

The tested workflow covers existing-model analysis and local AGAMA generation
through catalogue synchronization, verified transfer and a real short Slurm
isolation run. Full-scale, long-duration stability still needs scientific validation.

[IC model catalogue](IC_MODELS.csv) lists historical parameter summaries and
verified new builds. Local ICs live under each model’s `ic/` directory. Browse isolation data directly in
`<model>/isolation/`; transport receipts live separately in `workspace/receipts/isolation/`.
New particle exports use the model ID (e.g. `Fornax7003.ini`) with a local
`Fornax.ini` compatibility link. Isolation accepts multiple model IDs per command.
Server isolation stays separate under `isoFornax/<model>/runs/<run>/`;
every successful build automatically records and pushes its catalogue row.
See the workflow guide for layout, batch commands and synchronization retries.
