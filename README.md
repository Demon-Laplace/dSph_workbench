# dSph workbench

`IniFile/` contains the existing model parameter definitions; `repo/` retains
its upstream analysis and IC utilities. `production_analysis/` captures the
validated server analysis baseline without overwriting upstream scientific
implementations. `workflow/` coordinates local preparation, server Slurm work,
results retrieval, and Git-based analysis releases.

Read [the workflow guide](workflow/WORKFLOW.md) for setup and commands.
Machine-specific configuration, particle ICs, snapshots, logs, images, and
videos belong in the simulation workspace and must not enter this repository.
The `workflow/server-sync` branch receives scoped automatic agent updates;
changes are not automatically merged into `main`.

Install the `workflow/workflow*.py` modules together in the simulation data
root, copy and adapt `workflow/workflow.example.json` as `workflow.json`, and
run `python3 workflow.py doctor`. Keep the existing production source in the
configured `github.source` directory; no data migration is required.

Run the focused checks without scientific dependencies:

```bash
cd workflow
python3 -m unittest discover -s tests -v
```

The tested workflow covers existing-model analysis through validated tables,
images, video, and result collection. AGAMA construction and a new full GIZMO
simulation still require the existing scientific environment and validation.
