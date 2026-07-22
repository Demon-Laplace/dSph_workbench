# Central server workflow

The shared analysis and plotting programs live in
`/travail/xhuang/dSph_workbench`. Simulation directories contain model inputs,
snapshots, metadata, and outputs, but the processing code is not copied into
each model.

Run commands from the model directory so that `variable.py` resolves `output/`,
`output_img/`, `elinfo_<model>.csv`, and the simulation parameter files against
the selected model.

For a complete PlotFig rerun, submit the central launcher from the model:

```text
sbatch ../dSph_workbench/subchecknature.sh
```

Regenerating `elinfo` is independent and, when needed, uses:

```text
sbatch ../dSph_workbench/subgetinfo.sh
```

Plot frames remain in `<model>/output_img/`. The assembled video is written to
the model directory as `<model>_vN.mp4`, so deleting `output_img/` for a clean
frame rerun does not remove the previous video. Increment `PLOT_VIDEO_VERSION`
for every major plotting change that should produce a separately versioned
video.
