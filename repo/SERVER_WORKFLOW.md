# Central server workflow

The shared Python analysis and plotting programs live in
`/travail/xhuang/dSph_workbench`. Slurm and simulation launch templates remain
in `/travail/xhuang/Fornax`, so copying `Fornax/*` into a new simulation creates
the required launchers without duplicating the Python processing code.

Run commands from the model directory so that `variable.py` resolves `output/`,
`output_img/`, `elinfo_<model>.csv`, and the simulation parameter files against
the selected model.

For a complete PlotFig rerun, submit the copied launcher from the model:

```text
sbatch subchecknature.sh
```

Regenerating `elinfo` is independent and, when needed, uses:

```text
sbatch subgetinfo.sh
```

Plot frames remain in `<model>/output_img/`. The assembled video is written to
the model directory as `<model>_vN.mp4`, so deleting `output_img/` for a clean
frame rerun does not remove the previous video. Increment `PLOT_VIDEO_VERSION`
for every major plotting change that should produce a separately versioned
video.

PlotFig uses resume behavior by default: complete PNG files are not overwritten
unless `--replace` is supplied. To force a complete redraw without `--replace`,
first move any existing video out of `output_img/` into the model directory and
then delete `output_img/` before submitting `subchecknature.sh`.

If Slurm reports `(launch failed requeued held)`, treat it as a transient server
launch failure. Wait a few seconds, cancel the held job, and submit the same
launcher once more before investigating the plotting program itself.
