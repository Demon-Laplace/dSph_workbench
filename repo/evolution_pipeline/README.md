# Fornax physical-evolution pipeline

This directory contains a configuration-driven, resumable analysis of a GIZMO
dwarf-galaxy run.  It produces one CSV time series and publication-ready PDF
and PNG figures.  Plotting and re-deriving time derivatives use only the CSV;
they do not reopen the snapshots.

## Consistency with the observational analysis

For a native sightline the pipeline calls `prepare_snapshot_context`,
`compute_snapshot_summary`, and the old-star projected-kinematics functions in
the parent `dSph_workbench` directory.  The main `re_major_kpc` and
`sigma_los_kms` columns therefore use the paper definitions:

- the standard MW centring and outer stellar-density-peak dwarf selection;
- the standard projected old-star centre, ellipticity, and position angle;
- the semi-major-axis old-star half-light radius for `R_e`;
- old stars inside the directly measured circular half-light aperture for
  `sigma_los`, after fitting and removing a planar LOS velocity gradient.

The pipeline also supports a configured Euler sightline, but `native` is the
default and the appropriate mode for reproducing the paper figures.

The H I columns distinguish two definitions:

- `hi_mass_particle_msun`: direct sum of `mass * NeutralHydrogenAbundance`
  below the configured temperature threshold and inside the gas aperture;
- `hi_mass_contour_msun`: adaptive projected H I map integrated above the
  configured fixed column-density contour.

`hi.mass_definition` selects which one is copied to `hi_mass_msun` and plotted.
The Fornax4001 configuration uses the contour definition employed by the paper.

## Physical definitions

- `gas_mass_msun` is all gas inside the configured 3D dwarf aperture.
- The dwarf COM velocity is the mass-weighted stellar velocity inside either a
  fixed 3D aperture or a configured multiple of the instantaneous `R_e`.
- Local CGM properties are measured in a configurable hot-gas shell outside
  the dwarf aperture.  `knn_shell` uses the nearest configured number of valid
  particles; `fixed_shell` uses all particles in the shell.
- Enclosed stellar and gas masses are 3D spherical quantities.  The fixed
  radii and instantaneous projected `R_e` are encoded explicitly in column
  names.
- `tidal_proxy_gyr2` is `G M_MW(<R_GC) / R_GC^3`.
- `tau_gas_smoothed_gyr` is the absolute value of smoothed gas mass divided by
  its time derivative.  Raw mass, raw derivative, smoothed mass, smoothed
  derivative, and raw and smoothed timescales are all retained.
- `stellar_dynamical_time_gyr = 0.9777922217 R_e[kpc] / sigma_los[km/s]`.

The CGM shell estimate is deliberately simple.  It is a local environmental
proxy, not an SPH density reconstruction at the dwarf centre.  Its exclusion
radius, search radius, neighbour count, temperature cut, and ionisation cut
are all recorded in the configuration and metadata.

## Commands

Run from the simulation directory because the inherited analysis code has
historically used model-relative paths:

```bash
cd /path/to/run
python3 /path/to/dSph_workbench/evolution_pipeline/fornax_evolution_pipeline.py \
  run --config /path/to/config.json
```

Subcommands:

- `run`: extract missing snapshots, derive time-series quantities, and plot;
- `extract`: extract missing snapshots and update the CSV;
- `derive`: recompute smoothing, derivatives, timescales, and event markers
  from the existing CSV only;
- `plot`: regenerate the figure from the existing CSV only.

The extractor checkpoints atomically and skips snapshots already present in
the CSV.  This makes it safe to rerun after a simulation creates new outputs.
If an analysis choice changes, use `--overwrite` to rebuild the selected range;
otherwise a configuration-hash mismatch is treated as an error.
Smoothing and event-marker choices have a separate derivation hash, so they can
be changed with `derive` without rereading any snapshots.

To use Slurm:

```bash
export CONFIG_PATH=/path/to/config.json
sbatch /path/to/dSph_workbench/evolution_pipeline/run_evolution_pipeline.slurm
```

## Changing simulations

Copy `config.example.json`, change `paths.run_dir`, and adjust only choices that
really differ (snapshot range, particle types, apertures, sightline, H I map,
CGM shell, smoothing, or comparison epoch).  No model name, snapshot number,
particle ID, or comparison time is embedded in the Python script.

## Event markers

The interaction start is the configured time or the first processed snapshot.
The comparison epoch requires a real crossing of the configured heliocentric
distance within tolerance.  Pericentre is marked only after a resolved local
minimum followed by the configured number of post-minimum points and minimum
rise.  A still-decreasing orbit therefore has no pericentre marker.
