# Fornax paper-figure renderer

This directory contains a self-contained renderer for turning one GIZMO
snapshot into a black-background mock stellar observation suitable for an
A&A Letter figure panel.

## What the renderer does

- reads the input HDF5 snapshot without modifying it;
- follows the server analysis logic for shrinking-sphere Milky Way centring,
  k-nearest-neighbour dwarf identification, and the old-star definition
  (`PartType2/3`);
- transforms the dwarf to an ICRS tangent plane in degrees, with R.A.
  increasing to the left as in astronomical images;
- fixes angular and photometric scaling to the heliocentric distance defined
  by `d_today` in `variable.py` (currently 139.4 kpc), while recording
  simulated Galactocentric and transformed heliocentric distances separately;
- converts stellar mass to V-band light with `M/L_V = 2.6`;
- applies a two-scale PSF smoothing and a fixed `22--34 mag arcsec^-2` display
  range;
- reads candidate gas within 20 kpc of the dwarf stellar centre and uses the
  snapshot SPH smoothing lengths to construct the projected map;
- displays only NH-weighted cold H I (`T < 20000 K`), with no hot-gas layer
  and no additional density cut;
- draws fixed H I contours at `(0.5, 2, 5, 10) x 10^19 cm^-2`; the outer
  `5 x 10^18 cm^-2` level matches the sensitivity adopted for the paper
  comparison rather than scaling each panel to its own peak;
- annotates the selected dwarf stellar mass and the H I mass measured in a
  projected square 10 percent larger than the plotted field of view;
- restores the original normalized proper-motion direction arrow;
- displays the heliocentric distance as `D` and the lookback time in the
  lower-left corner of each paper panel;
- writes a compact PNG paper panel, a diagnostic PNG, and metadata.

No random sky background or artificial detector noise is added. This keeps
the morphology deterministic and prevents aesthetic noise from being mistaken
for a simulated tidal feature.

## Run

```text
python3 make_snapshot_mock.py <path-to-snapshot.hdf5> \
  --outdir output
```

The active Python environment must provide NumPy, SciPy, h5py, Matplotlib,
and Astropy. PNG and metadata outputs are written directly to the selected
output directory. By default, the renderer reads `elinfo` and defines
lookback-time zero at the row whose heliocentric distance is closest to the
adopted `d_today` distance unless overridden. Use `--elinfo` when the table
cannot be discovered next to the model or snapshot.

## Parameters to keep fixed for the final evolution montage

Use the same `--field-half-deg`, `--bright-limit`, `--faint-limit`, `--npix`,
`--smooth-pixels`, and `--gas-smooth-pixels` for every snapshot. A shared
surface-brightness/gas legend can then be added to the assembled multi-panel
figure. The current field is `+/-2.1 deg`, matching the existing diagnostic
pipeline.

## Six-panel evolution figure

`make_evolution_montage.py` reads both the model snapshots and `elinfo`, then
writes a two-row, six-panel PNG figure. By default it selects a 0.5 Gyr
lookback-time grid ending at the distance-matched present frame. If the model
starts slightly less than 2.5 Gyr before that frame, the first panel uses the
initial snapshot and the remaining five panels retain exact 0.5 Gyr spacing.
For Fornax2073 this gives snapshots 000, 043, 093, 143, 193, and 243, with
lookback times 2.43, 2.00, 1.50, 1.00, 0.50, and 0.00 Gyr. The short axis
labels are `RA` and `Dec`; tick values are in degrees.

```text
python3 make_evolution_montage.py /travail/xhuang/Fornax2073 \
  --elinfo /travail/xhuang/Fornax2073/elinfo_Fornax2073.csv \
  --outdir output
```
