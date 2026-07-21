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
- fixes angular and photometric scaling to the adopted heliocentric distance
  of 139.6 kpc (Li et al. 2021, distance modulus 20.72 mag), while recording
  simulated Galactocentric and transformed heliocentric distances separately;
- converts stellar mass to V-band light with `M/L_V = 2.6`;
- applies a two-scale PSF smoothing and a fixed `22--34 mag arcsec^-2` display
  range;
- reads candidate gas within 20 kpc of the dwarf stellar centre and uses the
  snapshot SPH smoothing lengths to construct the projected map;
- displays only NH-weighted cold H I (`T < 20000 K`), with no hot-gas layer
  and no additional density cut;
- annotates the selected dwarf stellar mass and the H I mass measured in a
  projected square 10 percent larger than the plotted field of view;
- restores the original normalized proper-motion direction arrow;
- displays the heliocentric distance as `D_sun` in each paper panel;
- writes a compact paper panel, a diagnostic panel, metadata, and vector PDF.

No random sky background or artificial detector noise is added. This keeps
the morphology deterministic and prevents aesthetic noise from being mistaken
for a simulated tidal feature.

## Run

```text
python3 make_snapshot_mock.py <path-to-snapshot.hdf5> \
  --distance-kpc 139.6 \
  --outdir output
```

The active Python environment must provide NumPy, SciPy, h5py, Matplotlib,
and Astropy. PNG and metadata outputs are written directly to the selected
output directory, while the vector version is placed in its `pdf` subdirectory.

## Parameters to keep fixed for the final evolution montage

Use the same `--field-half-deg`, `--bright-limit`, `--faint-limit`, `--npix`,
`--smooth-pixels`, and `--gas-smooth-pixels` for every snapshot. A shared
surface-brightness/gas legend can then be added to the assembled multi-panel
figure. The current field is `+/-2.1 deg`, matching the existing diagnostic
pipeline.

## Six-panel evolution figure

`make_evolution_montage.py` reads both the model snapshots and `elinfo`, then
writes a two-row, six-panel PNG/PDF figure. By default it uses snapshots 000,
050, 100, 150, 200, and 243. With the run cadence of 0.01 Gyr, the final panel
is therefore 2.43 Gyr. Snapshot 243 is the present-day match because its
`elinfo` heliocentric distance, 139.615 kpc, is closest to the adopted
139.6 kpc. The short axis labels are `RA` and `Dec`; tick values
are in degrees.

```text
python3 make_evolution_montage.py /travail/xhuang/Fornax2073 \
  --elinfo /travail/xhuang/Fornax2073/elinfo_Fornax2073.csv \
  --distance-kpc 139.6 --outdir output
```
