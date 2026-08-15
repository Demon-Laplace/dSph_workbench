# Current-epoch observable comparison

`make_current_observables.py` builds the second paper figure from the single
snapshot whose heliocentric distance is closest to `variable.d_today`.

The three panels compare the simulation with the adopted observations for:

- V-band surface brightness, using corrected contiguous annulus boundaries
  and the observed `r_mid_kpc` values;
- line-of-sight velocity dispersion of old stars after subtracting a fitted
  planar velocity gradient, sampled in the Walker et al. radial bins;
- the recent star-formation history, using only stars already formed in the
  selected snapshot, fixed 0.5-Gyr bins and step curves as in `PlotFig.py`,
  and the same elliptical aperture of `r_ell <= 0.8 deg`. Stellar masses are
  converted directly to the de Boer et al. plotting unit of
  `1e-4 Msun yr^-1`, without an additional normalization. Integrating the
  plotted simulation curve therefore recovers the stellar mass in the selected
  snapshot and aperture.

Run from the model directory so that the established pipeline resolves all
relative paths consistently:

```text
python /path/to/make_current_observables.py .
```

The script writes one PNG and one JSON metadata file to `paper_figure2/`.
