# Server workflow entry

The maintained entry is `workflow.py` in the simulation data root. See
`WORKFLOW.md` there for configuration, analysis submission, collection, GitHub
synchronization, and releases. Follow `AGENTS.md` in this source directory after
changing code. Historical per-model launchers remain for compatibility.

## H I statistics correction (2026-10-16)

NeutralHydrogenAbundance is n_HI/n_H. The solver uses HYDROGEN_MASSFRAC=0.76,
so each gas particle contributes 0.76 * mass * NeutralHydrogenAbundance to H I.
Both temperature phases contribute; the 20,000 K display split is not an H I cut.
GetInfo now writes hi_mass, hi_half_mass, hi_fraction, hi_half_fraction,
gas_mass and gas_half_mass. The user defines HI fraction as M_HI/M_star,
not M_HI/(M_HI+M_star). Numerator and denominator use identical apertures.
The legacy coldgas_mass/coldgas_half_mass columns alias the corrected HI masses.
Gas mass fields contain all gas, including helium and ionized material; the
baryonic theoretical dispersion uses total gas mass instead of summing HI and
legacy hotgas columns. hotgas_* retain their legacy ionized-hot-phase definition.
HI contours use both temperature phases and the hydrogen factor; contour HI
fractions use stellar mass in the same thresholded pixels.

Old elinfo files must be recomputed with GetInfo before using the new PlotFig;
old values cannot be corrected by changing labels. Existing frozen analyses and
immutable releases are preserved; publish/activate a new release for future jobs.
The build_production_ic report's legacy matched_gas_fraction key now reads the
explicit HI mass and returns M_HI/M_star; it is NaN for old CSVs without hi_mass.

Validation: all 16 analysis unit tests passed. Fornax7001 snapshots 100–102
completed GetInfo and PlotFig on astro_thin (job 305985). An independent
particle/aperture calculation for snapshot 100 (job 305987) matched the new
HI masses, total gas masses, ratios and CSV values; unaffected summaries
matched the previous immutable release exactly. HI mass changed from
82,045,848 to 62,355,117 solar masses; HI/stellar mass=1.01747367. The warm
phase contributes 268 solar masses in this aperture. This validates the
implementation on a representative case, not every historical simulation.
