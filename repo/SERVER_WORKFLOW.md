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

## Canonical heliocentric observations (2026-10-18)

The current observer is `astropy_native_direct_v1`. Native MW-centred axes are
explicit Astropy Galactocentric axes: R0=8.122 kpc, z_sun=0.0208 kpc, solar
velocity=(12.9,245.6,7.78) km/s, GC ICRS=(266.4051,-28.936175) degrees, roll=0.
The Sun lies on negative X. `canonical_observer.py` applies the complete six-D
transformation; `canonical_display.py` handles display geometry. Subtract the
same-snapshot MW position and bulk and apply units exactly once. Stars and both
gas temperature phases use the same conversion. Camera rotations preserve LOS.
This convention does not itself fit the model's absolute six-D state to Fornax.

Provide an explicit, unique uint64 `particle_id` array in `cohort_ids.npz` at
`<snapshotfolder>/analysis/astropy_native_direct_v1/`, or set
`DSPH_PREEXISTING_ID_PATH` to an audited frozen cohort. Its definition is
selected-then-frozen; it is not automatically all initial IC stars. No dynamic
density fallback is permitted. A new model needs its own documented cohort.
This frozen cohort identifies the old-star photometric reference, not an age
cut for velocity statistics. Current dwarf membership adds newly formed stars
using the existing spatial classification. Sigma (circular/elliptical/fixed
apertures and XYZ), Walker R5, the velocity profile/map, and R5 motion use all
current old+new members in their respective apertures. Light and half-light
radii retain the old-star reference. `sigma_re_nstar` counts all stars in the
circular half-light aperture; `sigma_re_noldstar` counts its old-star subset.
The observer config SHA includes this kinematic tracer definition, so outputs
from the old-only velocity analysis require a fresh analysis directory.

GetInfo and PlotFig write new results beside their snapshots under
`<snapshotfolder>/analysis/astropy_native_direct_v1/`. An explicit
`DSPH_CANONICAL_OUTPUT_ROOT` selects an isolated output directory. CSV identity
comments lock observer, method/config SHA and cohort-file SHA. Existing legacy
CSV, images and videos must remain unchanged; do not relabel or append to them.
Per-snapshot JSON records H/B, units, cohort/selected IDs, frame and helper SHA.
Its input path/size/mtime fingerprint is explicitly not a full snapshot hash.

The default sky display is north up/east left with raw wrapped DeltaRA and
DeltaDec each +/-2.1 degrees, 40x40 angular bins, occupied bins coloured and
velocity limits +/-12 km/s. Subtract the declared Walker v0 for display only.
The Heliocentric proper motion arrow starts at (0,0), with exact spherical
endpoint geometry. No gradient PA or angular-momentum-axis annotation is used.

Walker61 uses spherical angular distance in degrees and a forward projected
circle R=D0*tan(rho)<5 kpc, scanning 0,3,...,180 degrees, equal tracer weights,
Pi=1 and zero measurement errors. The display rectangle never defines this
fit. New CSV columns are walker_R5_k_kms_per_deg, walker_R5_v0_kms,
walker_R5_sigma_raw_mle_kms, walker_R5_sigma_after_mle_kms and walker_R5_nstar.
The historical half-light planar detrending statistic is a separate diagnostic.
Finite-sampling results require their observer, cohort, selection, seed and
draw IDs; old-observer probabilities do not transfer to this convention.

Workflow prepare pins the observer module and cohort; run/collect/fetch use the
versioned output directory and verify its identity. Previously prepared runs
retain their frozen code and legacy paths. Publish and activate a new immutable
release to select this method for subsequent workflow analyses. User shell
launchers, snapshots, ICs, simulation settings and manuscripts are preserved.
