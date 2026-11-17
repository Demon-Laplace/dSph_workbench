# Spatial structured parameters

The owner selected current physical/spatial membership on 2026-11-17.
Particle IDs, uniqueness, origin, parent-bit dictionaries and frozen ID rosters
are not membership gates. The stellar and gas loaders do not read ParticleIDs.
Repeated IDs preserve every physical particle row; they do not establish a
physical duplication or permit deduplication. Temporary row labels inside the
unchanged observer mathematics only align one in-memory catalogue.

Gas membership is strict three-dimensional r<5 kpc and loose r<15 kpc. Both
inventories and the 5<=r<15 annulus include all gas, including environment-origin
gas. These spatial samples are not certified bound samples. Geometry, phase
and each physical threshold are saved. r is 3D radius; R is circular projected
radius. The observer remains astropy_native_direct_v1 with its fixed constants.

Stars use physical stellar properties inside 3D15 kpc of the locally tracked
density centre. Pre-existing versus formed populations use the simulation's
stellar types2/3 versus4. Recorded stellar ages are handled independently, so
a positive formation time is not mistaken for proof of a different birthplace.
This population-based
old/new label does not claim certified birthplace. A local shrinking-sphere
centre is seeded from existing geometric trajectories or merger translations,
then refined from actual stellar positions. It does not use a frozen ID centre.
Stellar ages are read from the particle's age/formation-time properties,
independently of IDs and origin. A recorded formation time gives age as current
snapshot time minus formation time, including pre-existing stars when recorded.
Only a genuinely missing or undefined physical age field is marked unavailable.

The old30kpc gas radius is a preselection, not the final HI-fraction aperture.
Save fHI=MHI/(Mstar+MHI) separately inside3D5 and3D15kpc, matching numerator and
denominator and excluding any preselection-only outer gas. An r30gas inventory
must never be combined with an inner stellar denominator under that label.
Walker uses current spatial old stars with forward projected R<5 kpc, equal
stars, zero measurement error and 61 directions; the same fit is applied to
all aperture residuals. This is a newly named membership contract, not a
silent replacement of legacy elinfo or the old-ID comparison products.

The new primary products are raw/residual LOS moments, separate residual
mean/RMS, 3D covariance and sqrt(trace)/sqrt(trace/3), radial/tangential moments,
beta, age groups, shape/depth, mass profiles, gas/HI/cold/density proxies/SFR,
feedback field inventories, wind/environment pressures, instantaneous fluxes
and a spherical host GM/Rgc^3 proxy. No particle is removed just because its
ID is repeated. A missing optional physical field is null plus availability.
Unknown compiled star-formation density thresholds are not guessed from an
unused configuration parameter; fixed nH thresholds are named proxies.

Current production audit (tan-hi-20261026-v1, also matching editable source):

* snapshot_context.py sets dwarf_gas_radius=30, passed to
  particle_selection.find_dwarf_gas: every TP0 within3D30kpc of the current
  stellar centre, without gas-origin filtering. It has not implemented the
  new separate5/15kpc membership inventories.
  snapshot_metrics then applies the projected ellipse R_ell<8*Rhalf to
  stars and both gas phases for hi_fraction; the gas sample is the intersection
  of that ellipse and the30kpc preselection sphere. It is incorrect to describe
  the current fraction as the unfiltered total of all30kpc gas.
  hi_half_fraction uses a matched3Dold-star half-mass sphere for both masses.
* Dynamic stars are selected by density peak/position (about15.3278kpc for
  the present Fornax defaults), but old stars are replaced by a fixed ID cohort.
  canonical_observer._ids/frozen_cohort_mask reject duplicated IDs.
* The historical canonical_model_comparison constructor requires unique source
  and loaded IDs and authenticated gas parent encoding; its old+new stellar
  selection excludes environment-born stars. Its gas mass profiles nevertheless
  include all gas in spatial apertures, so these are not all origin-only data.
* The environment_timing tool contains all/dwarf_parent/other_origin inventories,
  while its old fluxes use own-origin gas. Preserve those distinctions when
  retaining historical references. Do not rename them to the new all-spatial flux.
* Existing sigma_fixed_500pc actually uses r_3d_cut=0.3; tsigma_xyz is the
  one-dimensional equivalent, not the three-dimensional total dispersion.

This batch rereads physical fields for all registered runs, including all runs
previously blocked solely by ID/parent failures. Existing completed results
provide geometric search seeds and separately labelled historical references,
not a permission to mix different membership estimators. Snapshot/config/file
hashes still identify physical inputs; file integrity is not particle-ID review.

One Slurm allocation processes models sequentially with16workers. Per-epoch
checkpoints stay beside each simulation in output/analysis/sigma_spatial_parameters.
Final run-prefixed CSV, JSONL, NPZ, context and metadata are published to the
atomic ready-deliveries.json. The local transfer helper only consumes delivery
manifests, validates file bytes/SHA and atomically writes files to each model
directory at F:/Fornax_data/Fornax处理/. It never queries Slurm. A terminal
transfer-final.json ends the helper after all declared deliveries are received.
No plots, PDF production, correlations or causal conclusions are performed here.
