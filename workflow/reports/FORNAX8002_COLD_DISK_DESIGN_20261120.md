# Fornax8002: finite-dispersion, rotation-dominated stellar disk candidate

User request, 2026-11-20: create **8002**, retaining 8001 and every old model.
This is an exploratory initial-condition and isolation test, not an observational
fit or an already proven stable cold disk.

## Physical recipe

| Quantity | 8001 | 8002 candidate | Reason |
|---|---:|---:|---|
| Stellar intrinsic density c/a | 0.55 | 0.12 | Low vertical dispersion requires a thinner equilibrium density |
| Stellar seed kappa | 0.70 | 1.0 | Stronger mean streaming, not a prescribed prograde percentage |
| Stellar seed beta_z | 0.0 | 0.65 | Retain radial pressure relative to the thin disk's smaller vertical pressure |
| Final spherical beta-second constraint | -0.1 | removed | This otherwise forces substantial radial second moments, contrary to a cold disk |
| Stellar mass / Sersic n / radial scale | 6e7 Msun / 0.6 / 0.54 kpc | unchanged | Keep mass and face-on radial-profile parameters |
| Gas density and seed recipe | parent | unchanged | Do not retune feedback or gas structure at the same time |
| Stellar/gas particle numbers | 600000 / 14500 | unchanged | Keep the first test at the parent mass resolution |

The stellar potential uses 40 radial and 35 vertical grid points, rmin=0.01 kpc,
zmin=0.005 kpc. Its density target uses a nonuniform 40-point radial grid. These
are explicitly recorded resolution changes to represent a compact thin disk;
force softening in the evolution is not silently changed. The gas and stellar
velocities are both generated in the new combined potential; unchanged gas
parameters do not mean identical gas phase-space particles.

No velocity is added to existing particles, no new dark halo is inserted, and
no old particle/snapshot data are needed. The Schwarzschild solver still fits
density and total mass with nonnegative orbit weights; a cold seed is not itself
an enforced final rotation curve. The optimized N-body IC must be measured.

## Why the choices are coupled

Vertical Jeans balance links thickness and vertical pressure:
`d(rho*sigma_z^2)/dz = -rho*dPhi/dz`; a sheet estimate gives
`sigma_z^2 ~ pi*G*Sigma*h`, with profile-dependent coefficients.
Thus forcing a thick c/a=0.55 population to have tiny sigma_z would not be an
equilibrium construction.

For a nearly planar system, the old target
`beta_second=1-(<v_theta^2>+<v_phi^2>)/(2*<v_r^2>)=-0.1`
requires tangential second moments approximately 2.2 times the radial moment.
Keeping it while raising streaming alone cannot produce arbitrarily small
radial random motion. It is deliberately removed, not renamed as beta_z.

The sampling beta_z is `1-sigma_z^2/sigma_R^2`; at 0.65 the seed has
sigma_R/sigma_z about 1.69 before optimization. This is a stability compromise,
not a claim that all velocity dispersions should be minimized.

Toomre's local stellar thin-disk indicator is
`Q_star=sigma_R*kappa_epicycle/(3.36*G*Sigma_star)`.
Here kappa_epicycle is an orbital frequency, **not** the seed rotation parameter.
Surface density and the total gravitational potential set a lower useful radial
dispersion. Q>1 alone does not guarantee gas/star, bending or global bar stability.
This baryonic disk is not given an additional stabilizing dark halo.

## Preliminary seed study (not a complete IC)

60000 points per recipe, same declared random seed 42, no N-body IC export:

- q=0.12, beta_z=0, kappa=0.98 gives main-region V/sigma_R about 4.20 and 99.6%
  prograde, but Q is only about 0.71--0.81 at R=0.25--0.85 kpc.
- q=0.12, beta_z=0.65, kappa=0.98 gives main V/sigma_R about 2.17,
  V/sigma_1d about 2.37, 96.4% prograde, with Q about 1.20--1.81 over
  R=0.25--1.2 kpc.

8002 starts from the latter compromise with kappa increased to 1.0.
These seed results are not substituted for the optimized IC or an evolution.

## Declared measurements before constructing the full IC

Use old physical stars with spatial membership, no IDs. Report the full aperture
and radial profiles, with a main region `0.25 <= R < 1.5 kpc, r<5 kpc` (roughly
84% of the seed stellar mass). The central rising rotation curve is not hidden:
its separate bins are also reported.

Exploratory indicators, not universal definitions of a cold disk:

- main V/sigma_R >= 2 and V/sigma_1d,cyl >= 2;
- main prograde fraction >=95%;
- median initial orbital circularity Lz/Lcirc(E) >=0.8;
- initial single-component thin-disk Q target floor about 1.2 over R=0.25--1.2 kpc;
- inspect measured thickness, all density residuals and their radial dependence.

`disk_diagnostics.py` computes Q and circularity in the **initial** smooth
potential, only for a completed new IC. Evolving snapshots must not be assigned
an instantaneous Q or binding energy from an unchanged initial potential.
The ordinary rotation diagnostics separately report V/sigma_R, V/sigma_1d,
vertical thickness, radial streaming and angular momentum.

If the optimized IC is a useful cold-disk candidate, run isolation to code time
0.5 with the same feedback=16.6 and PAF=4.0 settings as 8001. A lightweight
startup monitor and the post-run scientific measurements run inside that same
allocation, without another monitoring job or editing 8001's frozen driver.
An initially cold equilibrium orbit superposition can still develop bars,
spirals, clumps or vertical heating. That is a scientific outcome, not something
to mask by relabelling the model or changing the data.
