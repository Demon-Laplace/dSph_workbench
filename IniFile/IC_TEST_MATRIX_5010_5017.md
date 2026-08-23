# Fornax5010--5017 isolation test matrix

All models are derived from `IC_Fornax5008.ini`. Unless stated otherwise they retain
`M_gas=1.6e8 Msun`, `M_star=6e7 Msun`, stellar `n=0.6`, `R_star=0.55 kpc`,
`q=0.55`, and gas `R_gas=0.90 kpc`, `H_gas=0.45 kpc`.

| Model | Change relative to 5008 | Purpose | Priority |
|---|---|---|---|
| 5010 | beta=-0.20, kappa=0.75 | Compact-gas-only test; restores 5007 stellar kinematics | A |
| 5011 | R_gas=1.00 kpc, H_gas=0.60 kpc | Kinematics-only test; restores 5007 gas structure | A |
| 5012 | beta=-0.20 | Isolate beta at fixed kappa=0.70 | A |
| 5013 | kappa=0.75 | Isolate kappa at fixed beta=-0.30 | A |
| 5014 | q=0.60 | Moderate rounder stellar model | A |
| 5015 | q=0.65 | Upper-q bracket; run after 5014 | B |
| 5016 | R_gas=0.95 kpc, H_gas=0.525 kpc | Intermediate gas compactness between 5007 and 5008 | B |
| 5017 | q=0.60, beta=-0.25 | Combined candidate if 5012 and 5014 remain viable | B |

## Suggested execution order

Run 5010, 5011, 5012, 5013, and 5014 first. These five models separate gas
compactness, anisotropy, rotation, and intrinsic thickness. Models 5015--5017 are
bracketing or combined tests and may be deferred until the first group is inspected.

## Acceptance criteria at 3--6 Gyr

- No isolation gas instability: at least 90% of the extant gas remains within 20 kpc.
- Mean line-of-sight dispersion over 0.25--1 kpc remains at least 9.5 km/s in a
  Fornax-compatible projected shape.
- The 1--1.8 kpc dispersion does not fall materially below the 5008 result.
- The central dispersion does not develop a large peak or a persistent central dip.
- The stellar half-light radius retains enough headroom for modest expansion in the
  full Milky Way run.

The production-environment reference remains 5008. The isolation matrix is intended
to identify why it improved over 5007, not to delay the first full-orbit 5008 test.
