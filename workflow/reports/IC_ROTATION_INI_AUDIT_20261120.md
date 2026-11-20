# Historical IC rotation-argument inspection

2026-11-20. Small AGAMA **parameter INIs only**; no old particle INIs,
orbit arrays or simulation snapshots were opened for this inventory.
Historical files were not modified.

The current Schwarzschild wrapper forwards `icbeta` as `Density.sample(beta=...)`
and `ickappa` as `kappa=...`. The AGAMA sampler selects Eddington first if beta is
missing, so nonzero kappa alone does not enable rotating axisymmetric Jeans
sampling. A final `beta` constraint is not an `icbeta` argument.

## Findings before adding the new Fornax8001 test

| Source group | Stellar recipes with nonzero ickappa but missing icbeta | Other stellar recipes |
|---|---:|---|
| Current canonical IniFile library | 96 of 96 | none |
| Local frozen isolation parameters in the merged model workspace | 90 of 90 | none |
| Earlier local frozen parameter copies | 73 | 129 without ickappa |
| Local frozen IC builds | 8 | none |
| Legacy generated builds | 2 | none |
| Server frozen workflow runs | 13 | none |

The server canonical library independently has the same 96 affected stellar
recipes. These groups overlap in model IDs and include multiple historical
versions; the counts must not be summed as a number of distinct simulations.
No particle provenance was inferred from matching model numbers.

Canonical families: 2xxx=36, 3xxx=10, 4xxx=12, 5xxx=29, 6xxx=3, 7xxx=6.
All 96 canonical gas components supply both icbeta and ickappa, so this specific
missing-beta branch-selection issue does not apply to them. This is not a claim
that their hydrodynamic initial conditions are fully equilibrated.

Some earlier stellar recipes explicitly requested no rotation by omitting
ickappa; they should not be relabelled as suffering from an ignored requested
rotation. Three snapshot-export filenames were excluded from the parameter
inventory. ZENO/mkic5 Milky Way generators use a different path and are
not covered by this AGAMA keyword defect.

These are parameter-and-call-path findings, not a remeasurement of historical
N-body spins. Recipes generated with a different backend need separate
interpretation. Old INIs and frozen results stay unchanged.

## New controlled test: Fornax8001

Parent physical recipe: IC_Fornax5028.ini. The only parameter addition is
Component_stellar.icbeta=0.0; density, masses, particle/orbit counts, stellar
ickappa=0.70 and final spherical beta-second target=-0.1 remain unchanged.
Gas parameters are not changed. Seed 42 is newly declared, not claimed to match
an unknown historical random realization.

The first experiment tests the minimal branch-selection correction. It does
not silently introduce a new rotation-curve constraint or reinterpret the final
beta as a streaming-subtracted anisotropy. Final particle rotation must be
reported as a scientific measurement; a rotating seed is not a completed
evolution result.

The normal build/catalogue path uses the declared build ID rather than particle
hashes or certification rereads. Legacy digest fields, if present, remain
historical metadata. No old model is regenerated or globally edited.
