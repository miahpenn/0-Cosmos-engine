# Initial radial quadrature + CMC centre-row candidate — 13 October 2026

## Scope and provenance

This isolated candidate is stacked on the strict campaign/radiation integration source `9f0bd16d403ada22a4f21413048c93b37c9924a7` with vendor submodule `d6052d605673ce9d82cdc99b0df79855fa2ea215`. It changes only the numerical initial-data quadrature in `engine/v55_initial.py` and the first finite-difference row in `engine/cmc_gauge.py`, with focused regressions. It does not change the physical source equations, matter model, parameters, radiation cone, RHS or vendor source.

## Why both layers are included

Exact-source diagnostics separated two issues:

- With the default production initializer, the CMC alternatives regularize the even `alpha(r²)` fit and reduce the maximum `|K_rhs-target|` from about `1.58e-3` to `6.3e-6` at N=160, but the centre `Aa_rhs/r0²` remains about `1.67e-2` and the Hamiltonian witness remains about `4.0e-3`. The row-only change is therefore not a complete repair.
- With the diagnostic high-order radial integral, the same CMC row alternatives greatly suppress the centre `Aa_rhs/r0²` and improve the immediate M/connection witnesses. This indicates the initial constrained radial integration is a second, independent numerical layer.

## Candidate numerical changes

1. Replace the origin half-cell trapezoid with the analytic even-series weights `(401/720, -31/480, 11/1440)`, exact for a smooth even source through the (r^4) term.
2. Use integrated cubic-Lagrange interval weights, including an even-ghost first interval and a one-sided final interval, for cumulative radial integration.
3. Replace `alpha[0]=alpha[1]` in the CMC banded system with the discrete CMC PDE row obtained by even reflection at the origin. Interior rows and the outer `alpha[-1]=1` normalization remain unchanged.

## Admission boundary

This is a candidate, not a production promotion. Focused polynomial tests and the repository suite must pass, followed by a bounded full-machine step and review of centre-series/H/M/connection, radiation admissibility, Misner–Sharp conservation, matched spatial refinement, and the pinned archived-state radiation continuation. A short smoke cannot admit the full trajectory. No long campaign should be launched until the cross-lane gates are reconciled.


## Exact-head candidate validation

- Candidate unit + smoke workflow: [run 38328832767](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38328832767), candidate head `9a8432a76b4a169c5c68f59093d055d3bef8a924): 11 focused tests passed, 163 full-suite tests passed, and the N=40 / r_max=40 / t=0.03 production smoke completed with campaign admission true and strict JSON.
- Five-resolution real-kernel one-step sweep: [run 38328997773](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38328997773), final diagnostic commit `6e7479a66bfef51ae2aaa5dcb21ff9378c91614d`. All N=40/60/80/120/160 arms completed one real production step to `t=0.03*dr/dx` (dt=0.03*r_max/N). Production kernel/matter/vendor source pins are guarded against PR #16's integration head `9f0bd16d403ada22a4f21413048c93b37c9924a7`; vendor submodule `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Sweep artifact: [11781768503](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38328997773/artifacts/11781768503), ZIP SHA-256 `efc521559f8e2f01a0c9b932b3ac37940d09fc9c9ebfe31f142d610af8fdc570`.

### One-step refinement evidence for the combined candidate

| N | Initial full H max | Final full H max | Final full M max | Final full C_lambda max | Final central Aa_rhs/r0² | Final lapse even-series fit residual |
|---:|---:|---:|---:|---:|---:|---:|
| 40 | 3.865e-4 | 3.837e-4 | 2.068e-5 | 1.266e-6 | +1.043e-4 | 2.811e-4 |
| 60 | 1.003e-4 | 9.985e-5 | 3.367e-6 | 1.376e-7 | +2.256e-5 | 3.721e-5 |
| 80 | 3.677e-5 | 3.670e-5 | 4.191e-7 | 1.299e-8 | -1.466e-5 | 6.563e-6 |
| 120 | 1.057e-5 | 1.056e-5 | 3.340e-7 | 6.748e-9 | -4.452e-5 | 3.911e-7 |
| 160 | 6.881e-6 | 6.876e-6 | 2.965e-7 | 4.506e-9 | -5.568e-5 | 4.550e-8 |

The initial/post-step full-domain Hamiltonian maximum decreases monotonically over this refinement arm, and the immediate momentum and connection witnesses remain small. The centre `Aa_rhs/r0²` also stops growing like the default-initializer row-only control and appears to approach a finite negative value across the upper resolutions. The first version of this sweep's determinant field was misreported: the pinned vendor's `constraints["determinant"]` is already the residual `a*b² - 1`, but the diagnostic subtracted 1 again. Therefore its determinant column was invalid and has been removed from this summary; a corrected report is required before using determinant evidence from the sweep. The engine's separate production diagnostic and the full-machine smoke do report determinant residuals independently. These are encouraging one-step refinement indicators, not yet a proof of continuum convergence: each resolution has evolved for a different coordinate time because `dt=0.03*dr`, and a synchronized fixed-time comparison plus fixed-physical-radius profiles remains required.

The radiation face-jump metric remains about 0.50 in this diagnostic. It is a normalized difference between adjacent-face energy fluxes, not itself a conservation residual; inspect the signed outer/interior divergence and metric-source terms along with it. Same-step Misner–Sharp mass, energy flux and work terms were persisted for every resolution, but their history balance needs multi-step integration. The archived outer-boundary radiation continuation failure remains a separate blocker.

## Next gates before promotion

1. Run a short candidate trajectory with the same candidate source at matched physical horizon/timestep comparison for N=80/120/160, retaining per-step H/M/C_lambda/determinant, centre series, radiation cone and signed face divergence/source terms.
2. Complete matched fixed-physical-radius convergence checks; do not infer convergence from different-time endpoints alone.
3. Reconcile the archived radiation continuation replay and close Misner–Sharp flux/work conservation with the existing owners.
4. Only then request independent review of the candidate's source changes. No long synchronized t=12 or t=160 campaign, no merge to production `main`, and no claim of bounce/re-expansion/cycle admission follows from one-step tests.
