# Continuation handoff — discrete-consistent initialization and quadrature
**Date:** 2026-10-12  
**Repository:** `miahpenn/0-Cosmos-engine`  
**Diagnostic branch:** `diag/discrete-consistent-initial-production-test-2026-10-12`  
**Pinned production source:** `6d8b3f2790644591326df57c55908863f994cb0c`  
**Pinned vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`

## Scope and guardrails

This handoff records new centre/initial-data diagnostics added after the existing continuation handoffs. These are **diagnostic-only branches and workflows**. No production physics equation was edited, no diagnostic branch was merged to `main`, and no long evolution campaign was launched. Keep the existing rule: establish a reproducible numerical defect, alter only its owning layer, then rerun the same bounded full machine. Do not claim physical admission from initializer-only experiments.

## New workflows and evidence

- [Origin half-cell quadrature A/B](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290342829) — completed successfully. Earlier failed attempts are preserved in run history; final run explicitly restored the pinned initializer sequence before comparing.
- [Origin quadrature initial refinement, N=40–320](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38291486881) — succeeded; artifact `origin-quadrature-initial-refinement`, artifact ID `11767325140`, ZIP SHA-256 `42d898b6f195bdd823fba35a240cdb58bb557146c828c87f8b3a8b2a58636638`.
- [Uniform high-order radial quadrature initial profile, N=40–320](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38291938761) — succeeded; artifact `uniform-radial-quadrature-initial-profile`, artifact ID `11767500772`, ZIP SHA-256 `937a139ffd360ccae753f06572da34a400c26183dbdd0f308f1f8d67ef2bc918`.
- [Matched one-step uniform quadrature propagation](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38292249854) — failed before producing a report. The script compiled and dependency install passed. During the first attempted ordinary production step, `matter_system.evolve_species` rejected the radiation RHS at cell (i=0): (E=8.8188697393\times10^{-9}), (|S|=8.9975835306\times10^{-9}), ratio (1.02026493>1); reported source/tendency components were approximately (10^{-9}). The workflow uploaded no result artifact. Treat this as a failed diagnostic harness attempt, **not** evidence that a quadrature arm is physically invalid or that production evolution passed. Fix only the test protocol: capture each arm separately and record the exact first failing arm, initial state and RHS; do not add floors/clamps or relax production admissibility.

## Main result: the origin half-cell integral is a real contributor, but not the whole defect

The origin-only A/B comparison changed only the origin half-cell quadrature. On the fixed admitted region (i\ge2, r<5), the maximum absolute Hamiltonian residual decreased substantially at N=160/240/320:

| N | Pinned baseline | Origin half-cell only |
|---:|---:|---:|
| 160 | (6.6860\times10^{-5}) | (1.8530\times10^{-5}) |
| 240 | (7.9164\times10^{-5}) | (2.9550\times10^{-5}) |
| 320 | (8.3705\times10^{-5}) | (3.3877\times10^{-5}) |

The first-cell H residual falls from approximately (+4.04\times10^{-3}) to (+6.5\times10^{-5}). But the first-cell evolved (A_a) RHS remains finite rather than scaling away like (r^2): it tends toward about (-9.5\times10^{-5}) in the origin-only arm. Thus the origin half-cell error is real, but another discrete integration/operator mismatch remains.

## Stronger initializer-only result: uniform fourth-order interval integration

A separate frozen-(t=0) three-arm audit compared:
1. pinned baseline;
2. origin half-cell correction only;
3. origin correction plus cubic Lagrange integration on all centre-to-centre intervals.

The third arm changes only cumulative radial-integral quadrature. It leaves the physical source, six fixed-point iterations, matter initialization and CMC lapse solve untouched. Interior interval weights are analytic cubic-integral weights ((-1,13,13,-1)/24); the first interval uses an even-reflection ghost; the final interval uses one-sided cubic weights ((1,-5,19,9)/24). No fitted parameters, damping, clamps or projections are used.

Fine-grid values:

| N | Baseline H max, (i\ge2,r<5) | Origin-only | Uniform fourth-order |
|---:|---:|---:|---:|
| 160 | (6.6860\times10^{-5}) | (1.8530\times10^{-5}) | (3.8050\times10^{-7}) |
| 240 | (7.9164\times10^{-5}) | (2.9550\times10^{-5}) | (2.2079\times10^{-7}) |
| 320 | (8.3705\times10^{-5}) | (3.3877\times10^{-5}) | (1.5963\times10^{-7}) |

At N=320, the first-cell H residual is (+4.0543\times10^{-3}) baseline, (+6.6634\times10^{-5}) origin-only, and (+4.9875\times10^{-7}) uniform fourth-order. The fixed annulus (2\le r<5) H RMS at N=160/240/320 changes from baseline (3.7627\times10^{-6},1.7277\times10^{-6},9.8324\times10^{-7}) to uniform fourth-order (7.9046\times10^{-8},1.5753\times10^{-8},5.0056\times10^{-9}). This is a large, systematic initializer-level improvement in H consistency. The momentum residual remains at floating-point scale (~(10^{-18})) in these initial profiles.

However, the uniform arm's first-cell (A_a) RHS is still nonzero: (3.8843\times10^{-5}), (4.0613\times10^{-5}), (4.1368\times10^{-5}) at N=160/240/320, respectively. It is much smaller than baseline ((sim1.0874\times10^{-3})) but does not yet show the expected (O(r^2)) centre scaling. This is a promising numerical repair candidate, **not yet admitted**.

## Interpretation limits

- These are frozen-initializer constraint and regularity tests, not a physical evolution result.
- Improvement in H does not establish that the full discrete Einstein/scalar/matter operator set is consistent.
- Do not infer that the long-run t=12 H defect is entirely explained by this initializer quadrature. Need one ordinary step and then a bounded N=40/60/80 t=12 matched candidate-vs-baseline campaign only if the one-step harness is corrected and all gates pass.
- The failed one-step workflow points to a strict radiation admissibility issue at the centre in the attempted initial step. Preserve it exactly as evidence; do not hide it with a clamp or weakened criterion.
- Keep the scalar lapse-gradient inverse-metric lead separate. It is a distinct source-term consistency hypothesis and remains unproven unless independently compared.
- Do not merge this diagnostic branch into `main`. The production t=12 campaign remains execution/horizon admitted only; Hamiltonian/momentum and balance admission remain open. Do not launch t=160.

## Next bounded action

1. Correct the one-step diagnostic harness so each initializer arm is run independently and its initial state/RHS is captured before attempting the step. Report the first failing arm and cell, without changing production admissibility.
2. Check whether the radiation centre failure is present in the pinned baseline at the same (N), CFL and exact source SHA. If so, classify it as an existing one-step admissibility gate; if only candidate arms trigger it, compare the initial radiation moments and source terms.
3. Independently verify the cubic interval weights and endpoint treatment against polynomial moments and an analytic smooth radial integral; record convergence order and positivity for every (B) profile.
4. Only after the one-step test has useful, complete evidence, consider a bounded full (N=40/60/80), (t=12) matched candidate run. Keep all physics untouched and compare H/M/C\(_\Lambda\), determinant, radiation admissibility, scalar/matter source balance, proper time and Misner–Sharp profiles at identical times/radii.
5. No longer run and no claim of cycle/bounce/physical success before constraint and balance admission.

## Reproducibility files on the diagnostic branch

- `tools/audit_origin_quadrature_initial_refinement.py`
- `tools/audit_uniform_radial_quadrature_initial_profile.py`
- `tools/audit_uniform_radial_quadrature_one_step.py`
- `tools/origin_half_cell_quadrature_regularity_ab.py`
- `.github/workflows/origin-quadrature-initial-refinement.yml`
- `.github/workflows/uniform-radial-quadrature-initial-profile.yml`
- `.github/workflows/uniform-radial-quadrature-one-step.yml`
- `docs/DISCRETE_CONSISTENT_INITIAL_PRODUCTION_AB_2026-10-12.md` (previous handoff, updated in this branch)

