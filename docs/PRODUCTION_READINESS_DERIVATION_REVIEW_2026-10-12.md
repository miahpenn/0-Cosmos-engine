# 0-Cosmos Engine — Production-Readiness Derivation Review

**Date:** 12 October 2026 (UTC)  
**Production branch observed:** main, commit 2035acc1210db18da0cededea18750b7ed813260  
**Frozen source used in the centre/CMC numerical audit:** 6d8b3f2790644591326df57c55908863f994cb0c  
**Pinned geometry/PIRK vendor:** d6052d605673ce9d82cdc99b0df79855fa2ea215  
**Status:** Research/doc-only. No production physics or code changed, and no merge is requested.

## Executive result

The machine is implemented and operational at component/smoke level, but is **not yet scientifically admitted as a full production evolution**. Four findings must remain separate:

1. Radial quadrature materially changes the initial Hamiltonian residual, but the best tested initializer still fails the required regular-centre Aa_dot = O(r²) test.
2. Initial H/M/C_Lambda closure at roundoff does not stop a nonzero momentum-constraint tangent from developing under the unchanged production RHS. The tangent agrees closely with the one-step quotient.
3. A new frozen-state audit finds the CMC solver’s tridiagonal interior rows close at roundoff, while the production higher-order K operator differs from the solved CMC target on the uniform-quadrature initializer. The difference decreases about as expected at fixed radii for the mixed second-order-versus-higher-order discretizations, but the centre-adjacent maximum remains elevated.
4. Radiation admissibility remains a separate blocker: N=240/320 one-step cases violate |S|/E ≤ 1 in the first cell for all tested initializer arms, and a distinct long-run outer-cell trace shows the explicit predictor metric taking a near-saturated radiation state outside the cone.

The evidence does not support inserting a guessed global source Q, wiring a duplicate homogeneous COSMOS state into the V5.5 spatial step, or declaring a physical bounce/cycle.

## 1. Blueprint authority and architecture

The Library master archive identifies CYCLIC-COSMOS ENGINE-.docx as the preserved GEAR two-lane baseline: local GEAR-03 S/D modes and GEAR-63/65 barrier/ledgers, plus homogeneous COSMOS scalar, dark matter, baryons, radiation, shear, expansion/contraction and cycle coordinates. The archive leaves local proper time ↔ cosmic time, D energy → ordinary matter, local source → global Q, full-GR oscillaton endpoint and a first-principles bounce open. Preserve the source byte-for-byte.

The later V5.5 production target is a single spatially resolved spherical spacetime and one total stress-energy tensor. Local S/D, spatial COSMOS phi and conservative DM/baryon/radiation evolve on the shared metric. The homogeneous CosmosState in engine/cosmos.py is a reference/control, not a second live state to insert into ProductionState. A missing call to cosmos_rhs() in the production step is not, by itself, a missing-wire defect; the spatial COSMOS sector needs a matched homogeneous-limit regression.

REV16 remains authoritative: 0*_eff iff t ∧ b ∧ (w ≥ 1) and 0*_global iff 0*_eff ∧ G. The recorded score has t paid, w=1 saturated, b not payable in the matter lane, and G open, so 0*_eff remains FALSE. Numerical progress does not amend that result.

References: [GEAR dynamic-interface audit](https://github.com/miahpenn/0-Cosmos-engine/blob/research/gear-dynamic-interface-audit-2026-10-12/docs/GEAR_DYNAMIC_INTERFACE_BLUEPRINT_AUDIT_2026-10-12.md), [integration gates](https://github.com/miahpenn/0-Cosmos-engine/blob/research/gear-kernel-integration-gates/docs/GEAR_KERNEL_ARCHITECTURE_AND_INTEGRATION_GATES_2026-10-10.md), [canonical Library archive index](https://github.com/miahpenn/0-Cosmos-engine/tree/main).

## 2. Current-main test accounting

At review time main is 2035acc1210db18da0cededea18750b7ed813260. No combined commit status checks were returned for this exact SHA; the recently triggered Actions jobs on main were skipped because the workflow path/branch conditions did not match. The main STATUS.md records an earlier 79/79 test suite and a resolution-40, r_max=40, t=1 smoke run. That status explicitly says implemented/tested at unit and smoke level, not full-campaign-validated.

Do not combine test counts from independent branches as though they were all run on main:
- Historical architecture-gate branch: 110 tests passed; optional long campaigns skipped.
- Historical lapse-source candidate: 111 tests passed; no long physical trajectory.
- Current-main beta-source candidate at commit 88e4d2e809: 1 targeted test and 80 full-suite tests passed.
- Extended current-main beta-source candidate: [run 38302775872](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38302775872), commit 23462c7e2e433402f9fef4b5789a8ebfa9e7795c, 2 targeted source-stage tests and 81 full-suite tests passed in 7.01 s. This remained isolated from main and was not a production-admitted trajectory.

## 3. Initial data and constraint propagation

### Quadrature helps initial H, but is not the repair

The origin half-cell diagnostic replaced the baseline (dr/2)·source[0] term by an integral of a quadratic interpolant in r²; the uniform arm also replaced centre-to-centre trapezoids with cubic Lagrange interval integration. These were diagnostic initializers only; production field equations were unchanged.

The uniform arm reduces the initial first-cell |H| from roughly 4.05e-3 in the pinned initializer to 4.99e-7 at N=320; its formal-stencil H maximum reaches about 1.60e-7 and the fixed annulus 2 ≤ r < 5 reaches about 7.95e-9. The improvement survives one normal step through N=160.

Yet the uniform arm’s first-cell Aa_dot stays around +4.0 to +4.14e-5 from N=80 through N=320 instead of scaling as r². Aa_dot/r² increases from about 1.53e-4 at N=80 to 6.21e-4 at N=160 and about 2.65e-3 at N=320. Initial-H improvement is real but does not pass the centre-series gate.

### Constraint-consistent initial data still develops M

On the discrete-consistent candidate, initial H/M/C_Lambda close at roundoff by construction. The unchanged production RHS gives candidate admitted-inner M tangents of approximately -8.76e-5 (N=60), -1.24e-4 (N=80) and -1.91e-4 (N=120); the one-step Delta M/Delta t values agree closely. The first-admitted radius changes with resolution, so those values are not a fixed-radius convergence sequence. This localizes the leading problem in the semi-discrete RHS, beyond the initial Hamiltonian residual or report-join logic.

Smooth manufactured regular-centre Ricci/Aa controls pass at approximately fourth order on a fixed annulus. That disfavors a generic vendor parity/sign error but does not validate the actual initialized centre profile.

References: [latest centre-series/constraint handoff](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/discrete-consistent-initial-production-test-2026-10-12/docs/CONTINUATION_HANDOFF_LATEST_2026-10-12.md), [initial-metric A/B run](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38284401548), [matched one-step quadrature run](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38292762118), [constraint-propagation handoff](https://github.com/miahpenn/0-Cosmos-engine/blob/docs/constraint-propagation-handoff-2026-10-12/docs/CONTINUATION_HANDOFF_DISCRETE_CONSTRAINT_PROPAGATION_2026-10-12.md).

## 4. New uniform-profile CMC operator audit

**Run:** [38302763542](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38302763542)  
**Audit commit:** 9e8567f52a4b782f8003b24109707c0d0f330792  
**Artifact ID:** 11770503145  
**Artifact ZIP SHA-256:** 0fb5038ad589a6850d38d38c83ff9c4cf63c88c33f6b464349fb67b5cf3f9349  
**Internal report.json SHA-256:** dc79858d39ffc68d2ff4dffc651c30e7e2d13ffc4ccecb63a6f5e517b1f30bb8  
**Source/vendor:** 6d8b3f2790644591326df57c55908863f994cb0c / d6052d605673ce9d82cdc99b0df79855fa2ea215  
**Scope:** frozen initial profiles, N=40/60/80/120/160; no time advance, no production edit.

The CMC solve uses a tridiagonal centred three-point interior alpha row. The production K RHS uses the pinned grid’s higher-order derivative operators. The exact solver row closes to roundoff, while the production operator residual against the target is finite:

| N | dr | Max production-operator residual on PDE rows | Max CMC discrete-row residual |
|---:|---:|---:|---:|
| 40 | 2.000 | 8.328e-5 | 1.58e-16 |
| 60 | 1.333 | 6.779e-5 | 2.34e-16 |
| 80 | 1.000 | 5.646e-5 | 6.22e-16 |
| 120 | 0.667 | 4.577e-5 | 1.08e-15 |
| 160 | 0.500 | 4.143e-5 | 2.06e-15 |

The first PDE row moves inward with refinement (r=3, 2, 1.5, 1, 0.75), so the max-over-row sequence is not a fixed-radius convergence measure. Interpolating only interior PDE rows gives this signed production residual at fixed radii:

| N | r=2 | r=3 | r=4 | r=5 |
|---:|---:|---:|---:|---:|
| 40 | unavailable | -8.328e-5 | -4.783e-5 | -1.238e-5 |
| 60 | -6.779e-5 | -3.175e-5 | -1.375e-5 | -5.939e-6 |
| 80 | -3.635e-5 | -1.328e-5 | -7.725e-6 | -3.264e-6 |
| 120 | -8.559e-6 | -5.985e-6 | -3.426e-6 | -1.350e-6 |
| 160 | -4.854e-6 | -3.383e-6 | -1.925e-6 | -7.724e-7 |

At fixed radii, the gap decreases approximately as expected from a second-order CMC solve being compared with higher-order production derivatives. This is evidence of a mixed-discretization compatibility gap at finite resolution, not by itself proof that the continuum CMC equation is wrong. The centre-adjacent maximum remains more elevated than the fixed-radius sequence, so it deserves a same-state centre-series/RHS comparison with the signed Aa/K budget.

The audit column “production minus CMC lapse term” is exactly zero because it reconstructs the same high-order lapse term using the same derivative expression; that zero is an identity, not independent validation. The discriminating result is the complete production-operator residual versus the actual discrete row solved. The lapse was finite and positive at all resolutions, alpha[0]-alpha[1] was roundoff, and the outer normalization residual declined from about -2.74e-10 (N=40) to -3.02e-12 (N=160).

This result is specific to the uniform-quadrature initial profiles. The prior archived-state centre-row audit used different frozen profiles and reported roundoff-level operator closure; do not extrapolate either result to every evolving state.

## 5. Covariant scalar–DM source pair: candidate passes, trajectory gate remains open

The locked V5.5 covariant exchange law is Q_nu = beta_DM rho_DM partial_nu(phi), with the reciprocal source entering scalar and DM equations. With the active ADM coordinate-time convention, the candidate source pair requires alpha weighting in the scalar Pi source and DM radial conservative-momentum source.

On isolated candidate branch diag/current-main-beta-source-gate-2026-10-12, [run 38302775872](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38302775872) passed both source-pair stage tests and 81 full-suite tests. Its one small step was from an initialized, deliberately perturbed diagnostic state:
- Stage 0 energy-pair residual 3.98e-19, radial-current residual 2.49e-24.
- Stage 1 energy-pair residual 5.52e-19, radial-current residual 2.22e-24.
- Maximum lapse deviation from unity about 9.03e-3.

This is useful source-stage evidence, not full stress-energy conservation at every PIRK update, constraint admission, radiation admission, or a production trajectory. Main remains unchanged.

## 6. Radiation admissibility remains independent

The matched quadrature one-step campaign fails at N=240/320 in the first cell for every tested initialization arm, with |S|/E > 1. A distinct long-run trace at t≈46.98, central proper time tau≈22.0161, records a near-saturated outer cell crossing the cone during the explicit predictor-metric update. On the accepted metric the flux/source increments increased the margin; the predictor metric reduced it by more than the remaining margin. Expanding the domain moved the failure outward rather than curing it. These are captured diagnostic failures, not successful physical trajectories and not proof of a bad radiation law.

Do not clip conserved variables, introduce floors/damping, or relax |S|/E ≤ 1. Keep first-cell high-resolution failures separate from the long-run outer-domain / global-CMC / geometry-stage interaction.

References: [radiation predictor report](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/radiation-predictor-stage-ledger/docs/RADIATION_PREDICTOR_TRACE_2026-10-10.md), [P-minus audit run](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38028306184).

## 7. Other valid but non-dominant equation lead

The spatial metric has gamma_rr = a/X² and gamma^rr = X²/a. An isolated scalar candidate corrected a lapse-gradient contraction that used X^-2 instead of the declared inverse radial metric. Its matched N=40/60/80 campaign reached t=12, but H_eff/proper-time changes were tiny and centre-dominated constraint/Misner–Sharp defects remained open. Keep this source-correctness issue isolated; it is not the whole machine repair. See the [candidate assessment](https://github.com/miahpenn/0-Cosmos-engine/blob/fix/scalar-inverse-metric-lapse-gradient-strict-2026-10-12/docs/SCALAR_INVERSE_METRIC_CANDIDATE_ASSESSMENT_2026-10-12.md).

## 8. Recommended next gates

1. Reconcile CMC fixed-radius residuals, centre-cell lapse profile, signed Aa/K source split and M tangent on identical frozen states with the same source hash, initializer, grid and time level. Keep the centre boundary row, first interior PDE row, first formally admitted cell and fixed physical radii separate.
2. If a specific discrete identity is implicated, run one diagnostic-layer counterfactual changing only that layer. Compare H, M, C_Lambda, Aa_dot/r², determinant/Lambda identities, radiation admissibility and Misner–Sharp balance together. Earlier CMC counterfactuals improved one local K witness but worsened admitted-inner M by about 3e-6, so a local residual is not enough.
3. Keep the passed synthetic source-stage cancellation result; next use a production-relevant saved state with nonzero phi_r to confirm the same paired source terms at actual stage levels.
4. Reuse existing radiation predictor/P-minus artifacts first. If another diagnostic is needed, record raw metric RHS, unprojected predictor and post-regularity-projection geometry separately.
5. Complete the spatial homogeneous-limit regression and state-ownership guards. Keep CosmosState a reference/control; do not invent Q, D-to-matter, clock scaling, or a spherical-shear insertion.
6. Only after these gates pass: freeze an exact commit, pass exact-commit CI, preregister multi-resolution tolerances, and run one full controlled trajectory at a time.

## Standing guardrails

Preserve the GEAR archive and REV16. Do not merge diagnostic branches merely because a unit suite passes. No fitted coefficients, arbitrary damping, floors/clamps, centre projection, forced bounce/reset, D-to-matter law or phenomenological global Q. Distinguish workflow success, target-time completion, unit-test success and scientific admission. Do not infer convergence from centre-cell sequences whose physical radii change with N.

**Bottom line:** the next useful action is a matched frozen-state centre/RHS derivation, not another invented interface or long campaign. The source-stage beta pair has a successful isolated gate; the centre-series/RHS compatibility and radiation-admissibility gates remain open.
