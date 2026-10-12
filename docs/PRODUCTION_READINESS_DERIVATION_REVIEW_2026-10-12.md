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
3. A frozen-state audit first found a difference between the old CMC solver and production K operator on the uniform-quadrature initializer. A newer cross-agent derivation corrected the more fundamental problem: the cell-centred grid’s first value is at r=dr/2, but the old CMC solve imposed alpha[0]=alpha[1] as its centre row, which is not the even-parity boundary condition for samples at two distinct positive radii. An isolated candidate now solves the CMC PDE at the first cell using the same parity-aware fourth-order stencils used by the production K RHS. The frozen rows close at roundoff through N=320 and the even-series fit improves dramatically. A matched one-step A/B then showed substantial momentum/connection improvements at N=60–160, but a mixed Hamiltonian response and no demonstration of regular-centre convergence or radiation admission.
4. Radiation admissibility remains a separate blocker, but the two observed failures must be separated. The old N=240/320 combined CMC one-step capture was later traced to the artificial dt=1 RHS-extraction guard, not a real production-stage rejection; a bounded post-contract production-step check passes at both resolutions. The distinct long-run outer-cell trace still shows the explicit predictor metric taking a near-saturated state outside the cone.

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
- Extended current-main beta-source candidate: [run 38302775872](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38302775872), commit 23462c7e2e433402f9fef4b5789a8ebfa9e7795c, 2 targeted source-stage tests and 81 full-suite tests passed in 7.01 s.
- Latest extended current-main source candidate: [run 38303904514](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38303904514), commit 92bfab1e330c9754d9909f0b2cb066baef7e8944, 2 targeted stage-pair tests passed in 0.38 s and the full suite 82 passed in 7.38 s. This version also ports the spatial homogeneous-COSMOS limit regression; it remains isolated from main and is not a production-admitted trajectory.
- CMC centre-row candidate regression suite: [run 38303841261](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38303841261), isolated branch diag/cmc-lapse-full-ci-2026-10-12, 89 tests passed in 9.43 s. Comparison against the CMC candidate head shows this commit added the workflow only; it did not alter the tested production modules. This is regression coverage, not trajectory or science admission.
- The corrected CMC row is now also tracked in open draft [PR #7](https://github.com/miahpenn/0-Cosmos-engine/pull/7), branch fix/cmc-cellcenter-row-2026-10-12 at commit 72483464c6008b4529985c8767b501b11f57fcc9. Both [run 38304045777](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304045777) and [run 38304170442](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304170442) passed 2 focused row tests and the full 81-test repository suite; the latter was the latest code test on that branch. PR is still draft; no merge.
- The latest beta source-candidate A/B snapshot, [run 38304517568](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304517568), commit cd835df3d87f30f92d29afcd9a4b50a54b2877cf, passed 4 focused tests and 83 full-suite tests (11.38 s). An earlier compare attempt [run 38304363097](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304363097) failed during test collection because test_covariant_beta_source_one_step_ab.py had an unclosed parenthesis; that harness failure is superseded by the subsequent successful source A/B run and must not be reported as a physics failure.

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

## 5. Updated CMC centre-row finding — supersedes the earlier mixed-stencil lead

**Independent SI-agent candidate:** [cell-centred CMC row repair branch](https://github.com/miahpenn/0-Cosmos-engine/tree/diag/cmc-lapse-cellcenter-row-repair-2026-10-12).  
**Run:** [38302913119](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38302913119) — successful focused regression and frozen-state audit.  
**Artifact ID:** 11771098012; ZIP SHA-256: 2c64f2f5b8db074f40748e50f22cfa14da71c98398d52bb626282932a79f9d5c.  
**Candidate commit:** 9f7ebf1aab2076fe8925737591366b416427305e.  
**Scope:** N=40/60/80/120/160/240/320 frozen uniform-quadrature initial profiles; no time advance; main and archived GEAR untouched.

### Derivation

The radial mesh is cell-centred, (r_i=(i+1/2)\,\Delta r\). Thus the first stored lapse is at (r=\Delta r/2\), not the origin. The existing CMC matrix row `alpha[0] - alpha[1] = 0` forced equal lapse values at (r=\Delta r/2\) and (3\Delta r/2\). It was not the regular-origin parity condition; it also omitted the CMC PDE at the first cell.

The grid’s even-reflection fourth-order formulas at the first cell are:

`D_r alpha_0 = (-8 alpha_0 + 9 alpha_1 - alpha_2)/(12 dr)`

`D_r^2 alpha_0 = (-14 alpha_0 + 15 alpha_1 - alpha_2)/(12 dr^2)`

The candidate places that PDE row at cell 0 and uses parity-aware fourth-order rows in the interior, a one-sided fourth-order closure at cell N-2, and retains alpha[N-1]=1 as the outer normalization. The CMC target and physical sources are unchanged.

### Frozen-state results

| N | dr | alpha[0]-alpha[1] | first-eight-cell even-series max residual | max CMC PDE-row residual | max production Kdot-target residual |
|---:|---:|---:|---:|---:|---:|
| 40 | 2.000 | 1.528e-3 | 2.554e-4 | 1.77e-16 | 1.77e-16 |
| 60 | 1.333 | 7.175e-4 | 3.500e-5 | 5.62e-16 | 5.61e-16 |
| 80 | 1.000 | 4.122e-4 | 6.295e-6 | 9.81e-16 | 9.79e-16 |
| 120 | 0.667 | 1.862e-4 | 3.830e-7 | 2.21e-15 | 2.20e-15 |
| 160 | 0.500 | 1.053e-4 | 4.494e-8 | 5.12e-15 | 5.11e-15 |
| 240 | 0.333 | 4.700e-5 | 1.964e-9 | 8.73e-15 | 8.73e-15 |
| 320 | 0.250 | 2.648e-5 | 2.046e-10 | 2.00e-14 | 2.00e-14 |

The nonzero `alpha[0]-alpha[1]` is expected for even-profile samples at two different positive radii; forcing this difference to zero was the bad condition. The even-series fit residual falls by roughly six orders of magnitude between N=40 and N=320. All lapses were finite and positive.

**Critical limitation:** the frozen candidate still fails the required `Aa_dot = O(r²)` centre-series criterion. At N=80, `Aa_dot/r_0² ≈ -1.43e-5`; at N=160, `≈ -5.59e-5`; at N=320, `≈ -6.72e-5`. The first-cell `Aa_dot` tends towards zero more slowly than r², while fixed-radius Aa RHS values remain nonzero; this is not centre regularity admission. The row repair resolves a genuine gauge-discretization defect candidate, but it is not a full-machine repair. At the time of this frozen-state report, no one-step test had yet been run; it was subsequently completed and is summarized in the next section.

### Reconciliation with the earlier audit

My preceding frozen uniform-profile audit correctly excluded the old centre row from its interior PDE residual and compared fixed physical radii. Its mixed-order production residual was real for the old solver, but should no longer be treated as the leading explanation in isolation: the actual first-row closure was itself inconsistent with a cell-centred mesh. The new parity-aware candidate deliberately changes the numerical row so the CMC and production K RHS use compatible fourth-order derivative rows. Since the candidate's near-zero production Kdot-target residual is constructed from matching operators, it is a consistency result, not independent proof of constraint propagation or continuum correctness.

**No production patch or merge is justified yet.** The matched one-step A/B is now complete; its mixed results are summarized in the next section. The next test must extend the same A/B witness set to an explicit, separately instrumented radiation-cone budget and a closed Misner–Sharp balance residual before any trajectory or candidate promotion.

The parity-aware candidate also passed the full regression suite on isolated run [38303841261](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38303841261): 89 tests passed in 9.43 s. The comparison from candidate head 9f7ebf1aab2076fe8925737591366b416427305e to that branch contains only the newly added workflow, so the tested candidate sources were unchanged. This is a code-regression pass, not physics or trajectory admission.

## 6. Matched one-step A/B — centre-row candidate has mixed constraint results

**Run:** [38303550121](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38303550121) — completed successfully.  
**Branch:** `diag/cmc-lapse-cellcenter-one-step-ab-2026-10-12`  
**Commit:** `0de3fa09511c919941072390a24426d10681b1fd`  
**Artifact ID:** `11770494466`  
**ZIP SHA-256:** `9c4934ddafed682044aced6222d94caf7c8fcff2434ec00a72a423c4abccf19d`  
**Internal report.json SHA-256:** `1b32400c29ce7340d68b9db590942144fa1d5d57ac71fa6b170ccd28c6ea9bdd`  
**Vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`  
**Scope:** one ordinary production step at N=40/60/80/120/160, r_max=80, dt=0.03 dr. Same initial non-lapse geometry, scalars and conservative matter in each pair; old and candidate differ only in the CMC lapse row/operator. No physical source/parameter/boundary/damping/floor/clamp changed.

### Full-domain constraint maxima after one step

| N | dt | H old → candidate | M old → candidate | C_Lambda old → candidate |
|---:|---:|---:|---:|---:|
| 40 | 0.060 | 3.847e-4 → 3.837e-4 | 1.530e-5 → 2.062e-5 | 9.321e-7 → 1.261e-6 |
| 60 | 0.040 | 1.009e-4 → 9.985e-5 | 5.850e-6 → 3.359e-6 | 2.390e-7 → 1.372e-7 |
| 80 | 0.030 | 3.776e-5 → 3.670e-5 | 6.123e-6 → 4.175e-7 | 1.875e-7 → 1.294e-8 |
| 120 | 0.020 | 1.019e-5 → 1.056e-5 | 6.185e-6 → 3.342e-7 | 1.262e-7 → 6.750e-9 |
| 160 | 0.015 | 5.926e-6 → 6.876e-6 | 6.166e-6 → 2.965e-7 | 9.437e-8 → 4.506e-9 |

The candidate improves full-domain M and C_Lambda strongly at N=60–160, but worsens both full-domain maxima at N=40. H improves at N=40/60/80 and worsens at N=120/160. Because full-domain peaks include the centre cells, these results must be read alongside regional masks—not reduced to a single “pass”.

### Fixed physical annulus (2 ≤ r < 5)

| N | H old → candidate | M old → candidate | C_Lambda old → candidate |
|---:|---:|---:|---:|
| 40 | 2.682e-5 → 2.591e-5 (one cell only) | 4.101e-6 → 1.675e-6 | 2.532e-7 → 1.003e-7 |
| 60 | 9.816e-6 → 1.075e-5 | 5.850e-6 → 2.084e-7 | 2.390e-7 → 8.800e-9 |
| 80 | 3.754e-6 → 4.270e-6 | 3.832e-6 → 1.316e-7 | 1.167e-7 → 4.107e-9 |
| 120 | 3.613e-7 → 3.613e-7 | 9.732e-7 → 4.575e-8 | 1.978e-8 → 9.514e-10 |
| 160 | 1.261e-7 → 1.261e-7 | 8.525e-9 → 3.436e-9 | 1.242e-10 → 5.976e-11 |

For N=60–160, the candidate cuts the annulus M and C_Lambda witnesses materially—often by roughly 20–30×—while H is a little worse at N=60/80 and essentially unchanged at N=120/160. N=40 has no cells in the predeclared admitted-inner mask (i≥2, r<5), so its one-cell annulus value is not a meaningful regional convergence sample.

### Centre evolution witness

| N | old Aa_dot/r_0² | candidate Aa_dot/r_0² |
|---:|---:|---:|
| 40 | +1.446e-4 | +1.053e-4 |
| 60 | +1.166e-4 | +2.318e-5 |
| 80 | +1.543e-4 | -1.426e-5 |
| 120 | +3.385e-4 | -4.433e-5 |
| 160 | +6.270e-4 | -5.557e-5 |

The candidate sharply reduces the magnitude and removes the old positive sign at N≥80, but the candidate’s `Aa_dot/r_0²` remains finite/nonzero and drifts toward approximately `-6e-5` when compared with the frozen N=240/320 profile. That is not a demonstration of `Aa_dot=O(r²)` at the regular origin. This is a promising numerical-layer correction whose centre-series gate remains open.

### Additional controls and interpretation limits

- The CMC target computed from the common alpha=1 input slice is identical between arms at all five resolutions; the candidate’s matched production Kdot-target residual is roundoff, whereas the old row’s final max is approximately 1.42e-3–1.58e-3. The candidate row is built from the same derivative operators as production K RHS, so this residual is a discrete-consistency check—not an independent physics proof.
- At R_sigma≈10, the recorded Misner–Sharp mass, H_eff, flux and work observables differ only slightly between arms. This A/B report does not publish an explicit closed Misner–Sharp balance residual, so conservation admission cannot be inferred from the small observable differences.
- The workflow did not measure radiation cone margin/admissibility at every predictor stage. It used only N=40–160 and cannot close the separate N=240/320 first-cell radiation failure or the long-run outer-cell cone failure.
- The candidate is not authorized to merge. The next gate must join one-step H/M/C_Lambda, signed Aa/K series, radiation predictor-cone terms and a closed Misner–Sharp budget on the same exact states and stages. Only promote if the implicated witnesses improve together.

### Three-arm row-order separation (new, 12 October)

**Run:** [38304709682](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304709682) — successful diagnostic capture.  
**Branch:** `diag/cmc-lapse-three-arm-row-order-2026-10-12`; commit `1fca3f2f03626346720c95f4d4c114b884d1be96`.  
**Artifact:** `cmc-lapse-three-arm-row-order`, ID `11771896263`; ZIP SHA-256 `59918d53bc883fbc796bf79057666692a01c49a0a355bd26f89cbac22fbb08d8`.  
**Scope:** one-step matched A/B/C on N=80/120/160.

The third arm separates the *type of centre row* from its *derivative order*:

- **A:** historical `alpha[0] = alpha[1]` equality row, omitting the first-cell CMC PDE.
- **B:** second-order PDE row at the first cell, instead of the equality row.
- **C:** parity-aware fourth-order PDE row, matched to the production derivative operator.

This is an important discriminator. At N=80, changing A→B reduces admitted-inner momentum max by `3.70e-6`; B and C then differ in annulus M by only `2.64e-10`, annulus H by `1.27e-11`, annulus C_lambda by `8.30e-12`. The same pattern holds at N=120 and N=160: B and C constraint differences remain around (10^{-11}) or smaller for the fixed-annulus witnesses, while A→B changes M at approximately (10^{-6}). C therefore does not need its fourth-order row to reproduce the main one-step H/M/C improvement seen when the missing first-cell PDE replaces the invalid equality condition.

| N | A old-row annulus M | B second-order PDE-row annulus M | C fourth-order PDE-row annulus M | A full H | B full H | C full H |
|---:|---:|---:|---:|---:|---:|---:|
| 80 | (3.83225\times10^{-6}) | (1.31315\times10^{-7}) | (1.31579\times10^{-7}) | (3.77559\times10^{-5}) | (3.67008\times10^{-5}) | (3.67010\times10^{-5}) |
| 120 | (9.73239\times10^{-7}) | (4.57432\times10^{-8}) | (4.57517\times10^{-8}) | (1.01913\times10^{-5}) | (1.05613\times10^{-5}) | (1.05613\times10^{-5}) |
| 160 | (8.52531\times10^{-9}) | (3.43453\times10^{-9}) | (3.43645\times10^{-9}) | (5.92559\times10^{-6}) | (6.87585\times10^{-6}) | (6.87585\times10^{-6}) |

The row-order comparison corroborates that **the principal one-step benefit is from solving the omitted first-cell CMC equation, rather than from fourth-order versus second-order PDE stencil accuracy alone**. The C row remains preferable as the operator-matched candidate, but the B/C similarity on constraints does not mean their centre-series behavior is interchangeable: their cell-0 (\dot A_a/r_0^2) witnesses differ by about (3.96\times10^{-7}) at N=80, (1.93\times10^{-7}) at N=120, and (1.12\times10^{-7}) at N=160. Both still have a nonzero centre-series witness; neither is centre-regularity admission.

Most importantly, the three-arm test still lacks predictor-stage radiation margins and a closed Misner–Sharp interval residual. The new combined one-step radiation/MS diagnostic is running on isolated branch [diag/cmc-lapse-rad-ms-one-step-ab-2026-10-12](https://github.com/miahpenn/0-Cosmos-engine/tree/diag/cmc-lapse-rad-ms-one-step-ab-2026-10-12); no production file is changed.

### Three-arm row-order separation (12 October)

**Run:** [38304709682](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304709682), successful diagnostic capture; commit 1fca3f2f03626346720c95f4d4c114b884d1be96. **Artifact:** ID 11771896263, ZIP SHA-256 59918d53bc883fbc796bf79057666692a01c49a0a355bd26f89cbac22fbb08d8. Resolutions N=80/120/160, one matched production step per arm.

The third arm separates the *type of centre row* from derivative order:

- A: historical alpha[0]=alpha[1] equality row, omitting the first-cell CMC PDE.
- B: second-order CMC PDE row at the first cell instead of the equality row.
- C: parity-aware fourth-order PDE row, matched to the production derivative operator.

At N=80, changing A→B reduces the admitted-inner momentum maximum by 3.70e-6. B and C differ by only about 2.64e-10 in annulus M, 1.27e-11 in annulus H, and 8.30e-12 in annulus C_Lambda. The same pattern holds at N=120/160: the B/C constraint differences are around 1e-11 or smaller for the fixed-annulus witnesses, while A→B changes M at approximately 1e-6.

| N | A old-row annulus M | B second-order PDE-row annulus M | C fourth-order PDE-row annulus M | A full H | B full H | C full H |
|---:|---:|---:|---:|---:|---:|---:|
| 80 | 3.83225e-6 | 1.31315e-7 | 1.31579e-7 | 3.77559e-5 | 3.67008e-5 | 3.67010e-5 |
| 120 | 9.73239e-7 | 4.57432e-8 | 4.57517e-8 | 1.01913e-5 | 1.05613e-5 | 1.05613e-5 |
| 160 | 8.52531e-9 | 3.43453e-9 | 3.43645e-9 | 5.92559e-6 | 6.87585e-6 | 6.87585e-6 |

This supports the conclusion that the principal one-step M/C improvement comes from solving the omitted first-cell CMC equation, not primarily from fourth-order versus second-order stencil accuracy. C remains the operator-matched candidate. But the cell-0 Aa_rhs/r0^2 values differ between B and C by about 3.96e-7 at N=80, 1.93e-7 at N=120, and 1.12e-7 at N=160; neither row has passed centre-series admission. H worsens for B/C relative to A at N=120 and N=160. This is a diagnostic attribution, not permission to merge.

### Matched radiation-cone and one-step Misner–Sharp budget (new)

**Run:** [38304860664](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304860664), successful capture workflow on isolated branch diag/cmc-lapse-rad-ms-one-step-ab-2026-10-12. **Audit commit:** 759c5ead2d1b591667ddbe8ae82dbdc622aee15d. **Candidate engine source:** 9f7ebf1aab2076fe8925737591366b416427305e. **Vendor:** d6052d605673ce9d82cdc99b0df79855fa2ea215. **Artifact:** ID 11772240963, archive digest SHA-256 b63218c17b91a85d329441d9252124f96b6b7095602009cd1bf5e0d24e29e045; report.json SHA-256 0a938cfca543e838cf0415ecd37c450375c34f4c103faab9ba28f212b862d90f. The artifact is also present in the local runtime for this handoff. The helper tests passed 2/2 in 0.27 s; the complete repository suite passed 91 tests in 5.41 s.

The candidate check fail-closed gate verified there were no engine/vendor/dependency changes from the pinned CMC candidate source. The old run did **not** complete both arms at N=240/320: it captured a first-cell exception in both arms. Subsequent source inspection showed that this capture came from the legacy `species_rhs -> evolve_species(dt=1.0)` path and its unconditional cone check, so it is a valid capture of the API-contract defect but not a valid actual-step failure. The old workflow’s “success” still means expected diagnostic capture, not a seven-resolution physical-step pass. A bounded follow-up using the corrected PR #10 RHS/admissibility contract now completes one real kernel step at N=240 and N=320; see the correction section below.

At N=40–160, each arm remained within the radiation cone for this one-step sample. The largest recorded full-domain predictor |S|/E ratios were around 0.01013–0.01018 and the minimum cone margins stayed positive. This low-resolution one-step state is not representative of the long-run near-saturated outer-cell failure at t≈46.98 from the separate radiation trace.

At N=240 and N=320, both arms failed **inside the radiation RHS before the first explicit predictor’s CMC solve**, at cell 0. The recorded failure is therefore not caused by the changed CMC row or its lapse solve. Values:

| N | Old-row RHS failure ratio | Candidate-row RHS failure ratio | Candidate minus old |
|---:|---:|---:|---:|
| 240 | 1.0202778 | 1.0201269 | -0.0001509 |
| 320 | 1.3602899 | 1.3601766 | -0.0001134 |

At N=240, the *legacy unit-step surrogate output* was reported as E≈8.81887e-9 and |S|≈8.99770e-9 (old), versus E≈8.81887e-9 and |S|≈8.99636e-9 (candidate); at N=320 its respective ratios were 1.36029 and 1.36018. These are not the state reached after the requested finite timestep: the old guard rejects `U+L(U)` while `species_rhs` is subtracting the base state to estimate a rate. Do not use these ratios as evidence of a physical first-cell step failure.

The run also computed a one-step trapezoidal Misner–Sharp current residual using ΔM_MS minus dt times the averaged endpoint (flux_T+work_pR). It decreases in normalized magnitude with refinement but remains nonzero, and is almost unchanged by the CMC row:

| N | dt | Old scaled integrated MS residual | Candidate scaled integrated MS residual | Candidate − old integrated residual (absolute) |
|---:|---:|---:|---:|---:|
| 40 | 0.0600 | -5.9588e-2 | -5.9002e-2 | 8.11e-9 |
| 60 | 0.0400 | -1.9503e-2 | -1.9487e-2 | 1.31e-10 |
| 80 | 0.0300 | -6.7155e-3 | -6.7132e-3 | 1.34e-11 |
| 120 | 0.0200 | -9.8935e-4 | -9.8985e-4 | 2.21e-12 |
| 160 | 0.0150 | -3.8669e-4 | -3.8679e-4 | 3.20e-13 |

The unscaled integrated residuals range from approximately -8.48e-7 (N=40) to -1.22e-9 (N=160). The decreasing normalized mismatch is encouraging as a resolution trend, but no tolerance was preregistered and a one-step endpoint quadrature is not a multi-step conservation proof. It must not be reported as “Misner–Sharp conservation passed”. For this A/B it does show the centre-row switch does not materially account for the balance residual.

### Covariant source candidate — latest one-step A/B evidence

Run [38304517568](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304517568), isolated current-main source-candidate snapshot commit cd835df3d87f30f92d29afcd9a4b50a54b2877cf. The run passed 4 targeted tests and 83 full-suite tests. The source stage-pair checks at the two actual kernel RHS stages reported rest-source errors 3.96e-19 and 4.15e-19, energy-pair residuals 3.96e-19 and 6.73e-19, and radial-current residuals 2.51e-24 and 1.72e-24.

For its matched N=40, rmax=16, dt=0.012 one-step source A/B, the candidate-minus-baseline changes were extremely small in the derived global witnesses: max H changed by about -8.49e-11, admitted-inner M by about -6.0e-16, all-domain connection max by about -3.7e-17, Misner–Sharp mass at R≈10 by about -2.27e-13, and first-cell radiation margin by about -1.67e-22. Source terms themselves differ nontrivially and close pointwise, but the measured global step/constraint response in this test state is tiny. Therefore this supports the correctness of the reciprocal covariant source accounting on the tested stages; it does **not** show that this term is the dominant full-machine constraint defect or that the physical trajectory is materially repaired.

## 7. Covariant scalar–DM source pair: candidate passes, trajectory gate remains open

The locked V5.5 covariant exchange law is Q_nu = beta_DM rho_DM partial_nu(phi), with reciprocal scalar and dark-matter sources. On the ADM coordinate-time convention used by production, the source pair requires:

- Scalar momentum: Pi_dot += alpha * beta_DM * rho_DM.
- Conservative DM rest-current source: q_rest = -alpha*sqrt(gamma)*beta_DM*rho_DM*W*(Pi + v^r*phi_r).
- Conservative DM energy source: q_e = -sqrt(gamma)*(Q_t - beta^r*Q_r), where Q_t uses the coordinate-time phi_dot.
- Conservative DM radial-momentum source: q_s = alpha*sqrt(gamma)*Q_r.

The newest isolated current-main source candidate is [commit 92bfab1e330c9754d9909f0b2cb066baef7e8944](https://github.com/miahpenn/0-Cosmos-engine/commit/92bfab1e330c9754d9909f0b2cb066baef7e8944), run [38303904514](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38303904514):
- Targeted source-stage pair tests: 2 passed in 0.38 s.
- Full repository suite: 82 passed in 7.38 s.
- Stage 0: lapse deviation 9.029416e-3; rest-source error 3.984231e-19; energy-pair residual 3.983172e-19; radial-current residual 2.487600e-24.
- Stage 1: lapse deviation 9.026827e-3; rest-source error 4.206577e-19; energy-pair residual 6.680972e-19; radial-current residual 1.772703e-24.

The same current-main candidate also ports tests/test_spatial_homogeneous_cosmos_limit_current_main.py. That deterministic unadvanced FLRW-slice test compares the spatial scalar, DM rest/energy, baryons/radiation and the BSSN trace equation against the archived homogeneous COSMOS RHS and checks the initial Hamiltonian/momentum constraints. It excludes the archival Bianchi-I shear consistently with the exact-spherical V5.5 scope. This materially advances the homogeneous-limit/state-ownership gate while remaining a code/RHS regression, not an evolved cosmological trajectory.

These passes establish a source-level candidate and tightly closed synthetic/actual-stage identities on the tested states; they do not establish full stress-energy conservation through every PIRK stage, constraint propagation, radiation admissibility, or production campaign success. They also do not close the separate D-mode-to-matter or old two-lane local-source-to-global-Q bridges. Main remains unchanged.

## 8. Radiation admissibility remains independent

The previous combined one-step CMC A/B recorded N=240/320 first-cell ratios above unity, but these were generated by the obsolete dt=1 RHS-extraction guard. A direct bounded follow-up on the corrected RHS/admissibility candidate performs one normal production step at N=240 (`dt=0.01`) and N=320 (`dt=0.0075`); the initial, explicit-predictor, primary-predictor, and final cone checks all pass. See [high-resolution actual-step report](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/high-resolution-radiation-real-step-after-rhs-contract-2026-10-12/docs/HIGH_RESOLUTION_RADIATION_REAL_STEP_AFTER_RHS_CONTRACT_2026-10-12.md), run [38318958640](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318958640), report SHA-256 `21742fcefebb3155b61aab4963086a0e31875756fe17c1af03cce09146f07453`. A distinct long-run trace at t≈46.98 still records the outer near-saturated cell crossing during the explicit predictor-metric update; on the accepted metric, flux/source increments increased margin, while the predictor metric reduced it by more than the remaining margin. The N=240/320 one-step pass does not establish long-run high-resolution stability or identify a bad radiation law.

Do not clip conserved variables, introduce floors/damping, or relax |S|/E ≤ 1. Treat the obsolete N=240/320 RHS-extraction exception as an API-contract artifact, not a physical first-cell failure. Keep it separate from the still-open long-run outer-domain / global-CMC / geometry-stage interaction.

References: [radiation predictor report](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/radiation-predictor-stage-ledger/docs/RADIATION_PREDICTOR_TRACE_2026-10-10.md), [P-minus audit run](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38028306184).

## 9. Other valid but non-dominant equation lead

The spatial metric has gamma_rr = a/X² and gamma^rr = X²/a. An isolated scalar candidate corrected a lapse-gradient contraction that used X^-2 instead of the declared inverse radial metric. Its matched N=40/60/80 campaign reached t=12, but H_eff/proper-time changes were tiny and centre-dominated constraint/Misner–Sharp defects remained open. Keep this source-correctness issue isolated; it is not the whole machine repair. See the [candidate assessment](https://github.com/miahpenn/0-Cosmos-engine/blob/fix/scalar-inverse-metric-lapse-gradient-strict-2026-10-12/docs/SCALAR_INVERSE_METRIC_CANDIDATE_ASSESSMENT_2026-10-12.md).

## 10. Recommended next gates

1. Reconcile the completed one-step A/B against the signed Aa/K RHS split and discrete constraint-tangent result on identical states. The candidate improves M/C_Lambda at N=60–160, but full-domain H is slightly worse at N=120/160, and the first-cell Aa_dot/r² remains nonzero.
2. Keep the parity-aware centre row isolated; do not merge from this one-step result. Add a predictor-stage radiation-cone budget and an explicit Misner–Sharp balance residual to the same A/B instrumentation before making a promotion decision.
3. Keep the passed synthetic source-stage cancellation result; next use a production-relevant saved state with nonzero phi_r to confirm the same paired source terms at actual stage levels.
4. Reuse existing radiation predictor/P-minus artifacts first. If another diagnostic is needed, record raw metric RHS, unprojected predictor and post-regularity-projection geometry separately.
5. Complete the spatial homogeneous-limit regression and state-ownership guards. Keep CosmosState a reference/control; do not invent Q, D-to-matter, clock scaling, or a spherical-shear insertion.
6. Only after these gates pass: freeze an exact commit, pass exact-commit CI, preregister multi-resolution tolerances, and run one full controlled trajectory at a time.

## 11. Standing guardrails

Preserve the GEAR archive and REV16. Do not merge diagnostic branches merely because a unit suite passes. No fitted coefficients, arbitrary damping, floors/clamps, centre projection, forced bounce/reset, D-to-matter law or phenomenological global Q. Distinguish workflow success, target-time completion, unit-test success and scientific admission. Do not infer convergence from centre-cell sequences whose physical radii change with N.

**Bottom line:** the cell-centred CMC equality row remains a supported numerical-discretization defect candidate; replacing it with a first-cell PDE row improves one-step M/C, but H responds non-monotonically and `Aa_dot/r²` remains nonzero. The earlier N=240/320 first-cell radiation rejections have now been traced to the artificial unit-step RHS guard, and one-step production evolution on the corrected RHS/admissibility candidate passes at both resolutions. The separate t≈46.98 outer-cell predictor-metric cone crossing still remains, and the one-step Misner–Sharp balance residual is nonzero with a declining but unadmitted resolution trend. The covariant scalar–DM source pair passes source-stage tests plus a homogeneous-limit RHS regression, while its tested one-step global response is very small. Do not start a long campaign or merge from these results alone. Main and the canonical GEAR archive remain untouched.


### 2026-10-12 post-contract correction: high-resolution one-step artifact reclassified

The combined CMC/radiation/Misner–Sharp artifact was pinned to source commit `9f7ebf1aab2076fe8925737591366b416427305e`. At that SHA, `species_rhs` obtains the rate by calling `evolve_species(dt=1.0)`, and `evolve_species` unconditionally checks radiation admissibility on that output before subtracting the input state. Thus the N=240/320 `E,|S|` values from the old captured rejection describe the artificial `U+L(U)` state, not the normal `dt=0.01/0.0075` production step.

After the RHS contract fix, a bounded test ran one real production step at each resolution on the isolated PR #10 candidate, leaving the cone criterion and all physics unchanged. CI [38318958640](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38318958640) passed 4 radiation admissibility tests and 26 full-suite tests; all four actual-stage validators passed at both N=240 and N=320. The final max `|S|/E` was 0.00991926171186 (N=240) and 0.00991880271985 (N=320). The archived JSON is checksum-verified in the [actual-step report](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/high-resolution-radiation-real-step-after-rhs-contract-2026-10-12/docs/HIGH_RESOLUTION_RADIATION_REAL_STEP_AFTER_RHS_CONTRACT_2026-10-12.md); report SHA-256 `21742fcefebb3155b61aab4963086a0e31875756fe17c1af03cce09146f07453`.

This correction supersedes the earlier interpretation of those particular N=240/320 rejections as real finite-step radiation failures. It does not show multi-step stability at high resolution and does not change the status of the separate long-run outer-cell crossing near t≈46.98. The old artifact remains admissible evidence of the RHS API-contract defect, not as a physical-state failure capture.
