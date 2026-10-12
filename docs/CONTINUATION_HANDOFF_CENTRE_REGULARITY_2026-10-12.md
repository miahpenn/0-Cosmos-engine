# Continuation handoff — central regularity and constraint propagation
**Date:** 2026-10-12  
**Status:** diagnostic-only; production admission remains blocked  
**Scope:** initial radial quadrature, independent regular-centre curvature/RHS identity, and frozen-state constraint tangents.

## Executive finding

The investigation now separates three effects that must not be conflated:

1. **Initial Hamiltonian residual is strongly affected by radial quadrature.** Replacing only the origin half-cell integral substantially reduces the first-cell H residual, but it leaves a plateau near 6.5e-5–6.7e-5 at N=160–320. Replacing the remaining center-to-center trapezoids with a uniform fourth-order interval quadrature reduces the first-cell H residual further, to 4.99e-7 at N=320, and the formal-region maximum to 1.60e-7.
2. **The production trace-free curvature/Aa tendency still fails the regular-centre limit test.** Even with the uniform fourth-order quadrature, the first-cell raw Aa RHS remains about +4.14e-5 at N=320 rather than tending to zero; Aa_rhs/r² grows to about 2.65e-3 as the first cell moves inward. The improved H initial slice is therefore not, by itself, an admitted production repair.
3. **The semi-discrete momentum constraint has a real RHS tendency.** On the paired discrete-consistent initial states, finite-differenced production RHS tangents agree closely with the ordinary one-step change. The initial H Newton solve does not make the momentum RHS tendency vanish. This is not explained solely by the initial H residual or by timestep differencing.

No production physics equation, scalar/matter source, coefficient, gauge, or boundary condition was changed. No long trajectory was run in these latest audits. Main/production is not admitted by this work.

## Runs and artifacts

| Diagnostic | Run | Artifact | Artifact ZIP SHA-256 |
|---|---|---|---|
| Uniform fourth-order radial quadrature initial profile | [38291938761](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38291938761) | [11767500772](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38291938761/artifacts/11767500772) | `937a139ffd360ccae753f06572da34a400c26183dbdd0f308f1f8d67ef2bc918` |
| Origin quadrature refinement N=40–320 | [38291486881](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38291486881) | [11767325140](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38291486881/artifacts/11767325140) | `42d898b6f195bdd823fba35a240cdb58bb557146c828c87f8b3a8b2a58636638` |
| Signed initial H budget | [38291098188](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38291098188) | [11765137827](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38291098188/artifacts/11765137827) | `92bff9be02569772505794f3f4bad5853bd451c745c32d287db12951a51333c2` |
| Matched origin-quadrature one-step propagation | [38290654079](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290654079) | [11766655192](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290654079/artifacts/11766655192) | `a511f67db4d9a1a12ce27652c82a200dbb4b876266fa3045bf8f995ddc087987` |
| Independent regular-centre Ricci/Aa RHS identity | [38290099699](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290099699) | [11765776185](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290099699/artifacts/11765776185) | `345f91141bb19920051935b9c5c516bf1ddc798046012a36583426d9e5e34e97` |
| Discrete-consistent initial Lambda/Ricci sensitivity | [38289657512](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38289657512) | [11764557943](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38289657512/artifacts/11764557943) | `5cb556cc9a9311fb6ba1a86b1daf398e4c01c9ffeb7f458960ff45e767bb47a2` |
| Production initial trace-free Ricci cross-check | [38289275758](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38289275758) | [11764109840](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38289275758/artifacts/11764109840) | `d0525ab5494c7230c428bcaba1503aab3f12e7720316aefd4c50bbb604b82718` |
| Discrete-initial multi-constraint RHS tangent | [38290633971](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290633971) | [11766405316](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290633971/artifacts/11766405316) | `2011162ba822bea1b52d2b8c2099284f97f7a0ca5e76e2d2bdc1788c324a7186` |

Pinned production source for the frozen-state tangent was `6d8b3f2790644591326df57c55908863f994cb0c`; the vendor used in the independent Ricci cross-check is `d6052d605673ce9d82cdc99b0df79855fa2ea215`. The audit runs are diagnostic branches; their completion is not production-campaign admission.

## A. Initial radial quadrature: what changed

The A/B test modifies numerical quadrature only; physical source definitions, matter initialization, CMC solve and the six initialization iterations remain fixed.

### First-cell H residual

| N | Baseline center-to-center trapezoid | Origin half-cell only | Origin + uniform fourth-order intervals |
|---:|---:|---:|---:|
| 160 | +4.0416e-3 | +6.5061e-5 | +3.5598e-6 |
| 240 | +4.0510e-3 | +6.5985e-5 | +1.0699e-6 |
| 320 | +4.0543e-3 | +6.6634e-5 | +4.9875e-7 |

At N=320, the formal-region H maximum (i>=2) is 8.37e-5 baseline, 3.39e-5 origin-only, and 1.60e-7 with uniform fourth-order intervals. On the fixed annulus 2<=r<5, the maximum H residual in the uniform-fourth-order arm is 7.95e-9 at N=320; the corresponding baseline is 1.23e-6. This is strong evidence that the quadrature reconstruction contributes materially to the initial Hamiltonian defect.

However, the first-cell raw Aa RHS in the uniform-fourth-order arm remains about +4.14e-5 at N=320. Because the first cell is at r=dr/2, Aa_rhs/r² rises from 6.21e-4 at N=160 to 1.46e-3 at N=240 and 2.65e-3 at N=320. A regular even trace-free RHS should scale as O(r²), so the raw first-cell tendency should approach zero. The initial H improvement cannot be treated as proof of centre regularity.

The origin-half-cell-only arm illustrates the same distinction: its H0 plateaus near +6.5e-5 to +6.7e-5 at N=160–320, while raw first-cell Aa RHS tends toward roughly -9.5e-5 and its quotient Aa_rhs/r² grows in magnitude. Fixing the origin interval alone is insufficient.

### One-step propagation (origin-half-cell A/B, not the uniform-fourth-order arm)

The separate matched one-step test uses a CFL-scaled step and the exact production RHS. Starting M is roundoff in both arms. After one step:
- N=40: baseline max |M| 1.62e-4 vs origin candidate 9.63e-6; max |C_lambda| 9.86e-6 vs 5.86e-7.
- N=80: baseline max |M| 1.56e-4 vs candidate 2.08e-5; max |C_lambda| 4.75e-6 vs 6.36e-7.
- N=120: baseline max |M| 1.55e-4 vs candidate 2.41e-5; max |C_lambda| 3.15e-6 vs 4.92e-7.

The candidate substantially reduces one-step M and connection-constraint errors in this particular A/B, but the errors remain nonzero. This is a local diagnostic, not evidence of long-time constraint preservation. Do not apply these exact one-step numbers to the later uniform-fourth-order arm; a matched uniform-fourth-order one-step test still needs to be run/admitted.

## B. Independent regular-centre curvature/Aa RHS identity

An independent analytic metric-curvature reference was compared against production Ricci and trace-free Aa RHS on smooth manufactured regular-centre profiles, with N=40/60/80/120. The flat control gives maximum Aa RHS about 4.8e-19 and Ricci/lapse trace-free pieces about 5.9e-16, as expected. For non-flat smooth profiles, fixed-annulus errors show descriptive fourth-order behavior: approximately order 4.02 for the full Aa RHS and Ricci trace-free term, and 4.006 for the lapse-Hessian trace-free term. First-cell quotient errors decrease monotonically with refinement.

This establishes that the independent reference and production curvature pathway agree at the expected rate on these manufactured smooth profiles, including the tested analytic-Lambda and discrete-constraint-Lambda modes. It does **not** prove the evolved production initial profile satisfies the centre series. The production initial profile still has a nonzero first-cell Aa tendency, and the trace-free Ricci cross-check shows the vendor and independently reconstructed metric Ricci can differ in sign at i=0 under the current discrete treatment. This needs exact term-by-term attribution before any code repair.

A code-path discrepancy is also recorded in `docs/CENTRE_SERIES_CURVATURE_RHS_LEAD_2026-10-12.md`: the vendor offers `enforce_algebraic_regularity`, but the pinned `V55ProductionKernel.step()` path audited there does not invoke it, despite older generic status text implying a centre projector is reapplied. Treat this as an explicit audit target, not permission to add a projector.

## C. Frozen-state semi-discrete momentum tangent

The multi-constraint tangent uses six paired t=0 states (baseline/candidate × N=60/80/120), evaluates the exact production semi-discrete RHS, perturbs both directions with epsilons 2e-5 and 1e-5 while re-solving the CMC lapse, and compares the tangent to one ordinary production step. Signed H/M/C_lambda decompositions close to roundoff; this is accounting closure, not an independent proof of covariant stress-energy conservation.

At the first admitted cell:

| N | Baseline Mdot | Candidate Mdot | Candidate one-step ΔM/dt |
|---:|---:|---:|---:|
| 60 | +1.2113e-4 | -8.7553e-5 | -8.6967e-5 |
| 80 | +1.4873e-4 | -1.2354e-4 | -1.2295e-4 |
| 120 | +2.1261e-4 | -1.9063e-4 | -1.8993e-4 |

The initial H solve does not make Mdot vanish. Tangent and one-step rates agree closely, and the candidate flips sign rather than eliminating the tendency. At fixed radii, results are sign-changing/non-monotone near r=2 but approach near-cancellation at r=3–5 on the finest grid. The scalar source-side derivative is dominated by local S current; the main A/B response is in geometry. No single signed RHS term has yet been proved to be the unique faulty equation.

## Interpretation limits

- Quadrature is a demonstrated contributor to the initial H defect, but not a complete repair.
- The manufactured-profile Ricci/Aa test passes its intended convergence/regularity checks, but does not certify the actual evolved production profile.
- The RHS tangent demonstrates a nonzero semi-discrete M tendency on the tested initial states; it does not by itself distinguish a wrong continuum equation from a discrete compatibility error.
- The initial-metric candidate is not admitted as a production repair. Do not claim pointwise centre convergence, covariant source conservation, bounce, turnaround, ejection, or cyclic behavior from these results.
- No long production trajectory should be launched until the matched uniform-fourth-order one-step test and the signed centre-term audit pass.

## Next work — bounded and decisive

1. **Run the matched one-step propagation test on the uniform-fourth-order initializer** at N=40/60/80/120 with the same CFL, source state, and recorder as the baseline and origin-only arms. Record initial and post-step H, M, C_lambda, determinant, Aa/r², Lambda/r and fixed-radius residuals. Keep all arms paired; no long trajectory.
2. **Decompose production Aa RHS into signed components** at i=0,1,2 and fixed r=2,3,4,5: trace-free Ricci, trace-free lapse Hessian, anisotropic stress, K/A nonlinear terms and any explicit regularity terms. Compare against independent radial-series expressions on the exact frozen initial states.
3. Check the K RHS and matter-current continuity independently. Do not infer stress-energy conservation from a residual budget closure.
4. Verify exact parity/stencil identities for a,b,X,alpha,Aa,K,Lambda and the corresponding odd/even derivatives. Check the actual pinned production call graph for any centre regularity enforcement.
5. If and only if one discrete identity fails, make one isolated numerical-layer repair on a diagnostic branch, rerun the exact matched tests, and keep the physics lane untouched.
6. Resume the production campaign only after initial H/M/C_lambda admission, one-step propagation, centre regularity, and balance checks pass at matched resolution. No t=160 campaign before these gates close.

## Repository state

This report belongs on diagnostic branch `diag/discrete-consistent-initial-production-test-2026-10-12`. Do not merge it as a physics repair. `main` remains untouched by these diagnostic runs.
