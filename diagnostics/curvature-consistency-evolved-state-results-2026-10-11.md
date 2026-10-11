# Evolved-State BSSN Curvature Agreement and Constraint Drift — Diagnostic Results

**Date:** 11 October 2026  
**Repository:** `miahpenn/0-Cosmos-engine`  
**Review branches:** `diag/curvature-consistency-evolved`, cross-checked from `diag/independent-metric-ricci-audit`  
**Primary workflow:** [Run #38093029222 — success](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38093029222)  
**Later reruns:** [#38093222412 — downstream failure](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38093222412), [#38093319567 — downstream failure](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38093319567)

## Executive result

The evolved-state curvature comparison itself passed, along with the complete job in run #38093029222. On identical states from the discrete-consistent candidate at N=40, r_max=40, vendor `geometry_terms()["R"]` and the independently written reference `ricci_terms()` matched exactly at t=0, 1, 2, and 3. Every compared first- and second-derivative array also matched exactly at floating-point output precision. The initial comparison gate passed at tolerance (10^{-12}).

This supports **consistency between the two BSSN Ricci implementations on these states**. It does not validate the separate metric-derived/polar-areal Ricci identity, nor prove that the production strong-D initial-data mismatch is fixed.

The same successful run found that the central Hamiltonian residual on this *different discrete-consistent candidate* grew from numerical roundoff at t=0 to about (6.54	imes10^{-5}) by t=3. This is a second, distinct signal: agreeing curvature implementations do not by themselves guarantee that the evolved state remains constraint-consistent.

## Matched curvature and Hamiltonian checkpoints

| t | Max absolute H across grid | Vendor-reference Ricci difference across grid | Initial operator gate |
|---:|---:|---:|---|
| 0 | (1.974	imes10^{-15}) | exactly 0 in output | PASS |
| 1 | (8.758	imes10^{-6}) | exactly 0 in output | n/a |
| 2 | (2.925	imes10^{-5}) | exactly 0 in output | n/a |
| 3 | (6.539	imes10^{-5}) | exactly 0 in output | n/a |

At all four checkpoints, every logged derivative comparison had zero maximum absolute difference for the listed first derivatives (`a`, `b`, `X`, `Lambda`, `alpha`) and second derivatives (`a`, `b`, `X`, `alpha`). The first five cells of the vendor and reference curvature arrays were identical in the logs.

The discrete-consistent candidate's initial (|H|) is far below the strong-D production initializer's approximately (1.08	imes10^{-3}) central offset. These are different initial-data constructions and settings; do not directly equate their residual histories.

## One-step residual-rate probes

For each sampled evolved state, the script advanced separate copies by (Delta t/Delta r = 0.03, 0.015, 0.0075), then measured central (Delta H).

| Base t | (Delta t/Delta r=0.03) | 0.015 | 0.0075 |
|---:|---:|---:|---:|
| 0 | (-7.8922	imes10^{-9}) | (-1.7913	imes10^{-9}) | (-3.8194	imes10^{-10}) |
| 1 | (-4.3677	imes10^{-7}) | (-2.2029	imes10^{-7}) | (-1.1025	imes10^{-7}) |
| 2 | (-8.3699	imes10^{-7}) | (-4.1312	imes10^{-7}) | (-2.0568	imes10^{-7}) |
| 3 | (-1.2844	imes10^{-6}) | (-6.4192	imes10^{-7}) | (-3.2086	imes10^{-7}) |

At t=1–3, halving the one-step increment approximately halves (Delta H), while the per-unit-time rates remain similar. This is consistent with a nonzero local residual derivative on the evolved candidate state; it does **not** by itself establish whether the full-trajectory half-CFL residual growth is dominated by temporal discretization, spatial truncation, or structural constraint propagation.

The connection-RHS budget algebraic closure was tiny: (1.23	imes10^{-32}) initially and (O(10^{-20})) or lower at later checkpoints. This verifies the tested bookkeeping identity, not a vanishing physical connection constraint.

## Related matched spatial probes

The same successful workflow then ran `momentum_resolution_compare.py` at t=1. At fixed dimensionless probe step (Delta t/Delta r=0.015), the central H-rate measurements were approximately:

- N=40, (Delta r=1): (-1.4529	imes10^{-5})
- N=80, (Delta r=0.5): (-1.4986	imes10^{-6})
- N=160, (Delta r=0.25): (-2.3272	imes10^{-7}) (printed final estimate in the run summary: (-2.8352	imes10^{-7}) using the run's additional sampling/estimate)

The reported N=80/N=40 estimated absolute H-rate ratio was 0.1067; the N=160 estimate is smaller again. This indicates strong spatial-resolution sensitivity in this candidate-state local residual rate, but these measurements should not be transferred directly to the strong-D production initializer. The N=160 summary's distinction between probe-derived and aggregate estimated rate should be retained rather than collapsed into a single exact order estimate.

## What the later red runs mean

Runs #38093222412 and #38093319567 each passed dependency installation, diagnostic compilation, and **“Compare initial and evolved curvature.”** They failed only in the subsequent `momentum_resolution_compare.py` step while attempting to construct the discrete-consistent initial state:

`ConvergenceError: did not converge in 20 iterations; final max |H|=9.529304157828572e-13, tolerance=1e-13`

The final iterate was only about (9.53	imes10^{-13}) in max (|H|), but still exceeded its declared (10^{-13}) convergence tolerance. This is a narrow initialization/convergence-gate failure of that later spatial-refinement step, **not evidence that the curvature comparison failed**. Do not silently relax the tolerance or change the solver; diagnose the iteration history and the exact commit if that subtest is resumed.

## How this relates to the strong-D production/shadow test

The strong-D paired full evolution (#38100254183) used N=160, (r_{max}=80), (D=10^{-4}), production CFL 0.0075, and ran to t=22.5. Its cubic/Gauss shadow greatly reduced the central H offset, but max off-centre (|H|) grew to (4.638	imes10^{-4}) in both trajectories by t=22.5. Lumen correctly notes that the shadow is not the same initial-data construction as the discrete-consistent candidate used here. The two investigations answer different questions and should remain separate.

## Ordered next steps

1. Let half-CFL paired strong-D run #38103369774 finish; it was still running at the latest check. Compare matched actual times around t=4, 8, and 12, central H, off-centre max (|H|), and minimum lapse.
2. Use that comparison to judge temporal-step sensitivity in the strong-D production/shadow pair.
3. Then run spatial refinement of that same production/shadow off-centre growth at fixed CFL, with all other settings and initializer differences controlled.
4. Separately, if needed, revisit the later `momentum_resolution_compare.py` convergence failure by examining why its 20-iteration solution stops at (9.53	imes10^{-13}) instead of (10^{-13}). Don't combine that issue with the strong-D quadrature defect.

No production physics, source equation, gauge, projection, default, or `main` branch was changed by these measurements.
