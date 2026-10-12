# Discrete-consistent initial metric — matched production A/B result

**Date:** 2026-10-12  
**Workflow:** [38284401548](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38284401548)  
**Branch:** \`diag/discrete-consistent-initial-production-test-2026-10-12\`  
**Artifact:** \`discrete-consistent-initial-production-ab-t2p16\` (ID \`11762503161\`)

## What was tested

The baseline and candidate both used the same strict-admission production kernel, source equations, matter state, CMC gauge, PIRK stepping, outer domain and CFL. The candidate changed only the *initial metric construction*: its opt-in Newton solve varied the radial metric variable \(B\) until the discrete Hamiltonian residual used by the production constraint diagnostic was zero to tolerance. The baseline continued to use the existing cumulative-integral initializer.

Both arms ran at N=40,60,80,120, \(r_{\max}=80\), CFL \(=0.03\), radiation on, through \(t=2.16\). Each requested run reached its horizon and passed the strict execution gate. The test suite and explicit initializer regression passed. The baseline histories/checkpoints from this matched run agree with the archived exact-source t=2.16 profiles to roundoff (maximum norm differences typically \(10^{-14}\) or smaller).

No evolution equation, stress source, physical coefficient, gauge equation, or production default was changed.

## Initial-slice result

All four candidate solves converged in **two Newton iterations**. Maximum discrete Hamiltonian residual at \(t=0\) was:

| Resolution | Baseline max \(|H|\) | Candidate max \(|H|\) | Maximum Jacobian condition number |
|---:|---:|---:|---:|
| N=40 | \(3.77\times10^{-3}\) | \(4.27\times10^{-15}\) | \(2.73\times10^3\) |
| N=60 | \(3.93\times10^{-3}\) | \(3.34\times10^{-15}\) | \(6.15\times10^3\) |
| N=80 | \(3.99\times10^{-3}\) | \(2.31\times10^{-15}\) | \(1.09\times10^4\) |
| N=120 | \(4.03\times10^{-3}\) | \(6.11\times10^{-15}\) | \(2.46\times10^4\) |

The candidate retains \(K\), \(A_a\), all scalar fields and every conservative matter field exactly from the baseline initializer. The initial momentum residual is unchanged at roundoff; \(C_\Lambda\) remains at roundoff because Lambda is reconstructed from the metric; and \(a b^2=1\) remains at machine precision. Thus this is a genuinely constraint-satisfied discrete initial slice for H/M/C—not an H-only tweak that immediately injects an initial momentum error.

## Results at t=2.16

The large full-domain central residuals improve substantially at N=60 and N=80, but **the result is not uniform under refinement**. The finest tested grid develops a larger first-cell H and M residual than its own baseline.

### Full-domain maximum norms

| N | \(\max |H|\), baseline → candidate | \(\max |M|\), baseline → candidate | \(\max |C_\Lambda|\), baseline → candidate |
|---:|---:|---:|---:|
| 40 | \(7.41\!\times10^{-4}\to7.26\!\times10^{-4}\) | \(1.78\!\times10^{-3}\to9.17\!\times10^{-5}\) | \(1.64\!\times10^{-3}\to2.80\!\times10^{-4}\) |
| 60 | \(9.27\!\times10^{-3}\to9.93\!\times10^{-4}\) | \(1.11\!\times10^{-3}\to2.39\!\times10^{-4}\) | \(1.99\!\times10^{-3}\to2.67\!\times10^{-4}\) |
| 80 | \(6.97\!\times10^{-3}\to1.60\!\times10^{-3}\) | \(2.54\!\times10^{-3}\to5.04\!\times10^{-4}\) | \(9.73\!\times10^{-4}\to5.96\!\times10^{-4}\) |
| 120 | \(4.74\!\times10^{-3}\to9.85\!\times10^{-3}\) | \(2.43\!\times10^{-3}\to3.66\!\times10^{-3}\) | \(1.20\!\times10^{-3}\to1.40\!\times10^{-3}\) |

The coarse N=40 constraint peak is itself at the centre, and N=40 contains no cell satisfying both \(i\ge2\) and \(r<5\). Use the admitted-inner/fixed-radius measures below before interpreting any resolution trend.

### Admitted inner region (\(i\ge2,\ r<5\))

| N | \(\max |H|\), baseline → candidate | \(\max |M|\), baseline → candidate | \(\max |C_\Lambda|\), baseline → candidate |
|---:|---:|---:|---:|
| N=40 | No cells in this mask | No cells in this mask | No cells in this mask |
| N=60 | \(7.51\times10^{-4}\to6.19\times10^{-4}\) | \(5.00\times10^{-5}\to3.40\times10^{-5}\) | \(4.74\times10^{-4}\to2.61\times10^{-4}\) |
| N=80 | \(1.07\times10^{-3}\to1.02\times10^{-3}\) | \(3.04\times10^{-4}\to2.27\times10^{-4}\) | \(3.24\times10^{-4}\to2.38\times10^{-4}\) |
| N=120 | \(1.54\times10^{-3}\to1.67\times10^{-3}\) | \(9.36\times10^{-4}\to9.09\times10^{-4}\) | \(2.56\times10^{-4}\to1.84\times10^{-4}\) |

The fixed exterior \(r\ge5\) improves at every resolution. For instance, N=120 exterior H falls from \(5.46\times10^{-6}\) to \(1.50\times10^{-7}\), while N=80 falls from \(8.92\times10^{-6}\) to \(1.09\times10^{-6}\). The discrepancy remains concentrated at the centre/first retained cells.

The proper-time changes at t=2.16 are small (about \(-2.10\times10^{-5}\), \(-5.92\times10^{-6}\), \(-2.42\times10^{-6}\), and \(-5.99\times10^{-7}\) for N=40/60/80/120). Misner–Sharp current RMS also improves at N=40/60/80 but worsens at N=120, so the balance witness is not uniformly repaired either.

## Interpretation

This experiment separates two facts:

1. **The old initializer contributes a genuine discrete H-constraint defect.** The candidate removes it at t=0 without spoiling the initial M or connection constraint.
2. **That initial defect is not the whole production failure.** Starting from the constraint-satisfied slice, the existing evolution regenerates central H/M/C residuals by t=2.16. The residual response is strongly resolution-dependent: N=60/80 improves, whereas at N=120 the central H/M maxima are worse than baseline. That is consistent with a central discrete compatibility layer that an H-only initial solve does not cure.

The candidate's \(t=0\) H closure is built into its Newton solve, so it is a solver check, not independent validation. The short evolution and independent same-source baseline comparison are the meaningful tests. No long campaign was run.

## Decision and next discriminator

**Do not replace the production default initializer, merge this candidate, or authorize t=160 from this evidence.** Keep the solver opt-in.

The next test is a one-step comparison from the exact saved baseline and candidate \(t=0\) states at all four resolutions. Record \(\Delta H/\Delta t\), \(\Delta M/\Delta t\), and \(\Delta C_\Lambda/\Delta t\) after one unmodified production step, including the first five cells and fixed physical radii. This will locate how much of the regenerated defect is already in the first discrete update versus emerging over the first 2.16 units. If the constraint-satisfied initial data immediately develops the same signed defect, the next repair must be in evolution/stage regularity or derivative compatibility rather than further initial B fitting.

## Reference context

Constraint propagation studies emphasize that satisfying the constraints initially is necessary but does not guarantee a discrete evolution preserves them; the continuum BSSN constraint subsystem and its zero-speed modes make the evolution/constraint interface important. See the Living Reviews treatment of Einstein initial-boundary value problems and BSSN constraint propagation ([link](https://link.springer.com/article/10.12942/lrr-2012-9)) and Bernuzzi & Hilditch's comparison of BSSNOK and Z4 constraint propagation in spherical symmetry ([Phys. Rev. D 81, 084003 (2010)](https://doi.org/10.1103/PhysRevD.81.084003)). Neither paper validates this particular coupled CMC/PIRK kernel; they support the decision to test constraint propagation directly rather than infer success from initial residual closure.


---

## Continuation addendum — one-step propagation and matched checkpoint time series (2026-10-12)

### Frozen-input one-step propagation

- Workflow: [38285057887](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38285057887)
- Artifact ID: `11763037058`; GitHub ZIP SHA-256: `6e3321690510089c3e677e90cdcc4af3bac621eeb479c6033b8d50d5d13c5c3b`
- Exact source: `6d8b3f2790644591326df57c55908863f994cb0c`
- One unmodified production step was taken from each saved (t=0) baseline/candidate state, at N=40,60,80,120. No production equations or input states were modified by the audit.

The candidate begins with H at roundoff (maximum (2.3\times10^{-15}) to (6.1\times10^{-15})); M and (C_\Lambda) also begin at roundoff. After one normal step, H is already regenerated: the admitted-inner maximum is approximately (5.05\times10^{-7}), (5.14\times10^{-7}), and (5.20\times10^{-7}) for N=60,80,120. At those same resolutions, admitted-inner M reaches (3.48\times10^{-6}), (3.69\times10^{-6}), and (3.80\times10^{-6}), while (C_\Lambda) reaches (1.41\times10^{-7}), (1.12\times10^{-7}), and (7.72\times10^{-8}). (N=40 has no cells in the chosen mask (i\ge2, r<5); its full-domain H after one step is (9.82\times10^{-7}).)

The corresponding candidate constraint-rate maxima increase toward finer resolution: admitted-inner \\(\max|\dot H_{\rm discrete}|\\) rises from (8.47\times10^{-6}) (N=40 formal-stencil region) through (1.26\times10^{-5}) (N=60), (1.71\times10^{-5}) (N=80), to (2.60\times10^{-5}) (N=120); M-rate peaks rise from (6.30\times10^{-5}) (N=40 formal-stencil region) to (8.70\times10^{-5}), (1.23\times10^{-4}), and (1.90\times10^{-4}). This is direct evidence that the existing update regenerates constraints from the discrete-consistent slice on the first step; it does not by itself identify which RHS term is responsible.

### Completed saved-checkpoint time series

- Workflow: [38287772729](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38287772729)
- Artifact ID: `11764187318`; GitHub ZIP SHA-256: `c9981f1bdb2c5970a3736a8c3571148a6107bc25e26b5837da87472f6ab0ee7c`
- Script repair commit: [`77cf2e0c4fcdfa8ff33fb6a594921d7d345f1221`](https://github.com/miahpenn/0-Cosmos-engine/commit/77cf2e0c4fcdfa8ff33fb6a594921d7d345f1221). The prior attempt failed only because the audit paired scheduled labels against raw binary floating-point timestamps; the corrected script pairs on checkpoint labels and keeps the stored timestamps in the record.
- Post-processing only: 80 saved states (10 times from (t=0) through (2.16), both arms, N=40/60/80/120), 40 matched baseline/candidate comparisons; no evolution was performed by this audit.

The candidate’s initial H maximum is roundoff at every resolution. By (t=0.24), admitted-inner H maxima are (1.74\times10^{-5}) (N=60), (3.12\times10^{-5}) (N=80), and (7.01\times10^{-5}) (N=120); M is (1.79\times10^{-5}), (2.41\times10^{-5}), and (3.12\times10^{-5}). Thus the initial H defect is removed, but central H/M/connection residuals regenerate during the first 0.24 coordinate-time units.

At (t=2.16), admitted-inner baseline \(\to\) candidate maximum norms are:

| Resolution | H | M | (C_\Lambda) |
|---:|---:|---:|---:|
| N=60 | (7.51\times10^{-4}\to6.19\times10^{-4}) | (5.00\times10^{-5}\to3.40\times10^{-5}) | (4.74\times10^{-4}\to2.61\times10^{-4}) |
| N=80 | (1.07\times10^{-3}\to1.02\times10^{-3}) | (3.04\times10^{-4}\to2.27\times10^{-4}) | (3.24\times10^{-4}\to2.38\times10^{-4}) |
| N=120 | (1.54\times10^{-3}\to1.67\times10^{-3}) | (9.36\times10^{-4}\to9.09\times10^{-4}) | (2.56\times10^{-4}\to1.84\times10^{-4}) |

N=40 has no points in the admitted-inner (i\ge2, r<5) mask. Fixed-exterior H ((r\ge5)) improves strongly for the candidate at all resolutions (e.g. (5.46\times10^{-6}\to1.50\times10^{-7}) at N=120), while some exterior M/connection comparisons are mixed. The central response is resolution-dependent: N=60/80 improve some admitted-inner norms, whereas N=120 admitted-inner H is worse for the candidate at (t=2.16). The metric solve therefore addresses a real initial discrete H defect but does not establish a constraint-preserving evolution or justify changing the production initializer.

### Interpretation and next bounded test

The combination of (a) first-step regeneration from a discrete-consistent initial slice, (b) the resolution-dependent central (D_rK/r) layer already present in the centre-series audit, and (c) the earlier finding that changing the CMC center row or using a production-matched fourth-order CMC operator does not improve the full H/M/(C_\Lambda) set, narrows the next test to signed production (K)-RHS attribution.

Audit the exact lapse-Laplacian, nonlinear (K/A_a), and remaining signed terms separately on frozen (t=0) and (t\approx2.16) states at N=40/60/80/120; record the CMC elliptic residual and same-cell (D_rK/r) alongside them, plus fixed-radius witnesses. This is diagnostic-only: no field projection, CMC patch, production equation change, long trajectory, or physical interpretation is authorized by these results. The objective is to determine whether the central layer originates in the lapse solve, the (K)-RHS term balance, or their discrete interface before selecting any repair.


---

## Signed production K-RHS attribution — frozen states (2026-10-12)

**Workflow:** [38288404354](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38288404354)  
**Artifact:** discrete-consistent-initial-signed-k-rhs, ID 11764337761  
**Artifact ZIP SHA-256:** b68e0a0679674d2e4942003427ce8bbd0d66a58a041f7830035683cb0bdc521e  
**Source checkpoint artifact:** run [38284401548](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38284401548), exact production source 6d8b3f2790644591326df57c55908863f994cb0c, pinned vendor d6052d605673ce9d82cdc99b0df79855fa2ea215.

This is saved-state post-processing only: 16 frozen checkpoints covering N=40/60/80/120, baseline and discrete-consistent-initial candidate, t=0 and t=2.16. There was no evolution. The signed K RHS decomposition closes to a maximum absolute error of 5.42e-20 or less across the audited records.

### Exact term split

The production K RHS is decomposed into:
- lapse Laplacian contribution, -D²(alpha), from the pinned fourth-order geometry operator;
- nonlinear extrinsic-curvature contribution, alpha * (Aa² + 2 Ab² + K²/3);
- shift advection;
- KO dissipation;
- matter source, exactly as added by the production adapter, 4 pi alpha (rho + pr + 2 pt).

The sum reconstructs the exact adapter K RHS to roundoff. Shift advection and KO dissipation are zero in the reported inner slices for this zero-shift CMC branch.

### Main findings

1. **The large lapse and matter pieces cancel strongly near the centre.** At t=2.16 and the first admitted cell i=2:
   - N=60, r=3.333: lapse term -7.4415e-4, nonlinear A/K +9.4309e-5, matter +6.9525e-4; total K RHS +4.5407e-5 versus recomputed CMC target +6.2224e-5, leaving a signed difference -1.6818e-5.
   - N=80, r=2.5: lapse -8.9081e-4, nonlinear +9.4200e-5, matter +8.4497e-4; total +4.8349e-5 versus target +6.1858e-5, difference -1.3512e-5.
   - N=120, r=1.667: lapse -1.0109e-3, nonlinear +9.4144e-5, matter +9.7063e-4; total +5.3856e-5 versus target +6.1592e-5, difference -7.7358e-6.

2. **The discrepancy falls at a fixed physical radius with refinement, but the central K-gradient layer is not uniform.** At fixed r=3, production K RHS minus the recomputed CMC target is:
   - t=0: +8.305e-5 (N40), +3.168e-5 (N60), +1.325e-5 (N80), +5.967e-6 (N120).
   - t=2.16: -7.215e-5 (N40), -2.718e-5 (N60), -1.125e-5 (N80), -5.116e-6 (N120).

   The sign reverses during evolution; these are numerical/operator residuals, not evidence of a physical turnaround. On the first cell at t=2.16, D_r K/r changes from about -1.00e-3 (N40) to -2.35e-3 (N60), -4.22e-3 (N80), and -9.56e-3 (N120). At the first admitted cell i=2 it is +7.41e-5 (N60), +1.26e-4 (N80), and +2.78e-4 (N120). The near-centre derivative ratio therefore remains strongly resolution-sensitive.

3. **The stored lapse nearly satisfies its own second-order CMC stencil, while the fourth-order production K RHS does not match that target pointwise.** Re-evaluating the CMC target on the saved state (the historical target passed to the final lapse solve is not serialized) gives an interior second-order elliptic residual of about 1.43e-8, 4.35e-9, 1.85e-9, 5.52e-10 at t=0 for N40/60/80/120, and 1.20e-8, 3.67e-9, 1.57e-9, 4.69e-10 at t=2.16. This is a reproducible frozen-state residual, not an exact replay of the historical solver target. The production fourth-order K-RHS mismatch remains around 10^-5 at admitted inner cells.

4. **The discrete-consistent initial metric does not materially change this K-RHS pattern.** Baseline/candidate term values and target residuals are nearly indistinguishable at matched radii. The initializer removes the initial Hamiltonian defect, but the unchanged evolution and its lapse/K operator interface still regenerate central constraints.

### Decision and next gate

This audit confirms a concrete second-order-CMC versus fourth-order-production-K operator mismatch and a strong lapse/matter cancellation in the centre. It does **not** establish that the mismatch alone causes the H/M constraint drift: earlier one-step CMC-operator counterfactuals did not improve the full original H/M/connection constraint set. Do not patch the CMC solver, replace the K RHS stencil, use a field projection, or launch t=160 from this result alone.

The next decision should reconcile this term budget with the parallel manufactured-momentum-operator audit and any newer stage-level evidence. Only a term-level prediction that survives a frozen-input counterfactual with the original H/M/C diagnostics can justify a numerical-layer change.


---

## Frozen centre-series RHS-rate audit — completed (2026-10-12)

**Workflow:** [38289053751](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38289053751)  
**Artifact:** discrete-consistent-initial-centre-series-rhs, ID 11763734669  
**Artifact ZIP SHA-256:** c29969262e7bd79a8dde8ace357945685a5b5a418db3344f5aa08fda94b52021

This diagnostic evaluates the same 16 frozen baseline/candidate checkpoints at N=40/60/80/120 and t=0/2.16. It reproduces the production first-stage sequence: solve CMC lapse on a copy, zero shift and B as in the production stage, compute the actual scalar and geometry RHS, then evaluate rates of the centre-series ratios. It takes no evolution step.

### Main measurements

At t=0, the candidate's discrete H constraint is at roundoff, but the first-stage rate of the even-field ratio D_r K/r is almost unchanged by the initializer. Its first-cell rate grows from +5.107e-4 (N40) to +1.195e-3 (N60), +2.149e-3 (N80) and +4.869e-3 (N120); baseline and candidate agree closely at every resolution.

The first-cell rate of A_a/r² does respond to the B-solved candidate initializer. Baseline values grow from +1.104e-3 (N40), +2.456e-3 (N60), +4.357e-3 (N80) to +9.792e-3 (N120); candidate values are +3.949e-5, +7.431e-5, +1.463e-4 and +3.575e-4. Thus solving the initial H defect substantially changes this one tendency, but it does not remove the resolution-sensitive K-gradient tendency.

At t=2.16, the saved first-cell D_r K/r values themselves are -1.004e-3 (N40), -2.348e-3 (N60), -4.220e-3 (N80) and -9.558e-3 (N120). Their first-stage rates remain similar between baseline and candidate: -4.35e-4, -1.021e-3, -1.839e-3 and -4.172e-3, respectively. The first-cell Lambda/r RHS rate likewise grows with refinement, from about +1.218e-3 at N40 to +2.987e-3, +5.443e-3 and +1.247e-2 at N60/80/120, again almost unchanged by the candidate initializer.

For the first-cell A_a/r² tendency at t=2.16, the candidate is not uniformly better: it changes the N40/N60 values and reduces the N120 rate, but the N80 candidate rate remains +2.276e-3 against baseline +2.826e-3. The corresponding (a/b-1)/r² rate also changes sign/magnitude non-uniformly, especially at N120. These are centre-adjacent values at different physical radii, so they are a diagnostic flag rather than a fixed-radius convergence estimate. The separate fixed-r=3 witnesses remain sign-changing and resolution-sensitive for the K-gradient ratio rate.

### Interpretation

The initialized H defect is real, but it is not the whole evolution problem. The first-stage K-gradient tendency and the subsequent Lambda/r tendency are present at essentially the same scale on the baseline and the discrete-consistent initial slice; the A_a tendency is partly changed by the initializer but does not uniformly converge. Together with the signed K-RHS budget, this localizes the remaining issue to the central coupled geometry / connection / lapse-RHS interface more strongly than to initial H alone. It does not identify a single defective term, and the rates do not prove a physical centre singularity.

---

## Manufactured trace-free radial Ricci centre control — completed (2026-10-12)

**Workflow:** [38289053430](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38289053430)  
**Artifact:** manufactured-tracefree-ricci-center-audit, ID 11764423233  
**Artifact ZIP SHA-256:** 4e5ad8862499b8964025c3a94940b4a3e6bf772a6e1f7cd8f39a325a6a149c4f  
**Source:** production kernel 6d8b3f2790644591326df57c55908863f994cb0c, pinned vendor d6052d605673ce9d82cdc99b0df79855fa2ea215.

This independent no-evolution test compares the vendor radial Ricci and trace-free Ricci operators against analytic derivatives of smooth, regular spherical-metric profiles. It runs N=40/60/80/120 on r_max=40 and separates (a) the continuum analytic Lambda profile from (b) Lambda reconstructed using the production discrete metric-connection formula.

On the fixed physical annulus 3≤r≤6, the measured descriptive error order is 4.014 for continuum Lambda and 4.025 for production-discrete Lambda. That is a negative control against a generic low-order Ricci/stencil failure for smooth profiles.

The centre-adjacent comparison is more discriminating. The exact analytic R_TF/r² in the first cell is about +9.980e-6 across resolutions. With continuum Lambda, the vendor value changes from -4.248e-5 (N40) and -1.547e-5 (N60), through -4.789e-6 (N80), to +3.266e-6 (N120). When Lambda is instead reconstructed by the production-discrete metric-connection formula, the vendor R_TF/r² approaches the analytic value: -5.509e-6 (N40), +5.431e-6 (N60), +7.984e-6 (N80), +9.278e-6 (N120). The discrete-Lambda case has C_Lambda exactly zero by construction; that is a controlled operator comparison, not independent validation or a proposal to project the live state.

### Combined decision and next bounded gate

The smooth manufactured tests do not show a generic Ricci or momentum-operator defect. They do show that the trace-free Ricci centre response is sensitive to whether Lambda matches the same discrete metric-connection relation used by the production constraint. The evolved checkpoints have a much larger and sign-changing centre C_Lambda residual (e.g. baseline t=2.16 first cell: -1.990e-3 at N60, -9.728e-4 at N80, +1.197e-3 at N120), while the frozen centre-series audit finds a large, resolution-dependent Lambda/r RHS tendency.

The next bounded test should therefore freeze each production checkpoint and compare the current Lambda input to a counterfactual Lambda reconstructed from the exact production discrete metric-connection formula, recording the resulting trace-free Ricci and A_a RHS differences at the first cells and fixed radii. If this diagnostic materially changes the A_a tendency, the next step is one matched ordinary production step with the original H/M/C diagnostics retained. This is not authorization to insert a projection or replace Lambda evolution: any such intervention must earn support from a multi-constraint one-step comparison across resolutions. No t=160 campaign or physics change is authorized by these diagnostics.


---

## Frozen evolved-Lambda / metric-connection RHS sensitivity — completed (2026-10-12)

**Workflow:** [38289657512](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38289657512)  
**Artifact:** discrete-consistent-initial-lambda-ricci-sensitivity, ID 11764557943  
**Artifact ZIP SHA-256:** 5cb556cc9a9311fb6ba1a86b1daf398e4c01c9ffeb7f458960ff45e767bb47a2

The diagnostic reads the same 16 frozen baseline/candidate checkpoints. It re-solves the CMC lapse on a copy and compares the live saved Lambda with Lambda reconstructed from the exact production-discrete metric-connection expression. Only the Lambda input to the frozen RHS evaluation is changed; no saved field is altered, no step is advanced, and the original H/M/C_lambda constraints are evaluated on the untouched state.

### Main findings

At t=0, C_lambda is at roundoff and replacing Lambda by the same discrete metric-derived expression changes the trace-free Ricci and A_a RHS only at roundoff. Yet the independent same-slice Ricci-route control described below still disagrees strongly with the production trace-free curvature at the centre. Therefore evolved-Lambda drift cannot explain the entire initial centre discrepancy.

By t=2.16, the live Lambda constraint residual is O(10^-3) in the first few cells and the frozen A_a RHS is sensitive to Lambda replacement. At the first cell, baseline delta(A_a-dot/r²), defined as metric-Lambda RHS minus live-Lambda RHS, is -1.069e-4 (N40), +2.162e-3 (N60), +3.100e-3 (N80), and -7.396e-3 (N120). Candidate-arm values are -4.055e-5, +3.657e-4, +1.318e-3 and -1.008e-2. The sign and magnitude change with resolution, so the replacement is not a consistent repair.

The effect is localized: at the same frozen states the K RHS is unchanged to exact zero under the Lambda-only substitution, while the A_a RHS and Lambda RHS do change. At fixed r=3 the baseline delta(A_a-dot/r²) is +3.17e-5 (N40), -7.06e-5 (N60), -1.74e-5 (N80), +1.16e-5 (N120); the candidate values also change sign across resolutions. This demonstrates a strong centre sensitivity but fails the cross-resolution consistency test required for a correction.

**Decision:** do not project Lambda onto the metric connection or replace the evolved Lambda equation. The control proves that a Lambda mismatch can materially affect the trace-free Ricci / A_a tendency after evolution, but not that enforcing the metric connection would preserve the full H/M/connection system.

---

## Production-initial trace-free Ricci cross-check from parallel branch — completed (2026-10-12)

**Workflow:** [38289275758](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38289275758)  
**Artifact:** production-initial-tracefree-ricci-crosscheck, ID 11764109840  
**Artifact ZIP SHA-256:** d0525ab5494c7230c428bcaba1503aab3f12e7720316aefd4c50bbb604b82718  
**Source:** production initializer / kernel commit 6d8b3f2790644591326df57c55908863f994cb0c.

This is a frozen initial-slice audit only. It calls the pinned initializer and CMC lapse solve unchanged, compares the vendor trace-free Ricci with a separate metric-derived expression, and decomposes the actual A_a RHS. All split-closure checks pass at roundoff. The initial metric connection constraint is zero at every audited resolution.

Despite C_lambda=0, the first-cell production trace-free Ricci divided by r² is about +1.065e-3 (N40), +2.347e-3 (N60), +4.148e-3 (N80), and +9.297e-3 (N120). The independent metric-derived route gives about -5.915e-4, -1.286e-3, -2.262e-3 and -5.052e-3. Correspondingly, the live production A_a RHS at that cell is about +1.10e-3 before division by r² at each resolution, whereas the same-state metric-Ricci counterfactual yields about -5.7e-4 before division. The disagreement therefore exists on the initialized slice when the discrete connection constraint is already satisfied; it is not explained solely by evolved C_lambda error.

At the same time, the smooth manufactured-metric control shows near-fourth-order convergence of the pinned vendor Ricci operator. These findings must not be collapsed into a claim that one curvature formula is universally right: the production initializer’s centre-adjacent metric profiles may fail the discrete regular-centre series assumptions, and the independent metric route is a discriminator rather than an admitted replacement. The centre-series rate audit independently confirms that the production A_a tendency does not scale as r² in the first cell; its A_a-dot/r² witness grows roughly as dr^-2 with refinement, while D_r(K-dot)/r does likewise.

### Revised engineering conclusion

The evidence now separates three effects:
1. **Initial centre compatibility:** the pinned initializer satisfies C_lambda and the discrete H candidate can satisfy H at t=0, but the unmodified production A_a RHS still has a nonzero centre limit instead of an O(r²) tendency.
2. **Curvature-route sensitivity:** the independent metric Ricci route changes that initial A_a tendency substantially, but has not been established as the governing discrete evolution operator.
3. **Subsequent Lambda drift:** by t=2.16, C_lambda grows and a Lambda-only RHS counterfactual strongly changes A_a, but its effect is resolution-dependent and does not justify a projection.

The most discriminating next experiment remains the same: repair only the initialization/centre quadrature candidate’s reporting, compare its first-cell series and first-stage A_a/K/connection tendencies with the pinned baseline at matched physical radii, then carry only a justified numerical candidate through one ordinary bounded production step and re-evaluate H/M/C_lambda. No physics change, damping, clamp, centre projection, or long-horizon run is authorized at this point.


---

## Independent smooth regular-centre Ricci / full A_a RHS identity — completed (2026-10-12)

**Workflow:** [38290099699](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290099699)  
**Artifact:** regular-center-ricci-aa-rhs-identity, ID 11765776185  
**Artifact ZIP SHA-256:** 345f91141bb19920051935b9c5c516bf1ddc798046012a36583426d9e5e34e97

The diagnostic uses smooth analytic even metric and lapse profiles, with exact analytic Ricci and lapse-Hessian terms derived from proper radial distance and areal radius. It runs N=40/60/80/120, in both analytic-Lambda and production-discrete-connection-Lambda variants. No evolution is taken. The flat-space control returns Ricci trace-free, lapse-Hessian trace-free and total A_a RHS maxima at approximately 6e-16, 6e-16 and 5e-19.

On the fixed physical annulus 3≤r≤6, the measured descriptive error orders are:
- trace-free Ricci: 4.017 (analytic Lambda) and 4.031 (discrete-constraint Lambda);
- lapse-Hessian trace-free term: 4.006 in both variants;
- complete A_a RHS: 4.017 (analytic Lambda) and 4.030 (discrete-constraint Lambda).

The first-cell errors also decrease monotonically. For example, absolute error in (A_a RHS)/r² declines from 3.04e-4 to 4.06e-5 in the analytic-Lambda case, and from 1.45e-4 to 1.44e-5 in the discrete-constraint-Lambda case, from N40 to N120. The connection-matched form is more centre-consistent, but both variants approach the analytic answer as resolution increases.

**Interpretation:** this is a negative control against a generic vendor Ricci/parity-operator defect. It does not validate the production-initialized centre profile; that profile is the observed outlier and must be tested against its required regular-centre series.

---

## Regular-origin half-cell quadrature A/B — initial-slice result only

**Successful workflow:** [38289853406](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38289853406)  
**Artifact:** origin-half-cell-quadrature-regularity-ab, ID 11764434292  
**Artifact ZIP SHA-256:** bef576a15cc620803982e7e79c27009018617c3f63a24cb8e90ddb83973d9eca  
**Exact source:** 6d8b3f2790644591326df57c55908863f994cb0c. At this exact source commit, the production initializer does not call the separate centre-regularity helper present in a newer branch. The successful A/B test follows the pinned commit's actual initialize sequence and changes only the first origin half-cell integral.

The candidate replaces the baseline half-cell rule 0.5*dr*source[0] with a three-node polynomial-in-r² quadrature exact for source terms through r⁴. The pinned six-iteration initialization loop and all later centre-to-centre trapezoids remain unchanged. No evolution step is taken.

The first-cell Hamiltonian residual improves materially:
- N40: +3.769e-3 → +3.926e-4;
- N60: +3.932e-3 → +1.314e-4;
- N80: +3.989e-3 → +8.172e-5;
- N120: +4.028e-3 → +6.600e-5.

The raw first-cell production A_a RHS falls from approximately +1.09e-3 across N60/N80/N120 to +2.11e-5 (N40), -3.41e-5 (N60), -6.86e-5 (N80) and -8.65e-5 (N120). The first-cell a/b series ratio also falls by roughly a factor of three. Initial momentum remains roundoff and C_lambda remains zero by construction.

This is promising initial-data evidence, not an admitted fix. At the fixed physical annulus, the H RMS improves for N40/N60 but is slightly worse for N80/N120. The candidate's first-cell A_a RHS remains nonzero and does not yet display the required clear O(r²) scaling. A follow-up edit briefly tried to call a centre-regularity helper from a newer, different source version; that run failed with AttributeError and is invalid. It has been reverted. The result retained here is the earlier successful artifact tied to the exact campaign source commit.

### Next decision

Run one ordinary production step from the pinned baseline and the corrected origin-quadrature candidate at N=40/60/80/120, using the same dt=CFL*dr. Preserve the same H/M/C_lambda, determinant, and fixed-physical-radius witnesses before and after the step. Do not run a full trajectory yet. If the candidate retains its central improvement after the unmodified RHS step without worsening the rest of the constraint set, it earns a bounded longer validation; if the centre residuals immediately regenerate, the initializer is only a partial fix and the evolved source compatibility remains open. No physics change, fitted coefficient, damping, clamp, centre projection, or t=160 campaign is authorized.


---

## Matched ordinary one-step propagation of the origin half-cell quadrature A/B — completed

**Workflow:** [38290654079](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38290654079)  
**Artifact:** origin-quadrature-one-step-propagation, ID 11766655192  
**Artifact ZIP SHA-256:** a511f67db4d9a1a12ce27652c82a200dbb4b876266fa3045bf8f995ddc087987  
**Pinned source:** 6d8b3f2790644591326df57c55908863f994cb0c; vendor d6052d605673ce9d82cdc99b0df79855fa2ea215.

This campaign took one ordinary production kernel.step(dt) from the baseline and origin-quadrature candidate at N=40/60/80/120, with dt=0.03 dr: 0.06, 0.04, 0.03 and 0.02 respectively. Both initializers use the exact pinned CMC-lapse / shift / B sequence. No physical equation, source term, clamp, or projection was changed. All states remained finite and all eight cases reached the one-step target.

### What persists after one production step

The first-cell Hamiltonian defect remains much smaller for the candidate at every resolution:
- N=40: 3.769e-3 -> 3.747e-3 baseline versus 3.926e-4 -> 3.915e-4 candidate.
- N=60: 3.932e-3 -> 3.911e-3 baseline versus 1.314e-4 -> 1.331e-4 candidate.
- N=80: 3.989e-3 -> 3.969e-3 baseline versus 8.172e-5 -> 8.460e-5 candidate.
- N=120: 4.028e-3 -> 4.009e-3 baseline versus 6.600e-5 -> 6.929e-5 candidate.

Thus the initial central H improvement persists through one actual update; the candidate is not immediately reverting to the baseline error. The one-step update does, however, generate nonzero M and C_lambda in both arms.

### Why this is not yet an admitted initializer repair

The error is redistributed across the adjacent cells rather than reduced uniformly:
- At N=60, the first formally admitted cell i=2,r=3.333 has |H| 6.25e-6 -> 5.94e-6 baseline but 7.64e-5 -> 7.58e-5 candidate.
- At N=80, the same witness i=2,r=2.5 has |H| 1.77e-5 -> 1.81e-5 baseline versus 5.44e-5 -> 5.38e-5 candidate.
- At N=120, i=2,r=1.667 improves: 5.13e-5 -> 5.16e-5 baseline versus 2.19e-5 -> 2.13e-5 candidate.

The first-cell M residual magnitude falls by roughly a factor of 6–17 for the candidate, and the first-cell connection residual also falls. But at the first admitted cell, the post-step M residual is somewhat larger for the candidate at N=60/80/120. Across the fixed annulus 2<=r<5, H improves at N=40/60, worsens at N=80 and is nearly equal at N=120. This fails the preregistered decision condition that the full H/M/C_lambda set must improve without a new admitted-inner regression.

**Decision:** do not advance this quadrature candidate to a longer trajectory. Keep it as a bounded initial-data result and continue isolating the centre/first-admitted-cell compatibility issue.

## Signed initial H budget — completed

**Workflow:** [38291098188](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38291098188)  
**Artifact:** origin-quadrature-initial-h-budget, ID 11765137827  
**Artifact ZIP SHA-256:** 92bff9be02569772505794f3f4bad5853bd451c745c32d287db12951a51333c2

The audit evaluates the production geometric Hamiltonian term and 16*pi*rho_total on the same frozen t=0 baseline/candidate profiles. The exact source assembly closes H=H_geometry-16*pi*rho_total to zero residual in stored floating-point arithmetic. The independent decomposition R+(2/3)K²-A_ij A^ij also matches the production geometric term exactly for these states.

The half-cell quadrature change leaves the matter source effectively unchanged: at the centre, the absolute candidate-minus-baseline change in 16*pi*rho_total is only about 5.3e-9 (N=40), 1.8e-9 (N=60), 8.3e-10 (N=80), and 2.5e-10 (N=120). In contrast, the centre geometric term changes from roughly +6.38–6.72e-3 to +2.76–3.00e-3, accounting for essentially the entire H improvement. This is a geometry/initial-integral effect, not a matter-source correction.

The adjacent-cell witness changes sign and location: at N=60 the next cell H changes from -3.05e-4 to +2.05e-4; at N=80 it changes from -2.46e-4 to +2.78e-4; at N=120 it changes from -1.97e-4 to +3.35e-4. At fixed radius r=5, the H change between arms decreases strongly with refinement (about -6.23e-5, +2.87e-6, -1.63e-7, -1.99e-9 from N40 to N120), but that alone does not establish full-profile convergence.

**Interpretation:** the half-cell integral changes the radial metric reconstruction and its discrete curvature, while leaving the density source nearly fixed. It lowers the central H residual but relocates a substantial part of the mismatch into the next cells. The evidence justifies further centre-compatible initial-data analysis, not treating the quadrature alone as a production repair.

---

## Parallel invariant-diagnostic correction: areal-radius parity

**Workflow:** [38283346865](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38283346865)  
**Artifact:** areal-radius-odd-parity-audit, ID 11762750863  
**Artifact ZIP SHA-256:** 88625901862948db55519f5b888f476927ef8acc857e72b1fd18c1043e830743

A separate diagnostic branch found that R=r*sqrt(b)/X, the areal radius, was differentiated with even parity in two invariant-diagnostic paths, although regular-centre R is odd and dR/dr is even. The flat-space control showed the old first-cell derivative was 7/12, rather than 1, creating a false central trapped-surface indicator chi approximately 0.34. The correction uses odd parity and changes exactly cells i=0/1; all chi values at i>=2 are bitwise identical. This repairs trapping/Misner-Sharp diagnostic interpretation only; H/M/C_lambda constraints are unaffected. Old central min-chi approximately 0.356 claims should be treated as diagnostic artifacts, not physical trapped-surface evidence. Keep this independent of the unresolved evolution-constraint defect.


---

## Origin half-cell quadrature refinement through N=320 — completed

**Workflow:** [38291486881](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38291486881)  
**Artifact:** origin-quadrature-initial-refinement, ID 11767325140  
**Artifact ZIP SHA-256:** 42d898b6f195bdd823fba35a240cdb58bb557146c828c87f8b3a8b2a58636638

Frozen t=0 initializers were compared at N=40/60/80/120/160/240/320 using the exact pinned source and vendor. Only the origin half-cell integral changed; the six initialization iterations, later center-to-center trapezoids, physical source definitions, CMC solve, and matter initialization were untouched.

The central H residual is consistently reduced: by N=160/240/320 the baseline first-cell values are approximately +4.042e-3, +4.051e-3, +4.054e-3, while the candidate gives +6.506e-5, +6.599e-5, +6.663e-5. On the formally admitted region i>=2, r<5, the H maximum also falls from 6.69e-5/7.92e-5/8.37e-5 to 1.85e-5/2.96e-5/3.39e-5 for N=160/240/320. The peak moves inward as the grid is refined, so this is improvement, not convergence proof.

The effect is confined to the central numerical layer. On the fixed annulus 2<=r<5, baseline and candidate H RMS are indistinguishable to the displayed precision at N=160 (3.76e-6), N=240 (1.73e-6) and N=320 (9.83e-7). The fixed-radius H difference at r=3 decreases from -1.60e-8 at N=160 to -5.69e-10 at N=240 and -6.90e-11 at N=320. The fixed exterior r>=5 is likewise unchanged at fine resolution.

A crucial regularity witness remains unresolved: the candidate first-cell production Aa RHS divided by r^2 grows from about -1.46e-3 (N160) to -3.38e-3 (N240) and -6.07e-3 (N320). Multiplying back by r^2 shows the unscaled first-cell Aa RHS approaches roughly -9.5e-5 rather than zero. The baseline unscaled value approaches about +1.1e-3. Thus the origin correction greatly reduces but does not eliminate the nonzero centre-limit tendency. A regular evolution should give an even trace-free curvature / Aa tendency that vanishes as O(r^2) at the centre.

**Decision:** the origin half-cell error is a real contributor to the initial centre curvature residual, but the remaining finite first-cell Aa tendency shows that the rest of the radial integration/discrete-operator coupling still needs attention. Do not launch a longer trajectory. The next bounded candidate should test a uniformly fourth-order center-to-center quadrature, derived from the cubic interpolant on each interval and combined with the even-in-r origin half-cell rule; audit positivity, H/M/C_lambda and the Aa/r^2 scaling before any evolution.
