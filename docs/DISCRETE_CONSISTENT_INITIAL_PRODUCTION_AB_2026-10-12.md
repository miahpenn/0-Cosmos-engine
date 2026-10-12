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
