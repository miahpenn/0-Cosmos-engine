# Build Status — 2026-10-12

## Objective and current classification

Objective: get the unified V5.5 production machine running through a controlled, synchronized multi-resolution campaign without changing locked physics before numerical admission.

The canonical state graph is implemented. The first synchronized production campaign at N=40,60,80, r_max=80, CFL=0.03, radiation enabled, through t=12 has completed, and the fail-closed completion gate passed. This is **horizon completion, not physical/numerical admission**.

Current classification:

**IMPLEMENTED / 83-TESTED / SHORT-SMOKE-TESTED / t=12 MULTI-RESOLUTION RUNS COMPLETE / CONSTRAINT AND BALANCE ADMISSION STILL OPEN**

## Locked physics and scope

- Corrected archive COSMOS background equations.
- Local GEAR-03 S/D sector with GEAR-63/65 barrier and corrected COSMOS potential normalization V0=8.242522415500654e-5.
- Conservative spherical dark matter, baryons and radiation, with the locked beta DM/scalar exchange.
- One total Einstein stress-energy assembly and pinned reference-metric BSSN/PIRK numerical geometry.
- No fitted coefficient, imposed clock conversion, D-to-matter law, shell EOS, phenomenological interface source, artificial bounce/reset/branch flip, lapse floor, ejection threshold or physical stop.
- Bianchi-I shear remains outside the exact spherical source. Diagnostic-only forced CMC histories are not part of production physics.

## Engineering admission fix and checks

Draft PR #3: https://github.com/miahpenn/0-Cosmos-engine/pull/3

Branch: fix/strict-campaign-admission-2026-10-12

- The campaign runner now writes an admission result and raises if requested resolutions are missing, a resolution did not complete, or a trajectory fails to reach the requested horizon.
- Available campaign outputs are uploaded even when the run fails.
- The runner checks numeric geometry, scalar and matter state fields and stored diagnostic values for non-finite entries.
- The initial t=0 slice is now in the history, and an endpoint discrepancy within floating-point roundoff is recorded and canonicalized without a tiny extra evolution step.
- The strict smoke workflow passed. At the last smoke revision: full repository suite 83 passed; N=40, r_max=40, CFL=0.03, t=1 finished and produced its ledger. Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38207394186

## Synchronized production campaign through t=12

Baseline run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38207515326

- N=40: 200 steps, tau_final=11.9559993814.
- N=60: 300 effective steps after endpoint correction, tau_final=11.9284259877.
- N=80: 400 effective steps after endpoint correction, tau_final=11.9184114053.
- All requested resolutions reached t=12 and the state/diagnostic values were finite.
- Endpoint canonicalization in the rerun corrected only +7.9936057773e-14 for N=60 and +6.3948846218e-14 for N=80; no additional physical evolution step was inserted.
- Post-correction Misner-Sharp current residual max/RMS: N=40 2.3899e-3 / 1.2621e-3; N=60 5.5462e-4 / 2.3152e-4; N=80 6.0062e-4 / 2.1438e-4. The enormous old N=60/N=80 maxima were endpoint-sliver artifacts. The remaining current balance is smaller and improves strongly in RMS from N=40, but its resolution behavior is not yet a completed balance proof.

### Constraint-location finding

The earlier norm called hamiltonian_max excluded the first two *cells*, not the same physical region at each resolution. After correcting the recorded peak location, that maximum sat exactly on the first included cell: r=5.0 for N=40, r=3.3333 for N=60, r=2.5 for N=80. This means the former cross-resolution comparison used a moving near-center cutoff and is not a valid convergence norm.

At t=12, the full-domain Hamiltonian residual maxima were:
- N=40: 4.0417e-2 at r=1.0.
- N=60: 7.9698e-2 at r=0.6667.
- N=80: 1.3494e-1 at r=0.5.

Those maxima lie in the central cells and grow as the grid is refined; this remains a real unresolved numerical/regularity concern. It must not be hidden by an exterior-only norm. The new diagnostics compare like-for-like exterior norms over the fixed physical domain r >= 5.0, with separate inner (r < 5.0) and full-domain maxima and matching peak locations. These fixed-mask metrics await the next t=12 rerun.

The determinant constraint remains small at the reported end slices (roughly 2.8e-9, 4.4e-10, and 7.7e-10 for N=40,60,80), and CMC residuals are small, but that does not compensate for the central Hamiltonian growth or establish full constraint satisfaction.

## CMC target counterfactual and causal limits

R80 with jointly prescribed archived R160 outer-lapse boundary history and full-domain CMC Kdot target completed through t=12:
https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38206154450

It reduced the final proper-time discrepancy relative to archived R160 from about +1.1743e-2 to +6.30e-7 and maximum saved-slice discrepancy from about 1.5702e-2 to 1.2192e-5. This is a **joint gauge/boundary counterfactual**, not a self-consistent R160 solution on R80. It does not isolate Kdot.

The separate frozen-state Kdot target swap (https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38205934746) did not recover the R160 central lapse at t=12; its residual difference was about -3.793e-3. Therefore the near-match in the evolved run cannot be assigned to Kdot alone. If needed, separately evolved counterfactuals must vary outer-lapse boundary and Kdot independently.

The interpolation audit passed 24 held-out times with maximum absolute Kdot interpolation error 5.29e-8; the archived outer-lapse history has a larger reported held-out absolute interpolation discrepancy around 1.5e-4 over half-time-unit spacing.

## N=640 signed-pole audit

Run https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38183314088 was cancelled at the 240-minute limit near t=8.375, short of t=12. Its partial artifact uploaded. Early sign crossings disappeared by about t=5.125; a later residual ridge grew afterward. Term reconstruction closes to roundoff, but this is not independent constraint satisfaction.

## Next gates

1. Complete the fixed-radius-mask diagnostic rerun at N=40,60,80 through t=12 and review inner, fixed-exterior and full-domain constraint maxima with correct peak locations.
2. Independently audit the Misner-Sharp current definition and its flux/work terms at the same fixed-coordinate worldtube; don't infer physical ejection from residuals alone.
3. Check matter and radiation boundary flux direction and local conservation terms from raw signed terms, not ratios whose denominators may be near zero.
4. Preserve frozen N=320 reconstruction tolerance 1e-8 and report host/runtime variance; do not relax thresholds after seeing results.
5. Do not launch the t=160 campaign until the t=12 fixed-domain constraint behavior and balance residual are reviewed.
6. After the code and diagnostics are reviewed, enable same-workflow concurrency control so a sequence of pushes cannot launch overlapping expensive campaigns.

A passing unit suite or a completed run does not itself prove physical turnaround, ejection, or multi-cycle behavior. Numerical/gauge failures remain numerical/gauge results unless invariant, independent and resolution-consistent evidence supports a physical interpretation.
