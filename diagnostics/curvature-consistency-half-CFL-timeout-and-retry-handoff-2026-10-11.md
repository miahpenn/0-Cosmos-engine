# Strong-D curvature audit: half-CFL timeout and retry handoff

**Date:** 2026-10-11 UTC  
**Branch:** `diag/independent-metric-ricci-audit`  
**Scope:** diagnostic workflow/reporting only. No production physics, evolution equations, projection, gauge, boundary condition, or defaults changed.

## 1. First half-CFL attempt

Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38103369774  
Triggered commit: `f8be7b7ba0ad3ae2b1ab8dacc8c8caeda687d431`

Observed result:
- GitHub status: `completed`, conclusion: `cancelled`.
- Setup, checkout/submodule, Python, dependency installation, and compilation succeeded.
- The paired numerical evolution step was cancelled; the upload step failed because the expected final `report.json` did not exist.
- No artifact was produced, so there is **no valid half-CFL numerical comparison** from this attempt.
- The run began at 2026-10-11 01:53:57 UTC and ended at 02:54:17 UTC, almost exactly the workflow's configured 60-minute job limit. That strongly suggests the timeout is the cause, but the cancellation logs were not retrievable, so this is recorded as an evidence-based inference rather than a confirmed runner message.

The old script only wrote its report after both trajectories completed. Consequently, samples from a partially executed trajectory would be lost when the job was terminated.

## 2. Retry and diagnostic durability repair

Retry workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38106791775  
Retry commit: `ba118fff7e8a3a1353506019e928fb9d58d79186`

Changes committed on the audit branch:
- Raised the workflow job timeout from 60 to 240 minutes.
- Added per-mode progress checkpoints written after the initial sample and each requested checkpoint, so already-measured samples are persisted before the whole paired campaign ends.
- Added a heartbeat every 128 accepted steps.
- Changed artifact upload to include the run directory, containing each per-mode progress JSON and the final paired report if produced.
- Kept the report upload step unconditional; missing files warn rather than masking the true numerical/run status with an upload error.

The new workflow had started when this note was authored. Read the workflow's live status before treating it as complete; no numeric conclusion is assumed here.

## 3. Scientific question being tested

Reference successful baseline, production versus regular-F shadow, original CFL:
https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38100254183

At the original CFL, the central Hamiltonian residual in the production trajectory remains about (1.1\times10^{-3}), while the cubic/Gauss regular-F shadow starts near (-3.0\times10^{-7}) and is about (-8.0\times10^{-6}) at (t=22.5). However, the off-centre maximum residual grows in both trajectories and becomes nearly identical by (t\approx12). This already separates a large initializer/quadrature offset from a later evolution residual issue.

The half-CFL campaign asks whether changing the time step materially changes that later off-centre growth. Only after that answer is available should the next expensive campaign start.

## 4. Lumen's tightened domain-test requirements

Lumen's additions are adopted:
1. The first-cell lapse α(0) must be sampled together with the residual. Since the lapse is obtained from a domain-wide elliptic solve, it is not a passive control variable when (r_{max}) changes.
2. Record the maximum-|H| radius and the **full radial H profile** every 0.25 time units. This distinguishes a propagating peak from the argmax jumping between two separate peaks.
3. Save enough metadata to compare production and regular-F shadow trajectories at both domain sizes with identical physical parameters and sampling.
4. Keep Δr = 0.5 fixed for the domain comparison: (N=160, r_{max}=80) versus (N=320, r_{max}=160). Compare the feature's absolute radius, its trajectory, full H profile, onset time, central lapse, minimum lapse, and numerical health. Both domains must use the same code path apart from the explicit domain/resolution parameters and the already-defined shadow initializer.
5. A later spacing refinement is a separate third campaign, after the domain-size result: halve Δr and compare the front's radius-time relation. Domain-independence alone cannot distinguish a physical characteristic from a numerical feature propagating at a fixed grid spacing.

## 5. Guardrails and interpretation

- No physics patch and no production admission.
- No assumed shell, characteristic, or boundary explanation before data are compared.
- Do not interpret only the argmax. Preserve the full profile and identify whether the peak continuously moves or switches between local maxima.
- A stable residual trajectory with a domain-dependent central lapse would implicate coupling through the elliptic lapse solve; it would not, by itself, prove that mechanism.
- The archived GEAR finite-radius account is a hypothesis to test. The currently observed maximum around (r\approx10) for (r_{max}=80) is far inside the outer edge, weakening a simple edge-tracking explanation but not resolving it.

## 6. Lumen's production-arm reproduction record

Lumen reports a local production-arm rerun at (N=160, r_{max}=80, D=10^{-4}), CFL (=0.0075), to (t=12), matching the prior hourly-sampled trajectory to about 0.03%. Reported maximum-residual location:
- (t=0): cell 3, (r=1.75), |H| approximately (1.44\times10^{-5}).
- (t=1) through (8): cell 2, (r=1.25), |H| on the order of (10^{-4}).
- (t=9): cell 14, (r=7.25).
- (t=10): cell 16, (r=8.25).
- (t=11): cell 17, (r=8.75).
- (t=12): cell 19, (r=9.75), |H| approximately (2.83\times10^{-4}).
- The first-cell lapse falls sharply, to about 0.017 by (t=12).

These are **Lumen-reported local reproduction observations**, not yet a committed, independently retrievable artifact in the repository. They should be retained with Lumen's script and log, clearly labelled as reproduction evidence and not as the official half-CFL or domain campaign. The hourly samples are insufficient to decide whether the (t=8\rightarrow9) argmax jump is a moving front or a switch between peaks.

## 7. Next action

Wait for retry #38106791775 to finish and inspect its checkpoint artifact even if final status is not success. Compare production and shadow at matched times, including off-centre max-|H|, its radius, α(0), and lapse/geometry health. Do **not** launch the domain-size experiment until the half-CFL campaign's status and partial/final records have been assessed.
