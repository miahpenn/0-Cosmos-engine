# Experimental status — 0star central clock

**Branch:** `0star-central-clock`  
**Production:** `main`  
**PR:** #1, draft

## Completed evidence

- D-dependent interior evolution is reproducible.
- Rephasing by `D_active` produces strong trajectory collapse in D, momentum, and geometry proxies.
- The early-phase long campaign reached coordinate time 40 but accumulated only about 0.78 proper-time units after synchronization because the central lapse collapsed.
- The pre-repair resolution campaign showed an inward-moving Hamiltonian-constraint maximum tied to the first retained radial cell.
- CMC elliptic residuals were small while the full Hamiltonian constraint was not, so the old lapse-collapse result cannot yet be treated as clean physical evidence.

## Repair now under validation

The production kernel reapplies the existing spherical-center algebraic regularity projection during the evolution stages. This is a numerical repair, not a new physical law.

## Required gates

- Smoke/regression suite.
- N=80/120/160 center-regularity resolution audit.
- N=160 repaired Dhalf/Dbase/Ddouble phase campaign.
- CFL/domain audit.
- Only after these pass: extended turnaround/bounce campaign.

## Interpretation rule

A persistent lapse collapse after the repair must be tested against regional/normalized constraint diagnostics and resolution/CFL/domain convergence before being called a physical clock-stalling or singular behavior.

## Central proper-time campaign result

Run #107 on the repaired branch completed its CI job with all 83 repository tests passing, but the numerical campaign is quarantined. All three resolutions encountered genuine radiation admissibility failures near t≈21.2–21.5, accompanied by severe lapse/CMC deformation and growing outer/constraint defects. This is evidence against the central-proper-time CMC representation as the current strong-field continuation method, not evidence of a radiation-model failure.

## Moving-gauge boundary/centre consistency repair

The moving-gauge kernel was audited against the pinned reference implementation. Two representation inconsistencies were isolated and repaired without changing the V5.5 field equations:

- The legacy frozen R0 light-boundary reconstruction is derived for the nonadvective-lapse/algebraic-B characteristic system. The moving gauge instead uses advective 1+log lapse and independently evolved B, so the moving PIRK step no longer applies that boundary reconstruction.
- The vendor centre projector sets B=3/4 Lambda as an algebraic convenience. Because B is independent in the moving gauge, the repaired kernel preserves the evolved B and imposes the same odd spherical centre regularity on B directly.

Both changes are numerical boundary/regularity consistency repairs only. The next moving-gauge campaign will determine whether the strong-field failure persists.

## Moving-gauge repair now under validation

A separate kernel, `engine/moving_pirk_kernel.py`, now uses the pinned reference moving-puncture/1+log PIRK2 ordering for the strong-field representation while retaining the existing V5.5 scalar, matter, and stress-energy equations. The vendor center-regularity projection is retained, but its algebraic `B=3/4 Lambda` assignment is not allowed to overwrite the independently evolved moving-gauge B field.

This kernel is **unvalidated** until the smoke test and controlled strong-field resolution/CFL gates pass. No bounce, horizon crossing, or physical singularity is claimed from the implementation alone.

## True-CMC stage-aware candidate

Because the moving-gauge campaign reached a common late outer-gauge pathology before the strong-field solution could remain clean, an archive-backed true-CMC candidate has been added as an independent representation test. It is stage-aware: CMC is solved at the beginning of the step, on the explicit predictor, and on the completed state. The spatial shift is held at zero for this isolated foliation gate.

This candidate retains the current conservative radiation/DM/baryon state and validates physical matter admissibility on real RK stages. It does not reuse the legacy nonadvective/algebraic-B outer characteristic boundary.

The archived CMC gate is the reference witness for this path: at t=24 and CFL=0.06 it reported H maxima of 1.95e-3, 1.66e-3, and 1.02e-3 for N=40/60/80, with minimum lapse about 6.78e-2, 6.70e-2, and 6.66e-2. Those values are archive evidence, not yet a result of the new repository kernel.

## Campaign launch plumbing repair

The moving-gauge strong-field and resolution workflows retain `workflow_dispatch` for compatibility and also accept an explicit branch-push launch token. They run only when the commit message contains `[run-moving-gauge-gates]`, so ordinary experimental commits cannot launch the expensive campaigns. The campaign jobs now fail CI when `run_campaign` reports a numerical failure instead of allowing a numerical failure to appear as a successful workflow.

## Campaign launch

Moving-gauge smoke/strong-field/resolution gates are now launched from the experimental branch through the explicit commit token `[run-moving-gauge-gates]`. This keeps `main` untouched and prevents ordinary branch pushes from starting expensive campaigns.

## Current next gate

1. Run the moving-gauge smoke/regression path.
2. Run the manual moving-gauge strong-field gate.
3. Run the manual N=40/60/80 resolution audit at controlled CFL.
4. Compare invariant witnesses, constraint localization, radiation admissibility, trapped roots, and Misner-Sharp current against the archived reference.
5. Only then decide whether to extend through the former ~28 lapse-collapse region and reconnect/extend long bidirectional runs.


## Latest moving-gauge diagnostic run

The first controlled moving-gauge campaign reached the numerical campaign step at all requested resolutions, but the jobs failed and discarded their ledgers because artifact upload was conditional on success. The workflow is now repaired to preserve ledgers on numerical failure. The next controlled run will capture the exact first-failure state before any physics-layer repair.

## True-CMC campaign launch
The N=40/60/80 stage-aware true-CMC gate is now authorized to run at Rmax=40, CFL=0.06, final_time=24.0. The campaign is diagnostic/acceptance gating only; no physical interpretation is attached to a successful numerical continuation.

## Archived moving-gauge failure state

The captured moving-puncture runs confirm a common late numerical boundary in both resolutions before the current boundary/centre repair was applied:

- N=40, CFL=0.015: failure at t≈27.18; lapse minimum ≈1.76e-21 at the first retained cell r=0.5; lapse maximum ≈3.42e-9 at r=34.5; Hamiltonian maximum ≈7.34 with outer H L2 ≈2.85; determinant constraint remained at machine precision.
- N=60, CFL=0.010: failure at t≈27.36; lapse minimum ≈7.24e-22 at r=1/3; lapse maximum ≈4.41e-9 at r=35.67; Hamiltonian maximum ≈6.15 with outer H L2 ≈1.97; determinant constraint remained at machine precision.

The dominant late errors are therefore in the outer gauge/geometry region while the determinant identity stays exact. Radiation remains physically admissible in the reported final ledgers. This points to the moving-gauge outer continuation layer, not to the radiation source model or centre determinant repair.


## Controlled repair campaigns launched
The repaired moving-gauge gate and the stage-aware true-CMC gate are now launched from the current branch. Both preserve the V5.5 matter/source model and are acceptance campaigns only.
