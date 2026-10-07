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

## Moving-gauge repair now under validation

A separate kernel, `engine/moving_pirk_kernel.py`, now uses the pinned reference moving-puncture/1+log PIRK2 ordering for the strong-field representation while retaining the existing V5.5 scalar, matter, and stress-energy equations. The vendor center-regularity projection is retained, but its algebraic `B=3/4 Lambda` assignment is not allowed to overwrite the independently evolved moving-gauge B field.

This kernel is **unvalidated** until the smoke test and controlled strong-field resolution/CFL gates pass. No bounce, horizon crossing, or physical singularity is claimed from the implementation alone.

## Campaign launch plumbing repair

The moving-gauge strong-field and resolution workflows retain `workflow_dispatch` for compatibility and also accept an explicit branch-push launch token. They run only when the commit message contains `[run-moving-gauge-gates]`, so ordinary experimental commits cannot launch the expensive campaigns. The campaign jobs now fail CI when `run_campaign` reports a numerical failure instead of allowing a numerical failure to appear as a successful workflow.

## Current next gate

1. Run the moving-gauge smoke/regression path.
2. Run the manual moving-gauge strong-field gate.
3. Run the manual N=40/60/80 resolution audit at controlled CFL.
4. Compare invariant witnesses, constraint localization, radiation admissibility, trapped roots, and Misner-Sharp current against the archived reference.
5. Only then decide whether to extend through the former ~28 lapse-collapse region and reconnect/extend long bidirectional runs.
