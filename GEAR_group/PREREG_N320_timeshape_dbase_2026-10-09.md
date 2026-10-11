# Preregistration: N = 320 resolution run, Dbase, r_max = 80 (Lumen, 2026-10-09)

Written before any run. Only the resolution changes from the admitted N = 160 Dbase case.

## Fixed configuration (unchanged from the admitted Dbase r_max = 80 case)
- Case: Dbase, initial amplitude 1e-10. Radiation ON. Physics, gauge, boundary conditions, output definitions unchanged.
- Domain r_max = 80.0. CFL = 0.0075. Final time 22.5. Checkpoint interval 0.25.
- Cases run: one (Dbase). D0 at N = 320 is NOT part of this preregistration.

## Grid convention (verified from code at 1179e08, vendor r0_operators.py)
- dr = r_max / N. Cell faces at i * dr. Cell centers at (i + 0.5) * dr, i = 0..N-1.
- N = 160: dr = 0.5; first center 0.25; index-148 center 74.25 (matches admitted data).
- N = 320: dr = 0.25; first center 0.125; last center 79.875.
- Step: dt = CFL * dr = 0.001875. Steps to t = 22.5: exactly 12000.
- Expected native ledger rows: 12000. Checkpoints: 90 (0.25 interval, unchanged).

## Acceptance gates (integrity; fixed now)
1. Run completed; no failure; no cycle events, handoffs, turnarounds, re-expansions, trapped roots (as at N = 160).
2. Rows: 12000. Allowed roundoff exception, defined now: 12001 rows is accepted only if the final row is at
   t = 22.5 and the last step size is < 1e-9. Any other count fails.
3. 90 checkpoints at 0.25 intervals.
4. All five required fields finite in every row.
5. Artifact ZIP, summary, and ledger SHA-256 recorded and pinned.
6. Final tau reported with full precision.

## Convergence criteria (numerical tolerances fixed now, justified below)
- Final tau: |tau(N=320) - tau(N=160)| <= 1.0e-3 absolute (about 5e-5 relative at tau ~ 20.86).
  Justification: the admitted N = 160 -> N = 320 change for the r_max = 80 D1e-4 control was 1.9e-5 relative;
  this tolerance is about 2.5 times that, a deliberately looser bound, fixed in advance.
- Onset time for each of the 8 quantities at f = 0.25: |onset(320) - onset(160)| <= 0.05 time units
  (about 13 steps at N = 320; about 27 steps at N = 160; the sampling of onsets is resolution-limited).
- Lag for each ordering quantity at f = 0.25: |lag(320) - lag(160)| <= 0.05 time units.
- Decision rule (timing structure, from the frozen primary rule): every ordering |lag| <= 1.0 at N = 320, f = 0.25.

## Failure handling
- A failed gate is reported as a failure of THIS configuration. No reruns with changed parameters, no tuning,
  no physics change, no change of marker.
- A failed convergence tolerance is reported as NON-CONVERGED at the stated tolerance. It is not reinterpreted
  after the fact.
- Any diagnosis of a failure is a separate, separately registered run.

## Interpretation limits
- Passing the gates supports resolution stability of the Dbase r_max = 80 result at N = 160 and N = 320.
  It does not establish the physical interpretation, universal convergence, or the time map to cosmic time.
- The marker sensitivity at f = 0.10 remains open and is reported unchanged.

## Launch
- One run, Dbase, via manual dispatch with case = Dbase on a new workflow commit with its own tag or a
  dedicated N = 320 input. Launch only after this preregistration is reviewed and the workflow change is verified
  as workflow-only. One expensive run at a time.
