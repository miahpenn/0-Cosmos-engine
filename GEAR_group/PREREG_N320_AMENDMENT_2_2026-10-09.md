# Amendment 2 to PREREG_N320_timeshape_dbase_2026-10-09 (Lumen, 2026-10-09)
Written BEFORE any N = 320 launch. No run has occurred under this preregistration.

## 1. Row definition (verified against the admitted r_max = 80 native ledger)
- The native ledger has ONE row PER STEP. There is no t = 0 row. Verified: the N = 160 Dbase and D0 ledgers
  each have 6000 rows, first row at t = dt = 0.00375, last row at t = 22.5.
- N = 320: dt = 0.001875. Expected native rows: 12000. First row at t = 0.001875. Last row at t = 22.5.
- Checkpoints are separate from rows: 90 checkpoints at the 0.25 interval, unchanged.
- Acceptance: 12000 rows. The roundoff exception in the original preregistration (12001 rows with the final step < 1e-9)
  is retained as written.

## 2. Origin of the final-tau tolerance (stated accurately)
- Tolerance: |tau(N=320) - tau(N=160)| <= 1.0e-3 absolute, at tau ~ 20.86 (about 4.8e-5 relative).
- This tolerance is a CHOSEN acceptance threshold. It is not derived from the convergence order.
- Rationale: the admitted r_max = 80 resolution comparison for the D = 1e-4 repair-ON control (N = 160 -> N = 320)
  changed central tau by 1.93e-5 relative. The chosen tolerance is about 2.5 times that observed change.
- This is a RESOLUTION comparison (same domain, same case family, N changed only). It is not a domain comparison.
  The observed change comes from a different case (D = 1e-4, not Dbase). It supports the threshold as a rationale;
  it does not measure the Dbase resolution change directly.

## 3. Undefined and edge cases (fixed now)
- UNDEFINED onset: a quantity whose maximum excursion |q(t) - q(0)| is <= 1e-12 has an undefined onset.
  Its lag is undefined.
- Gate effect: an UNDEFINED value at N = 320 is a gate FAILURE only if that quantity was DEFINED at N = 160.
  If it was undefined at both resolutions, it is excluded from the comparison and reported as such.
- Edge case: if the onset falls at the final sample (t = 22.5), it is flagged EDGE and reported; the value is kept.
  An edge flag is not, by itself, a failure.
- Always defined otherwise: when max excursion > 1e-12, the first crossing exists by construction and lies inside
  the run. "Onset moves outside the sampled interval" cannot occur under this definition.
- Multiple crossings: the first crossing is used (unchanged rule).

## 4. What did NOT change
Physics, configuration, marker rules (f = 0.25 primary; 0.10 and 0.50 sensitivity), lag definition,
decision rule (every ordering |lag| <= 1.0), failure handling (no reruns, no tuning), interpretation limits.

## 5. Launch status
Not launched. Launch waits for: this amendment reviewed; the preregistration committed to GEAR_group/ with hashes;
workflow diff verified as workflow-only; the workflow's expected row and checkpoint counts matching section 1.

## 6. Sampling source for onset and lag (stated explicitly; method unchanged)
- Onsets and lags are computed from the NATIVE PER-STEP LEDGER history rows (one row per step), NOT from the 90 checkpoints.
- Onset resolution equals the step: dt = 0.001875 at N = 320 and dt = 0.00375 at N = 160. Both are far below the 0.05 tolerance.
- The 90 checkpoints (0.25 interval) are used only for the checkpoint-count gate, not for onset or lag values.
- This matches the analysis applied to the admitted r_max = 80 and r_max = 160 ledgers (D_TIME_SHAPE_ANALYSIS and T9/T10 both used the full history arrays).