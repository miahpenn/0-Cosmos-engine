# Amendment 3 to PREREG_N320_timeshape_dbase_2026-10-09 (Lumen, 2026-10-09)
Written BEFORE any N = 320 launch. No run has occurred under this preregistration.

## Checkpoint timestamp tolerance (new acceptance criterion; not present in the original preregistration)
- Checkpoints are written at the nearest completed integration step. Their embedded time t may differ from the
  nominal 0.25 target.
- At N = 320, dt = 0.001875. Multiples of 0.25 fall between steps (0.25 / 0.001875 = 133.33), so an exact match
  is impossible for most targets.
- Accepted tolerance: |t_checkpoint - t_target| <= dt/2 + 1e-9 = 0.000938 + 1e-9, for each of the 90 checkpoints,
  matched in sorted order to targets 0.25, 0.50, ..., 22.50.
- Rationale: the tolerance is derived from the step grid, not chosen to fit an outcome.
- Checkpoint count (90), checkpoint grid (N = 320 cells, first center 0.125, last 79.875, spacing 0.25),
  and finiteness requirements are unchanged.
- The native ledger row gate, onset rule, lag rule, and tolerances are unchanged. Checkpoints are still not used
  for onset or lag values.

## Launch status
Not launched. Launch requires: this amendment committed under GEAR_group/ with its hash recorded; workflow approval
input set explicitly to "approved" by the owner; one run at a time.
