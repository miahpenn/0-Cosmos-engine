# f_ring A6: where ring M sits in r, and its per-step shape (reanalysis of #591)

Stored A6 data only; no time step. `research/w1_fring_a6/m_profile.py` re-evaluates M on each saved state with the runs' own `snapshot()`; it matched the stored nine sample cells exactly at every save of every arm. Output: `m_profile.json`.

The A6 record doesn't define a ring edge. For position, the initial field's extent (static slice, dr 0.05) is: S above half its peak deviation out to r 2.08, above 10% to 4.18, above 1% to 7.03; D to 1.78, 3.43, 5.18. Origin is r 0; the grid runs to r 60.

## Peak |M| in r, ring arms
| save t | lambda .03 (A1) | lambda .015 (A2) |
|---|---|---|
| 0.015 | r 0.025 (cell 0), -3.7e-11 | r 0.025, -5.5e-12 |
| 0.8235 | r 0.025, 5.0e-10 | r 1.075, 9.1e-11 (half-peak span r 0.025-1.53) |
| 1.6455 | r 0.075 (cell 1), 4.3e-10 | r 0.075, 6.6e-10 |
| 3.291 | r 0.025, -1.7e-9 | r 0.025, -1.5e-9 |

Away from the peak, max |M| at t 3.291 (A1/A2): r 0.5-2 5.3e-11/5.5e-11; r 2-5 2.4e-11/1.2e-11; r 5-10 ~6e-12; r 10-30 ~2e-12; r 30-60 3.0e-11/2.5e-11 (a rise toward the open outer edge). Control peaks sit at the outer edge or the origin at 1e-14 to 3e-11.

## Time series
Only M at cells 0-5 (r 0.025-0.275) was kept per step; no full-grid max|M| or norm series. Every save but A2 at t 0.8235 has its peak inside those cells. The runs span one ring period (t 0.015 to 3.291, 0.996 x 3.29), so cycle-to-cycle steps can't be seen.

Max |M| over cells 0-5 at quarter periods (t 0.015, 0.82, 1.65, 2.47, 3.29):
- A1: 3.7e-11, 5.0e-10, 4.3e-10, 8.4e-10, 1.7e-9 (rise, dip near mid-period, rise).
- A2: 5.5e-12, 8.4e-11, 6.6e-10, 1.1e-9, 1.5e-9 (rise at every quarter).

Both arms carry a fast wiggle: 21-step-smoothed maxima recur about every 0.08 in t at both step sizes (A1 0.33, 0.41, 0.49, 0.57, 0.65, 0.73; A2 0.32, 0.40, 0.49, 0.56, 0.65). Each arm also has a stretch with no turning points: A1 0.73-1.02, A2 0.65-1.38.
