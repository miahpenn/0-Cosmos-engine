# Band M against time (stored data only)

Question (Jeremiah, Nov 24 23:24 UTC): does band M's growth peak and then fade? Band R 0.05-0.25, read on the A6 ring arms (#591/#592, dr 0.05, A1 dt 1.5e-3, A2 dt 7.5e-4, to t 3.291) and the dr 0.0125 replays (t 0.015). Script `research/w1_fring_a6/band_time.py`, numbers in `band_time.json`. Cell centres every step (cells 1-4); Aster's 17-radius cubic exactly at the four saves and approximately every step (spline through cells 0-5; within 0.6% of exact at t 3.291).

## A6, dr 0.05
- No single rise-then-fade. Band M jumps by about 30x in the first 0.25-0.5, dips (A1 lowest t 1.75-2.0, A2 t 0.75-1.0), then climbs again.
- Cell centres: A1 peaks at the last step, 8.18e-10 at t 3.291 (peak/last 1.00); A2 peaks at 9.91e-10 at t 2.555, last 8.67e-10 (1.14).
- Aster's radii (approx.): both peak at t of about 3.0 (A1 1.50e-9, A2 1.83e-9) and fall to 9.25e-10 and 9.93e-10 at t 3.291 (peak/last 1.62 and 1.84). Exact at the saves: A1 2.69e-11, 4.80e-10, 4.47e-10, 9.30e-10; A2 4.83e-12, 8.75e-11, 6.95e-10, 1.00e-9.
- The late fall is in the cubic reading only (it draws on cells 0-1 near R 0.05) and the run ends 0.3 after it, so it cannot be told apart from the earlier dip.

## dr 0.0125, t <= 0.015
No turnover: Aster-input samples rise at every sample both ways and at both dt (17 radii 9.7e-13 to 4.69e-12 and 5.5e-13 to 4.67e-12). Our per-step series peaks 1.5% above its last value at t 0.01425 (dt 0.0001875), a step-level wobble.
