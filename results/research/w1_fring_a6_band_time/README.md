# Band M against time (stored data only)

Question (Jeremiah, Nov 24 23:24 UTC): does band M's growth peak and then fade? Band R 0.05-0.25, read on the A6 ring arms (#591/#592, dr 0.05, A1 dt 1.5e-3, A2 dt 7.5e-4, to t 3.291) and the dr 0.0125 replays (t 0.015). Script `research/w1_fring_a6/band_time.py`, numbers in `band_time.json`. Cell centres every step (cells 1-4); Aster's 17-radius cubic exactly at the four saves and approximately every step (spline through cells 0-5; within 0.6% of exact at t 3.291).

## A6, dr 0.05
- No single rise-then-fade. Band M rises 20-35x between t 0.015 and t 0.5, dips (A1 lowest t 1.75-2.0, A2 t 0.75-1.0), then climbs again.
- Cell centres: A1 peaks at the last step, 8.18e-10 at t 3.291 (peak/last 1.00); A2 peaks at 9.91e-10 at t 2.555, last 8.67e-10 (1.14).
- Aster's radii (approx.): both peak at t of about 3.0 (A1 1.50e-9, A2 1.83e-9) and fall to 9.25e-10 and 9.93e-10 at t 3.291 (peak/last 1.62 and 1.84). Exact at the saves: A1 2.69e-11, 4.80e-10, 4.47e-10, 9.30e-10; A2 4.83e-12, 8.75e-11, 6.95e-10, 1.00e-9.
- The late fall is in the cubic reading only (it draws on cells 0-1 near R 0.05) and the run ends 0.3 after it, so it cannot be told apart from the earlier dip.

## dr 0.0125, t <= 0.015
No turnover: Aster-input samples rise at every sample both ways and at both dt (17 radii 9.7e-13 to 4.69e-12 and 5.5e-13 to 4.67e-12). Our per-step series peaks 1.5% above its last value at t 0.01425 (dt 0.0001875), a step-level wobble.

## Lumen²'s rounding prediction (centre-component Amendment 3), scored as written
Relayed Nov 24 23:28 UTC, after the time series above was read. If the band is rounding: (a) size grows roughly steadily with step count, with no peak then fall; (b) the sign flips cell to cell; (c) the finer-dt arm is larger. A clear peak then fall, smooth across cells, would contradict rounding. Numbers in `rounding_prediction.json`.
- **(a) Not met.** The growth is not steady. It rises 20-35x by t 0.5, then dips, then rises again. There is no sustained peak then fall: the cell-centre band ends at or near its maximum.
- **(b) Not met at dr 0.05.** Band cells 1-4 share one sign at t 0.015, 0.8235 and 1.6455 in both arms (negative, then positive). At t 3.291 there is one flip (cell 3 negative), the same in both arms. **Met at dr 0.0125, t 0.015:** 9 and 10 sign flips across the 15 gaps between band cells 4-19.
- **(c) Mixed.** The finer/coarser dt ratio (A2/A1) is 0.19 and 0.18 at t 0.015 and 0.8235, so the finer arm is smaller. It is 1.55 at t 1.6455 and 1.06-1.08 at t 3.291, so the finer arm is larger. Both band readings agree. At dr 0.0125, t 0.015, the ratio is 1.00 on the 17 radii and 1.13 on cell centres.
- **The contradicting pattern is not seen either.** At dr 0.05 the band is smooth across cells, but there is no clear peak then fall.
- **Reading.** At dr 0.05 the band does not look like rounding. At dr 0.0125, over t 0.015, its cell-to-cell signs do.
