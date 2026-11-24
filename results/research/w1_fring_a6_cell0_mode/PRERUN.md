# Cell-0 fast mode in the A6 data: method fixed before any reading

Ask: Lumen²'s derivation (Drive 1eWfUeGlmyh6uD4BjHdETcwtfZLGrQ6JW-8GbEXw03kk). Jeremiah approved it in the Lumen² channel at 20:50 UTC Nov 24 2026. Stored A6 data only (#591, #592): no run and no engine change.

## Vendor functions
`grep -rn --include=*.py` over the repo (results excluded) on main after #592: `enforce_algebraic_regularity` (vendor/bb-palatini-unified-r0/vacuum_bssn.py:152) and `regularized_aa_prime` (:117) are defined and have no callers.

## What is stored
- Per step, ring arms A1 (lambda .03, dt 1.5e-3) and A2 (lambda .015, dt 7.5e-4): H and M at cells 0-5 only.
- Not stored per step: Aa and a/b. They exist only at the four saves (t 0.015, 0.8235, 1.6455, 3.291), and re-evaluating saved states gives no time series. So prediction (1) is read on M only; for Aa and a/b it is reported as not stored.
- #592 had already shown smoothed M maxima near the origin recurring about every 0.08 in t.
- The stated period is 0.0843. Note that 3.685 x dr = 0.184 at dr 0.05, not 0.0843. Both numbers are reported against the measured spectrum. Pass/fail uses 0.0843 as stated.

## Methods (fixed now)
1. **Series.** Signed M at cell c, every accepted step with t >= 0.015, per arm.
   - Detrend by subtracting a centered running mean of width 0.5 in t (333 steps at A1, 667 at A2).
   - Keep only points with the full window (t 0.265 to 3.041).
   - Multiply by a Hann window and FFT, zero-padded to 2^18 points.
   - The dominant fast period is 1/f at the largest |spectrum| over f from 2 to 50 per unit t (periods 0.02-0.5), found on cell 0.
   - (1) passes if the cell-0 period is within 5% of 0.0843 in both A1 and A2.
2. **Amplitude per cell.** |spectrum| of each cell 0-5 at that arm's cell-0 peak frequency (nearest bin), each series treated as in 1.
   - (2) passes if, in both arms, cell 0 is the largest of cells 0-5 and cells 1 and 2 are both below cell 0.
3. **Aa/r^2 near the origin.** Stored Aa/r^2 at cells 0, 1, 2 (r 0.025, 0.075, 0.125), each save, both arms.
   - Smooth even Aa makes Aa/r^2 linear in r^2 near the origin. Extrapolating from cells 1-2 predicts cell 0 = v1 - 0.5 (v2 - v1).
   - Cell 0 "departs" if |v0 - predicted| > |v2 - v1|.
   - (3) passes if it departs at every save in both arms. Per-save results are reported too.
