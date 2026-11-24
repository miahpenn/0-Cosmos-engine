# Frozen-slice D_r Aa comparison: declarations before any reading

Reanalysis of stored data only; no time step. Asked by the coordinator (Nov 24 21:36 UTC) on Lumen²'s readings R1-R4 in AMENDMENT 1 to the vendor audit (Drive 1sP1TjeIq_Fr_a4C4DN8h2JAWNY_bRjFbsONJ6h0BDm8), taken as written. Engine and vendor as on main; nothing written to either. Team log not used. Script `research/w1_fring_a6/frozen_aa_prime.py`.

## Data
- Primary: the 16 saved A6 states (A1, A2, B1, B2 at t 0.015, 0.8235, 1.6455, 3.291), dr 0.05, 1200 cells (`results/research/w1_fring_a6/arms/`).
- Also (R3, R4 at the finer grids): the 8 t 0.015 split files of #595 (`results/research/w1_fring_spatial_rerun/runs/`). They hold 3Aa/r, not Aa, so Aa is recovered as r * (3Aa/r) / 3, a few ulp from the evolved Aa. The (B) files are the reconstructed slices for R4.

## Definitions
- dM = regularized_aa_prime(grid, Aa) - _d(grid, Aa, 1) per cell (vendor functions). Aa' enters momentum_constraint once, so this is reg-minus-current signed M; checked on the A6 states against M recomputed with momentum_constraint patched (as in spatial_replay), gap reported.
- delta = Aa[0] - r0^2 Q, Q = vendor _quadratic_extrapolate_first(grid, Aa/r^2) (the reconstruction's own formula). Aa_reg = Aa with Aa[0] = r0^2 Q. dM_reg = dM on Aa_reg. kernel k = (4, -14/3, 2) / dr at cells 0, 1, 2, zero beyond.
- Roundoff scale per cell rnd = eps * (|d|(|Aa|) + (4/3) r |A_l| + (2/3) r^2 |d|(|A_l|)), A_l = 1.5 Aa/r^2, |d| = the same five-point stencil with absolute weights on absolute values (even ghosts). "At roundoff" = |x| <= 10 rnd.
- Band: cell centres with 0.05 <= r <= 0.25 (dr 0.05: cells 1-4; dr 0.025: cells 2-9; dr 0.0125: cells 4-19).

## Input checks first
1. Each A6 state, loaded, reproduces its stored production M at cells 0-5 (arms/*.json save rows) bitwise, else STOP.
2. Synthetic discrimination on the dr 0.05 grid: Aa = r^2 and r^4 give dM at roundoff (R1 passes); r^6 gives the predicted +24 A6 dr^4 r; a cell-0 departure gives exactly delta*k (R2 passes); a departure at cell 1 fails R2; dM's linear split dM - dM_reg - delta*k is at roundoff.

## Readings (pass/fail as written)
- R1 (regular part): every cell, |dM_reg| <= 10 rnd. Reported: max |dM_reg|/rnd overall and on the band; the cell of the max.
- R2 (cell-0 departure): resolved if |4 delta/dr| > 10 rnd[0]; else "unresolved". When resolved, pass if dM at cells 0, 1, 2 each equals delta*k within max(1% of |delta*k|, 10 rnd) AND cells >= 3 are at roundoff. Reported: delta, dM[0..5], ratios dM1/dM0 (vs -7/6) and dM2/dM0 (vs 1/2), and the share of the stored current M[0] that (16/3) delta/dr accounts for. Lumen²'s sizes (delta ~4.5e-13, changes ~+3.6e-11, -4.2e-11, +1.8e-11 at dr 0.05) are compared, not gated.
- R3 (band): pass if every band cell with index >= 3 is at roundoff. Also the cubic (4-point Lagrange) interpolation of dM at R = 0.05 against 2 delta/dr times cell 2's weight, reported.
- R4 (reconstructed slices, the (B) files): pass if every cell is at roundoff. Also: the A6 states after replacing Aa[0] by r0^2 Q (that is dM_reg, so R4 there equals R1).
- Ring and control reported separately; the control is expected at roundoff throughout and nothing is chosen from it.
