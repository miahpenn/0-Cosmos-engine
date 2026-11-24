# Frozen-slice D_r Aa comparison: result

Readings R1-R4 of AMENDMENT 1 to the vendor audit (Drive 1sP1TjeIq_Fr_a4C4DN8h2JAWNY_bRjFbsONJ6h0BDm8), method fixed in PRERUN.md before reading. Stored data only, no time step. Numbers in `frozen.json`.

## Checks
- Synthetic cases discriminate: Aa = r^2, r^4 give dM at 0.7 and 1.2 rnd; r^6 gives +24 A6 dr^4 r to 7e-8 relative; a cell-0 departure passes R2; a cell-1 departure fails it.
- All 16 A6 states reproduce their stored M at cells 0-5 bitwise.
- M recomputed with momentum_constraint patched minus current M equals dM to 1e-22 (ring); dM - dM_reg - delta*k is at most 1.1 rnd everywhere, so the R2 kernel (4, -14/3, 2)/dr is exact on the data.

## Readings (dr 0.05, ring A1/A2 at the four saves)
- **R1 FAIL at every save.** The regular part is not at roundoff: |dM_reg| reaches 3.5e-12 (t 0.015) to 2.2e-10 (t 3.291) at cells 0-5 and on the band, 1e10-1e13 x rnd, up to 3% of the current Aa' on the band at t 3.291. At cells 0 and 1 it is in ratio 1:3 (proportional to r), the shape of +24 A6 dr^4 r; read that way at A1 t 0.015 it implies A6 of about -2.8e-7 (inferred), the same size as Aa/r^2 itself (-1.85e-7), not the small A6 R1 assumed.
- **R2 FAIL at every save.** delta = -2.2e-13 (A1 t 0.015), +2.5e-12, +1.2e-12, -4.1e-12 later; A2 -2.8e-14 to -2.9e-12. The delta kernel holds exactly, but the regular part adds a comparable signed change at cells 0-2, so dM misses delta*k by 5-32% at cell 0, 12-82% at cell 1, 7-177% at cell 2; ratios dM1/dM0 -0.33 to -3.1 (vs -7/6), dM2/dM0 -0.18 to 2.0 (vs 1/2). (16/3) delta/dr accounts for 21-65% of the current M[0]. Lumen²'s size (delta ~4.5e-13, changes ~+3.6e-11, -4.2e-11, +1.8e-11) compared: at A1 t 0.015 ours is -2.2e-13, giving -1.8e-11, +2.1e-11, -8.9e-12 (half size, opposite sign).
- **R3 FAIL.** Band cells 3-4 change by 3e10-9e12 rnd (the regular part). At dr 0.025 and 0.0125 (#595 files) also FAIL; the cubic value at R = 0.05 is -1.6e-12 vs predicted -1.8e-12 at dr 0.0125 dt 0.0001875, otherwise off in size or sign.
- **R4 FAIL.** On the reconstructed (B) slices delta is 0 to 1e-24 but dM still reaches 4e-13 to 3.8e-12 (1e10-1.6e11 rnd); same on the A6 states with Aa[0] projected.
- Control B1/B2: Aa is itself at roundoff level, dM 1e-14 to 1e-11, reported only.

## Added after reading, not gated
The regular-part change does not fall with refinement as a dr^4 truncation would (16x per halving): band max |dM_reg| at t 0.015 goes from 5.6e-13 (dr 0.025) to 1.4e-12 (dr 0.0125) at dt 0.000375, and 7.0e-13 to 4.0e-12 at dt 0.0001875. Open.
