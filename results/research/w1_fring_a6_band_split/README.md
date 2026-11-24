# Band M split into Aa, K' and j parts: readings (stored data only)

Method and amendments: `PRERUN.md`. Script: `research/w1_fring_a6/band_split.py`. Numbers: `band_split.json`.

**Band data.** No dr 0.0125 band data is in the repo, and Aster's two-grid data (fed74597) isn't either. The two-grid reading is one halving (dr 0.05 to 0.025) at t 0.015, on the stored overlap r 0.05-0.15. The new reruns were not made (Amendment 2).

**Size of the parts.** On the band at t 0.015, the parts are about 1e4 times larger than M and cancel to it:
- Aa part 9.1e-9, j part 8.1e-9, K' part 9.8e-10 (G, integrated |part| dr);
- M itself: 2.6e-13 (dr 0.05) and 3.7e-13 (dr 0.025) at matched dt 7.5e-4.

**Readings.**
1. j's share shrinks with refinement: PASS as declared, by a hair. 0.446074 to 0.446053 at matched dt, the same at matched Courant. The share is effectively unchanged.
2. Aa and K' parts grow about 8x per halving and carry the band growth: FAIL. Under the halving at matched dt, Aa grows 1.00005x, K' 1.0009x and Aa+K' 0.99994x, while band M grows 1.45x.
3. K' carries the outer-edge rise: FAIL. At the cell of max |M| in r 30-60, K' is 0.16-0.54 of the parts' magnitudes; it reaches 0.5 only at lambda .015, t 0.015. The Aa part, mostly Aa', is the largest elsewhere: Aa' 3.6e-11 and K' -7e-12 at lambda .03, t 3.291.

**Aster's (A) target at dr 0.025, dt 3.75e-4.** The archived run (results/research/w1_fring_origin_trace/A-0.025-half.json) reads peak |M| at t 0.015 = 3.55737e-12 against 3.5574e-12: reproduced. Nothing else was reproduced, because no rerun was made.
