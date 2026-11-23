# Frozen-S threshold follow-up (research only)

`python research/frozen_s_threshold_followup_v1.py` (61 runs, T = 10000, ~49 s on 4 cores); numbers in `summary.txt`.
Design: Drive DESIGN-frozen-S-threshold-followup-Nov23-2026 (Lumen², predeclared); pre-run entry in the coordination log 2026-11-23 (Projects). Run on Jeremiah's go, 2026-11-23 01:34 UTC. `run()` from `frozen_s_threshold_scan_v1.py` is imported unchanged; the first scan's files are untouched.

**By the rule as declared:**
- **R 24: PASS.** Bracket (6.0, 6.1), prediction 6.0 (interpolation 6.04, Account). P1-P7 and B1 pass; witness 6.72 reads 0.0062 (0.007 +/- 30%).
- **R 12: MIXED.** Bracket (5.7, 5.8), prediction 5.78 (interpolation 5.79). Everything passes except P7: dr 0.05 moves the 5.8 row from 0.0022 (growth) to 0.0020 (no growth, just under the 2e-3 line), and 5.7 from 0.0008 to unresolved. By the union rule the recorded bracket is (5.7, 5.9); the 5.9 edge assumes 5.9 grows at dr 0.05 (primary 0.0035; not run at dr 0.05). 5.78 is inside either bracket.
- **Item 3 (edge width, R 24): Delta = -0.2, OPPOSITE SIGN; Lumen's ~6.3 prediction FAILS.** First growth row wd 0.3 is 5.9 vs 6.1 at wd 1.0; the sharper edge grows faster at every row (e.g. 6.5: 0.0062 vs 0.0048). dr 0.05 leaves both wd 0.3 bracket rows' class unchanged. Delta sits on the -0.2 boundary of the rule.
- **Item 4 (R 12, 6.5): 4a PASS** (H 0.01 reads no growth); **4b PASS** (H 0.004: 0.0058 vs 0.0114 - 0.006 = 0.0054).

Disclosures: all six P5 rows reached the cap inside their earlier T, so P5 replayed the first scan's fits exactly (a reproduction check, as declared before the run). Measured floors are 0.99-1.33x the analytic estimate (none below 0.5x). Not tested: seed independence at T = 10000, S back-reaction, non-uniform stone, expansion, l > 0. A pass licenses only the linear, frozen, uniform-pump l = 0 threshold for x < 16/15.

**5.9 ruler check (after #565, on Jeremiah's go 2026-11-23 02:08 UTC):** R 12 x_A R 5.9 at dr 0.05 reads 0.0033 (growth, 3.4 e-folds; primary 0.0035). The union bracket (5.7, 5.9) is now measured on both grids. R 12 stays MIXED as declared. `ruler_5p9.txt`, `research/frozen_s_threshold_followup_5p9_v1.py`.
