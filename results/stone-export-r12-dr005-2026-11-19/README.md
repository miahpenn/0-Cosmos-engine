# Stone-export R 12 at 8g, dr 0.1 vs 0.05 (research only)

`N2N_FINE=R12 python research/stone_export_n2n_v1.py results/stone-export-r12-dr005-2026-11-19`; rule fixed in the log first (N vs 2N rule).
In `summary.txt` the first value of each pair is dr 0.1, the second dr 0.05 (the E3 header line is labelled correctly in this version).
**Verdict: f 0 PASS, f 1 FAIL.**
- f 0: mid 66-72% at dr 0.2, 0.1 and 0.05 (three resolutions), lines identical (mid 1.728, near-gap 1.414/1.461, broad 2.01), w_tot 0.083 / 0.098 (1.5R / 3R at dr 0.05).
- f 1: mid 69-71% -> 48-52%, broad 17% -> 30%, broad peak 2.07 -> 2.36; near-gap 10-12% -> 15-19%. Not converged; drift 1.2e-2 (the largest of the runs).
- The summary files of the dr 0.05 R 24 run were regenerated with correct dr labels; numbers unchanged.
