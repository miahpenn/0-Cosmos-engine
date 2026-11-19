# Stone-export dr 0.1 vs 0.05 (R 24, research only)

`N2N_FINE=1 python research/stone_export_n2n_v1.py results/stone-export-dr005-2026-11-19`; rule fixed in the log first (same as N vs 2N).
**Verdict: R 24 f 0 FAIL, R 24 f 1 FAIL** (still moving 0.1 -> 0.05): near-gap share 8.5-24%, mid 28-40%, mid peak 1.57-1.73, near-gap peak 1.41-1.55.
**Stable across 0.1 -> 0.05 (read off the table, not part of the rule):**
- near-gap + mid combined (1.35-1.9): f 0 56-59% -> 57-58%; f 1 49-52% -> 48-49%. Their split is what moves.
- broad 35-48% (f 0 35-39 -> 38-39; f 1 44-47 -> 47-48), high 4-6%.
- total w at 3R: f 0 0.105 / 0.131 / 0.124 (dr 0.2 / 0.1 / 0.05), f 1 0.105 / 0.135 / 0.141; at 1.5R 0.091 -> 0.113, 0.110 -> 0.121.
- drift 1.7e-3 / 3.2e-3 at dr 0.05 (closure still above 1e-3, falling with dr). Energy inside 3R by T: f 0 -61 / -58%; f 1 -74 / -83% (0.1 / 0.05).
