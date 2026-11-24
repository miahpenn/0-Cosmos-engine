# W1 A_a rate term split: stored replay arrays (#589)

Written unchanged by `python -m research.w1_aa_split.run <dir> --divisor N` during the #589 run (Nov 24 2026, 06:46-06:55 UTC; Python 3.12.13 / NumPy 2.3.5 / SciPy 1.17.0). Input to that PR's `split.json`. No regrouping has been computed from them.

| file | replay | steps | dt |
|---|---|---|---|
| 1.npz | nominal, divisor 1 | 20 | 7.5e-4 |
| 2.npz | half, divisor 2 | 40 | 3.75e-4 |
| 4.npz | quarter, divisor 4 | 80 | 1.875e-4 |

Each file holds, over all 2400 cells, in extended precision (longdouble):
- `tag_<name>`: that tag's accumulated share of the ledger's tracefree tangent at t 0.015, for the 38 rate tags in `research/w1_aa_split/tags.py`, plus `tag_lapse`, `tag_initial` and `tag_remainder`.
- `ledger`: the conformal_ledger tangent (determinant, tracefree) with its own arithmetic.
- `endpoint` (a, b, X at t 0.015), `initial` (a, b at t 0).
- `metadata` (JSON): attempts, end time, source sha256s, closure-gate statistics, non-tag pieces, per-step tag contributions at cells 0-5 for the first and last steps.

Per-step rate values are not stored. A grouping's response is the frozen-background response of the summed tag arrays (`combine` in `run.py`). sha256 in `SHA256SUMS`.
