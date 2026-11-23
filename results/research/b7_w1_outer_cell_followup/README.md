# B7 W1 outer-cell trace, follow-up (research only)

Go: Jeremiah, 2026-11-23 23:07 UTC ("Go on the follow-up"), on the gates in Lumen²'s Amendment 3
(Drive 1FD7pQV5-gKMMPn5AwQjjg_NGFDXQPI4eUBFw6BI5WGU) as amended by Amendment 4 (1gDXlmh3H9SwrisSV1VojFGT65k8p1dY6NPtX_9cyoLU).
Follows PR #570. No engine or production change. Gates were not changed after the run.

## How it was run

`research/b7_w1_outer_cell_followup.py` drives `research/b7_w1_outer_cell_trace.py` (extended with step count,
G0 window and the per-operation split; the PR #570 defaults are unchanged). Like PR #570 it needs Aster's draft
initializer on `aster/w1-profile` at cfae43c:

    OPENBLAS_NUM_THREADS=1 python -m research.b7_w1_outer_cell_followup results/research/b7_w1_outer_cell_followup

Eight runs: W1 and PC2, 1200 and 2400 cells, dt 0.00075 and 0.000375, 64 steps each (t = 0.048 and 0.024).
Bootstrap: moving blocks of 8 stage values, 4000 resamples, fixed seeds.

## Results against the fixed gates

| Gate | Result |
|---|---|
| G0 | PASS. All eight peaks over the first 0.003 match Aster's to 3 digits. |
| T1 | PASS. Worst stored-vs-rebuilt residual 0.50 of tolerance; the write changes only Aa[−1]. |
| F1 | PASS. 128 stage values per run. |
| F2 | PASS, all four bind (relative 1-sigma 0.08–0.11). 2400/1200 RMS ratio of the arithmetic-only control: W1 2.26 (dt .00075), 1.73 (.000375); PC2 2.01, 1.55. Band [1.3, 3.0]. |
| F4 | PASS. Control at most 0.18F on every run (limit 2F). |

Report-only (F3, outer cell, RMS of per-stage values):
- (a) stencil sum carries essentially all of the arithmetic-only control: 4.1e-15 to 5.0e-15 at 1200 cells,
  6.7e-15 to 1.0e-14 at 2400.
- (b) 2(...)/r about 1e-23, (c) 1/b − 1/a with reciprocals about 2.6e-18, (e) other about 1e-23. All negligible.
- Closure: (a)+(b)+(c)+(e) equals the control to 3e-30.
- (d) departure difference minus control: 6.3e-15 to 8.0e-15 at 1200 cells, 1.3e-14 to 1.5e-14 at 2400.
  It does not fall with dt (W1 1200 rises 6.3 → 8.0e-15; the others change by under 10%). The Amendment 3/4
  prediction that (d) falls with dt is not supported.
- Lumen²'s prediction that (c) is roughly flat in dr holds, but (c) is about 3 orders below (a).
- All four F2 ratios are lower at the smaller dt (2.26 → 1.73, 2.01 → 1.55), as in PR #570. They stay in band.
- Departure-copy drift over 64 steps: 0.6 to 6.1 ulp of 1 in the outer window.

Full per-stage values are in `report.json`.

## Post-hoc: equal-time ratios (report-only, Lumen² Amendment 5, 1CfJRD-KduL5XlbopdgIjdlRiYcYP3C2yz61zfSuTZoA)

Recomputed from `report.json` only, no new simulation. Window t in (0, 0.024]: dt .00075 first 32 steps
(64 stage values), dt .000375 all 64 steps (128 values). Same bootstrap (block 8, 4000 resamples). Values in
`equal_time_ratios.json`.

| Fixture | dt .00075 | dt .000375 |
|---|---|---|
| W1 | 2.07 (rel 1-sigma 0.11) | 1.73 (0.11) |
| PC2 | 2.06 (0.14) | 1.55 (0.09) |

At equal elapsed time the ratio still drops when dt halves, so elapsed time does not explain the trend. Each drop
is about 1.1 sigma (W1) and 1.6 sigma (PC2) on its own, so the trend is not established either.
