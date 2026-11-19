# Spatial pumping with back-reaction, toy S (Nov 19 2026)

Script `research/pumping_spatial_toyS_v1.py` (= `pumping_spatial_v1.py`, S mass term 4s -> 1s, nothing else), table `summary.txt` (10 s). Same design, readings and rows as results/pumping-spatial-2026-11-18; validation first. Research only. Supersedes the numbers of that run.

## Validation (passed)
At d0 = 1e-8, G_k matches the toy-S single-mode table (results/toy-s-rerun-2026-11-19) within 25% in 5 of 6 rows: 679 vs 678 (H 0.0057, x0 2.0), 1.2e5 vs 1.23e5 (2.67), 40.2 vs 40.1 (H 0.01, 2.0), 848 vs 863 (2.67), 1.95e5 vs 1.98e5 (H 0.003, 2.0). The sixth (H 0.003, x0 2.67: 1.4e8 vs 3.1e9) is excluded: the linear gain times d0 = 1e-8 is already ~30, so the run is not in the linear regime.

## Result (Measured): peak-mode gain G_k against the starting D amplitude d0 (toy S)
| H, x0 | 1e-8 | 1e-5 | 1e-3 | 1e-2 | 1e-1 | S's ring change (x/x_S-only - 1) at 1e-2 / 1e-1 |
|---|---|---|---|---|---|---|
| 0.0057, 2.0 | 679 | 679 | 646 | 162 | 18.9 | -20% / -68% |
| 0.0057, 2.67 | 1.2e5 | 8.2e4 | 1569 | 253 | 35.9 | -92% / -97% |
| 0.01, 2.0 | 40.2 | 40.2 | 40.2 | 39.3 | 15.4 | -0.5% / -24% |
| 0.01, 2.67 | 849 | 849 | 698 | 129 | 17.6 | -15% / -67% |
| 0.003, 2.0 | 2.0e5 | 1.5e5 | 2352 | 308 | 31.2 | -82% / -95% |

- **Saturation sets in far earlier than with the core's S.** The gain departs from the linear value once D's final rms (d0 times G_rms) reaches roughly 0.1-1 (G_rms is in summary.txt): at production H, x0 = 2.0 the gain is linear to d0 ~ 1e-3 (5% off) and 76% down at 1e-2; at x0 = 2.67 it is already 32% down at 1e-5 and 98% down at 1e-3. With the core's S the same rows stayed within 1-12% to 1e-2.
- **Back-reaction is large:** S's ring energy at t_exit is cut by 20-97% (production H) once saturation sets in; at x0 = 2.67 it is nearly all drained (-92% at d0 = 1e-2). The pump can convert most of S's ring into D.
- **G = 10 survives amplitude in every row** (at d0 = 0.1 the gain is still 15-36); H for G = 10 does not move to a lower H across these rows at production-order H.
- **Attribution:** back-reaction on S (S changes by >1% wherever the gain moves by >40%); no frozen-S lattice control, so detuning is not separated.

## Limits
Planar not spherical, flat noise k <= 4, one seed, constant H, no gravity, uniform S; the row H 0.003 / 2.67 is not validated. The number to read is the saturation scale (final D amplitude ~ 0.1-1), not the gain at a fixed d0.
