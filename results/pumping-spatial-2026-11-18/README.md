# Spatial pumping under expansion with back-reaction (Nov 18 2026)

Script `research/pumping_spatial_v1.py` (repo root), table in `summary.txt`. Design, validation and pass/fail were fixed in the coordination log (Nov 18) before running. Planar periodic comoving 1D lattice (N=512, dx 0.25), S uniform at x0, D flat-spectrum noise (k<=4) at rms d0, full nonlinear V, friction 3H on both, k_phys = k/a. Research only.

## Validation (passed)
At d0 = 1e-8 the spectral-peak gain G_k matches the single-mode table within 25% in all six (H, x0) rows: 260 vs 255, 2744 vs 2950, 23.9 vs 24.2, 87.3 vs 87.1, 6.35 vs 6.4, 16.3 vs 16.6; t_exit agrees (28.1 vs 28.1 at 0.0057/2.0).

## Result (Measured)
G_k (peak-mode gain) against the starting D amplitude d0; relative to d0 = 1e-8:
| H, x0 | 1e-8 | 1e-3 | 1e-2 | 1e-1 | back-reaction on S at 1e-1 |
|---|---|---|---|---|---|
| 0.003, 2.0 | 260 | 259 | 163 | 29 | -31% |
| 0.003, 2.67 | 2744 | 1314 | 213 | 20 | -90% |
| 0.0057, 2.0 | 23.9 | 23.9 | 23.7 | 19.4 | -11% |
| 0.0057, 2.67 | 87.3 | 87.1 | 76.6 | 28.2 | -25% |
| 0.01, 2.0 | 6.35 | 6.35 | 6.35 | 6.25 | -1.6% |
| 0.01, 2.67 | 16.3 | 16.3 | 16.3 | 13.2 | -9% |

- **Q2, answered for this setup: yes, the gain depends on D's amplitude, and in one direction.** Large D lowers the gain. The 0.5-e-fold rule (G down by more than ~40%) is met at 0.003/2.0 (d0 = 0.1), 0.003/2.67 (d0 >= 1e-3), 0.0057/2.67 (d0 = 0.1). The rest stay within it.
- **Attribution (the rule in the log):** every row that moves also shows S's amplitude changing by more than 1% (R3 non-negligible), so the cause is back-reaction on S, not detuning alone. No separate frozen-S lattice control was run, so a detuning share is not excluded.
- **Where it bites:** when D's final rms (d0 times G_rms) reaches about 0.1-1. At production H (0.0057) and d0 up to 1e-3 the linear gain stands to within 1%; at d0 = 1e-2 it is within 12%. A stone's pump is therefore the linear one while D is small, and saturates (gain falls 20-70%) once D grows to O(0.1).
- **Edge in H:** G(0.0057, 2.0) stays above 10 up to d0 = 0.1 (19); G(0.01, 2.0) stays below 10 for every d0. The H = 10 edge does not move across the production-H rows for x0 = 2.0 or 2.67; at the pocket H (0.003) the gain falls from 260-2700 to 20-30 at d0 = 0.1 and remains above 10.

## Limits
- Planar, not spherical: a stone's waves also spread geometrically and leave the pocket, which lowers a real stone's D further than here.
- Flat noise on k <= 4, not a stone's actual D content; one noise seed; constant H; no gravity; S uniform (no pocket profile).
- Attribution to back-reaction rests on the 1% rule; detuning is not separately tested.
- Does not touch the contraction side or the clock-stop band.
