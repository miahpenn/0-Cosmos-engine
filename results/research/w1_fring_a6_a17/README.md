# A17 section 5: zero-step readback of the archived W1 states

Jeremiah's go (06:07 UTC, via the coordinator): "Have Wright do it". This is a read of stored states only: no step taken, no engine import, no engine file changed. Stage 1 of the control tests stays on hold.

Script: `research/w1_fring_a6/a17_readback.py`. Numbers: `readback.json`.

## Inputs

The inputs are our earlier local copies of the W1 spatial replays: dr 0.025 and 0.0125, dt 3.75e-4 (primary) and 1.875e-4 (supplementary), with five saves each (t = 0 to 0.015), as signed native-cell M and H. The SHA-256 of each file is in `readback.json`.

They match Aster's late audit on every count compared:
- t 0.015 native sign changes: 1/7 and 9/15;
- 8-shell signs: `-------+` and `-+--+-+-`.

## Definitions (A17 section 1)

- **Band cells:** native cells with r in [0.05, 0.25].
- **Sign rule:** at least 40% of adjacent pairs change sign (at least 3 of 7, at least 6 of 15).
- **Shell shares:** squares of the r^2-weighted means on the four coarse shells, as shares of their sum.

## Result, dt 3.75e-4 (M)

| step (t) | .025 changes | .025 coarse shares % | .0125 changes | .0125 coarse shares % | band RMS .025 / .0125 |
|---|---|---|---|---|---|
| 0 (0) | 3/7 | 0.6 51 19 29 | 11/15 | 27 52 7 14 | 4.0e-17 / 7.6e-17 |
| 10 (.00375) | 1/7 | 66 22 9 3 | 12/15 | 23 66 1 10 | 1.7e-13 / 5.8e-13 |
| 20 (.0075) | 1/7 | 63 28 6 3 | 9/15 | 36 63 0 1 | 3.3e-13 / 1.2e-12 |
| 30 (.01125) | 1/7 | 60 35 3 2 | 9/15 | 34 66 0 0 | 5.2e-13 / 1.4e-12 |
| 40 (.015) | 1/7 | 48 46 3 4 | 9/15 | 38 62 1 0 | 6.8e-13 / 2.1e-12 |

- **Growth:** the band RMS of M − M(0) grows with step-count exponents 0.97, 1.15, 0.92 at .025 and 1.02, 0.53, 1.23 at .0125. That is linear (systematic), not the square root a random walk would give.
- **M at t = 0:** it is at the rounding floor, about 7,600x below the step-10 value at .0125.
- **H at t = 0:** it alternates on the fine grid (12/15, RMS 4.3e-12) but not on the coarse grid (0/7, RMS 1.7e-12).
- **M against initial H:** M's band signs follow H(0)'s on 88% of cells at .025, and have the opposite sign on 75–88% of cells at .0125.

## Scored against A17's declared signatures (by Wright; Lumen² to confirm)

- **G7 (initial data, fine ≥ 40% at t = 0):** met as written (11/15). It does not discriminate, though. M(0) is rounding noise, and it can't seed an M that is 7,600x larger ten steps later. The fine initial H does carry grid-scale alternation. Reading that as the source is Account, and stays Open.
- **G5 (evolution rounding):** not met. Alternation does not build from zero: fine-grid M is at 12/15 by step 10, then falls to 9/15. The growth is linear, not a square root. The shell rule fails at every checkpoint on both grids (a coarse shell is above 40% every time).
- **Centre source:** not met. Fine-grid changes span pairs 0–13 already at step 10. Coarse-grid M has its single change at the outermost pair throughout.
- **Outcome:** MIXED by the declared rules. The G7 hit is a floor artefact in M. The fine-grid H(0) alternation is the lead.
