# The reset-to-lane handoff contract, run as a routine (research only)

Next step after PRs #490–491: the machine derives what it hands over. Stored states (five cells: x_c 2.0 at r_max 16,
24, 36; x_c 2.0 closure (i); x_c 1.5). No new runs; nothing in `engine/` changed.

- `research/handoff_contract_v1.py` → `contract.txt` (pieces, handoff time, lane from the derived handoff)
- `w_test()` in the same file → `w_test.txt`

Marks: **Derived**, **Read**, **Open**, **SET BY HAND**.

## Summary

- **Handoff time:** the routine finds t_h = 13 in all four x_c 2.0 cells, the same under all four separation rulers
  (13.5 under one in two cells). In x_c 1.5 three rulers give 10.5 and the Hann ruler gives 6, inside the release:
  **the time is not yet ruler-safe at x_c 1.5.** The threshold (2% of the S–D energy) is **set by hand**. The planned
  threshold, the rulers' own spread, cannot be used: the held residual stays above it (below).
- **Held from the sea's rest:** subtracting the sea's displaced rest value removes about half of the residual held
  step (raw 1–5% → 0.3–2.6% after the release), not all of it. **Open:** either the local x estimate is too coarse at
  the stored cadence, or PR #490's identification of the residual is incomplete.
- **Period-averaged handoff:** it removes the snapshot jitter where the jitter was large (x_c 1.5: lanes handed over at
  t = 10.5–24 now agree to 0.9% instead of 2.4%). It also exposes a **steady +0.5–1.5% offset** with the lane's dust
  (w = 0) ahead of the cell, shrinking for later handoffs.
- **That offset is the S–D species' small pressure (Read).** With the measured mean w (+0.004 to +0.008) instead of 0,
  the sign flips (−0.1% to −1.2%). w is not constant: it follows S's dilution, as in the EOS note's in-sea
  w = (2/3)g − (3/28)x. **The contract must pass the EOS w(x), not a number and not dust.**
- **D_r (drained D, PR #491) cannot be measured from the stored cells.** The estimator swings between 2% and 12% of
  the S–D energy from one snapshot to the next: output every 0.5 gives about 3 samples per forced-ring period. This
  needs **one short run with output every ~0.1**, the only new run proposed.

## 1. Pieces of the S–D energy (x_c 2.0, r_max 24; other cells alike)

| t | held from the sea's rest | held, raw (PR #490) | D_r estimate | rulers' spread of held |
|---|---|---|---|---|
| 5 | 0.107 | 0.176 | 0.04 | 0.001 |
| 8 | 0.013 | 0.049 | 0.11 | 0.012 |
| 12 | 0.026 | 0.054 | 0.08 | 0.006 |
| 24 | 0.009 | 0.018 | 0.10 | 0.001 |
| 54 | 0.006 | 0.012 | 0.11 | 0.002 |

The D_r column is not reliable (see summary). The held column is computed from the sea's rest:
h(D_slow) − h(d_sea), with d_sea² = 2 − (3/2)x and x = 2⟨s²⟩ over one S period (Derived form).

## 2. Lane from the handoff: ρ_SD against the cell (period-averaged), t = 30 / 40 / 54

| cell | handed over at | period-averaged state | snapshot state |
|---|---|---|---|
| x_c 2.0, r_max 24 | derived t_h = 13 | +1.2% / +1.5% / +1.2% | −0.1% / +0.2% / +0.1% |
| | 24 | +0.5% / +0.9% / +0.7% | −0.2% / +0.2% / +0.03% |
| x_c 1.5 | derived t_h = 10.5 | +1.1% / +1.3% / +1.3% | +0.9% / +1.2% / +1.3% |
| | 18 | +0.3% / +0.6% / +0.6% | −1.6% / −1.2% / −1.1% |

H agrees to 2–7e-4 of H(0) for the period-averaged handoff in all cases. Full table: `contract.txt`.

- In x_c 2.0 the snapshot's small residual is partly luck: it sits on a jitter of about 0.5–0.8% (the parenthesised
  column in `contract.txt`, against the unaveraged cell).

## 3. The contract, updated

| item | definition | status |
|---|---|---|
| **when** | first time after the release peak with held (from the sea's rest) below 2% for one D period | Derived form; threshold **SET BY HAND**; ruler-safe at x_c 2.0, **not at x_c 1.5** |
| held D | h(D_slow) − h(d_sea(x)) | Derived form; residual 0.3–2.6% **Open** |
| state handed over | the lane's eight averages, each averaged over one S period (π) | Read: removes jitter |
| S–D species' w | **the EOS w(x) as S dilutes**, not 0 and not a constant | Read (this note); form Derived (EOS note) |
| D_r | free D ring, w = 0 | **Open**: not measurable at the stored cadence |
| ρ_dm | passed through from the imported fluid | **dust still imported** (flag unchanged) |

## Proposed next

1. One short cell run (x_c 2.0 and 1.5, r_max 24, to t ≈ 30) with output every 0.1. It measures D_r, firms up the
   held residual, and checks the Hann ruler's early t_h at x_c 1.5.
2. Give the lane the S–D species' w(x) from the EOS note, and re-run the handoff ruler. Done in research only.

Unpaid inputs: the 2% threshold; the stored cadence (0.5); the D_r estimator's sinc correction (0.805), which assumes
D's own line sits at √2.
