# Amendment 1 to the handoff contract note: the fine-cadence run and the lane's w(a)

Both items of the note's "Proposed next", done. Research only; nothing in `engine/` changed.

- New short runs: x_c 2.0 and 1.5, r_max 24, N 240, closure (ii), to t = 30, **output every 0.1**
  (`research/lattice_cell_run_v1.py full X 240 24 30 OUT ii cold 0.2 0.1`; scratch, not committed: 2 × 13 MB).
- `research/handoff_contract_v1.py`: `line_energy()` (D_r by demodulation), `eos_terms()`, `lane_w()`, `w_lanes()`.
- Outputs: `fine_cadence.txt`, `dr_demod.txt`, `w_lanes_stored.txt`.

## Corrections to the note

1. **The D_r swings were not caused by the output cadence.** At output every 0.1 the boxcar estimator swings just as
   much (3–14%), with S's period. The cause is the estimator: a boxcar cannot separate D's own line at √2 from the
   lines at 2 and 4 (gaps 0.59 and 2.6). Replaced by demodulation at √2 (below).
2. **The held residual after the release is not a cadence effect either.** Held from the sea's rest is the same at both
   cadences to ≤ 0.003 (e.g. t = 12: 0.026 vs 0.026; t = 24: 0.009 vs 0.009). It stays **Open**: PR #490's in-sea
   identification accounts for about half of it.
3. **The x_c 1.5 Hann ruler's early handoff (t = 6) was a cadence effect.** At output every 0.1 all four rulers agree:
   t_h = 10.2–10.4 (x_c 1.5) and 12.9–13.1 (x_c 2.0). **The handoff time is ruler-safe in both cells** (threshold 2%
   still set by hand).

## D_r, measured (Read)

Demodulation of D_fast at √2 in local proper time, Hann window 10.7 (one beat period against the line at 2):

| t | x_c 2.0: D_r / S–D energy | line at 2 | forced line at 4 | x_c 1.5: D_r | line at 2 | at 4 |
|---|---|---|---|---|---|---|
| 9 | 0.056 | 0.020 | 0.004 | 0.062 | 0.011 | 0.003 |
| 13 | 0.067 | 0.013 | 0.003 | 0.075 | 0.022 | 0.004 |
| 17 | 0.078 | 0.017 | 0.003 | 0.077 | 0.028 | 0.001 |
| 21 | 0.082 | 0.027 | 0.003 | 0.079 | 0.032 | 0.002 |

- **D_r carries 6–8% of the S–D energy after the release, smooth and slowly rising** (about +0.2% per unit time),
  in both cells. The slow rise matches a drain that is on but far from complete, as PR #491's thresholds expect for
  these cells. That PR's numbers alone don't separate drain from what the release left ringing.
- The window resolves √2 from 2 only to its main lobe (cross-leak ≈ 25% in energy). Correcting for it changes D_r by
  ≤ 0.01. A window of 21.5 would resolve it cleanly but needs a record to t ≈ 45.
- D_r is free D ring: w = 0, matter-like. It belongs to the S–D species the lane already carries.

## The lane with w(a) (Read; form Derived)

Lane from the period-averaged handoff, with the S–D species' w given as the EOS note's in-sea form with its terms
diluting: w(a) = (2/3)g_h (a_h/a)² − (3/28) x_h (a_h/a)³. Here g_h and x_h are read at the handoff (energy-weighted
x). The dilution laws are assumed: gradients of free waves go as a^-2, and S's amplitude² goes as a^-3.
ρ_SD against the cell (period-averaged), t = 30 / 40 / 54:

| cell, t_h | w = 0 | w(a) | cell's own w(t) (ruler) |
|---|---|---|---|
| x_c 2.0 r_max 24, 13 | +1.2% / +1.5% / +1.2% | +0.2% / +0.3% / −0.1% | −0.2% / −0.2% / −0.1% |
| x_c 2.0 r_max 24, 18 | +1.0% / +1.3% / +1.1% | +0.2% / +0.2% / −0.2% | −0.2% / −0.1% / −0.0% |
| x_c 2.0 r_max 36, 13 | +1.1% / +1.4% / +1.2% | +0.1% / +0.2% / −0.1% | −0.3% / −0.3% / −0.1% |
| x_c 1.5 r_max 24, 13 | +0.3% / +0.6% / +0.6% | −0.1% / +0.1% / +0.0% | −0.0% / −0.3% / −1.0% |
| x_c 1.5 r_max 24, 18 | +0.3% / +0.6% / +0.6% | −0.1% / +0.1% / +0.0% | −0.0% / −0.3% / −1.0% |

- **The lane with w(a) follows the cell's S–D energy to ≤ 0.3% from t_h to t = 54**, against 0.5–1.5% with dust.
  Handing over at t = 24 drifts to −0.5% by t = 54 (the a^-2 law for g is the likely weak point; its value at the
  handoff is 0.019, read).
- The ruler (the cell's own measured w fed back) is no better than the Derived form at x_c 1.5 (−1% by t = 54).
  So what remains there is not pressure but exchange: S–D energy moving to or from something the lane doesn't carry.
  **Open.**

## The contract, after this amendment

| item | status now |
|---|---|
| when | t_h by the derived criterion; ruler-safe in both cells at fine cadence; **threshold 2% SET BY HAND** |
| held D | from the sea's rest; residual 0.3–2.6% **Open** (half explained) |
| state | period-averaged over one S period |
| S–D species' w | **w(a) = (2/3)g_h (a_h/a)² − (3/28)x_h (a_h/a)³**, g_h and x_h read at handoff; ≤ 0.3% (Read) |
| D_r | 6–8% of the S–D energy, w = 0, inside the species (Read, demodulation) |
| ρ_dm | imported (flag unchanged) |

Unpaid inputs: the 2% threshold; the dilution laws for g and x (assumed); the demodulation window (10.7); the
record length (t = 30) of the fine runs.
