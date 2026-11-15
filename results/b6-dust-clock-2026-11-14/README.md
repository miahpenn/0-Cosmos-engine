# B6: the frozen-state test of the local-to-cosmos channel, on the dust clock

HANDOFF Nov 14 2026 (5), Task B6, with the decisions of Nov 14 (Jeremiah) adopted. RESEARCH ONLY. No time stepping,
no new lattice runs, nothing in `engine/` changed; the cosmos lane's own `cosmos_rho_p` is called as it is.
`research/b6_frozen_test_v1.py`, `b6_output.txt`.

Marks: **Derived**, **Read**, **Open**.

## Decisions adopted (Jeremiah, Nov 14)

1. The cosmos lane's t is the proper time of its own fluids, weighted by rest mass. Derived from the lane's own
   conservation laws. When S–D matter becomes a cosmos-lane species it joins the weighting.
2. The pocket's inner redshift (S's clock against the dust's) is carried as its own labelled quantity. Its size is
   not quoted until read on a slicing-independent clock.
3. B6 runs on the dust clock.

Flag carried: the dust is an imported species (core-sheet gap), so t rests on an imported species for now. A label,
to be replaced when S–D matter enters the cosmos lane.

## The test

    cell:     dH/dτ = −(dK/dt)/(3N),  N = ⟨α/W⟩ weighted by dust rest mass,  dK/dt from the cell's lapse solve
    cosmos:   dH/dτ = −H² − (ρ + 3p)/6     Raychaudhuri form, ρ, p from cosmos_rho_p at the proper-volume averages
    channel:  −C/3,  C built from the six parts of B3

Two kinds of frozen state:

- **Fresh initial slices.** Every part of C is built from the fields and the fluid primitives (ρ, v), not from the
  engine's assembled ρ + S. Fully independent.
- **Stored evolved snapshots** (Nov 14 (4) runs). They hold the fields but not the baryon and radiation primitives. So
  the fluid part is (assembled total) − (scalar parts). Independent for the scalars, the shear, the lapse weighting
  and the clock; **not** for the fluids. A fully independent evolved test needs snapshots that store the fluid
  primitives.

## Result

**Fresh slices** (S at its turning point, so the S–D part is its potential alone and C is negative at this instant):

| state | miss without channel | C/F | miss with channel |
|---|---|---|---|
| x_c 2.0, r_max 16 | −9.5% | −0.40 | 2.6e-6 |
| x_c 2.0, r_max 24 | −3.2% | −0.13 | 8e-7 |
| x_c 2.0, r_max 36 | −1.0% | −0.044 | 3e-7 |
| x_c 1.5, r_max 24 | −2.5% | −0.10 | 6e-7 |
| x_c 2.0, r_max 24, hot | −3.1% | −0.13 | 8e-7 |
| no D, r_max 24 | −1.8% | −0.074 | 5e-7 |
| no pocket (uniform), r_max 12 | 1e-13 | 0 | 1e-13 |

Misses are (cosmos − cell)/|cell| in dH/dτ. F = ⟨½ρ_dm + ½ρ_b + ρ_r⟩.

**Stored states, means over windows of 6** (x_c 2.0 cold unless said):

| run | miss without channel, t = 3 / 15 / 39 / 63 | with channel | C/F, same times |
|---|---|---|---|
| r_max 24, closure (ii) | +0.45% / +0.29% / +0.085% / +0.059% | ≤ 9e-7 | +0.036 / +0.086 / +0.14 / +0.28 |
| r_max 24, closure (i) | +0.46% / +0.29% / +0.085% / +0.059% | ≤ 9e-7 | +0.036 / +0.087 / +0.14 / +0.28 |
| r_max 36 | +0.15% / +0.094% / +0.029% / +0.020% | ≤ 7e-7 | +0.011 / +0.028 / +0.047 / +0.096 |
| r_max 16 | +1.2% / +0.84% / +0.22% / +0.14% | ≤ 1.5e-6 | +0.11 / +0.26 / +0.37 / +0.71 |
| x_c 1.5, r_max 24 | +0.46% / +0.22% / +0.059% / +0.036% | ≤ 2e-8 | +0.036 / +0.066 / +0.097 / +0.17 |
| no D, r_max 24 | +0.51% / +0.13% / +0.026% / +0.016% | ≤ 2e-7 | +0.032 / +0.036 / +0.044 / +0.075 |

- **Read: the test passes.** With the derived channel the cosmos lane's Raychaudhuri equation reproduces the cell's
  dH/dτ to 1e-6 or better on every state. Without it the lane misses by the amounts above.
- **Read: sign and size as predicted for the cycle mean.** Positive C: the cosmos lane decelerates too little. At
  r_max = 24, x_c = 2.0, t = 15: C = 0.086 F (predicted 0.08). The instantaneous sign alternates within S's cycle; at
  the turning point it is negative (fresh slices).
- **Read: closure-free on the dust clock.** Closures (i) and (ii) agree to the digits shown.
- **Read: cell size.** The miss scales with the pocket's share: 1.2% / 0.45% / 0.15% at t = 3 for r_max 16 / 24 / 36.
- On the dust clock C/F grows with time (0.036 → 0.28 at r_max 24). The S–D part's ratio to F is steady; the growth
  is the rate difference between the dust clock and the focusing clock (0.04%) multiplying ⟨αQ⟩, which is hundreds of
  times F late in the run. In units of dH/dτ itself the channel shrinks, 0.45% → 0.06%.

## Negative controls

| control | expected | found | |
|---|---|---|---|
| 1. S's action as the averaging weight for the cosmos lane's fields, same channel | fail | misses by −0.03% at t = 3, then −2.2%, −5.9%, −5.8% (r_max 24): worse than no channel at all from t = 15 | as expected |
| 2. the lane's −(ρ + p)/2 form with −C/3 added | fail, by the Friedmann gap H² − ρ/3 | misses by +2.3% / +0.51% / +0.086% / +0.028% (r_max 24) | as expected |
| 3. pocket without D: C falls to 0.003 F and below | pass | **C = +0.03 to +0.07 F; the miss is +0.51% → +0.016%** | **prediction wrong** |

**Control 3, corrected.** The prediction came from a wrong number in the channel note: its no-D row left the restored
hand floor out of φ's potential, which inflated F for that row only. Corrected (Amendment 2 to the channel note): a
pocket of S without D carries a channel about 40% the size of the full one. S alone is a species the cosmos lane does
not carry; D is not needed for the gap. The nonlinear parts (shear, spread of Π) are about ten times smaller without D.

Control 1 at t = 3 nearly passes by accident: early on the J-weighted and volume-weighted averages of φ, Π and the
fluids almost coincide, because nothing has yet moved.

## What this does and does not establish

- **Established (Read):** the derived channel, on the dust clock, in Raychaudhuri form, closes the cosmos lane's dH
  against the cell to 1e-6. The measure and the equation's form matter (controls 1 and 2).
- **Not established:** that a cosmos lane *evolved* with this channel follows the cell. That is the few-step test,
  and it needs the companions in the φ and ρ equations (channel note, B4), not yet evaluated.
- **Not established:** anything at production λ = 0. All states here are one-budget, λ = 5000.
- By construction "with channel" closes whenever C is defined as the difference. The content of the pass is that C
  *built from its six parts and the fluid primitives* equals that difference (fresh slices: fully; stored: for the
  non-fluid parts).

Unpaid inputs: x_c, the pocket's shape and width, the cell sizes, the seed, the one-budget constants, the dust as the
clock's carrier. Forced: spherical symmetry, CMC slicing, the wall (identical neighbours in step), cold start.

## Next (not done; production, for review)

Writing `cosmos_rhs`'s dH in Raychaudhuri form and giving `exchange_Q` this meaning is a production change. Before it:
evaluate the φ and ρ companions on the stored states; choose a production (λ = 0) state for the first test.
