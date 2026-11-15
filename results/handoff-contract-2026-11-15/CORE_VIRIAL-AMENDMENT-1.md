# CORE_VIRIAL, Amendment 1: chat's paper cancellation tested, and the exact boundary terms

Chat's Amendment 1 to the derivation (Drive) cancels the two large terms on paper. Chat also withdrew two claims: the
(R·H)² scaling across r_max, and the guess that the leftover held step is the sea's pressure. Stored states only.
`research/core_virial_v1.py`: `analyse4()`, `report4()`. Output: `core_virial_form4.txt`.

## Summary

- **One term was missing from the paper formula, and it is the largest one.** Chat's formula gives the core's
  rest-frame pressure. The stored stresses are measured by observers moving with the Hubble flow, and the core moves at
  ν H ℓ against their H ℓ. So the measured (slice-frame) w carries the relative flow:
  **w_slice = w_rest + (1 + w)(ν − 1)² (R_c H)² / 3.** With ν ≈ 0.3–0.8 this is +0.003 to +0.012, most of the measured w.
- **With that term, the prediction passes after the core settles, and fails while it is still filling.**

  | window | measured − predicted, x_c 2.0 (r_b 6 / 8 / 10) | x_c 1.5 (r_b 6 / 8 / 10) |
  |---|---|---|
  | t = 30–54 | +0.0004 / +0.0004 / +0.0006 (r_max 36: +0.0010 / +0.0011 / +0.0013) | +0.0009 / +0.0006 / +0.0007 |
  | t = 13–30 | −0.0042 / −0.0064 / −0.0074 | −0.0008 / −0.0028 / −0.0040 |
  | t = 13–54 | −0.0016 / −0.0025 / −0.0027 | +0.0002 / −0.0008 / −0.0013 |

  Pass condition (chat): residual below about 0.002. **t = 30–54: passes everywhere (≤ 0.0013).** t = 13–30 fails at
  x_c 2.0 by about w itself. There the core is still gaining energy, so the formula's steady-energy, self-similar
  assumption doesn't hold. Nothing is fitted.
- **The core's rest-frame pressure is close to zero (Read).** The expansion's tug (−0.010 to −0.024) and the support
  against its own gravity (+0.010 to +0.019) nearly cancel. With the small positive boundary term, w_rest comes out
  −0.0007 to +0.0024. **The bound stone is near von Laue: its measured w is mostly its motion against the Hubble flow,
  not pressure.**
- **The exact boundary terms work, and r_b = 6 is now a second ruler.** Keeping the momentum and energy carried across
  the comoving boundary, the identity closes to −0.002 to −0.006 at r_b = 6 (it was −0.010 to −0.015 without them). At
  r_b = 8–10 it closes to −0.004 to −0.008. One more refinement was needed: the boundary stress must be the
  central-frame T^rr, whose boost part largely cancels the momentum term (−B1 net: ℓ_b² A H j).
- The identity's own residual (−0.002 to −0.008) is larger than the paper form's at late times. Taking Ï from data
  is noisier than evaluating the cancellation on paper, which was chat's point.

## Two cautions

- At t = 30–54 the window runs into the turnaround (H → 0 at t ≈ 55–60), so the column means of ν and ν′ are not
  meaningful there. The tug term is finite (R_c² H² ν′ = R_c² H dν/dτ), and it is what's tested.
- The lane uses the slice-frame w. In this cell, that w is mostly the relative-flow term of one bound core. What
  the right comoving w is for a population of bound cores (many stones, each falling behind the flow about its own
  centre) is a question this one-core cell can't settle. **Open; for chat.**

## The contract, after this

| item | status |
|---|---|
| S–D species' w after settling | **Derived and tested to ≤ 0.0013 (t = 30–54)**: w = w_rest + (ν − 1)²(R_c H)²/3, w_rest = tug − W_x/3E + S_b/3E ≈ 0 |
| while the core is still filling (t ≲ 30) | Open: the formula's steady-energy assumption fails; w(a) fit stays for this phase |
| outside piece | in-sea formula, a^-2 gradients (Amendment 2) |
| leftover held step | Open; not the sea's pressure (withdrawn) |
