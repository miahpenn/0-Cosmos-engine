# The stone split: rest-frame w in the species, the lag elsewhere (research only)

Chat's proposal (Nov 15): the S–D species carries its rest-frame w (≈ 0, von Laue: dust). The bound stone's lag behind
the Hubble flow goes into the channel as a spread-of-expansion (kinematical backreaction) term. Check that the two are
not counted twice. Also extend the paper formula with the core's energy gain for the early window. Stored states only;
core r_b = 8. Script `research/stone_split_v1.py`. Output: `stone_split.txt`.

## Summary

- **T1. The lag is one quantity (Read).** Measured from the core's size (ν = d ln R_c/d ln a) and, independently, from
  the field's momentum relative to the slice observers (ν_j), it agrees to 0.01–0.02 in every window and cell (e.g.
  x_c 2.0, t = 20–30: 0.704 vs 0.689). The slice-frame momentum carries exactly the lag. Its kinetic part is 2–4% of the
  cell's S–D energy (x_c 2.0) and 1–2.5% (x_c 1.5), falling with time.
- **No double counting is possible in the machine's own channel, because the slot chat names doesn't exist there.**
  The channel is CMC (CHANNEL_NOTE): K is uniform, so there is no spread-of-expansion term. The lag lives in one place
  only, inside the S–D sector's slice-frame ρ + 3p. Moving it to a backreaction term means changing the slicing or the
  averaging frame, and then it lives there instead. It is never in both.
- **T2. But the split doesn't hold at the 1–2% level: settled stones are not quite dust in their own frame.** A lane
  with the species as dust tracks the cell's rest-frame S–D energy (slice energy minus the lag's kinetic part) to
  −1.3% to −1.9% (handoff t = 13), −0.6% to −1.0% (t = 18), −0.1% to −0.5% (t = 24). The rest-frame energy grows
  against a^-3 while the slice-frame energy falls against it (+0.5% to +1.5%). **About half of the lag's kinetic energy
  ends up inside the stone** as it decays (the core's E′/(HE) ≈ +0.02 to +0.03), and the rest leaves. This is infall
  energy being absorbed, as in virialisation.
- **T3. The energy-gain extension does not pass.** Adding the self-similar gain terms (½E″R² + 2E′RṘ)/3E fixes x_c 1.5's
  early window (−0.0028 → −0.0006). It barely helps x_c 2.0 (−0.0064 → −0.0053) and spoils the late window
  (+0.0004 → −0.0024). The gain arrives at the core's edge (infall), not spread through it self-similarly. The paper
  formula without it stays the tested form, valid after settling.

## A caution on the name (for chat)

In the Newtonian limit, Buchert's kinematical backreaction Q = (2/3)(⟨θ²⟩ − ⟨θ⟩²) − 2⟨σ²⟩ vanishes identically on a
closed or periodic domain (Buchert & Ehlers 1997, as recalled here, not re-checked). The stone's lag is peculiar motion. Its
cosmic bookkeeping in that limit is the cosmic energy (Layzer–Irvine) equation: peculiar kinetic energy against the
binding energy, d(K + W)/dt + H(2K + W) = 0. That matches T2: kinetic energy of the lag passing into the bound
stone. **So I'd file the lag as the peculiar-kinetic/binding exchange, not as expansion backreaction (Read; borrowed
name, not derived here).** Chat's scale argument is unaffected: everything here goes as (R_c H)², so it fades for
stones small against the horizon.

## The contract, after this

| item | status |
|---|---|
| species' w, settled stones | rest frame ≈ 0 (tested ≤ 0.0013, CORE_VIRIAL Amendment 1); in the slice frame plus the lag term |
| the lag | one quantity (T1). In the CMC lane it sits in the S–D sector's slice-frame stress, where the w(a) fit and the oracle lanes put it. Not a separate channel term. |
| species as dust in its own frame | **fails at 1–2% here**: the stone absorbs about half of the lag's decaying kinetic energy. ∝ (R_c H)²: a toy's-shadow term, small for realistic stones (Read) |
| early window (core filling) | Open; self-similar gain extension fails |
| leftover held step | Open |

Unpaid inputs: r_b = 8; smoothing (two boxcars of π); the first-order frame boost; a one-stone cell standing in for a
population.
