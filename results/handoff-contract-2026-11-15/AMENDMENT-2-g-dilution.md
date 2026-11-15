# Amendment 2 to the handoff contract note: how the gradient piece dilutes, and what that does to w(a)

Research only; stored states (five cells); nothing in `engine/` changed.
`research/handoff_contract_v1.py`: `g_split()`, `w_lanes_split()`. Outputs: `g_dilution.txt`, `w_lanes_split.txt`,
`w_core_outer.txt`.

## Summary

- **Outside the pocket the a^-2 law holds (Derived; Read).** Free waves' gradient share falls as a^-1.87 to a^-1.91
  (t = 13–54, all cells).
- **The pocket core does not follow it.** The core (r < 6) holds 81–91% of the S–D energy after the release, and the
  share is rising. Its gradient share stays nearly flat: a^-0.2 to a^-0.3 at x_c 2.0, a^-0.8 at x_c 1.5. The pocket
  is bound: its physical size grows about ×1.9 while a grows ×2.8 (t = 13–54).
- **So Amendment 1's w(a) was right in effect but not in form.** Feeding in the measured core law (a^-0 to a^-0.5)
  makes the lane worse (−1.6% to −3.6% by t = 54). The all-a^-2 form keeps ≤ 0.3%. **The core's pressure is not its
  gradient piece:**

| cell, t = 13–54 | core w | core gradient piece (2/3)g | outside w | cell w |
|---|---|---|---|---|
| x_c 2.0, r_max 24 | +0.008 | +0.018 | +0.002 | +0.007 |
| x_c 2.0, r_max 36 | +0.008 | +0.019 | +0.002 | +0.007 |
| x_c 1.5, r_max 24 | +0.002 | +0.008 | −0.001 | +0.001 |

  The core carries about half (x_c 2.0) to a quarter (x_c 1.5) of its gradient pressure. The rest is cancelled by its
  binding: the S–D interaction and the non-virial remainder. That is what a bound object does: von Laue, ∫p dV → 0 for
  a stationary bound object. Here it is not fully stationary: the core still contracts in comoving terms.
- **The in-sea formula w = (2/3)g − (3/28)x was derived for a homogeneous sea.** It describes the outside, not the bound
  pocket where most of the S–D energy now sits.

## What changes in the contract

| item | before | now |
|---|---|---|
| S–D species' w | w(a) = (2/3)g_h (a_h/a)² − (3/28)x_h (a_h/a)³, "form Derived" | the same expression, **demoted to a fit** (≤ 0.3%, five cells); not derived for the core |
| outside piece | — | Derived: in-sea formula, g ~ a^-2 (measured a^-1.9) |
| core piece | — | **Open**: a bound pocket's pressure. Von Laue gives 0 if stationary; measured +0.002 to +0.008 |

## Proposed next (derivation, chat or here)

The bound pocket's pressure as it relaxes: virial with the time-dependent term. The volume-integrated trace for a
non-stationary bound object is ½ d²I/dt², where I is the energy-weighted ∫ r² ρ dV. The core's w should then be
(1/6)(d²I/dt²)/E in proper time. That can be tested on the stored cells with no new run. If it holds, the core's w
goes to 0 as the pocket settles, and the contract gets a derived w for both pieces.

Unpaid inputs: the core radius (r = 6, from the EOS note); period averaging over π for core w, which is noisy (±0.02
snapshot to snapshot), so only long-time means are quoted.
