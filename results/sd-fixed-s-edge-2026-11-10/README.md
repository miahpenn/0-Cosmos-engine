# Fixed-S, no-gravity edge test (HANDOFF-to-projects Nov 10, Task 1)

Research only; no production change. `research/sd_fixed_s_edge_v1.py`. Unpaid input: **λ = 5000** (GEAR-67 3(e)).
No gravity, so no vacuum constant enters (neutral on V0 vs V0\*); G_N plays no role; G_reset is not involved.

**Set-up.** S is prescribed from the x_c = 2.0 and 1.5 oscillaton profiles (`osc_l5000_x*.json`, from
`research/sd_oscillaton_v1.py`): `osc` = A(r) cos(ω(r) t) with ω = √(1 + 3x/8); controls `avg` (static S² = A²/2) and
`peak` (static S² = A²). Only D evolves, flat space, GEAR-67 quartic, seed 1e-10 inside only (width 6), r_max = 120,
Neumann outer edge. With no gravity, coordinate time is local proper time (one clock).

## 1. Energy-conserving runs: the ejected region floods the core (`edge_hot_no_friction.txt`)

With no dissipation, the energy released where D falls into its well (ΔV = 1/(2λ) per volume) stays in D as waves.
They run inward and focus at the origin: |D(0)| reaches 0.03–0.18 (1.5–9 × D\*) by t ≈ 50 for both pockets, for the
oscillating S **and** for the static averaged control (`avg`, where D = 0 is the only minimum in the core). So the
inward move is not the hold failing; it is the released energy having nowhere to go. The same happens in the engine
runs (no dissipation there either), which is a reading of the t = 80 inside-seed runs to keep in mind.

## 2. With the released energy removed: the energetic prediction holds (`edge_cold_friction.txt`, `settled_profiles.txt`)

A diagnostic friction −γ Ḋ (not physics; γ = 0.3, checked at γ = 0.1, dr = 0.05 and 0.025, final profile unchanged to
~1%) lets D settle. The settled D follows the cold static profile D² = (2/λ)(1 − 3x/4) outside and falls continuously
toward zero at **x = 4/3** (r = 9.75 for x_c = 2.0, 6.0 for 1.5), with a gradient tail inside:

| | x_c = 2.0 | x_c = 1.5 |
|---|---|---|
| D / D\* at r = 16 (x ≈ 0.63): run / prediction | 0.727 / 0.722 | 0.737 / 0.735 |
| D / D\* at r = 12 (x ≈ 1.08 / 0.93) | 0.466 / 0.437 | 0.563 / 0.553 |
| D / D\* at the x = 4/3 radius (osc / avg) | 0.29 / 0.28 (tail) | 0.35 / 0.27 (tail) |
| D / D\* at the centre | **0.015** | **0.14–0.16** |
| oscillating S vs static ⟨S²⟩ | same to ~2% | same to ~10% in the core |

- **There is no wall.** The transition is continuous (second order): D's amplitude goes to zero at x = 4/3 and leaks
  inward as a tail. "Edge = |D| ≥ half the cap" sits at x = 2/3 by construction for a settled profile
  (D² = (2/λ)(1 − 3x/4) ≥ 1/λ ⇔ x ≤ 2/3), so that edge measure must not be read as the hold boundary.
- **D feels the time average of S².** Oscillating S and static ⟨S²⟩ = A²/2 give the same settled profile (the
  Floquet-band result of the engine pockets, here in its nonlinear, settled form).
- **The lighter pocket does not collapse.** Its core holds D at 14–16% of D\*, against 1.5% for x_c = 2.0: weaker,
  not lost. No wall tension acts, because no wall forms.

## Limits

S is prescribed (no back-reaction: S cannot lose energy to D, and D cannot push S); flat space; friction is a
device to remove the released energy, not a mechanism. λ = 5000 only.
