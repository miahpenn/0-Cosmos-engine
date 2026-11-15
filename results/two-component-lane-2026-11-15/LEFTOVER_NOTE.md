# The dark-energy leftover: three routes, the sea's answer, and some pokes

HANDOFF Nov 15 2026 (2), Task 2 (paper only), with a third route added (agreed Nov 15) and some cheap pokes
(Jeremiah: "poke at it and observe"). Algebra and a homogeneous toy; `research/poke_saddle_v1.py` → `pokes.txt`.

Marks: **Derived**, **Read**, **Picture**, **Open**. Scaled units, μ = S's bare mass = 1, grip g = κ/λ.

## The question, and three routes

Today's dark energy must be a small leftover of D's step, not the whole step (POKE-1, C6 note). Three ways it could be:

| route | keeps the one-budget rule? | the small number comes from |
|---|---|---|
| A. a finite region size at ejection, picked by the dynamics | yes | the region's size in matter units (makes κ/λ small) |
| B. something parks φ at its potential's node | yes | φ held where its potential cancels D's step |
| C. the floor is only a fraction of the step (POKE-1) | no | part of the step never released |

## Route A: a finite region size (Derived; one poke)

- **In C6's form, size cannot enter.** The balance there is breathing against matter's clock, and a uniform region's
  breathing is fixed by its density alone: H² = κρ/(3λ), whatever its size. Size would only matter for a
  compactness criterion (ejection at a given 2m/R), which nothing yet proposes. **Not supplied** for C6 as stated.
- **The machine does make large sizes (poke P3).** In POKE-1's homogeneous toy the saddle vacuum inflates before D
  lets go:

| g = κ/λ | e-folds, seed 1e-6 | seed 1e-10 | H at release |
|---|---|---|---|
| 0.5 | 7.4 | 11.4 | 0.21 |
| 2 | 19.9 | 31.6 | 0.37 |
| 3 | 27.9 | 44.3 | 0.43 |
| 30 | 237 | 378 | 1.16 |

  **N ≈ (0.54–1.0) g ln(√2/seed)**: hilltop inflation off the saddle (D's tachyonic mass there is 1). POKE-1's 27–234
  reproduced.
- **So a region's size e^N is a power of D's seed,** and the seed is an unpaid input. Route A's large number is
  borrowed from the seed.
- **Even then, it does not make dark energy small.** Expansion dilutes matter, not vacuum.
- **Verdict:** supplied in form (a size generator exists), unpaid (the seed), and not a route to a small leftover.

**Read, a side observation:** H at release depends on g alone, not on the seed. For g = 0.5–2 it lands at 0.21–0.37,
inside the clock-stop band 0.14–0.37. In the toy, "H at ejection meets matter's clock" is automatic for any grip up to
about 2. It does not select g.

## Route B: parking φ at the node (Derived)

- **The lane's potential cannot park φ.** Under the one-budget rule the floor is D's step, which grows with φ through
  GEAR-91's factor: V = V0 e^(−φ) − S0 e^(cφ), c = −4β_D = 0.16. Then V′ < 0 everywhere: no minimum, φ rolls through
  the node.
- **The dark-matter coupling does not park it either.** φ̈ = V0 e^(−φ) + β ρ_dm has a moving balance point at
  V0 e^(−φ) = |β| ρ_dm. It tracks the diluting dark matter and passes the node once, when ρ_dm = f V0/|β|.
- **A pair of exponentials can have a minimum** (for a step falling with φ) only at negative V, unless the two are
  tuned. A minimum at V = 0 is the cancellation itself, put in by hand.
- **Verdict: not supplied by the machine.**

## Route C: a held fraction that never releases (Derived; poke P1)

The only thing in the machine that holds part of the step is matter in the hold band (GEAR-104's matched pair).

- **Cosmically it is not dark energy.** Held vacuum inside bound pockets dilutes with the pockets' number, as a⁻³:
  w = 0 on cosmic scales, matter-like.
- **Held regions that expand with the cosmos release.** Their S dilutes out of the band (Task 1).
- **No part of the machine keeps a held fraction constant per volume.** **Verdict: not supplied.**

**Poke P1, flagged:** inside a held pocket the held step's negative pressure lowers the pocket's active gravitational
mass, ρ + 3p = 2B + 6U − 1:

| x_A | 1.20 | 1.27 (band edge) | 1.5 | 2.0 | 2.67 (band edge) |
|---|---|---|---|---|---|
| (ρ + 3p)/ρ | −0.02 | +0.04 | +0.22 | +0.52 | +0.80 |

- Held matter is nearly weightless at the band's lower edge.
- Its active mass would turn negative at x_A = 1.222, 3.5% below the edge, where the step can no longer be held.
- Coincidence-class: the two conditions are different equations.
- It does say pockets near the lower edge barely gravitate, which bears on the knot and on the two-stones test.

## 2c. The sea's answer to the breathing (Derived)

**A slow breath gets no answer.**

- The sea (D_s at rest) carries waves only above its gaps: √2 for D, 2 for S. A uniform breathing slower than that
  is answered adiabatically.
- **The machine is classical (no ħ),** so there is no vacuum to amplify. Only waves already present can be pumped.
  A slow breath leaves the sea with no waves at all, not merely few.

**Parametric pumping needs the breath to oscillate at twice a gap or faster:** 2√2 = 2.83 for D, 4 for S.

- The smooth breathing at the reset (H ≈ 0.2–0.45) is far below.
- S's ringing modulates the breathing at 2ω_S: 2.42–2.81 across the band, just below D's threshold.

**Poke P2, flagged:** at the band's upper edge 2ω_S = 2.814, 0.5% below 2√2. Checked: the Floquet trace at
ω_S = √2 exactly is 2.28, not 2, so the edge is not where S reaches D's gap. Coincidence-class.

**Verdict:** the answer removes no energy from the step (it is the vacuum the sea sits in; answering does not change
its value), and it does not bear on Route A. **Not supplied.** It requires the outer sea to be the machine's own D_s
(GEAR-107's open edge, gap #3, parked).

## Verdicts

| route | verdict |
|---|---|
| A. region size | size generator exists (saddle inflation, N ∝ g ln(1/seed)); borrowed from the seed; does not shrink the leftover |
| B. φ parked at the node | not supplied by the machine (no minimum at V = 0 without tuning) |
| C. floor a fraction of the step | not supplied (held vacuum dilutes as matter or releases) |
| 2c. the sea's answer | not supplied (classical, below the gaps: adiabatic) |

**The leftover remains unpaid.** None of the machine's parts makes today's dark energy a small leftover of D's step.
The hold on promoting λ > 0 stands. What would change this, as far as I can see: a mechanism that keeps a held
fraction constant per volume as the cosmos expands, or a dynamical reason for φ to stop at the node. Neither exists in
the machine yet.
