# What λ is: the rescaling, the one-budget relation, and H/μ as a ruler

Task 1 of the Nov 14 order. Algebra and 0* only. RESEARCH ONLY; no runs; nothing in `engine/` changed.
Pinned by `tests/test_research_sd_rescaling.py`.

Marks: **Derived**, **Read**, **Borrowed**, **Forced**, **Open**. μ = S's mass = 1 (the machine's time unit).
ε ≡ H_today/μ.

## Summary

- **Derived:** λ is not a parameter of S–D matter's own dynamics. It is how strongly that matter gravitates, κ/λ.
- **Derived:** the one-budget rule ties λ to ε: **λ ε² = 4π/(f·v) = 49.0**, with f = 0.10 the floor fraction and
  v = V0/H_today² = 2.565.
- So the machine has one unpaid scale number, ε = H_today/μ (map row 14), plus f and v. λ is not a fourth.
- λ is **not settled** by this: f is put in by hand and v is borrowed through the imported "today" (map row 10).
- ε is a ruler. λ = 5000 is ε = 0.099; production densities are ε = 0.0057. Both are toy universes.
- Proposed second value: **λ = 20 000 (ε = 0.0495)**, same grid, twice the run length. Not run.

## 1. The rescaling (Derived)

V_SD = ½S² − ½D² + (λ/8)(S⁴ + 6S²D² + D⁴). Put S = s/√λ, D = d/√λ:

    V_SD = (1/λ) [ ½s² − ½d² + ⅛(s⁴ + 6s²d² + d⁴) ],   kinetic and gradient terms likewise × 1/λ.

- The S–D action is 1/λ times a λ-free action. The field equations for s and d contain no λ.
- D's resting value: d* = √2. Its vacuum energy: −1/(2λ). S's mass with D at rest: 2. D's own mass there: √2.
  All at every λ.
- x = λA² is the scaled amplitude squared. Results stated in x are already in scaled form.
- λ survives in one place: Einstein's equations see κ·T_SD = (κ/λ)·t_sd. **κ/λ is the sector's gravitational
  strength.**

A classical field has no ħ, so there is no "mass in Planck units" here. The amplitude scale 1/√λ plays that part:
κ/λ is (field amplitude / Planck scale)².

## 2. The one-budget relation (Derived from the rule; the rule's inputs are not)

The rule (AMENDMENT-REV20-2 §4, as coded in `research/dm_one_budget_v1.py:65`): D's ejected vacuum is the floor.

    κ/(2λ) = f · V0′,    V0′ = the cosmic potential scale in μ² units.

With V0′ = v·ε² (v = V0/H_today² = 2.565 from the present-day state):

    λ ε² = 4π/(f v) = 49.0,     κ/λ = 2 f v ε² = 0.513 ε².

| universe | ε = H_today/μ | λ | density scale c² | D* |
|---|---|---|---|---|
| research runs | 0.0990 | 5 000 | 304.92 | 0.020 |
| proposed second value | 0.0495 | 20 000 | 76.23 | 0.010 |
| production densities | 0.00567 | 1.52 × 10⁶ | 1 | 0.00115 |
| S read as ordinary matter (order of magnitude) | ~10⁻⁴⁰ | ~5 × 10⁸¹ | — | ~10⁻⁴¹ |

- Without the rule there would be two unpaid numbers: ε and κ/λ. The rule makes the second a function of the first.
- The last row is conditional. Production asserts ε = 5.67e-3, i.e. μ ≈ 2.55e-31 eV (map row 14). 10⁻⁴⁰ is what ε
  would be if S's clock were a particle's. Which reading is meant is the question row 14 already holds open.
- The rule does not explain why ε is small. It moves the smallness into λ. In that sense this is the
  cosmological-constant problem restated in the machine's variables, not solved.

## 3. 0* on the inputs

| input | label | why |
|---|---|---|
| λ | **not independent** under the one-budget rule | = 49.0/ε² |
| the one-budget rule itself | Derived (Nov 13, chat) | not re-derived here |
| f = 0.10, the floor fraction | **Forced** (put in by hand) | map row 1: "written by hand". No derivation in the repo. It sets λ ε² directly: f → 2f halves λ. |
| V0 | **Borrowed through "today"** | `construct_present_day` hard-codes the state from Ω_de = 0.69 (map row 10); V0 is the potential normalisation reconstructed to fit it |
| v = V0/H_today² = 2.565 | Borrowed, same route | a pure number of the imported present-day state |
| ε = H_today/μ | **Unpaid; asserted by a shared unit** | map row 14: one length unit shared between lanes with μ = 1. "V0/μ² not fixed by the equations." |
| μ = 1, G = 1 | units | COORDINATE |

**Verdict.** λ is RULER-DEPENDENT through ε and is not a derived value at any of the rows above. The single unpaid
number is ε; f and v ride on it. λ ≈ 1.5 million is "borrowed through gap #6 (today) and the hand floor", as chat
put it. Agreed.

## 4. ε and the reset clock

- The reset, in Jeremiah's reading, is set by matter's master clock: S at the saddle. ε is the expansion rate
  measured in ticks of that clock. Any statement "the reset happens when the expansion meets matter's clock" is a
  statement that ε-like quantity reaches order one.
- The repo's only number for it (map row 9): the analytic clock-stop candidate **|H| = 0.14–0.37 μ**. Not in code.
- **Read:** the λ = 5000 universe starts at H = 0.100 μ. Its whole expansion runs within a factor 1.4–3.7 of that
  band; it turns round at t = 55 and by t = 70 is contracting at |H| = 0.044 μ, heading back toward it.
  The production universe sits a factor 25–65 below the band.
- So the research runs are close to the regime where the candidate reset condition lives. Results from them may be
  near-reset physics and not generic. A reason to treat ε as a ruler, independent of cost.
- "C6" is chat's label; the repo does not carry the core sheet, so I connect to rows 9, 10 and 14 and leave the
  mapping onto C6 to chat. **Open.**

## 5. What can depend on ε (Derived classification; none of it checked at a second ε yet)

Held fixed when ε changes: the scaled fields (x_c), the pocket's size and the cell's size in 1/μ, the one-budget rule.
Then every gravitational source scales as ε² (cosmic densities and κ/λ alike), K as ε, and the number of S
oscillations per e-fold as 1/ε.

| class | expected | results we hold |
|---|---|---|
| **A. ε-free** (ratios of two ε² quantities; identities; field-sector facts) | carries over | D* and S's mass 2 at rest; the pocket's share of the cell's density (hence the channel's leading S–D part as a fraction of the fluids' source); e-folds to turnaround (1.75–1.79); the κ conversion; dK/dt = ⟨αQ⟩; β companion = 0 on the dust clock; the forms of the averaged equations |
| **B. ∝ ε²** (one pocket's own gravity, at fixed size in 1/μ) | **the toy's shadow**: shrinks 4× at the second value, 305× at production densities | the pocket's inner redshift (5%); every difference between clocks (wall, volume, dust, S); lapse dips and 2m/R; shear and spread-of-Π parts of the channel; companion terms; T5's label difference |
| **C. ∝ ε or better** (expansion per oscillation; adiabaticity) | improves as ε falls | S's action drift under expansion (2–3%); "D costs matter no action" (1–2%); the S–D species' departure from dust after settling (w ≈ 0.01); the share of an e-fold spent settling (w = −0.21 lasts 0.6–1.2 Hubble times at ε = 0.1); lane-follow residuals |
| **D. cannot say from scaling** | must be run | no LOCK; D swinging through zero at the centre and the dust cusp (focusing is field-sector, the dust's response is gravity); no trapped surface in cold starts (B suggests safer at small ε, not proved); hot-start collapse; anything through the turnaround in oscillation counts |

- In the research runs the pocket (radius 8/μ) is 0.8 Hubble radii across. At production densities it is 0.045.
  The toy's pockets are horizon-sized; most of class B follows from that.
- **A consequence to flag (Open):** at fixed x the pocket's density is a fixed multiple of the cosmic density at every
  ε (about 0.2x here). A pocket with x of order one is never dense matter. Reading S as tangible matter needs either
  very large x or a different relation than the one-budget rule gives. For chat.

## 6. The second value

**λ = 20 000, ε = 0.0495** (c² = 76.23).

- Same grid: r_max = 24, 36; dr = 0.1; x_c = 2.0; cold. The pocket's share is unchanged by construction.
- Run length doubles for the same e-folds: turnaround near t = 110, so T = 140. About twice the T = 70 runs, same
  cells. Short enough to count as a test, not a campaign; four runs (two cells, with and without D).
- Why not production densities (ε = 0.0057): 17 times longer, and class B falls to 1e-4, below what dr = 0.1
  resolves cleanly. A factor of two in ε already separates ε⁰, ε¹ and ε² (1, 2 and 4).
- Why not larger ε (λ = 1250, ε = 0.2): cheap, but it starts inside the clock-stop band. Useful later as a third
  point for exactly that reason; not as the control.

**Predictions to be checked, stated before the run:**

| quantity | at ε = 0.099 | predicted at ε = 0.0495 |
|---|---|---|
| S–D density / dark-matter density, r_max 24 | 0.084–0.093 | the same |
| channel's S–D part / fluids' source | 0.073–0.077 | the same |
| e-folds to turnaround | 1.75–1.79 | the same |
| S's clock rate / dust's clock rate − 1 | −4 to −5% | −1 to −1.3% |
| wall clock − volume clock | 0.25% | 0.06% |
| shear part of the channel / fluids' source | 0.002–0.010 | ÷ 4 |
| S's action drift over the expansion | 2–3% | about half, or better |
| settling as a share of the first e-fold | 0.6–1.2 Hubble times | half |
| LOCK, central knot, cusp | none / present / present | not predicted |

A result that lands on its predicted row is classified. One that does not is the more interesting.

## Not done

The second-value runs. The equation of state (next, in scaled variables). Any production change.
