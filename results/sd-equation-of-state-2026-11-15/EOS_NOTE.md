# S–D matter's equation of state, and the spectrum check

> **Corrected since (Nov 17):** the hold band's edges are 1.270 and 2.670 (#470); the 1.28 and 2.68 rows are the
> nearest samples. See [`AMENDMENT-1-band-edge.md`](AMENDMENT-1-band-edge.md).

HANDOFF Nov 15 2026 (1), Tasks 1 and 2, with the adjustments agreed that day. RESEARCH ONLY. Algebra and stored states;
**no new runs**; nothing in `engine/` changed.

- `research/sd_eos_v1.py` → `eos_algebra.txt`, `eos_stored_states.txt` (plain windows of 6), `eos_stored_hann.txt` (smooth windows of 12)
- `research/sd_spectrum_v1.py` → `spectrum.txt`
- Pinned by `tests/test_research_sd_eos.py`

Marks: **Derived**, **Read**, **Open**. Scaled variables throughout: s = √λ S, d = √λ D, energies in units of 1/λ, μ = 1.
x_A = λA² at S's turning point (the convention of the hold band, PR #470). x = 2λ⟨S²⟩ where a mean is meant.

## Summary

- **Derived:** over a quasi-period, the virial identity fixes w exactly from four averages: kinetic aside, the bond
  B, each element's own quartic U, and the gradient energy G. No fitting.
- **Derived:** three regimes.
  - **Held** (D_s kept out, the band): w = −0.32 to −0.07. It is S's own ring (+0.09 to +0.14) plus the held vacuum
    step (w = −1).
  - **In the sea** (below the band, D_s present): **w = (2/3)g − (3/28)x**. Both terms are small; after settling
    w = +0.006 to +0.015 on the stored states, matching the formula to 0.003.
  - **Dense** (above the band): **w → 1/3 exactly** as soon as the bond is negligible. Dense S–D matter is hot.
- **Jeremiah's prediction:** holds exactly at the saddle (but as an identity); fails in the sea, where the
  coupled/uncoupled split stops meaning anything; the −0.21 is the held vacuum step, as he read it, with one
  refinement.
- **Closure:** no pressure-as-function-of-density exists. The lane needs two components: S–D matter (dust to 1.5%
  after settling) and the held vacuum step (w = −1, decaying as D_s refills).
- **Spectrum:** after settling, one line at 2.00 ± 0.02 in local proper time everywhere except the knot. At the knot,
  a two-frequency quasi-period (1.84 and 1.24). While D_s settles, no clean clock.
- **GEAR-104 confirmed exactly:** released vacuum / matter energy = 1/[x_A(1 + x_A/4)].

## 1. The identity (Derived)

In the elements x₁, x₂ = (s ± d)/√2 the S–D potential is exactly

    V = B + U,   B = x₁x₂ = (s² − d²)/2   (the bond: coupled stiffness)
                 U = (x₁⁴ + x₂⁴)/4 = (s⁴ + 6s²d² + d⁴)/8   (each element's own quartic: uncoupled stiffness)

For bounded motion with slow expansion, the time derivative of ⟨s ṡ + d ḋ⟩ averages to zero over a quasi-period, so
⟨kinetic⟩ = G + B + 2U. Then, above D_s's vacuum (ρ_vac = −½, p_vac = +½):

    p = (2/3)G + U − ½,    ρ = 2G + 2B + 3U + ½,    w = [(2/3)G + U − ½] / [2G + 2B + 3U + ½]

- **Read:** the identity holds on the stored states with D to a virial residual of 0.4–3% (smooth windows of 12).
  Plain windows of 6 are not whole periods (S's period is π with D_s at rest, 2π without D) and leave non-virial
  pieces up to ±0.2 in w. All numbers below use the smooth windows.
- Gradients enter only as (2/3)G in p: a gradient-dominated field has w = 1/3, as does a quartic-dominated one.

## 2. The three regimes

### Held: D_s kept out (Derived, exact by quadrature)

With d = 0, S rings in s²/2 + s⁴/8. Its own w_S = U/(2B + 3U) = 1/(2B/U + 3).

| x_A | S's clock ω_S | w_S | ρ_S/ΔV = x_A(1 + x_A/4) | ΔV/ρ_S | GEAR-104 1/[x(1 + 9x/32)] | w with the held step | uniform D_s mode |
|---|---|---|---|---|---|---|---|
| 1.00 | 1.171 | +0.073 | 1.25 | 0.800 | 0.803 | −0.40 | grows |
| 1.28 | 1.214 | +0.088 | 1.69 | 0.59 | — | −0.32 | held (band edge) |
| 1.50 | 1.247 | +0.098 | 2.06 | 0.485 | 0.488 | −0.26 | held |
| 2.00 | 1.318 | +0.119 | 3.00 | 0.333 | 0.336 | −0.16 | held |
| 2.68 | 1.408 | +0.142 | 4.45 | 0.225 | 0.227 | −0.07 | grows (band edge) |
| 8 | 1.98 | +0.229 | 24 | 0.042 | 0.043 | +0.18 | grows |
| 100 | 6.08 | +0.322 | 2600 | 4e-4 | 4e-4 | +0.32 | grows |

- The hold band comes out at x_A = 1.28 to 2.68 (Floquet, exact ringing): PR #470 reproduced.
- **GEAR-104, exact:** ΔV/ρ_S = 1/[x_A(1 + x_A/4)] with x_A the turning point. The filed 1/[x(1 + 9x/32)] is the same to
  ≤ 1.5% with x = 2λ⟨S²⟩. Scale-free and λ-free, as filed.
- **w with the held step** is exactly (w_S ρ_S − ΔV)/(ρ_S + ΔV): S's ring plus a vacuum component holding the step.

### In the sea: D_s present, below the band (Derived to second order; Read)

S at small amplitude pushes D_s down: the static shift is d² = 2 − (3/2)x. Expanding to second order in amplitude and
using the virial identity for s and for d separately (D's forced ring at twice S's frequency included):

    p = ⟨s⁴⟩/8 + (c/2)⟨δ s²⟩ = −(6/7)⟨s²⟩²,   ρ = 4⟨s²⟩   ⇒   w = −(3/28) x  (+ (2/3)g with gradients)

c = (3/2)√2; δ = d − √2. The static part of δ alone would give −(15/128)x; the forced ring at frequency 4 against D's
mass √2 brings it to −(3/28)x.

- **Read (homogeneous pair, long averages):** D's shift sqrt(1 − 3x/4) confirmed to 0.1%; w = −0.00214 at x = 0.020
  against −0.00215; −0.0114 at 0.101 against −0.0108 (higher order).
- **Read (stored states, after settling):** w − (2/3)g against −(3/28)x, with x weighted by S's energy (S sits in the
  pocket; the cell mean is 30–100 times smaller):

| run, zone | t = 15 | 27 | 39 | 51 | 63 |
|---|---|---|---|---|---|
| x_c 2.0, r_max 24, cell: w | −0.001 | +0.010 | +0.007 | +0.009 | +0.014 |
| gradient piece (2/3)g | +0.021 | +0.018 | +0.015 | +0.015 | +0.023 |
| w − (2/3)g | −0.022 | −0.008 | −0.008 | −0.006 | −0.010 |
| −(3/28)x | −0.016 | −0.007 | −0.005 | −0.005 | −0.010 |
| x_c 2.0, r_max 36, cell: w − (2/3)g / −(3/28)x | −0.023 / −0.016 | −0.008 / −0.008 | −0.008 / −0.005 | −0.007 / −0.005 | −0.010 / −0.011 |
| x_c 1.5, r_max 24, cell: w − (2/3)g / −(3/28)x | −0.010 / −0.010 | −0.004 / −0.004 | −0.000 / −0.002 | −0.001 / −0.002 | −0.002 / −0.002 |

  Agreement to 0.003 from t ≈ 25; to 0.006–0.010 while the virial is still settling (t = 9–21).
- **In units of ρ, both terms are small:** the gradient share g is 2–4% (pocket width against S's Compton length) and
  x after settling is 0.02–0.15.

### Dense: above the band (Derived; homogeneous pair Read)

When the bond is negligible against U and G, the identity gives w = (2G/3 + U)/(2G + 3U) = **1/3 exactly**, however
the energy divides between gradients and the quartic.

| x_A (homogeneous pair, D seeded) | 3 | 4 | 10 | 30 | 100 | 1000 |
|---|---|---|---|---|---|---|
| w above the vacuum | −0.02 | +0.16 | +0.28 | +0.32 | +0.325 | +0.332 |

The homogeneous pair leaves out the Lamé conversion of D_s into D_r at finite wavelength. That changes how the energy
is shared, not the conclusion: any quartic- or gradient-dominated ringing has w = 1/3. **Dense S–D matter is
radiation-like.** If S includes tangible matter, it is not dense S–D matter in this sense. That bears on the open flag
of PR #486.

## 3. Jeremiah's prediction (verdicts)

**(a) "w tracks the coupled-to-uncoupled ratio as x changes."**

- At the saddle (D_s held out): **holds exactly**, w_S = 1/(2B/U + 3), monotone from 0 (bond-dominated: dust) to 1/3
  (quartic-dominated: radiation). But it is the virial identity, true of any quadratic-plus-quartic ringing. It
  confirms the potential's form, not the picture beyond it.
- In the sea: **does not hold.** With D_s at rest, B and U each carry D_s's rest terms (−1 and +½ per unit volume),
  which cancel against the vacuum constant. On the stored states the "bond share" of ρ is −1 in the core and −10 to
  −2400 over the whole cell. The split stops meaning anything, and w is negative, set by how far S pushes D_s down.
- Above the band: holds at the quartic end (w → 1/3).

**(b) "The w = −0.21 phase while D settles is the bond's downhill (saddle) side acting as vacuum."**

- **Holds in substance (Read).** It is the held-out vacuum step. In the pocket core during t = 0–3 the identity gives
  −0.18 (x_c 2.0) and −0.22 (x_c 1.5); the static held formula gives −0.16 and −0.26 at those x_A. The plain-window
  values (−0.30 to −0.39) include non-virial transients of −0.13: the first windows are not stationary.
- **Refinement (Derived):** the step is not the bond alone. Going from the saddle to D_s's rest, the bond's downhill
  direction contributes −d²/2 = −1 and D's own quartic +d⁴/8 = +½. The step ΔV = ½ is their sum: the uncoupled
  stiffness sets where the downhill stops.

**One mechanism at every x (Derived):** matter displacing D_s from rest acts as vacuum. Fully inside the band (the held
step, w = −1 component); partly in the sea (D_s pushed to d² = 2 − (3/2)x, giving −(3/28)x). Not at all once the pair
mixes above the band.

## 4. The closure the cosmos lane needs

There is **no barotropic closure**: w depends on x where S's energy sits (not on the mean density), on the gradient
share g, and, during settling, on how much of the region still holds D_s out.

Minimal closure, two components:

| component | w | variable it needs | status |
|---|---|---|---|
| S–D matter in the sea | (2/3)g − (3/28)x_eff: 0 ± 0.015 after settling | none at the 1.5% level (dust); x_eff and g for corrections | Derived; Read |
| held vacuum step | −1 | the held volume fraction; decays as D_s refills (settling time 6–12 here) | Derived (static); decay Read only |
| dense S–D matter (if ever present) | 1/3 | none | Derived |

- This accounts for the lane-follow results: started after settling, dust was enough to 1e-3 (PR #485). Started at
  t = 0, dust failed by 22% in ρ_SD: the held step is missing.
- **0\*:** w(x, g) is λ-free and ε-free at fixed x and g: it carries over (class A of PR #486). How x_eff and the held
  fraction evolve depends on expansion against S's clock (class C). g depends on the pocket's width in S's Compton
  lengths (unpaid shape input).

**Correction (kept):** the companions note (PR #484) gave the S-alone control as "w = +0.13 at t = 9, +0.09 at 15".
Those were plain windows of 6, about one period without D. Smooth windows give +0.058 and +0.021 over the cell. The
control still carries a non-virial excess of +0.003 to +0.006 late in the run, larger than the expansion correction
(9/8)(H/m)² (≤ 0.0007 after t = 33). **Open:** S alone disperses (x falls 50-fold), which the identity does not cover.

## 5. Spectrum check (Task 2; Read)

Native clock: 2, S's mass with D_s at rest, in **local** proper time (lapse taken out).

| signal (t = 18–54, after settling) | main line / 2 | power in it | broadband floor | phase wander about the mean rate |
|---|---|---|---|---|
| S on its initial profile, wall proper time | 1.002 (r_max 24), 1.004 (36) | 95–96% | 2–3% | 0.46 rad |
| S at r = 3, local proper time | 1.008–1.010 | 99% | 1% | 0.15 rad |
| S at r = 8 | 0.997–0.999 | 100% | 0% | 0.06–0.10 rad |
| S at r = 15 | 1.017–1.019 | 100% | 0% | 0.18 rad |
| S at the centre, x_c 1.5 (no knot) | 1.003 | 99% | 1% | 0.12 rad |
| **S at the centre, x_c 2.0 (the knot)** | **0.917–0.920, with sidebands at 0.62 and 3.07** | **70% + 20% + 9%** | **1%** | **0.54 rad** |

- **Outside the knot: one line at the native clock, to 0.2–2%.** The phase wander, 0.06–0.18 rad over 35 time units,
  bounds the mean frequency to about 0.5%, tighter than the spectral resolution (±18%). The +1.7% at r = 15 fits
  travelling S waves, ω² = 4 + k² with k ≈ 0.37.
- **At the knot: two incommensurate frequencies.** A carrier at 1.835 (8% below native) and sidebands at carrier
  ± 1.24. Quasi-periodic, not broadband (floor 1%), the same in both cells. **Reading, not derived:** 1.24 is D's
  large swing through zero at the knot modulating S's mass; D's small-amplitude mass is √2 = 1.41 and a swing over
  both wells is slower.
- **While D_s settles (t = 4–18):** record 13 units, resolution ±45%; broadband floor 18–34%; the main line moves.
  No clean clock.
- **Without D (control):** one line at 1.00–1.02. Matter's clock is set by D's state: 1 with no D, 1.21–1.41 at the
  saddle with D_s held out (the band), 2 in the sea.
  > Pointer (Nov 22): in-band omega_S is 1.22-1.41, not 1.21-1.41 (Lumen's AUDIT-6 reply, C5; core sheet v3's omega_S^2 = 1 + 0.375x at x = 4/3).
- **Limits:** lines closer than about 18% are not resolved; a Lamé gap structure finer than that is untested. The
  phase wander's own spectrum has its strongest line at 1.36–1.44 near the pocket (D's mass √2: D_r ringing modulates
  S's phase, weakly) and near 1.9 farther out.

**What Task 1 averages over:** after settling, any smooth window of at least four of S's periods (π each); used 12. At
the knot, at least three beats of 2π/1.24 ≈ 5. During settling no window is clean, which is another reason to carry
the held step as its own component.

## Not done

The Lamé conversion at finite wavelength (needs the D_r return-path derivation, queued). The held fraction's decay as a
derived rate (only read here). Any change to `cosmos_rhs`.
