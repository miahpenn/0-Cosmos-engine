# The channel's companion terms, the S–D species, and whether the lane then follows the cell

Follow-up to the channel note and B6 (proposed handoff of Nov 14, items 1–3). RESEARCH ONLY. Stored frozen states and
stored trajectories of the Nov 14 (4) runs; **no new runs**; nothing in `engine/` changed.
`research/channel_companions_v1.py` → `companions_frozen_states.txt`; `research/cosmos_lane_follow_v1.py` →
`lane_follow.txt`.

Marks: **Derived**, **Read**, **Open**. All numbers: one-budget rule, λ = 5000, cold start, dust clock.

## Summary

- **Derived:** the averaged equation for every cosmos-lane variable, and for the S–D species the lane does not carry.
- **Read:** on the dust clock the companions in a, φ, π and ρ_dm are small: 1e-5 to 3e-3 of each equation's own
  right-hand side at r_max = 24. The β term's companion in the π equation is exactly zero on the dust clock.
- **Read:** what the lane lacks is the species. S–D matter is 8–9% of the dark-matter density at r_max = 24 and,
  once D has settled, dilutes as dust (w = +0.004 to +0.010).
- **Read:** a lane carrying that species as dust, started after settling and told nothing more, follows the cell to
  1e-4 – 1e-3 for 54 time units, through the turnaround. Without the species it is off by 5e-3 – 1e-2.
- **Contradicts the channel note (kept):** with the species carried, the lane's own −(ρ + p)/2 form does about as
  well as the Raychaudhuri form. "dH must move to Raychaudhuri form" holds for a correction *term*, not for a carried
  *species*.
- **At production λ = 0 (Derived):** S alone is a clean dust-like species. D is not a species at all: it has active
  gravitational mass with no energy density, and no floor.

## 1. The companion terms (Derived)

Averaging rule on a CMC slice with zero shift, proper volume: d⟨f⟩/dt = ⟨∂t f⟩ − K[⟨αf⟩ − ⟨α⟩⟨f⟩].
Dust clock: dτ = N dt, N = ⟨αρ_dm⟩/⟨ρ_dm⟩. Bars are proper-volume averages. K = −3H.

| lane equation | the cell's averaged equation, per unit dust time | companion (cell − lane) | fed by |
|---|---|---|---|
| da/dτ = aH | aH⟨α⟩/N | aH(⟨α⟩/N − 1) | clock: volume against dust |
| dφ/dτ = π | [⟨αΠ⟩ − K Cov(α, φ)]/N | (⟨αΠ⟩/N − Π̄) − K Cov(α, φ)/N | lapse where Π and φ differ from their means |
| dπ/dτ = −3Hπ − V′(φ) + βρ_dm | [K⟨α⟩Π̄ − ⟨αV′(φ)⟩ + β⟨αρ_dm⟩]/N | KΠ̄(⟨α⟩/N − 1) − [⟨αV′⟩/N − V′(φ̄)] + **0** | clock; Jensen gap of V′; the β part vanishes |
| dρ_dm/dτ = −3Hρ − βρπ | [K⟨α⟩ρ̄ − β⟨αρ_dm Π⟩]/N | Kρ̄(⟨α⟩/N − 1) − βρ̄(Π_dust − Π̄) | clock; Π where the dust is, not the mean Π |
| dρ_b/dτ = −3Hρ_b | K⟨α⟩ρ̄_b/N | Kρ̄_b(⟨α⟩/N − 1) | clock |
| dρ_r/dτ = −4Hρ_r | K[⟨α⟩ρ̄_r + ⟨αρ_r⟩/3]/N | clock and lapse weight | clock |
| dH/dτ | channel note, B3 and B4 | C | — |

- The gradient terms of the π equation are a divergence, ∇·(α∇φ), and integrate to zero at the wall.
- **The β companion in the π equation is identically zero on the dust clock** (dust at rest): ⟨αρ_dm⟩/N = ρ̄_dm is the
  definition of N. The dust clock is exactly the clock in which φ's coupling to dark matter needs no correction.
  A second, independent reason for decision 1.
- Π_dust = ⟨αρ_dm Π⟩/⟨αρ_dm⟩: the dark matter exchanges with the φ it sits in.
- Velocities enter at second order and are left out (the stored states do not hold them).
- N uses the dark matter's rest mass only. The decision says all fluids; the baryons are not in the stored states.
  **Open**, small (the fluids are near-uniform).

**Read** (identities checked against time differences of the stored averages: φ 1e-6, Π 5e-6, ρ_dm 2e-5).
Companions as fractions of each lane equation's own right-hand side, window means, t = 15 / 39 / 63:

| r_max | a | φ (both parts) | π: clock + potential | ρ_dm: clock | ρ_dm: β |
|---|---|---|---|---|---|
| 16 | 5e-4 / 9e-4 / 1e-3 | 5e-4 / 8e-4 / 9e-4 | −8e-4 / −1.4e-3 / +1e-3 | −5e-4 / −9e-4 / +9e-4 | 5e-6 / 3e-5 / 4e-5 |
| 24 | 2e-4 / 4e-4 / 4e-4 | 3e-4 / 4e-4 / 4e-4 | −4e-4 / −7e-4 / +4e-4 | −2e-4 / −4e-4 / +4e-4 | 2e-6 / 1e-5 / 2e-5 |
| 36 | 8e-5 / 1.4e-4 / 1.5e-4 | 9e-5 / 1.3e-4 / 1.4e-4 | −1.5e-4 / −2.5e-4 / +1.5e-4 | −8e-5 / −1.4e-4 / +1.4e-4 | 9e-7 / 5e-6 / 7e-6 |

- In the π equation the clock part and the potential part are each up to 9e-3 and cancel to the figures shown.
- Closure (i) gives the same numbers as (ii). All scale with the pocket's share of the cell.
- **0\*:** derived from the machine's own equations; closure-free on the dust clock; cell-size dependent as the
  pocket's share; CMC and wall forced, as for C. Small enough here that their sign is the only robust reading.

## 2. The S–D species (Derived, then Read)

Energy conservation of the S–D sector (it exchanges with nothing but gravity when β_S = β_D = 0), averaged:

    dρ_SD/dτ = [ K⟨α⟩ρ̄_SD + (K/3)⟨α S_SD⟩ + ⟨α Aa (p_r − p_t)⟩ − ⟨j·∇α⟩ ] / N
                 dilution     pressure work     shear work          lapse work

with ρ_SD its energy density above D's vacuum, S_SD = p_r + 2p_t, p_r − p_t = |∇S|² + |∇D|², j = −P_S∇S − P_D∇D.
In lane form: dρ_SD/dτ = −3H(1 + w)ρ_SD + shear work + lapse work, w = ⟨S_SD⟩/(3ρ̄_SD).

- **It is not closed.** w is not a function of ρ_SD. Cycle mean, from the virial relation:
  ⟨p⟩ = ⅓⟨|∇S|² + |∇D|²⟩ + ⟨½ f·∂V/∂f − V⟩. Zero for a harmonic field without gradients; positive from gradients
  and from the quartic. Within a cycle w swings at twice S's frequency.

**Read** (identity against time differences: 4e-3; window means):

| | t = 3 | 9 | 15 | 27 | 39 | 63 |
|---|---|---|---|---|---|---|
| ρ_SD,E / ρ_dm, r_max 16 / 24 / 36 | 0.26 / 0.084 / 0.027 | 0.29 / 0.093 / 0.030 | 0.29 / 0.093 / 0.030 | 0.28 / 0.090 / 0.029 | 0.27 / 0.088 / 0.028 | 0.26 / 0.084 / 0.027 |
| w (r_max 24) | −0.21 | +0.007 | +0.010 | +0.009 | +0.008 | +0.016 |
| lapse work ÷ dilution | +0.001 | +0.005 | +0.005 | +0.005 | +0.005 | +0.011 |
| shear work ÷ dilution | −6e-4 | −9e-4 | −6e-4 | −4e-4 | −4e-4 | −1.3e-3 |
| gradient share of the energy | 0.047 | 0.037 | 0.031 | 0.027 | 0.022 | 0.033 |

- **While D settles (t < 6) the species is not dust:** w = −0.21. D is falling into its vacuum.
- **After that it is dust to 1%:** w = +0.007 to +0.016, and the lapse work (energy gained by flux moving to lower
  lapse) is +0.5 to +1%, the opposite sign. They nearly cancel. A lane species with the measured w but no lapse work
  does *worse* on ρ_SD (−2.8% by t = 66) than plain w = 0 (−0.6%). **Read.**
- Without D (control) the S pocket alone is stiffer early: w = +0.13 at t = 9, +0.09 at 15, falling to ±0.04. Its
  quartic is half its mass term at this amplitude; with D present S's mass doubles and the quartic matters less.
- x_c = 1.5: w = −0.02 to +0.02.

## 3. Does the lane then follow the cell? (Read, on stored trajectories)

The lane is integrated alone in dust time from the cell's averages at t0 and compared with the cell afterwards. The
companions of §1 are **not** fed in; their size is what the residual should show.

Started at t0 = 12 (after settling), residuals at t = 66, x_c = 2.0:

| lane | r_max | ΔH / H(t0) | Δ ln a | Δρ_dm / ρ_dm | Δρ_SD / ρ_SD |
|---|---|---|---|---|---|
| as it is (no S–D species) | 16 / 24 / 36 | +1.4e-2 / +5.1e-3 / +1.7e-3 | +9.5e-3 / +3.4e-3 / +1.1e-3 | −2.8e-2 / −1.0e-2 / −3.4e-3 | — |
| + species, w = 0 | 16 / 24 / 36 | +2.4e-3 / +1.2e-3 / +4.2e-4 | +8e-4 / +4e-4 / +1.5e-4 | −2.3e-3 / −1.2e-3 / −4e-4 | −3e-4 / −6e-3 / −3e-3 |
| + species, w = 0, the lane's −(ρ+p)/2 form | 16 / 24 / 36 | −1.0e-3 / −9e-4 / −6e-4 | −1.3e-3 / −9e-4 / −4e-4 | +4.2e-3 / +2.8e-3 / +1.4e-3 | +6e-3 / −2e-3 / −2e-3 |

- **Read:** carrying the species as dust improves the lane four- to tenfold in every variable. The remainder is
  1e-4 – 1e-3 to t = 54 and up to 2.4e-3 at t = 66 (the contraction amplifies it). That is the size of the companions.
- **Read:** closures (i) and (ii) agree; residuals scale with the pocket's share.
- **Contradiction with the channel note, B4.** There the −(ρ + p)/2 form failed (B6 control 2). Here, with the
  species inside ρ, it does as well in H and two to three times worse in a and ρ_dm. Reason (**Derived**): the
  Friedmann gap is H² − ρ/3, and most of it was the missing S–D density. Once ρ includes it, what is left is
  ⟨Aa²⟩/4 − ⟨R⟩/6, small. So: a correction *term* needs the Raychaudhuri form; a carried *species* does not strictly.
  The Raychaudhuri form remains the safer one (it does not assume the constraint), and it is clearly better from
  t0 = 0 (below).
- **Started at t0 = 0 (through the settling), r_max = 24:** with the species at w = 0, H stays within 2e-4 and ρ_dm
  within 0.5%, but ρ_SD is 22% low throughout. H comes out right partly by compensation (the species' early
  negative pressure against its too-low density). Not a pass. The −(ρ + p)/2 form is off by 1.2e-2 in H here.
- **Without D (control):** ρ_SD is 9% off, because that pocket is not dust-like early (§2).

**What this establishes:** after settling, in these cells, "cosmos lane + dust-like S–D species" is a good one-way
description at the 1e-3 level. **What it does not:** that it is two-way. The lane is not feeding anything back into
the cell here; the cell is a wall-bounded copy of itself. And ρ_b, ρ_r at t0 were solved from stored totals, not read
independently.

## 4. What the channel means at production λ = 0 (Derived; no runs)

V_SD = ½S² − ½D². There is no D vacuum and no one-budget constant.

- **S:** a harmonic field of mass 1. ρ + 3p = 2(P_S² − ½S²), cycle mean ½ρ_S,E·2 = its density: it gravitates as
  dust. w = ⅓⟨|∇S|²⟩/ρ_S, the gradient share over three, 1% here. A well-defined, closed, dust-like species.
  The no-D control above is *not* this case: it has the λ = 5000 quartic on S.
- **D:** tachyonic. The growing mode has P_D = D, D ∝ e^t. Then
  ρ_D = ½P_D² − ½D² + ½|∇D|² ≈ 0 (gradients only), while ρ + 3p = 2P_D² + D² = 3D².
  **Active gravitational mass without energy density.** No (ρ, w) can represent it: w is unbounded. To carry D the
  lane would need D̄ and P̄_D themselves, plus the spread.
- **No floor.** From a seed of 1e-10, D reaches 0.02 at t ≈ 19 and order one at t ≈ 23. The channel's D part grows
  without limit; nothing in the λ = 0 equations stops it. What stops it in the existing production runs is not
  derived (**Open**; those runs are PROVISIONAL, AUDIT §C).
- **Consequences:**
  - A production-λ test of the channel is clean for **S alone, for t ≲ 15** after seeding, while D is negligible.
  - Any longer statement needs the floor, i.e. λ > 0 under the one-budget rule. That is gap #4's implementation
    step. Until λ > 0 is promoted, the long-time channel is a research result, labelled as such.
  - "The cosmos lane carries S–D matter" means, at λ = 0, "carries S". D's place in the lane depends on λ.

## 0* on the whole

| item | verdict |
|---|---|
| companion terms | Derived; closure-free on the dust clock; ∝ pocket's share; ≤ 3e-3 here; CMC and wall forced |
| β companion = 0 on the dust clock | Derived (identity) |
| S–D species equation | Derived; **not closed** (w, shear work, lapse work are not functions of ρ_SD) |
| "dust-like" | Read, after settling, with D present; not during settling; not for S alone at λ = 5000 |
| lane + species follows the cell | Read, 1e-3 over 54 time units, one-way only |
| Raychaudhuri form required | Derived for a correction term; **not** required once the species is carried (revises B4) |
| λ = 0 | S a species; D not a species; long-time channel needs λ > 0 |

Unpaid inputs: x_c, pocket shape, cell sizes, seed, one-budget constants, the dust as clock carrier, the choice
w = 0. Forced: spherical symmetry, CMC, the wall, cold start, velocities neglected, dark-matter-only clock weight.

## Not done

- The evolved test with a *fresh* run that stores the fluid primitives (independent ρ_b, ρ_r, velocities). Needs a
  short cell run; handoff (5) says no new lattice runs. Waiting for Jeremiah's word.
- Any two-way test (the lane's state changing what the cell does). Needs a design: the wall cell has no outside.
- Any change to `cosmos_rhs`. Production; for review.
