# Gap #1: the local-to-cosmos channel. Derivation note with 0* verdicts

HANDOFF Nov 14 2026 (5), Task B. Algebra and 0* only. No production change. No new lattice runs: the numbers below
are read from stored frozen states of the Nov 14 (4) runs (`research/cosmos_channel_v1.py`,
`channel_frozen_states.txt`).

Marks: **Derived**, **Read**, **Open** as before. Normalisation: Task A's (`KAPPA_CONVERSION_DECISION.md`); ρ_E = 8πρ.

## Summary

- **Derived:** in a cell with zero lapse gradient at the wall, on CMC slices, dK/dt = ⟨αQ⟩ exactly, with ⟨ ⟩ the
  proper-volume average. It holds for either closure. **Read:** to 2e-6 in every stored run.
- **Derived:** the candidate channel C = ⟨αQ⟩ − Q(averaged fields) has six parts: shear, the S–D sector, the spread
  of Π, the Jensen gap of φ's potential, fluid velocities, and a lapse (clock) term.
- **Read:** its leading part is not subtle. It is the S–D pocket's own ρ + 3p, which the cosmos lane does not carry
  at all. The genuinely nonlinear parts are 10–100 times smaller in these cells.
- **0* verdict:** derived, not borrowed. Conditional on: proper-volume average (forced by the derivation), CMC
  slicing, the wall (identical neighbours in step), and a named clock. It scales with the pocket's share of the cell.
- B6 (the minimal test) is specified, not run.

## B1. What Q is, and what exchange_Q was

**Q in closure (ii)** (`research/lattice_cell_v1.py:292`, same expression as `engine/cmc_gauge.py:192`):

    Q = K_ij K^ij + 4π(ρ + S) = (3/2)Aa² + K²/3 + 4π(ρ + p_r + 2p_t)

It is the whole source of K's equation, ∂tK = −∇²α + αQ: a focusing rate, units 1/time². Most of it is what the
cosmos lane already has.

**exchange_Q in `cosmos_rhs`:** an argument that is accepted and ignored. The repository holds no definition of it:
the docstring calls it a "phenomenological reciprocal source", and it is already ignored in the earliest commit that
contains it. By its name and place it was one extra scalar source, an exchange rate.

**They are not the same quantity.** Q is the total; exchange_Q was meant to be an addition. What the derivation offers
in exchange_Q's place is the difference C of B3, entering dH/dt as −C/3.

**Correction to the Nov 14 (4) README:** it described closure (ii) as "K̇ = ⟨Q⟩". Stated properly: K̇ = ⟨αQ⟩ is
forced for any closure (below); closure (ii) is the choice of clock in which ⟨αQ⟩ = ⟨Q⟩, i.e. the Q-weighted mean
lapse is 1. Closure (i) is the clock with α = 1 at the wall. The handoff's lapse-weighted form is the right one.

## The identity (Derived)

CMC: K depends on t only. The lapse equation is ∇²α = αQ − dK/dt. Integrate over the cell with the proper volume
element. The left side is the flux of ∇α through the boundary, which is zero (regular centre; zero lapse gradient at
the wall, W4). So

    dK/dt = ⟨αQ⟩,    ⟨f⟩ = ∫ f √γ d³x / ∫ √γ d³x.

- The proper-volume measure is not a choice here. It is the only weight for which the lapse term integrates away.
- The closure does not enter. It sets the overall scale of α, that is, which clock t is.

## B2. Homogeneous limit (Derived)

Uniform cell, α = 1, Aa = 0, K = −3H, S = 3p:

    −3 dH/dt = 3H² + 4π(ρ + 3p) = 3H² + (ρ_E + 3p_E)/2
    dH/dt = −H² − (ρ_E + 3p_E)/6 = −(ρ_E + p_E)/2      using H² = ρ_E/3.

This is `cosmos_rhs`'s dH. **Read** (Task A): the production lapse solve gives K̇ = −3 dH/dt to 2e-12 on the frozen
chart states.

One point matters later: the cosmos lane's form −(ρ + p)/2 has already used the Friedmann constraint. In an
inhomogeneous cell that constraint changes too (B4).

## B3. The inhomogeneous cell (Derived)

For each scalar the gradient energy cancels in ρ + S: ρ + S = 2Π² − 2V. For dust ρ + S = ρ W²(1 + v²). So

    4π(ρ + S) = [Π² − V_c(φ)]                           φ (Einstein-normalised already)
              + 8π[P_S² + P_D² − V_SD(S, D)]              S and D
              + ½ρ_dm W²(1+v²) + ½ρ_b W²(1+v²) + ρ_r(…)   fluids, ρ_E densities

The cosmos lane evolves averages f̄ = (φ̄, Π̄, ρ̄_dm, ρ̄_b, ρ̄_r) and has no S, D, shear or velocities. Its source is
Q(f̄) = K²/3 + Π̄² − V_c(φ̄) + ½ρ̄_dm + ½ρ̄_b + ρ̄_r. The difference C = ⟨αQ⟩ − Q(f̄) is

| part | form | sign | fed by |
|---|---|---|---|
| shear | (3/2)⟨αAa²⟩ | ≥ 0 | any inhomogeneity, through gravity |
| S–D sector | 8π⟨α(P_S² + P_D² − V_SD)⟩ (less D's vacuum constant when λ > 0) | cycle mean ≥ 0 for oscillation about a minimum (then = half the sector's ρ_E: it gravitates as dust); instantaneous sign alternates at twice S's frequency | S and D themselves |
| spread of Π | ⟨Π²⟩ − ⟨Π⟩² | ≥ 0 | φ made non-uniform: β ρ_dm where dust gathers; the lapse |
| Jensen gap | −[⟨V_c(φ)⟩ − V_c(⟨φ⟩)] | ≤ 0 (V_c is convex) | the same |
| fluid velocities | ½⟨ρ_dm[W²(1+v²) − 1]⟩ + … | ≥ 0, second order in v | infall |
| lapse (clock) | ⟨(α − 1)·(cosmos-carried terms)⟩ | either; depends on the clock | lapse differences across the cell |

- Dust density itself adds nothing beyond its mean: ρ + S is linear in ρ. Clumping enters only through the lapse
  term and the velocities.
- Over the whole of Q, ⟨(α − 1)Q⟩ = 0 in closure (ii)'s clock, by construction.
- **Sign of the whole:** positive in every stored run with D, in closure (ii)'s clock. More focusing: dH/dt more
  negative by C/3.

**Read** (cold start, window means, in units of the fluids' own source F = ⟨½ρ_dm + ½ρ_b + ρ_r⟩):

| run | C/F, t = 15 … 63 | S–D sector | shear | spread of Π | Jensen | lapse on cosmos terms |
|---|---|---|---|---|---|---|
| x_c 2.0, r_max 16 | +0.256 … +0.263 | +0.240 … +0.224 | 0.005 … 0.015 | 0.002 … 0.013 | −0.0005 … −0.001 | +0.010 … +0.012 |
| x_c 2.0, r_max 24 | +0.084 … +0.093 | +0.077 … +0.073 | 0.002 … 0.010 | 0.001 … 0.006 | −0.0004 … −0.0007 | +0.004 … +0.005 |
| x_c 2.0, r_max 36 | +0.027 … +0.030 | +0.025 … +0.023 | 0.001 … 0.004 | 0.0003 … 0.002 | −0.0002 … −0.0003 | +0.001 … +0.002 |
| x_c 1.5, r_max 24 | +0.064 … +0.063 | +0.060 … +0.053 | 0.001 … 0.005 | 0.0006 … 0.003 | −0.0002 … −0.0004 | +0.002 … +0.003 |
| no D, r_max 24 | +0.003 … +0.0002 | +0.003 … +0.0002 | 1e-5 | 2e-5 … 5e-6 | −8e-6 | 9e-5 … 2e-6 |
| x_c 2.0, r_max 24, closure (i) | +0.001 … −0.89 | same as (ii) | same | same | same | −0.08 … −0.97 |

- The S–D part holds a fixed ratio to F through expansion and turnaround: the pocket's sector redshifts as dust does.
- Shear and the spread of Π grow relative to F (to 1% and 0.6% at r_max = 24). They are the parts no added species
  would capture.
- In the raw file S and D are listed apart. Apart they are misleading: S's extra mass energy (¾λS²D²) is booked under
  D. Only their sum is meaningful.

## B4. In the cosmos lane's variables

With H = −K/3 and Task A's conversion (everything already Einstein-normalised, factor 1):

    dH/dt = −H² − (ρ̄ + 3p̄)/6 − C/3                              (Raychaudhuri with the channel)
    H² = ⟨ρ_E⟩/3 + ⟨Aa²⟩/4 − ⟨R⟩/6                               (averaged Hamiltonian constraint, α-free)

- ⟨ρ_E⟩ here is the *full* local density: it includes S, D, gradient energies and W² factors, not only ρ(f̄).
- So the cosmos lane's dH = −(ρ + p)/2 cannot simply take "−C/3" added to it. Either dH is written in the
  Raychaudhuri form above, or a second difference (the constraint's) is carried with it. **Derived.** This is the main
  practical finding for any production change.
- The other cosmos equations get companions from the same averaging. For any f, with zero shift,
  d⟨f⟩/dt = ⟨∂t f⟩ − K[⟨αf⟩ − ⟨α⟩⟨f⟩]. Applied to Π (the gradient terms integrate away at the wall):

      d⟨Π⟩/dt = K⟨α⟩⟨Π⟩ − ⟨αV_c′(φ)⟩ + β⟨αρ_dm⟩

  against the cosmos lane's −3Hπ − V′(φ̄) + βρ̄_dm. **Derived; not yet evaluated on the stored states.**

## B5. 0* on the candidate

| question | answer |
|---|---|
| derived or borrowed? | **Derived**, from the machine's own K equation and the wall. Its shape is recognisable as a spatial-averaging (backreaction) identity; nothing was taken from that literature. |
| cell size? | **Depends, as 1/r_max³**: C/F = 0.26, 0.085, 0.028 at r_max 16, 24, 36 (ratios 3.0 and 3.1; (3/2)³ = 3.4). It is the pocket's share of the cell. The number of pockets per volume is an unpaid input. |
| closure? | The identity: no. The value of C per unit time: **yes, through the clock**. Lapse part +0.004 F in closure (ii), −0.08 to −0.97 F in closure (i) (r_max 24), because the large cosmos-carried terms are multiplied by ⟨α⟩ ≠ 1. The clock the cosmos lane's t stands for must be named. In closure (ii)'s clock the lapse term on the whole of Q vanishes. |
| averaging measure? | Proper volume is **forced** by the derivation. With S's action J as the weight for the averaged fields the split is unrecognisable: C/F = −0.6 to −27 against +0.08 to +0.09. The total dK/dt does not move; the division into "cosmos part" and "channel" does. C does **not** survive two measures as a number. It survives as an identity with the measure fixed by derivation. |
| slicing? | **CMC only** in this form. K is uniform, so there is no spread-of-expansion term. Another slicing adds one. Forced component (AUDIT §D). |
| other conditions | spherical symmetry; identical neighbours in step (the wall); the one-budget constants in the stored runs (λ = 5000). |

**Verdict: passes as a derived identity; RULER-DEPENDENT as a number.** Specifically:

- dK/dt = ⟨αQ⟩ with the proper-volume average: derived, verified to 2e-6, closure-free.
- "The channel" as one scalar handed to the cosmos lane is well defined only after three things are fixed: the measure
  (proper volume, derived), the clock (to be chosen and labelled), the slicing (CMC, forced).
- Its leading content, the S–D sector's ρ + 3p, is robust: same sign and a steady ratio to F in every run, both
  closures, three cell sizes. That part is better described as a missing species than as an exchange.
- The nonlinear remainder (shear, spread of Π, Jensen) is real, derived, and small here: up to 1.6% of F at
  r_max = 24, 3% at 16.

**Open:** which clock the cosmos lane's t is. Whether "identical neighbours in step" may stand in for a universe of
pockets. The companions in the φ and ρ equations (B4), not yet evaluated. The hot-start and knot behaviour (parked).

## B6. Minimal test (specified; not run; waits for Task A to merge and for review)

Frozen state, no evolution:

1. Take one stored cell state (or a fresh initial slice; no time stepping needed).
2. Build C from the six forms of B3 directly from the fields: Aa; S, P_S, D, P_D; φ, Π; the dust and radiation
   primitives (ρ, v); α. Do not use the engine's assembled ρ + S for this.
3. Convert the averaged fields with Task A's `local_to_cosmos` rule and evaluate the cosmos lane's right-hand side in
   Raychaudhuri form, dH = −H² − (ρ̄ + 3p̄)/6, with and without −C/3.
4. Compare with −K̇/3 from the cell's lapse solve.

Expected: without the channel the cosmos lane misses dH/dt by C/3, positive C (too little deceleration); with it the
two agree to the discretisation of the averages (1e-6 on the stored states). Size at r_max = 24, x_c = 2.0: C = 0.08 F.
Negative controls: J-weighted averages (must fail, by B5); the −(ρ + p)/2 form with −C/3 added (must fail by the
constraint's difference, B4); no-D pocket (C falls to 0.003 F and below).

A few-step version (C supplied step by step, H followed) only after the frozen test. Any change to `cosmos_rhs` is
production and needs Jeremiah's review. At production λ = 0 there is no D vacuum to rest in, so the first production
test needs a state of its own; that choice is open.
