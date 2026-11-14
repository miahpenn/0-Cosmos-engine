# HANDOFF Nov 14 2026 (3): GEAR-105 and GEAR-103 field checks, the settling clock, the static collapse, the spike

RESEARCH ONLY. No production code changed. No new expanding gravity runs (decision of Nov 14: no interim edge).
The dust-absorbing layer of PR #477 was not used for any result here; it stays a labelled research diagnostic.

Scripts: `research/d_ring_box_v1.py`, `d_ring_box_report_v1.py`, `d_ring_floquet_fine_v1.py` (Task A);
`sd_settling_clock_v1.py` (B); `static_collapse_analysis_v1.py` with `dm_one_budget_v1.py` cases `static_full` /
`static_noD` (C); `gear103_ring_coupling_v1.py` (D); `spike_from_disk_v1.py` (E). Outputs: `task*.txt`.

## Unpaid inputs and forced components (AUDIT §D)

| task | unpaid inputs | forced / imposed |
|---|---|---|
| A | λ (scales out; run at 5000), ring size ε, seed | the uniform ring placed by hand; the seed's size and spectrum |
| B | λ = 5000, x_in = 2.0, 3.27, seed 1e-10 | the pocket placed by hand; the cold start (artificial settle); the absorbing sponge at r > 40 |
| C | λ = 5000, x_c = 2.0, seed 1e-10, β_D = β_S = 0 | φ frozen and D's vacuum cancelled by a constant (no φ / no vacuum); a 1e-12 dust trace; cold start; maximal slicing with zero shift; lapse pinned to 1 at the outer cell |
| D | λ (scales out), β_D = ∓0.04 (value anchored by observation), ε | φ prescribed as a uniform ramp 0 → 5 (not evolved) |
| E | as PR #475–#477 (one-budget, ×304.92) | everything in those runs, and the production outer edge: PROVISIONAL |

Flat runs are "no φ / no vacuum". G_N = 1 in Task C; G_reset not involved.

## Task A. GEAR-105 field check

Full nonlinear S, D equations, 1D periodic box; D = D*(1 + ε) at rest, seed 1e-8 D* per box mode (k < 4) with random
phases. Rulers: (L = 800, dr = 0.1), (L = 800, dr = 0.05), (L = 400, dr = 0.1, different random seed).

**(a) Matter-wave seed: no growth.** Largest amplification of any S mode over the fit window (envelope, last quarter
over first quarter):

| ε | L 800, dr 0.1 | L 800, dr 0.05 | L 400, dr 0.1 |
|---|---|---|---|
| 0 (control) | 1.002 | 1.000 | 1.002 |
| 0.05 | 1.003 | 1.002 | 1.001 |
| 0.10 | 1.011 | 1.010 | 1.011 |
| 0.20 | 1.037 | 1.036 | 1.037 |

Fitted rates are ≤ 3e-4 (ε = 0.2), ≤ 1.2e-5 (0.1), ≤ 1e-6 (0.05), against D-wave rates 170× larger at the same ε.
The few-per-cent amplification is the envelope's bounded wobble (it scales as ε² and does not accumulate).

**(b) D-wave seed: one long band.**

| ε | peak rate (field; three rulers) | rate/ε² | Floquet, fine k grid | GEAR-105 table | peak k²/ε² (field) | band top k (field / Floquet) |
|---|---|---|---|---|---|---|
| 0.05 | 2.782e-3, 2.782e-3, 2.756e-3 | 1.113, 1.113, 1.102 | 1.114 | 1.07 | 2.99–3.55 | 0.110–0.118 / 0.125 |
| 0.10 | 1.167e-2, 1.167e-2, 1.165e-2 | 1.167, 1.167, 1.165 | 1.167 | 1.16 | 2.99–3.26 | 0.251–0.259 / 0.253 |
| 0.20 | 5.129e-2, 5.129e-2, 5.125e-2 | 1.282, 1.282, 1.281 | 1.283 | 1.28 | 3.12–3.26 | 0.518 / 0.506 |

(The L = 400 values differ in the third digit because the box has half as many modes in the band. The field band top
is where the rate falls to 10% of the peak, so it sits slightly inside the Floquet edge at ε = 0.05.)

- **Correction to GEAR-105 (for chat's amendment file):** the table's 1.07 at ε = 0.05 comes from scanning k² in steps
  of 0.01, as wide as the band there. On a fine grid the Floquet value is 1.114, and the field runs give 1.113. The
  small-ring coefficient is ≈ 1.11, rising to 1.17 and 1.28.
- Peak at k² ≈ 3.1–3.2 ε²; the band ends at k² ≈ 6.3–6.4 ε² (Floquet, fine grid).

**(c) Nonlinear end state** (seed 1e-6 D* in both S and D; run to 30/(1.1 ε²)). Shares of the energy above the vacuum,
mean over the last 20% of each run (three rulers):

| ε | uniform ring | long D waves (0 < k < 0.5) | short D waves (k ≥ 0.5) | S | mean k of the D waves |
|---|---|---|---|---|---|
| 0.05 | 0.22, 0.22, 0.34 | 0.78, 0.78, 0.67 | 0.000 | 0.000 | 1.1–1.3 ε |
| 0.10 | 0.14, 0.14, 0.14 | 0.86, 0.86, 0.86 | 0.003–0.004 | 0.000 | 1.0–1.05 ε |
| 0.20 | 0.09, 0.09, 0.11 | 0.86, 0.86, 0.84 | 0.073–0.076 | 0.000 | 1.06–1.08 ε |

- **Read (invariant):** the ring's energy ends in long D waves; nothing reaches S; almost nothing reaches short waves
  for small rings.
- **Read (range only):** the uniform ring does not vanish. It drops to 6–10% at the first breakup and then partly
  returns and goes again (e.g. 0.10 → 0.66 → 0.24 at ε = 0.05, L = 400). The late share depends on box and time.
- The waves that carry the energy late are longer (k ≈ 1.0–1.3 ε) than the fastest-growing ones (k ≈ 1.8 ε).

## Task B. The settling clock at field level

Spherical pocket (R = 8), flat, absorbing sponge; hot start (D released, flood present) and cold start; T = 80.
S's saddle crossings = zero crossings of S projected on its initial profile. Three pumping measures, pocket-integrated
(the whole-domain values, printed beside them in `sd_settling_clock_v1.py`'s report, differ by less than 0.1 in all but a few cycles):
P1 = rate of energy into D's own sector; P2 = rate of energy out of S's own sector; P3 = d/dt of D's kinetic energy.
Lead = t(crossing) − t(peak) per half-cycle (half-period ≈ 1.9–2.1). Both grid spacings (0.04, 0.02) give the same
P1 and P2 leads to ±0.01; P3 agrees within 0.12 except one cycle (x_in = 3.27, hot, cycle 2).

| run | measure | leads, cycles 1–10 |
|---|---|---|
| x_in 2.0, hot | P1 | +0.64 +0.58 +0.51 +0.48 +0.50 +0.84 +0.79 +0.68 −0.46 +0.51 |
| x_in 2.0, hot | P3 | +0.71 +0.68 −0.54 +0.55 +0.50 +0.93 +0.08 +0.45 −0.43 +0.51 |
| x_in 2.0, cold | P1 | +0.80 +0.77 +0.82 +0.71 −0.93 −0.87 +0.85 +0.88 +0.54 +1.06 |
| x_in 2.0, cold | P3 | −0.03 +0.00 −0.02 +0.07 −0.93 −0.87 +0.85 +0.88 +0.93 +1.02 |
| x_in 2.0, cold | P2 | −0.13 +0.10 −0.24 −0.15 −0.58 −0.62 −0.27 −0.85 +0.83 −1.06 |
| x_in 3.27, hot | P1 | +0.74 −1.00 +0.89 +0.69 +0.75 −0.93 −0.99 −0.91 +0.98 −0.82 |
| x_in 3.27, cold | P1 | +0.75 +0.80 −0.99 −0.96 −0.91 +0.73 +0.73 −0.42 −1.00 +0.81 |

- **Read (two grid spacings):** hot start at x_in = 2.0: the pumping (P1) leads the crossings and the lead shrinks over
  the first four cycles, 0.64 → 0.58 → 0.51 → 0.48, as GEAR 02b predicts for the violent phase. It then stops
  shrinking and scatters.
- **Read:** cold start at x_in = 2.0: D's kinetic-energy rate (P3) is locked on the crossings for the first four
  cycles (|lead| ≤ 0.07) while D is still out of the pocket, and the lock is lost when D floods in (cycle 5 on).
- **Contradicts the prediction:** there is no late lock onto the crossings in any run. From cycle 5 the peaks of P1 and
  P3 sit mostly at S's turning points (lead ≈ ±0.9; P1: 45–91% of cycles 5–15, P3: 36–91%), in hot and cold starts alike. P1 carries a
  factor S², which is largest at the turning points.
- **Contradicts the prediction:** at x_in = 3.27 there is neither a shrinking lead nor a lock at any stage.
- **Convention:** the answer depends on which rate is called the pumping. P2 (energy leaving S) trails the crossings
  by 0.1–0.9 throughout.
- One pocket, one geometry (sphere, absorbing edge): the row and the reflecting edge were not run.

## Task C. The static-gravity collapse with D

Pocket x_c = 2.0 with cold-start D on a static background (no φ / no vacuum), maximal slicing, zero shift. Four runs:
dr = 0.1 at r_max = 36 and 24 (to t = 40); dr = 0.075 at r_max = 24 (to t = 33.7) and at r_max = 36 (reached t = 27.7
inside the 2-hour run limit, before the horizon). Control: no D, dr = 0.1, r_max = 36. Read-outs are
slicing-independent: χ = (∇R)² = 1 − 2m/R, the Kodama clock N_K = √|χ|, the null expansions θ± = (2/R)(U ± Γ).

**It is a collapse (read, invariant under grid spacing and box size).** A trapped surface (θ+ < 0 and θ− < 0) forms:

| run | first trapped surface, t (centre's proper time) | apparent horizon: areal radius at first sight → later | mass inside (R/2) | max 2m/R at t = 0 |
|---|---|---|---|---|
| dr 0.1, r_max 36 | 30.5 (15.75) | 1.90 → 2.17 (t = 40) | 0.95 → 1.08 | 0.228 |
| dr 0.1, r_max 24 | 29.5 (15.74) | 1.78 → 2.17 (t = 40) | 0.89 → 1.08 | 0.228 |
| dr 0.075, r_max 24 | 29.7 (15.73) | 1.83 → 2.03 (t = 32.2) | 0.91 → 1.02 | 0.228 |
| dr 0.075, r_max 36 (to t = 27.7) | not reached; min χ = 0.25 at t = 26.2, as the others at that time (0.21–0.26) | — | — | 0.228 |
| no D, dr 0.1, r_max 36 | none; max 2m/R falls 0.144 → 0.045 (S disperses) | — | — | 0.144 |

- The horizon takes in nearly all of the pocket's mass (1.08 of ≈ 1.02 at the start, in G = 1 units per the
  Misner–Sharp mass).
- With D the pocket starts 60% more compact than without (2m/R 0.228 against 0.144): the interior, where D is held
  out, carries the vacuum step. In this static set-up that step is positive energy relative to the outside, because
  D's vacuum is cancelled by a constant (FORCED).
- **Clocks.** Before the horizon the Kodama clock's minimum falls 0.88 → 0.10–0.24 (t = 0 → 28–30) while the central
  lapse falls 1.4 → 0.03–0.05. The lapse goes on to 0.0006–0.002 by t = 40; that further fall is the slicing halting
  inside the horizon (the centre's proper time stops at 15.83), not an invariant statement.
- **The edge here.** The lapse is pinned to 1 at the outer cell. Two units inside the edge it drifts by +0.7% (r_max 36)
  and +1.8% (r_max 24) over the run; K at the edge reaches 7e-3 and 2.8e-2; S arriving at the edge is ≤ 7e-5. There is
  no dust to pile up. The horizon time differs by 1.0 between the two boxes, so the edge is felt at the few-per-cent
  level but does not make the collapse.

**(a) Is the core clear of D at rest and full of rung D? No (read, three runs agree).** Core r < 3, windows of 4:

| window | D at rest ⟨D⟩/D* | ⟨D²⟩/D*² | ring energy: slow / travelling |
|---|---|---|---|
| 0–8 | 0.01–0.02 | 0.001 | 1e-6 / 1e-6 |
| 12–16 | 0.05–0.06 | 0.05 | 5–8e-6 / 1.2–1.6e-4 |
| 16–20 | 0.10–0.17 | 0.43–0.57 | 1.4–1.8e-5 / 1.1–1.5e-3 |
| 20–24 | 0.44–0.49 | 1.4–1.5 | 7–8e-5 / 4.0–4.8e-3 |
| 24–28 | 0.73–0.78 | 2.3–2.6 | 5–6e-5 / 7.1–7.9e-3 |
| 28–32 | 0.92–1.06 | 2.3 | 2–3e-5 / 6.6–6.7e-3 |

The core starts clear of D altogether. D first arrives as ring (rung D with almost no rest value, t = 12–20), and
then D at rest fills in to ≈ D* by the time the horizon forms, with about as much again in ring (⟨D²⟩ ≈ 2.3 D*²).
The ring in the core is all at k ≥ 0.5 (the core is smaller than a slow-ring wavelength, so the slow/travelling split
says nothing here).

**(b) The hammer: yes (read).** Across the pocket (r < 10) the ring's strength ⟨δD²⟩(r) follows ⟨S²⟩(r) once D is in:
correlation of the logs +0.82 to +0.95 and slope 0.94–1.35 over t = 12–32, in all runs. Before D enters (t < 8) there
is no relation (−0.2), and after the horizon it is lost (the interior stops evolving in this slicing).

**(c) Trapped, not made in place (read).** D's energy flux through r = 2 is inward throughout the collapse and grows
(−1e-5 at t = 8–12 to −3e-3 at t = 20–28, per unit coordinate time); through r = 4 it is inward too (−2e-5 to −1e-4).
The core's D is carried in from outside. At r = 8 the flux turns weakly outward (+1e-6 to +8e-6) after t = 12. As the
lapse falls the inward flux at r = 2 drops per unit coordinate time (to −3e-5 by t = 36–40) but not per unit local
proper time (−3e-3 to −8e-3): the shut-off is the slicing's. Nothing flows out through r = 2 after t = 12.

**(d) The freeze in local proper time (read).** At the centre, S's proper-time rate √(⟨S_τ²⟩/⟨S²⟩) rises from 1.3
(t < 8) to 5–6 (t = 20–24) and 8–12 (t = 24–32): per unit of its own time the core rings faster as it collapses; it
does not freeze. The coordinate freeze is the lapse. The centre's proper time stops advancing at 15.83 in this
slicing, so the run cannot say whether the core releases later (GEAR-94 §4): maximal slicing does not reach that
part of the interior.

**Status.** Collapse with a trapped surface: read, two grid spacings and two boxes. Still one pocket (x_c = 2.0), one
start, and a static background that is itself forced (vacuum cancelled by hand, φ frozen). What happens inside the
horizon, and whether a pocket in a self-consistent background collapses, are open.

## Task D. GEAR-103 check

Uniform D ring (one ODE), S = 0, V = −e^{−2β_D φ} D²/2 + (λ/8) D⁴, φ = 0 → 5 ramped at speed v. Rulers: v = 0.005 and
0.0025; integrator tolerance 1e-11 and 1e-8 (all three agree to the digits shown).

| β_D | ε | slope d ln ω / dφ (prediction −β_D) | (ring part of ∂V/∂φ) / ρ_ring at start / end (prediction −β_D) | ρ_ring end/start |
|---|---|---|---|---|
| −0.04 | 0.05 | +0.04018 | +0.04012 / +0.04007 | 1.220 |
| −0.04 | 0.10 | +0.04076 | +0.04051 / +0.04028 | 1.222 |
| −0.04 | 0.20 | +0.04371 | +0.04244 / +0.04128 | 1.231 |
| +0.04 (control) | 0.05 | −0.04032 | −0.04011 / −0.04021 | 0.819 |
| +0.04 | 0.10 | −0.04141 | −0.04051 / −0.04094 | 0.817 |
| +0.04 | 0.20 | −0.04751 | −0.04244 / −0.04483 | 0.806 |
| 0 (control) | all | 0.00000 | 0.00000 | 1.000 |

- **Read (invariant):** the ring's frequency follows √2·e^{−β_D φ}; the departures are the ring's own anharmonic shift
  (∝ ε²: 0.5%, 2%, 9–19% at ε = 0.05, 0.1, 0.2), which changes as the ring's relative size changes along the ramp.
- **Read (invariant):** the ring's part of ∂V/∂φ equals −β_D ρ_ring within 0.3% at ε = 0.05. With GEAR-91's
  φ'' = … − 8π ∂V/∂φ, the ring's source in φ's equation is +8π β_D ρ_ring (sign convention: source on the right-hand
  side; in the handoff's wording, ∂V/∂φ|ring = −β_D ρ_ring).
- **Read:** the ring's energy over the ramp scales as e^{−β_D Δφ} (1.221 predicted for β_D = −0.04, Δφ = 5; measured
  1.220): energy follows the mass at fixed action.
- Not tested: φ's own response (φ is prescribed), gravity, a non-uniform ring.

## Task E. The spike, from existing output only (PROVISIONAL)

No new runs. ρ_dm(innermost cell)/ρ_dm(r = 2):

| run | t = 24 | 28 | 32 | 36 | 40 | spike half-width | S's core (half of ⟨S²⟩) | ⟨D²⟩(0)/D*² |
|---|---|---|---|---|---|---|---|---|
| x_c 2.0, D, dr 0.1, r_max 36 | 1.02 | 1.07 | 1.20 | 1.82 | 3.50 | 0.15 | 0.35 | 1.21 |
| x_c 2.0, D, dr 0.075, r_max 36 | 1.02 | 1.07 | 1.21 | 1.89 | — | 0.19 | 0.34 | 1.18 |
| x_c 2.0, D, dr 0.05, r_max 36 (to t = 32) | 1.02 | 1.07 | 1.24 | — | — | 0.28 | 0.43 | 0.82 |
| x_c 2.0, D, dr 0.1, r_max 24 | 1.02 | 1.07 | 1.21 | 1.92 | 4.20 | 0.15 | 0.35 | 1.11 |
| x_c 2.0, no D (r_max 36 and 24) | 1.03 | 1.04 | 1.04 | 1.04 | 1.04 | none | 6.3–6.8 | 0 |
| x_c 1.5, D, r_max 36 | 1.01 | 1.02 | 1.02 | 1.03 | 1.04 | none | 1.05 | 0.96 |
| x_c 1.5, no D, r_max 36 | 1.02 | 1.02 | 1.02 | 1.03 | 1.03 | none | 7.2 | 0 |

- **Read:** no spike without D, and none at x_c = 1.5 with D (to t = 40).
- **Read:** where the spike appears, S itself has contracted to a core of radius ≈ 0.35–0.43 (it stays spread over
  r ≈ 6–7 without D), and the spike's half-width (0.15–0.28) lies inside that core. The spike is dark matter falling
  into S's own central concentration.
- At x_c = 1.5 with D, S's core has already narrowed to r ≈ 1.05 by t = 40 without a spike yet.
- **Note for the map:** the dust species has no pressure, so once the centre deepens nothing in the model stops dust
  focusing there; central caustics are the model's own behaviour, not a numerical artefact.
- Half-widths at the last snapshot are 1.5–5.5 cells; the spike's growth is the resolved quantity (PR #477).

## Waiting on the derived edge (not started)

Rerun of the expanding one-budget results to T = 60–80 with whole-domain action; the spike after t = 40; GEAR-105's
ring on an expanding background; the closed row of pockets with gravity and expansion.

**Design note (closed row of pockets, no runs).** The spherical machine cannot hold a row, but it can hold one cell of
a lattice: a pocket at the centre of a spherical cell whose outer wall is a symmetry surface (every field even across
it, every radial flux zero, including the dust's). That wall is closed (nothing leaves, so whole-domain action is
well defined), it has no free boundary data, and on an expanding background it moves with the dust. It may be the
same object as the background-following edge chat is deriving; if so, one spec would serve both the edge and the
closed box. Offered as a question for the derivation, not as a spec.
