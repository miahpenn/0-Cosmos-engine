# The reset-to-lane handoff: separation rule, signed accounting, handoff-time ruler, and a draft contract

Next step after PR #488, as agreed Nov 15: the reset must derive what it hands over. RESEARCH ONLY. Stored states (five
cells: x_c 2.0 at r_max 16, 24, 36; x_c 1.5; closure (i)); **no new runs**; nothing in `engine/` changed.

- `research/reset_handoff_v1.py` → `handoff.txt`
- 2c, derived: `SEA_ANSWER.md`, `research/sea_answer_v1.py` → `sea_answer.txt`; correction filed as
  `results/two-component-lane-2026-11-15/AMENDMENT-2-sea-answer.md`

Marks: **Derived**, **Read**, **Open**.

## Summary

- **Separation rule (Derived):** D's slow part is D averaged over one period of D's ring (2π/√2) and one period of
  the forced ring (π/2), each in local proper time. The rule is ruler-invariant to 2–7% before and after the release,
  but **only to 20–25% during it**, because the release lasts 3–4 time units, about one ring period.
- **What is left "held" after the release is not held step:** it is S pushing D_s down in the pocket (the in-sea
  piece of the EOS note). The contract must measure the held step from the sea's displaced rest value.
- **Signed accounting:** the held → ringing transfer is one-way during the release (positive 91–95% of the time). It
  delivers 1.50–1.59 times the held energy present when the release begins. The extra is expansion's work on a w = −1
  component, paid by gravity.
- **Handoff-time ruler: passes once the step is released.** Lanes handed over at t = 12, 18 and 24 agree to about 1%
  in S–D energy and to a few 1e-5 in H. Handing over at t = 8, with 4–5% still held, misses by about 3%. The
  remaining 1% is the jitter of a single snapshot: the contract must hand over period-averaged quantities.

## 1. The separation rule (Derived; ruler Read)

D's fast part has two lines in local proper time: its own ring at √2 (period 4.44) and the ring forced by S at twice
S's clock (4, period π/2). A boxcar of one period puts a zero exactly on its line.

**Rule:** D_slow = boxcar(4.44/α) ∘ boxcar(1.57/α) applied to D, α the radius's mean lapse.

**Rulers:** a single boxcar; a Hann window of two D periods; no lapse stretching.

Held step / S–D energy, by rule (x_c 2.0, r_max 24; the other cells within 0.01):

| t | derived rule | single boxcar | Hann ×2 | no lapse | spread / value |
|---|---|---|---|---|---|
| 5 | 0.176 | 0.175 | 0.179 | 0.176 | 2% |
| 7 | 0.077 | 0.071 | 0.086 | 0.074 | 20% |
| 8 | 0.049 | 0.045 | 0.057 | 0.046 | 23% |
| 12 | 0.054 | 0.056 | 0.051 | 0.057 | 11% |
| 24 | 0.018 | 0.018 | 0.018 | 0.017 | 6% |

- **Read:** during the release (t ≈ 6–9) the split depends on the rule by up to 25%. The slow part changes as fast as
  D rings, so "slow" and "ringing" are not separable there. That is a property of the release, not of the rule. **The
  reset can hand over cleanly only after its release, not during it.**
- **The residual after the release (1–5%, decaying) sits in the pocket** (r = 1–8; the knot r < 1 adds ≤ 0.5%). It
  is S pushing D_s below its rest value: d² = 2 − (3/2)x gives h/ρ_S ≈ 0.14x, about 4% at t = 12 and 1.4% at t = 24,
  against 3.8% and 1.3% read. This is the vacuum-like piece already inside the in-sea w = (2/3)g − (3/28)x, not held
  step.
- **Refinement for the contract:** measure the held step from the sea's displaced rest value,
  h_held = h(D_slow) − h(d_sea(x)). After the release it then goes to zero instead of 1–5%.

## 2. Signed accounting across the release (Read)

Per proper-volume element, CMC with zero shift (d ln dV/dt = −αK):

    Q = −dE_h/dt − Σ ρ_h α K dV      (held → ringing; the second term is the expansion's work on a w = −1 component)

| cell | Q > 0, t = 4.4–15 | Q > 0, t = 15–60 | ∫Q dt over 4.4–15 / held energy at t = 4.4 | held left at t = 15 |
|---|---|---|---|---|
| x_c 2.0, r_max 16 / 24 / 36 | 95% | 69% / 66% / 69% | 1.58 / 1.59 / 1.59 | 0.12–0.13 |
| x_c 2.0, closure (i) | 95% | 67% | 1.59 | 0.13 |
| x_c 1.5 | 91% | 64% | 1.50 | 0.21 |

- **One-way during the release.** Q is positive 91–95% of the time from the start of the release until the
  held step is spent.
- **The release delivers more than was held: 1.5–1.6 times.** A w = −1 component gains energy as its volume grows, so
  the expansion keeps adding to the held step while it is being released. That extra is paid by gravity. Same in every
  cell and both closures.
- **After the release Q flickers both ways** (a third of the time negative), at about 1e-3 of the S–D energy per unit
  time. That is the in-sea displacement following S's dilution and ringing (§1), not a held step.
- In the S–D sector overall, the energy identity of the companions note holds to 4e-3 (PR #484).

## 3. The handoff-time ruler (Read)

Lane: S–D matter as dust, Raychaudhuri form, the fluids' proper time, started from the cell at t_h. Residual against
the cell; ρ_SD shown.

| cell | t_h = 8 (held 4–5%) | 12 | 18 | 24 | largest difference from t_h = 12, t = 30–66 |
|---|---|---|---|---|---|
| x_c 2.0, r_max 24 | −1.8% to −3.1% | +0.8% to −0.6% | +0.7% to −0.7% | +0.5% to −0.9% | 2.6% (t_h 8); 0.12% (18); 0.37% (24) |
| r_max 16 | −1.6% to −2.7% | +1.1% to −0.03% | +1.0% to −0.1% | +0.5% to −0.6% | 2.7%; 0.14%; 0.64% |
| r_max 36 | −1.8% to −2.8% | +0.75% to −0.33% | +0.64% to −0.43% | +0.48% to −0.58% | 2.5%; 0.11%; 0.27% |
| closure (i) | −1.7% to −2.9% | +0.8% to −0.4% | +0.7% to −0.5% | +0.5% to −0.7% | 2.5%; 0.14%; 0.33% |
| x_c 1.5 | −2.4% to −3.4% | +0.5% to −0.6% | −0.6% to −1.6% | +0.7% to −0.3% | 2.9%; 1.1%; 0.26% |

- **H:** all handoff times agree to a few 1e-5 of H(0) at t = 30–54 (1e-4 at t = 66, through the turnaround).
- **Passes from t_h = 12 on,** to 0.1–1.1% in S–D energy.
- **Fails at t_h = 8 by 2.5–2.9%,** where 4–5% of the S–D energy is still held step: handed over as dust, it is
  missing its negative pressure.
- **The ~1% that remains is snapshot jitter.** ρ_SD in one snapshot oscillates with S's ringing (the x_c 1.5, t_h = 18
  case is the largest). The contract must hand over period-averaged quantities.

## 4. The contract (draft, for Jeremiah and chat)

What the reset hands to the lane, and when:

| item | definition | status |
|---|---|---|
| **when** | once the held step, measured from the sea's displaced rest value, is below the separation rule's own spread (a few %) | handoff ruler passes after that (§3); Derived condition, threshold Read |
| held D | h(D_slow) − h(d_sea(x)), D_slow by the rule of §1 | Derived; clean only after the release |
| ringing D and S | the rest of the S–D energy | Derived |
| S–D energy and pressure | period-averaged (whole S periods, π), proper-volume mean; w = (2/3)g − (3/28)x (0 ± 0.015) | Derived (EOS note) |
| φ, π, ρ_dm, ρ_b, ρ_r | proper-volume means; fluids from their primitives | Read (PR #485) |
| H, a, the clock | H = −K/3; a from proper volume; the fluids' proper time (N, rest-mass weighted) | Derived (PR #481) |
| pocket-local information | not passed: after the release the lane needs none (§3 passes without it) | Read |
| signed transfer | Q over the release, delivered into S–D matter; ∫Q ≈ 1.5–1.6 × the step held at the release's start | Read; one-way |

**Not in the contract yet (Open):**

- **The reset's own derivation of the release.** Its timing and the held region's expansion are pocket-scale (PR #488).
  In a cell the release starts from D's tail out of the surrounding sea. In a cosmos-wide reset there is no
  surrounding sea, and the uniform toy says the seed sets the timing (e-folds ∝ ln(1/seed)): an unpaid input.
- **The two-way test** (parked).

## 0* verdicts

| item | verdict |
|---|---|
| separation rule | Derived (period zeros on both of D's fast lines); RULER-DEPENDENT during the release (to 25%), invariant to ≤ 7% outside it |
| held residual after the release | identified: the in-sea displacement (Derived), not a held step |
| one-way transfer during the release | Read, all cells and both closures |
| release delivers 1.5–1.6 × the held energy | Read; the extra is expansion's work on a w = −1 component (Derived) |
| handoff time | ruler-invariant to ~1% from t = 12 on (Read); limited by snapshot jitter |
| the reset's output as a first-principles initial condition | **not yet**: the release's timing is pocket-scale here and seed-set in a uniform reset |
