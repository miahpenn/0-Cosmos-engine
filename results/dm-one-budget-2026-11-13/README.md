# Dark matter under the one-budget rule, λ = 5000 (HANDOFF to projects Nov 11 (2), Task 4)

RESEARCH ONLY. No production code changed. Harness `research/dm_one_budget_v1.py`, analysis
`research/dm_one_budget_analysis_v1.py`; full output in `analysis_output.txt`, run logs in `log_*.txt`.

## Unpaid inputs (labelled)

| Input | Value |
|---|---|
| λ | 5000 |
| V0 | production V0 × c², c² = 8π/(2λ·8.2425e-6) = 304.92 → V0 = 2.5133e-2; all present-day densities × c² (H0 = 0.0990, Hubble radius ≈ 10) |
| Floor | D's ejected vacuum at the run's own λ (hand floor off, D_POT = 0) — AMENDMENT-REV20-2 §4 |
| φ_today | 0.179055 (π scaled by √c²) |
| β_D | 0 (default); −0.04 (labelled variant, GEAR-91 coupling via `research/gear91_coupling_v1.py`); β_S = 0 |
| D's seed | 1e-10 inside the pocket |
| Start | cold (D relaxed with friction against the time-averaged S before t = 0) |
| Pocket | S ball, x_c = λA² = 1.5 and 2.0, R = 8, tanh edge width 1 |
| Grid | N = 240, r_max = 24, dr = 0.1, T = 24 (≈ 2.4 e-folds); one resolution only |
| G | G_N = 1 throughout (S at G = 1; dark matter's archive energy divided by 8π). G_reset not involved. |

Controls: **no-D** (D off, hand floor restored at 0.10·c²V0, equal to D's step 8π/(2λ) by the identity, so the
background is the same) and **no-matter** (no S; D uniform at D*·e^{−β_D φ}).
Clock: coordinate time t of the CMC slicing (production gauge); all rates are per unit t, all compared on that clock.

## (a) Dark matter around matter

ρ_dm / ρ_dm(no-matter control), same t:

| run | centre at t = 8 | centre at t = 24 | r = 4 | r = 8 | r = 10 | r = 12 | r ≥ 16 |
|---|---|---|---|---|---|---|---|
| x_c = 1.5, β_D = 0 | 1.063 | 1.213 | 1.250 | 1.097 | 1.022 | 1.004 | 0.995 |
| x_c = 2.0, β_D = 0 | 1.048 | 1.256 | 1.295 | 1.110 | 1.023 | 1.003 | 0.993 |
| x_c = 1.5, β_D = −0.04 | 1.057 | 1.191 | 1.237 | 1.090 | 1.017 | 1.002 | 0.995 |
| x_c = 2.0, β_D = −0.04 | 1.039 | 1.215 | 1.275 | 1.100 | 1.018 | 1.001 | 0.993 |

Dark matter gathers into the pocket: +20–30% inside R by t = 24, still rising, with a slight deficit (−0.5 to −0.7%)
at r ≥ 16. The no-D runs gather about as much (centre/far field at t = 24: 1.173 at x_c = 1.5, 1.266 at x_c = 2.0,
against 1.209 and 1.252 with D), so the gathering is matter's gravity, and D's effect on it is a few per cent with
either sign. β_D = −0.04 lowers it by 1–3 points.

## (b) Rates: dark-matter energy into the pocket vs S's energy exchange with D

Pocket = coordinate r < 11, proper volume. Dark matter: conserved energy minus the no-matter control (inflow beyond the
background). S: E_S(full) − E_S(no-D) (S's own energy, excluding the S–D interaction term). Per unit t.

| run | window | DM energy in (beyond background) | d/dt [E_S(full) − E_S(noD)] | ratio |
|---|---|---|---|---|
| x_c = 1.5, β_D = 0 | 0–4 | +6.4e-3 | +1.9e-3 | 3.4 |
| | 4–12 | +8.1e-3 | +6.9e-4 | 12 |
| | 12–24 | +3.9e-3 | +6.3e-4 | 6.1 |
| x_c = 2.0, β_D = 0 | 0–4 | +7.2e-3 | +1.7e-3 | 4.1 |
| | 4–12 | +1.0e-2 | +2.0e-3 | 5.0 |
| | 12–24 | +4.3e-3 | +1.8e-3 | 2.4 |
| x_c = 1.5, β_D = −0.04 | 0–4 | +6.4e-3 | +2.0e-3 | 3.1 |
| | 4–12 | +7.7e-3 | +8.0e-4 | 9.6 |
| | 12–24 | +3.5e-3 | +7.9e-4 | 4.5 |
| x_c = 2.0, β_D = −0.04 | 0–4 | +7.1e-3 | +1.9e-3 | 3.8 |
| | 4–12 | +9.5e-3 | +2.1e-3 | 4.6 |
| | 12–24 | +3.8e-3 | +1.9e-3 | 2.0 |

E_S at start: 3.98e-2 (x_c = 1.5), 5.77e-2 (x_c = 2.0).

**The sign is the reverse of the question.** In this window S does not lose energy to D; it gains. With D, E_S rises to
1.7–1.9× its start by t = 24; without D it falls to 0.85×. Including the S–D interaction energy (3λ/4)S²D², S's share
rises from 1.3–1.5× to 2.5–2.7× its start. D refills the pocket (D(0)/D* goes from 0.03–0.17 at the cold start to 0.6–1.0 by
t ≈ 8, then oscillating between 0.2 and 1.3), and that raises S's mass. S's energy goes up with its frequency, at the expense of D's ejected vacuum, which
releases 1/(2λ) per unit volume where D falls into its well (≈ 0.2 inside R = 8, several times E_S).

Dark matter's energy enters the pocket 2–12× faster than S's energy changes through D, at both x_c and both β_D. This
is not a balance: the two flows run in the same direction (both into the pocket), and nothing in the runs connects
them beyond gravity.

## Claim status

- **Read:** dark matter gathers +20–30% in the pocket by t = 24; D barely changes this; β_D = −0.04 changes it by 1–3 points.
- **Read:** D refills the pocket under the one-budget rule from a cold start, as in the flat-space bag controls (PR #470).
- **Read:** S gains energy from D in t = 0–24 (both x_c, both β_D); dark-matter inflow exceeds that exchange 2–12×.
- **Derived (mechanism, not computed separately):** S's energy rising with its mass as D enters (the energy/frequency
  ratio of a slowly re-tuned oscillator), paid for from D's released vacuum energy.
- **Open:** whether S turns to losing energy once D has settled (T = 24 covers the refill, not what follows); one
  resolution only (dr = 0.1, no convergence check); the gradient read-out neglects the γ_rr correction; the split of
  the interaction energy between S and D is a convention; the rates use the CMC coordinate clock.
