# One energy budget: the cosmos lane's floor as D's ejected vacuum (HANDOFF v2 to projects, Nov 10 2026; T1, T2)

Research only; no production change. Scripts: `research/one_budget_cosmos_v1.py` (T1), `research/pocket_with_vacuum_tov_v1.py` (T2).
**Unpaid inputs: λ, V0 (production value in T1; "no V0" under the identity), φ_today, β_D, and D's seed (T1b).**
G_N = 1; G_reset not involved. "Today" is the production present-day state.

## T1a. Static D (D = D\*(φ)) — `T1a_static_D.txt`

| Control | e-folds to turnaround | a_turn/a_today | H0·t_turn | φ_turn |
|---|---|---|---|---|
| engine/cosmos.py (production) | 1.9878761623 | 7.3000 | 5.418340 | 3.100105 |
| (i) F on, D vacuum off | 1.9878761623 | 7.3000 | 5.418340 | 3.100105 |
| (ii) F off, β_D = 0, λ = 1.52458e6 (8πΔV = F) | 1.9878761623 | 7.3000 | 5.418340 | 3.100105 |
| (iii) F on AND D vacuum on, β_D = 0 (double count) | 1.3706635 | 3.9380 | 3.217011 | 2.408665 |
| (iii) the same with β_D = −0.04 | 1.2779004 | 3.5891 | 2.868750 | 2.284549 |
| (iv) F off, β_D = −0.04, λ = 1.56889e6 (matched today) | 1.7386400697 | 5.6896 | 4.260418 | 2.784959 |
| (iv−) F off, β_D = +0.04 (negative control) | 2.3318944112 | 10.2974 | 7.364018 | 3.528811 |

(i) and (ii) reproduce production to all printed digits. The double count (iii) removes a second floor from today's
budget (ρ_today 8.817e-5 instead of 9.641e-5, Ω_m 0.339 instead of 0.310) and brings the turnaround in to 1.37 e-folds.

## T1b. Dynamic D from a seed at the saddle, hand floor off — `T1b_dynamic_D.txt`

- **D falls in ≈ 8 to 27 time units (seed 1e-6 to 1e-14): 2.3 units per decade of seed, logarithmic as expected.**
- **The energy that is not the floor behaves as dust.** Averaged over oscillations w = 0.00 ± 0.02 once settled
  (−0.05 to −0.08 in the first half e-fold after release, while D still visits the top of the well). Its ratio to
  the dark matter stays constant afterwards (it dilutes as a⁻³ and, at β_D = −0.04, picks up the same φ factor).
- **Its size today depends only on when D was released:**

| D released at ln(a/a_today) | H/μ then | D energy above the floor today / ρ_total | / dark matter |
|---|---|---|---|
| 0 (today; seed 1e-6 … 1e-14) | 0.0057 | 0.04–0.06 half an e-fold later | 0.39–0.55 |
| −1.0 | 0.015 | 9.1e-3 | 0.034 |
| −1.5 | 0.030 | 3.3e-3 | 0.012 |
| −2.3 | 0.099 | 1.2e-3 | 0.0044 |

  (β_D = 0; β_D = −0.04 within 3%.) D can first move when H falls to about μ (ln a ≈ −3.8 in production units), which
  extrapolates to a few 1e-5 of the dark matter. At the real scale ratio it is nothing.
- **Effect on the turnaround:** released today, +0.04 to +0.05 e-folds (2.027–2.039 against 1.988); released at
  ln a = −2.3, +0.001 (1.98891). At β_D = −0.04: 1.7334 against the static 1.7386 (the early start uses the
  production lane's past, so 0.005 is the set-up's own mismatch).
- Fragmentation (not run here): inhomogeneous D carries the same energy as travelling waves (GEAR-95, PR #467), which
  are stiffer than dust while relativistic and so dilute faster at first. The homogeneous number is then an upper
  bound on what is left today. Reasoned, not measured.

## T2. Pocket with the one-budget vacuum — `T2_pocket_with_vacuum.txt`

In units of D's step ΔV (so for every λ): matter's pressure is P_S = 3x²/32, and the ejected vacuum outside the
matter has pressure +1.

- **The pockets studied so far cannot be in static balance with the ejected vacuum.** In the hold band (x = 4/3 to
  2.45) P_S = 0.17 to 0.56, against 1 outside. A smooth edge does not help: the total pressure has a minimum (0.143 at
  x = 1.143) and would have to rise outward to 1.
- **A static surface needs P_S = 1, i.e. x = √(32/3) = 3.27 at the surface and more inside: a bag with bag constant
  ΔV.** Such pockets exist (x_c = 4 to 100: R·H_today ≈ 0.46 to 0.75, 2m/R ≈ 0.13 to 0.5, with φ's vacuum from v = 0 up
  to φ ≈ 1). With today's φ vacuum (v = 8.36) most have no surface.
- **But x ≥ 3.27 is where the uniform D mode is not held** (single-harmonic estimate: unstable for 2.45 < x < ~50), so
  whether any D-free pocket survives in the ejected vacuum is open. The engine's inside-seed run is consistent with
  being squeezed once D saturates outside (the x = 4/3 radius moves from 9.6 to 16.6 between t = 45 and 50), but it
  ended at t = 57 before settling this.
- This replaces the earlier flat-space reading (pocket ≈ 2.5 Hubble radii, HANDOFF-REV20-ready 1.6 / R11).

Limits: Thomas–Fermi, static, φ uniform and frozen, β_D = 0, sharp surface (no gradient energy); band statement is
single-harmonic and for the uniform mode only.
