# Contraction approach and the fire at H = 0 (Nov 18 2026)

> **SUPERSEDED (Nov 19 2026) as numbers:** this run used #523's core S equation (S mass^2 = 4 at d = 0). The stone toy has S mass^2 = 1 + 1.5 d^2; with it the gains are larger and the 'needs a seed' contraction result changes. See results/toy-s-rerun-2026-11-19/. Verdicts that survive are listed there.


Script `research/contraction_crossing_v1.py` (repo root; ~6 min), full table `summary.txt`. Design and readings: coordination log, "DRAFT DESIGN: contraction approach" with chat's edits (A)-(E). Planar periodic comoving lattice (as the Nov 18 spatial run), H(t) = Hdot t through zero, friction 3H on S and D, a(0) = 1, k_phys = k/a. S starts at x = 0.3 with D at rest (d^2 = 2 - 1.5x) plus noise; t_i is solved (homogeneous 2-DOF) so that x(0) = x_c. Fire rule (set before running): e_D/e_S first >= 0.10 at t10; "way up" if t10 <= -20, "at crossing" if -20 < t10 <= 0, "after" if t10 > 0, "none" otherwise. Research only.

## Result (Measured)
Outcome with the only named seed that needs no luck, d0 = 1e-8 (x_c = 1.4 / 2.0 / 2.6):
| Hdot | 3e-3 | 1e-3 | 3e-4 | 1e-4 | 1e-5 | 4e-7 (machine order) |
|---|---|---|---|---|---|---|
| outcome | none / none / none | none / none / after | none / none / none | none / none / none | way up (t10 -22 / -109 / -81) | way up (t10 -811 / -1056 / -1265) |

- **At the machine's own turning rate the pump, not H = 0, is the trigger.** For Hdot <= 1e-5, every row fires on the approach for every seed, including 1e-8 (gains 1e8 saturate): t10 from -22 (Hdot 1e-5) to about -1300 (4e-7). By H = 0 the energy in D is 0.3-1.7 of S's remaining ring energy (the ratio exceeds 1 because S is drained).
- **At bounce-rate crossings (Hdot 1e-4 to 3e-3) the answer depends on the seed.** d0 = 1e-8: no fire in 17 of 18 rows (one "after" at 1e-3, x_c 2.6). The required seed d0* = 0.3/G_k(0) is: 1e-4/2.0: 3e-4; 3e-4/2.0: 6e-3; 1e-3/2.0: 0.02; 3e-3/2.0: 0.04; so a crossing at Hdot >= 3e-4 needs a seed of 1e-3-0.1, many orders above 1e-8. At d0 = 1e-2: "at crossing" (1e-4, 2.0), "way up" (1e-4, 2.6), "after" for the others.
- **Estimate check (chat's E_b = 1.75 Gamma sqrt(Gamma/(1.23 Hdot))):** the measured crossing exponent ln G over +-2 tau_w matches in 11 of the 12 non-saturating rows within the tolerance max(1, 0.3E): e.g. 3.59 vs 3.6 (1e-3, 2.0), 6.20 vs 6.5 (3e-4, 2.0), 11.35 vs 11.3 (1e-4, 2.0). The miss is 1e-4/x_c 1.4 (9.9 vs 6.6; its Gamma = 0.12 is an extrapolation).
- **Band entry:** x enters 4/3 at H = -0.0125, -0.027, -0.032 (Hdot 3e-3), -0.0076, -0.0155, -0.0207 (1e-3), -0.0041, -0.0095, -0.0118 (3e-4), -0.0026, -0.0055, -0.0070 (1e-4) for x_c 1.4/2.0/2.6. For Hdot <= 1e-5 the fire comes first: x is below 4/3 (about 0.8-0.9 at 4e-7) when D is already carrying the energy, and S is drained before it reaches the band (no x entry is recorded for 4e-7).

## What this says and does not say
- It says: if the machine's turning rate is as slow as chat's estimate (4e-7), the pump fires hundreds to a thousand time units BEFORE the bounce, from a state below the hold band, for any seed. The bounce is then not where S fires. If the turn is as fast as 1e-4 or more, it fires only from a seed far above the sea's noise.
- It does not test what brings H to 0 (prescribed H), nor stones (planar, uniform S).

## Limits and one open question (read first)
- **Below-band regime rests on the core's S equation.** #523's core has S mass^2 = 4 at d = 0; the stone toy has S mass^2 = 1 + 1.5 d^2 (4 only when D is at the sea). For x < 4/3 and the approach, these give different S frequencies (core 2.5 at x = 0.3, toy ~1.8), so the early fire at x ~ 0.8-0.9 may shift. I did not run the toy's S equation. The in-band numbers are the core's and agree with E_b.
- x is measured from mean(s^2); entry/exit times are at the 4/3 proxy and a smoothed window, not exact.
- Flat noise on k <= 4; one seed; planar; no gravity; constant dH/dt through zero is a stand-in for the real bounce. The ratio e_D/e_S uses total D energy (2 x mean kinetic) over S's current ring energy, so it can exceed 1.
- Hdot = 4e-7 comes from chat's estimate (stand-ins unchecked); the rows are only as good as it.
