# Dynamic bag test (HANDOFF Nov 11, Task 2) and chat's two controls

Research only. **Unpaid input: λ = 5000. No V0, φ_today or β_D enters (flat runs: no gravity, no vacuum constant).**
G_N plays no role in the flat runs; G_reset not involved. Engine run with gravity: `research/sd_bag_engine_v1.py` (pending).

## 1. First flat run, sharp edge, closed box (`flat_two_field.txt`, `research/sd_bag_flat_v1.py`)

Both fields evolved. With D ejected outside, the S ball keeps its radius through t = 35 (7.1–7.7 from 8); with no D
outside it has dispersed by t = 20. Under-pressured (x = 2) shrinks to 5.6, over-pressured (x = 5) swells to 9. But D at
the centre reaches D\* within 10–20 units in every case with D outside. The sharp starting edge sent a wave to the centre
at t ≈ 8, and the released energy could not leave: this run cannot separate the causes.

## 2. The two controls (`controls_cold_start.txt`, `research/sd_bag_controls_v1.py`)

Smooth edge (width 1, R = 8); open domain (sponge beyond r = 40); **cold start**: D first relaxed with friction against
the time-averaged S, so the energy released where D falls into its well is removed before the test (GEAR-95's settled
state with its gradient tail). Then, with no friction inside:

| x_in | S time-averaged (no pumping) | S oscillating, prescribed | S and D both evolved |
|---|---|---|---|
| 2.0 (uniform mode held) | D(0) stays 2.5e-2 D\*, flux 0 | D(0) 1–4.5 D\*, outgoing flux 0.2 | D(0) 1.4–2.9 D\*; S energy −70% by t = 79 |
| 3.266 (bag surface) | stays 6.5e-4 D\* | 3–4.7 D\*, flux 0.4 | 1.6–2.7 D\*; S energy −50% |
| 10 | stays 1.6e-8 D\* | 2–4.3 D\*, flux 1.7 | 2.5–4.7 D\*; S energy −95% by t = 40 |

(D(0) is the envelope per 10-unit window, in D\* = √(2/λ). Flux is D's mean outgoing energy flux at R + 3, t ≥ 50.)

1. **Once the released energy is removed, there is no flood:** with time-averaged S, D stays exactly at its settled
   level for 80 units at every x. The hot-start filling (section 1, GEAR-95) was the flood.
2. **With oscillating S, D fills the pocket anyway, at every x, including x = 2 inside the uniform hold band.** That is
   the finite-wavelength parametric growth (research/sd_exact_floquet_v1.py: 0.28 at x = 2, k² ≈ 1.2), emitted
   outward (large outgoing flux when S is held fixed).
3. **With back-reaction, the matter is drained:** S's energy inside the pocket falls 50–95%, converted into D and
   radiated. D in the core does not settle below D\*; it sits at 1.4–4.7 D\*.
4. **Chat's reading 3 (a steady flow of D through intact matter, D inside well below D\*) is not what happens here.**
   The flow exists, but the matter does not stay intact and D inside is not below D\*.

Limits: flat space, no gravity; λ = 5000; one radius (R = 8 ≈ 1.1/√ΔV); no quantum or thermal seeds (the 1e-10 seed and
the stage-1 profile seed D). The cold start is a diagnostic device, not a mechanism.
