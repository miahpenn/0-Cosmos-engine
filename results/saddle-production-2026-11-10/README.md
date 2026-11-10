# Does matter make D inside itself? (Nov 10, after AMENDMENT-2-GEAR-98)

`research/sd_saddle_production_v1.py`. Research only. Unpaid input: λ = 5000. Fixed-S checks are flat space with no
gravity (no V0, no G_N role); engine check reads the saved pocket runs (CMC, λ = 5000, V0 = 0). G_reset not involved.

## 1. With D confined to the matter, oscillating S makes D; static ⟨S²⟩ does not (`confined.txt`)

S prescribed from the oscillaton profiles; D = 0 beyond r_cut (Dirichlet), so nothing can flow in from outside; seed
1e-10. Core growth rate of D (free rate 1):

| x_c | r_cut | oscillating S, one ω | oscillating S, local ω(r) | static ⟨S²⟩ = A²/2 |
|---|---|---|---|---|
| 1.0 | 6.0 | +0.006 | +0.005 | +0.001 |
| 1.0 | 9.75 | +0.367 | +0.367 | +0.359 |
| 1.5 | 6.0 | **+0.137** | +0.104 | −0.007 |
| 1.5 | 9.75 | **+0.179** | +0.079 | +0.007 |
| 2.0 | 6.0 | **+0.258** | +0.176 | −0.002 |
| 2.0 | 9.75 | **+0.202** | +0.142 | +0.009 |

- For x_c = 1.5 and 2.0 the matter makes D only when S oscillates: 0.10–0.26 e-folds per unit, zero with static ⟨S²⟩.
  It is the oscillation (the crossings), not the average, that makes it.
- For x_c = 1.0 the core is not held at all (⟨S²⟩ alone lets D grow at 0.36 when the box is wide enough): that is the
  ordinary instability below the hold band, not production by the crossings.

## 2. Mechanism: parametric resonance at S's own frequency (`floquet_finite_k.txt`)

The uniform mode (k = 0) is held for x = 1.5–2.45 (rate 0). Shorter-wavelength modes are not: for k² ≈ 0.6–1.5 the
Floquet rate at ω = ω_loc is 0.16–0.28 (x = 1.5: 0.16 at k² = 1, 0.22 at 1.5; x = 2.0: 0.18–0.28 for k² = 0.6–1.5). These are
the D modes whose own frequency is close to S's: S² pumps at 2ω_S, the main parametric band. The rates match the
confined runs. Produced D oscillates at about S's frequency, above its own mass in the core, so it propagates: it is
emitted, not stored.

## 3. In the engine runs, the matter-made D is soon swamped from outside (`engine_core_frequency.txt`)

Zero-crossing frequency of D at the centre against S at the centre (centre proper time; rough, samples every 0.1):
early on the core D oscillates with S (x_c = 2.0 inside-seed run: D 1.349 = S 1.349 for t = 15–30; broad-seed run:
D = S = 1.024 for t = 2–15), at amplitudes near the 1e-10 seed. Once D outside has grown (t ≳ 30), the core D stops
oscillating and grows monotonically: the tail of outside D. The outside grows at ≈ 1 per unit, five or more times
faster than the matter's own production (0.1–0.26), so in these set-ups the core level is set by inflow.

## Consequences for AMENDMENT-2-GEAR-98

- Its 0.45 e-folds per unit in the x_c = 1.5 core, rising monotonically outward, is the inflow regime, not
  saddle-made D.
- Task 1's 14–16% of D\* in the x_c = 1.5 core is mostly the tail of outside D: static ⟨S²⟩, with no crossings at all,
  gives 0.14 of D\*; oscillation adds about 0.02.
- So a longer run of that set-up would test inflow, not the saddle reading. A test of matter-made D needs the
  outside D not to dominate (D seeded in the resonant band inside, or no growing outside region).

Limits: fixed S (no back-reaction); flat space; engine frequency check is coarse; one λ.
