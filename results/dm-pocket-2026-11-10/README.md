# Dark-matter pocket experiment (HANDOFF-to-projects Nov 10, Task 3)

`research/dm_pocket_v1.py` (runs), `research/dm_pocket_analysis_v1.py` (table: `dm_density_differences.txt`).
Research only. **Unpaid inputs: λ = 5000 (GEAR-67 3(e)); V0 and φ_today at production values ("at production V0",
GEAR-92); β_D = β = −0.04, β_S = 0 (the handoff's defaults; the physics choice is Jeremiah's).** G_N = 1; G_reset not
involved. Densities in machine units: direction and mechanism only, not size.

**Set-up.** Production kernel (CMC, β = 0), no wall, N = 240, r_max = 60, T = 45. Uniform dark matter at
ρ_DM,present (archive-normalised 2.59e-5), φ uniform at φ_today with Π_today, V0 on; a self-held S oscillaton
(x_c = 1.0, 1.5, 2.0; S at its turning point) with a D seed 1e-10 inside (width 6). Time-symmetric slice (K = 0, all
sources in the Hamiltonian constraint): **not today's expanding background**. Dark matter couples to φ (β) in every run.

| Run | GEAR-91 | Isolates |
|---|---|---|
| off | off | gravity only |
| on | β_S = 0, β_D = β | total effect of the coupling |
| nosource | on, but without φ's 8π ∂V_SD/∂φ | D's φ-dependent mass, with nothing flowing to φ |

So on − nosource is the chain matter → D → φ → dark matter (matter's frequency acting through D on φ), and
nosource − off is the change in D's own dynamics (its tachyonic mass is e^(−2β_D φ) ≈ 1.014 at φ_today, so D grows
slightly faster) acting on dark matter through gravity.

## Result at t = 44.6 (relative to the off run; centre r = 0.1 / r = 20 / r = 30)

| x_c | gravity: ρ_DM(centre) / initial | total (on − off) | φ channel (on − nosource) | D-mass channel (nosource − off) |
|---|---|---|---|---|
| 1.0 | × 5.0 | −2.5% / +0.5% / +0.7% | −0.08% / +0.01% / +0.03% | −2.5% / +0.5% / +0.6% |
| 1.5 | × 17 | −5.1% / +0.5% / +0.7% | −0.15% / +0.01% / +0.03% | −4.9% / +0.5% / +0.6% |
| 2.0 | × 85 | −10.5% / +0.5% / +0.7% | −0.32% / +0.01% / +0.04% | −10.2% / +0.5% / +0.7% |

1. **Gravity gathers dark matter into the pocket**, more for heavier pockets (× 5, 17, 85 at the centre by t = 45).
   Part of the late growth is the whole domain contracting once the ejected D region forms (its negative vacuum).
2. **The coupling hollows the centre slightly and gathers dark matter outside**, at every mass. It switches on only
   once D has grown (differences are round-off, ≤ 1e-11, before t ≈ 30).
3. **About 97% of the coupling's effect is the D-mass channel**, not the φ channel: shifting D's mass changes when the
   ejected region forms, and its gravity moves the dark matter.
4. **The φ channel (matter's frequency, through D, on φ, on dark matter) has the same sign**: the centre hollows
   (φ is lowered there by ≈ 8e-5 and raised outside by ≈ 8e-3, where D sits in its well), outside gathers.
5. **No mass-dependent sign flip between x_c = 1.0 and 2.0.** Both channels keep their sign; their size grows with the
   pocket's mass (centre −0.08%, −0.15%, −0.32% for the φ channel).

**Logged, not interpreted (AMENDMENT-3-GEAR-91 §2):** g/κ at the centre with β_S = 0 follows φ: −139.6 at t = 0
(φ = 0.179), −101 at t = 19.8, −72 at t = 39.6 (φ ≈ 0.36). κ < 0 throughout.

## Limits

Time-symmetric start, not an expanding background; production V0 (39–40 orders off if S is tangible matter, GEAR-92);
the ejected D region and its negative vacuum dominate after t ≈ 35, so only the direction of the late effects is read;
one resolution; x_c = 1.0 is outside the stable hold band (its core does not hold D), 1.5 and 2.0 inside.
