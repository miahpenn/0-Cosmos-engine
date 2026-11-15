# The lane-follow test and the evolved B6 build, made independent

Order item 1 (Nov 14): a short fresh cell run that stores the fluid primitives. RESEARCH ONLY; nothing in `engine/`
changed. Approved as a short run, not a campaign.

- Runs: `research/lattice_cell_run_primitives_v1.py`, x_c = 2.0, cold, closure (ii), r_max = 24 and 36, dr = 0.1,
  **t = 0–24**. (Longer than "a few units" because the lane must start after D settles, at t = 12. About five minutes.)
- Analysis: `research/lane_follow_independent_v1.py` → `output.txt`.

Marks: **Read** throughout; one-budget rule, λ = 5000.

## What the stored primitives change

| | r_max 24 | r_max 36 |
|---|---|---|
| ρ_b/ρ_dm at t = 12: read / solved from totals earlier | 0.151154 / 0.151216 | 0.151145 / 0.151165 |
| ρ_r/ρ_dm at t = 12: read / solved | 1.35e-4 / 1.88e-4 | 1.36e-4 / 1.53e-4 |
| largest dust speed, t = 12 / 24 | 0.031 / 0.043 | 0.031 / 0.043 |
| fluid velocities' part of ρ + 3p, t = 24 | +2.1e-4 | +7e-5 |
| clock rate, all rest mass (α/W) − dark matter only, t = 24 | −5e-5 | −2e-5 |

- The earlier solved baryon density was right to 4e-4. The solved radiation density was wrong by 39% (13% at
  r_max 36), but radiation is 1.4e-4 of the dark matter here, so it moved H by less than 1e-6.
- Velocities, left out before, are a 1e-4 effect.

## B6 on evolved states, fully independent

Every part of C built from the fields and the fluid primitives; none from the assembled ρ + S.

| t | independent build, \|⟨αQ⟩ − K̇\|/K̇ (24 / 36) | miss without channel (24 / 36) |
|---|---|---|
| 3 | 2e-7 / 1e-7 | +1.2% / +0.39% |
| 12 | 4e-8 / 6e-8 | +0.40% / +0.13% |
| 24 | 9e-8 / 7e-9 | +0.17% / +0.055% |

B6's limit ("stored states are only partly independent") is lifted for these states.

## Lane integrated alone from t = 12, initial densities read

Residuals at t = 24, clock = all rest mass:

| lane | r_max | ΔH/H(t0) | Δ ln a | Δρ_dm/ρ_dm | Δρ_b/ρ_b | Δρ_SD/ρ_SD |
|---|---|---|---|---|---|---|
| as it is | 24 / 36 | +6.2e-4 / +2.0e-4 | +1.4e-4 / +4e-5 | −4.0e-4 / −1.2e-4 | −3.9e-4 / −1.2e-4 | — |
| + S–D species, w = 0 | 24 / 36 | +8.5e-5 / +3.0e-5 | −1.3e-4 / −4e-5 | +4.0e-4 / +1.4e-4 | +4.0e-4 / +1.4e-4 | +3.8e-3 / +2.7e-3 |
| + species, the lane's −(ρ+p)/2 form | 24 / 36 | −2.1e-5 / −1.9e-5 | −1.7e-4 / −6e-5 | +5.2e-4 / +1.9e-4 | +5.2e-4 / +1.9e-4 | +4.0e-3 / +2.8e-3 |

- **Confirms the stored-trajectory result.** With the species, H improves sevenfold; at this early time a and the
  densities are no better (their residual is the companion terms either way). The stored run gave 7.1e-5 in H at
  t = 24 with the dark-matter-only clock; here 7.1e-5 with that clock, 8.5e-5 with all rest mass.
- **The two clocks do not differ enough to choose.** All rest mass against dark matter only: 5e-5 in rate, and
  residuals within 20% of each other. Decision 1 (all fluids) stands; the earlier shortcut was harmless.
- **The −(ρ + p)/2 form again does as well as the Raychaudhuri form** once the species is carried (better in H,
  slightly worse in a and the densities).
- Scales with the pocket's share: r_max 36 is about one third of r_max 24 throughout.

Still one-way: the lane feeds nothing back. Still λ = 5000, cold start, one seed.

Unpaid inputs and forced components: as in `results/channel-companions-2026-11-14/COMPANIONS_NOTE.md`.
