# Two-component closure in a research copy of the cosmos lane

HANDOFF Nov 15 2026 (2), Task 1, run on the basis agreed that day: the held step's fading is the reset's job; the lane
starts from the reset's output, after settling. RESEARCH ONLY. Stored states and the fresh runs of PR #485 (they
store the fluid primitives); **no new runs**; nothing in `engine/` changed.

- `research/two_component_lane_v1.py`: the lane, the release law, the controls → `release_law.txt`, `lane_follow.txt`
- Task 2 and the pokes: `LEFTOVER_NOTE.md`, `research/poke_saddle_v1.py` → `pokes.txt`

Marks: **Derived**, **Read**, **Open**. One-budget rule, λ = 5000, cold start, x_c = 2.0 unless stated.

## Summary

- **After settling, the closure works.** S–D matter carried as dust follows the cell to 1e-4 in H by t = 24 and to
  1e-3 by t = 66, through the turnaround (2e-3 in the smallest cell). Residuals scale with the pocket's share. Two rulers: stored cells with
  densities solved, fresh cells with densities read.
- **The held step's fading is derived in form, but it is pocket-scale.** The law gets the step's share right (19%
  of the S–D energy) but releases 15–50% late. The held region also gains energy at about half the w = −1 rate,
  because it sits inside the overdense pocket. Both need the pocket's profile, which the lane does not carry.
- **From t = 0 every lane fails on the S–D energy**, in the way each control was predicted to fail. H and a stay at
  1e-3 regardless, so they do not discriminate.
- **0\* verdict:** the held step is reset-scale (sub-grid) physics. It belongs with the reset (C6, POKE-1), not in the
  lane. The lane's initial state is the reset's output.

## 1a. The two components, and the fading law

Lane: Raychaudhuri form, the fluids' proper time, with

    rho + 3p = 2π² − 2V + ρ_dm + ρ_b + 2ρ_r − 2ρ_vac + ρ_m − 2ρ_h
    ρ_m: S–D matter, w = 0;   ρ_h: the held step, w = −1;   dρ_m/dτ = −3Hρ_m + Q,   Q = −dρ_h/dτ ≥ 0

**The fading law (Derived).** Where S holds D_s out, D's slow part obeys, averaged over S's period (EOS note, §2):

    d̈ + 3H ḋ = d (1 − (3/4)x) − d³/2,     h(d) = ½ (1 − d²/2)²   (scaled; the held energy)

- **Holding.** While S's local concentration x is above 4/3 (the averaged band edge; Floquet gives 1.27), D's
  uniform mode is held.
- **Release.** Below it, D grows at γ = √(1 − 3x/4) and the step's energy goes into S–D matter, one way.
- **What drives it** is S's concentration falling out of the band: the local expansion dilutes it as a⁻³.
- **Where D starts.** D grows from the initial slice's tail inside the pocket (2.5% of D* at the centre), not from
  the 1e-10 seed.
- **Inputs**, read from the cell at t0 shell by shell: D, its rate, x, and the volume weights. Nothing fitted.

**Read: what the cell does.**

- **x leaves the band early.** S's local concentration leaves the band everywhere by t ≈ 2.5 (x_c 2.0) and before
  t = 2 (x_c 1.5).
- **D lets go about five time units later.** At the centre it grows from 0.02 to 0.7 D* between t = 2 and 8, doubling
  each unit (γ ≈ 0.69; the law's √(1 − 3x/4) gives 0.71–0.79 there).
- **The held edge sweeps inward,** reaching the centre at t ≈ 7.5 (x_c 2.0) and t ≈ 5 (x_c 1.5).

**Read: law against cell.**

| cell | step's share of S–D energy, t = 0.5: cell / law | half-release time: cell / law |
|---|---|---|
| x_c 2.0, r_max 24 (fresh) | 0.197 / 0.189 | 3.85 / 5.85 |
| x_c 2.0, r_max 36 (fresh) | 0.198 / 0.189 | 3.86 / 5.88 |
| x_c 2.0, r_max 16 | 0.197 / 0.189 | 3.84 / 5.75 |
| x_c 1.5, r_max 24 | 0.196 / 0.187 | 3.19 / 3.66 |

- **Without D's gradient the law is late,** by 15% (x_c 1.5) and 50% (x_c 2.0).
- **With the gradient it is early, then unstable.** At t = 5 it leaves 14% of the step against the cell's 30%, then
  goes unstable at late times: D's invasion from the outer sea is a real part of the release.
- **The held region expands more slowly than the cell (Read).** Comoving, the cell's held energy grows from 0.19 to
  0.29 of the initial S–D energy by t = 3. A w = −1 component in a region expanding with the cell would reach about
  0.41. The pocket is overdense and its clocks run slow.
- **Both misses are pocket-scale.** One is how D_s invades; the other is how fast the held region itself expands.
  Neither is a function of the lane's variables. **Open:** a derived law for either.

## 1b. Lane against cell

**From t = 0 (the case that failed before).** Residuals at t = 24, r_max 24, fresh cell (stored cell the same to two
digits):

| lane | ΔH/H(t0) | Δ ln a | Δρ_dm/ρ_dm | Δρ_SD/ρ_SD |
|---|---|---|---|---|
| two components, derived release | +2.5e-4 | +5.7e-4 | −1.7e-3 | **+27%** |
| no S–D species (the lane as it is) | +1.2e-3 | +2.4e-3 | −7.1e-3 | — |

The +27% is the release law's lateness combined with the held region's slower expansion: the lane holds the step
undiluted too long, in too much volume.

**After settling (from t = 12).** S–D matter as dust:

| cell | ΔH/H(t0), t = 24 / 66 | Δ ln a, t = 24 / 66 | Δρ_SD/ρ_SD, t = 24 / 66 | lane as it is: ΔH/H(t0), t = 24 / 66 |
|---|---|---|---|---|
| r_max 16 (stored) | +1.5e-4 / +2.4e-3 | −2.3e-4 / +8.5e-4 | +0.6% / −0.03% | +1.7e-3 / +1.4e-2 |
| r_max 24 (stored) | +7.1e-5 / +1.2e-3 | −1.1e-4 / +4.1e-4 | +0.4% / −0.6% | +6.0e-4 / +5.1e-3 |
| r_max 24 (fresh) | +8.5e-5 / — | −1.3e-4 / — | +0.4% / — | +6.2e-4 / — |
| r_max 36 (stored) | +2.5e-5 / +4.2e-4 | −3.8e-5 / +1.5e-4 | +0.3% / −0.3% | +2.0e-4 / +1.7e-3 |
| r_max 36 (fresh) | +3.0e-5 / — | −4.5e-5 / — | +0.3% / — | +2.0e-4 / — |
| x_c 1.5, r_max 24 | +4.1e-5 / +6.8e-4 | −7.0e-5 / +2.3e-4 | −0.2% / −0.6% | +4.3e-4 / +3.6e-3 |

(The fresh runs end at t = 24.)

- **Read: pass.** Residuals are at the size of the companion terms (PR #484). Carrying the species improves H 8–11-fold at
  t = 24 and about fourfold at t = 66.
- **Scaling:** H's residual at t = 24 goes 1.5e-4 : 7e-5 : 2.5e-5 for r_max 16 : 24 : 36, about the pocket's share
  (0.26 : 0.09 : 0.03). Stored and fresh agree to 20%.
- The two-component lane started at t = 12 does slightly worse than dust (+2.4% in ρ_SD). At t = 12 the snapshot's D
  includes D_r's ringing, which the law reads as held step. The held measure needs the slow part of D, which a single
  snapshot does not give.

## 1c. Controls (each fails as stated; from t = 0, t = 24, r_max 24)

| control | predicted failure | Δρ_SD/ρ_SD | ΔH/H(t0) |
|---|---|---|---|
| without the held step (its energy left out) | S–D energy too low by about the step | −36% | +1.0e-4 |
| step held constant, no fading | S–D energy runs away as the cell expands | +1 060% (×11) | +4.7e-3 |
| S–D matter as dust throughout (the step carried as dust) | too low: misses the step's negative pressure | −22% | −1.5e-4 |

- The truth lies between "dust throughout" and the two-component lane. The release is faster, and the held region
  smaller, than the lane can know.
- With the step held constant, H fails at 1% by t = 24 and 6% by t = 66.
- The other controls stay at 1e-3 in H: S–D matter is only 9% of the density at r_max 24.
- x_c 1.5 and the other cell sizes give the same pattern.

## 1d. 0* verdicts

| item | verdict |
|---|---|
| S–D matter as dust after settling | **Read, passes**, at the companion terms' size, ∝ pocket share; two rulers (solved and read densities), three cells |
| in-sea correction (2/3)g − (3/28)x | not needed at this accuracy (≤ 1.5%); needs x_eff and g, which the lane does not carry |
| the held step's fading law | **Derived in form** (band exit, then γ = √(1 − 3x/4)); its timing and the held region's own expansion are **pocket-scale**: not functions of the lane's variables |
| where the held step belongs | **with the reset**: it is D falling off the saddle, the ejection of POKE-1 and C6. The lane starts from its output |
| the one-way-cosmos gap | **narrowed, not closed.** The channel's leading part is closed by carrying the species, validated on independent rulers. The lane still only receives: no two-way test exists yet (the two-kinds-of-cell design is parked) |

**Before a production change (for Jeremiah's review; nothing proposed yet):**

1. The floor question (λ > 0 promotion is held).
2. At λ = 0 the species is S alone (D is not a species there, PR #484).
3. dH in Raychaudhuri form.
4. The lane's initial state as the reset's output: needs the reset derived, not borrowed from "today".
5. The two-way test.
