# Bound stones on main (expanding) vs the frozen box (research only, Nov 17–18 2026)

Follow-up to #520, on chat's review and Jeremiah's go. #520's stone (R 30, x(0) 0.3) was carried off by the Hubble
flow. Question: with a bound stone, does the frozen box's large outflow survive on main, or is it a shadow of the
frozen background? Script: `research/bound_stone_main_v1.py` (the #520 fixture, conversion, speed shim and readings,
with R and x(0) as inputs, plus extra readings). `research/quiet_stone_main_v1.py` gains a `WINDOWS` setting; #520's
windows are unchanged. Readings: `summary.txt` (main runs + toys per configuration), `criterion.txt`. Nothing in
`engine/` or `vendor/` changed. No fitting.

## Binding criterion (Derived)

- **The draw.** The bag's draw (#518, from V_eff = 2n − (15/16)n²) is d²⟨r²⟩/dt² = −A, with A = (45/64) x_w. As the
  stone swells, x_w ∝ ⟨r²⟩^(−3/2), so A = A₀ (y₀/y)^(3/2), where y = ⟨r²⟩ (energy-weighted).
- **The flow.** A stone at rest in the expanding chart starts with the flow, v = H r, so ẏ₀ = 2H y₀.
- **The energy integral.** The draw holds the stone (ẏ turns back before y → ∞) iff ½ẏ₀² < ∫ A dy = 2A₀y₀:
  **B ≡ A₀ / (H² y₀) > 1.**
- **The flow's own push.** The virial term 2⟨v²⟩ = 2H²y₀ opposes the draw at the start. Initial deceleration needs
  A₀ > 2H²y₀, i.e. **B > 2** (the threshold used).
- **#520's stone:** B = 1.6, below both thresholds. It was seen to be carried off. The criterion is consistent with
  that (Read).

| R | x(0) | x_w | ⟨r²⟩ | B (H = 0.00567) | naive collapse √(2y/A) |
|---|---|---|---|---|---|
| 30 | 0.3 | 0.092 | 1255 | 1.6 (#520) | 197 |
| **12** | **0.3** | 0.092 | 201 | **10.1** | 79 |
| 30 | 1.2 | 0.370 | 1255 | 6.4 | 98 |
| **20** | **1.0** | 0.308 | 558 | **12.1** | 72 |

**Pick:** C1 = R 12 at #518's x(0) = 0.3 (the smaller-R route) and C2 = R 20 at x(0) = 1.0 (the denser route). Both
have B ≥ 10, five times the threshold. R = 30 cannot reach that margin: x(0) < 4/3 caps it at B ≈ 7. Both starts are
quiet: x(0) < 4/3, D at d² = 2u − 1.5x, fields at rest.

## Runs

- **Main:** production kernel, expanding lattice wall, production λ and φ, dt = 0.2 dr, r_max = 240, T = 320.
  C1 on grid A (dr 0.2) and grid B (dr 0.15); C2 on grid A only.
- **Toy:** #518/#519's `stone_toy_v1.run` with the same starts at dr 0.2/0.1/0.05 with gravity off, and dr 0.1 with
  gravity on. Same box and sponge, T = 320.
- **Windows:** t_b + 6 … t_b + 200, where t_b is the first 1.5R flux peak. No twin (not cheap).

## Readings (window t_b + 6 … t_b + 200; Read)

| | C1 main A | C1 main B | C1 toy (0.05 / 0.1 / 0.2 / 0.1 grav) | C2 main A | C2 toy (0.05 / 0.1 / 0.2 / 0.1 grav) |
|---|---|---|---|---|---|
| first burst t_b | 68.2 | 68.3 | 90.4 (all) | 78.2 | 62.4 (all) |
| ω (t; proper time at 1.5R) | 1.77; 1.77 | 1.77; 1.77 | 1.76–1.77 | 2.16; 2.16 | 1.96 |
| D flux out 1.5R | 1.13 | 1.38 | 1.98 / 2.11 / 1.91 / 1.90 | 44.8 | 97 / 111 / 96 / 97 |
| D flux out 3R | 1.20 | 1.41 | 3.95 / 4.21 / 3.81 / 3.88 | 38.1 | 104 / 105 / 95 / 105 |
| energy in areal 3R, change | −42.0% | −48.0% | −61 to −64% | −51.3% | −61 to −65% |
| energy in comoving 3R, change | −3.7% | −8.8% | = areal | −9.6% | = areal |
| D_rms(0) | 1.28 | 1.28 | 0.98–1.05 | 1.41 | 1.37–1.56 |
| halo energy (areal 1.5R–4R) | 38 | 40 | 88 / 79 / 77 / 70 | 1312 | 1075 / 888 / 1311 / 1420 |
| halo growth / unit | +0.002 | +0.003 | −0.005 to −0.011 (decays) | +0.008 | +0.010 to +0.014 |
| halo drift | +0.02 | +0.04 | +0.02 to +0.05 | +0.07 | +0.08 to +0.14 |

- **Bound on main (Read).** The S energy half-radius in areal (physical) radius against the Hubble factor:
  - C1: 12.1 at t = 0, 12.0 at t_b, 20.7 at t_b + 128, 47 at t_b + 200. The flow alone would give 12, 17, 28, 35.
    Toy: 12.0, 13.3, 30.6, 44.9.
  - C2: 20.1, 26.5, 54, 55 against the toy's 19.8, 17.9, 35.8, 59.4. The flow alone would give 20, 30, 48, 60.
  - Before the burst neither stone follows the flow: the coordinate half-radius shrinks (C1 12.1 → 8.5, C2 20.1 → 17.9).
  - After the burst both swell physically as the toy's stones do. C2 swells more before the burst (26.5 vs 17.9).
- **Lapse:** α at 1.5R and 3R stays within 0.997–1.002, so ω on proper time equals ω on t to 3 digits.
- **Kernel diagnostics (worst; C1 B: exterior |H| 1.8e-7, |M| 1.9e-7, origin 3.7e-4) (C1 A / C2 A):** exterior |H| 3.7e-7 / 9.2e-6, |M| 3.5e-7 / 1.3e-5. The origin cell
  (inside the exclusion radius) reaches 3.0e-4 / 1.4e-1. C2's origin residual is large (Open, as in #520).
  No lapse-solve failures.

## Answer (0*)

- **The outflow partly survives; it is not a shadow (Read, two configurations).** With bound stones the near-gap
  halo still forms on main:
  - C2's halo is inside the toy's own grid scatter (1312 vs 888–1420, a 60% spread), with similar growth and drift; no closer match is claimed.
  - C1's halo is about half the toy's (38 vs 70–88).
- **The radiated D flux is about 1.4–3.5× smaller on main (mostly 2–3×; C1 at 1.5R only 1.4–1.9×)** (C2 45/38 vs 96–111; C1 1.1/1.2 vs 1.9–4.2). The areal-3R
  energy loss is smaller too (42–51% vs 61–65%).
- **The burst moves:** C1 earlier (68 vs 90, same ω), C2 later (78 vs 62, ω 2.16 vs 1.96).
- **So the frozen box overstates the outflow by about ×1.4–3.5 (mostly ×2–3), but the burst-plus-halo picture is not an artefact of the
  frozen background.**
- **Comoving 3R (Read, not comparable to the toy's 3R):** its loss of 4–10% shows that most of main's areal-3R loss
  stays within the stone's comoving region. But comoving 3R is ~3× larger physically by the end, so it is no
  like-for-like ruler against the toy.

## Unpaid inputs

The criterion's model (energy-weighted y, x_w ∝ y^(−3/2), threshold 2 from the start-time virial terms), H at the
start only, the picks (R, x(0)), T = 320, dt, the windows (t_b + 6 … + 200, sub-windows in summary.txt), halo band
and shell (scaled with R), the S half-radius rule (r < 0.75 r_max), plus #518/#519/#520's inputs.

## Limits

Two configurations. One main grid for C2. No twin. The toy has a sponge; main has the wall. Spherical. Run data
(npz) are not committed; the script regenerates them.

## Amendment (18 Nov 2026; chat's review of #522 in the coordination log; docs only, no numbers changed)

- **Answer lines corrected** (above): flux ratio 1.4–3.5× (not 2–3.5×); C2's halo "matches" → inside the toy's grid scatter.
- **C1 is not converged:** main's 1.5R flux moves 22% between grids A and B (1.13 → 1.38); the toy's moves about 5%. The 1.4–1.9× at 1.5R is therefore uncertain.
- **The binding criterion is mixed:** the energy-integral threshold B > 1 has no virial term; B > 2 adds it by hand. Both are unpaid. Tested at one failure point (B = 1.6) and two bound points (≥ 10); nothing near B = 2.
- **Absent readings:** no Kodama clock beside the lapse reading (α only), and no D three-part split (at rest, slow ring, travelling). Not read here; label Open.
- **Centre cell, Open:** C2 origin residual 0.14; D_rms(0) main 1.28 / 1.41 vs the toy's about 1.0 / 1.4–1.6.
- Heavy follow-ups (C1 grid convergence, a stone with B near 2, the C2 centre residual) are not started; they wait for a decision that hinges on them.
