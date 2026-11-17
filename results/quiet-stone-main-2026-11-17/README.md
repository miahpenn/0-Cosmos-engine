# Chaos step 2: the quiet-start stone on main (expanding) vs the frozen box (research only, Nov 17 2026)

Step 2 of chat's brief (Jeremiah's go). Question: is the post-burst behaviour of #518's quiet stone changed by
expansion? The post-burst state is chaotic (#519, λ ≈ 0.058/unit), so only statistics are compared.
Script: `research/quiet_stone_main_v1.py` (run / toy / analyse). Readings: `summary.txt`. Nothing in `engine/` or
`vendor/` changed. `research/stone_toy_v1.py` gains two optional snapshot keys (`snap_r`, `snap_dt`); defaults unchanged.
No fitting.

## Set-up (labelled research fixture, not an engine default)

- **Main:** `V55ProductionKernel.initialize()` (today, background alone, lattice wall by `lattice_wall.default_edge`,
  wall clock "i", production λ = `SD_QUARTIC_LAMBDA`, φ = φ_today with production's coupling, production fluids incl.
  radiation). S, PS, D, PD are then overwritten and the slice rebuilt with `lattice_wall.conformal_wall_initial_data`
  (the rescaled-pocket fixture's own path) and the kernel's lapse solve.
- **Conversion (Derived):** engine units are μ = 1, G = 1, canonical fields; λ V_SD(S, D) = V_toy(√λ S, √λ D) − u²/2,
  so s = √λ S, d = √λ D, same lengths and times. Start: S = √(x(r)/λ), D = √(max(2u − 1.5x, 0)/λ), PS = PD = 0,
  x = `stone_toy_v1.quiet_profile(r, 30, 0.3, "sp")`. u = e^(−2β_D φ_today) = 1.0144, so main's sea is d = 1.4244
  (toy: √2 at φ = 0). Energies quoted ×λ (toy units), relative to the sea D_s(φ).
- **Expansion (Read):** H/μ = 0.00567 at the start, 0.00187 at T = 420; 1.414 e-folds (scale factor ×4.11).
- **Grids:** A dr 0.2 (N 1200), B dr 0.15 (N 1600); r_max = 240 = 8R; dt = 0.2 dr. Twin: A with S × (1 + 1e-10).
  T = 420 so that t_b + 250 fits (t_b ≈ 151). Toy: #519's three base configs rerun deterministically (gravity off,
  φ = 0, sponge, dr 0.05/0.1/0.2) with D snapshots every 0.4 out to 4.5R, needed for the halo (not in #519).
- **Early probe, amended:** a first grid-A run at dt = 0.4 dr (the toy's step) blew up at t ≈ 15 ("lattice-wall lapse
  solve has no positive solution"). All final runs use dt = 0.2 dr; none hit that error.
- **Speed shim (bit-identical, labelled in the script):** `matter_system._reconstructed_primitive` is O(N²) per step;
  the driver caches its arrays. Checked: max |difference| = 0.0 after 5 steps.

## Readings (windows in μ⁻¹ after the first burst t_b; Read unless marked)

Window t_b + 6 … t_b + 250 (sub-windows +6…+128 / +128…+250 in `summary.txt` are the ruler on the window).

| | main A (dr 0.2) | main B (dr 0.15) | toy 0.05 | toy 0.1 | toy 0.2 |
|---|---|---|---|---|---|
| first burst t_b | 150.7 | 151.2 | 145.6 | 144.0 | 144.0 |
| burst line ω (D at 1.5R, W 32) | 1.77 | 1.77 | 1.96 | 1.96 | 1.96 (#519: 2.01, its window) |
| mean D flux out, areal 1.5R | 11.6 | 9.1 | 41.2 | 44.3 | 53.3 |
| through 3R | 9.8 | 6.6 | 45.7 | 45.0 | 54.2 |
| energy inside areal 3R, change over window | −54.5% | −53.9% | −40.5% | −43.8% | −40.4% |
| E3(t_b)/E3(0) | 0.881 | 0.863 | 1.000 | 0.999 | 0.997 |
| D_rms(0) | 1.30 | 1.25 | 1.42 | 1.46 | 1.43 |
| halo: band 1.4–1.7 energy, areal 1.5R–4R (mean) | 329 | 281 | 797 | 927 | 1351 |
| halo centroid | 70.1 | 69.0 | 70.7 | 70.0 | 72.7 |
| halo centroid drift (units/μ⁻¹) | +0.145 | +0.106 | +0.092 | +0.117 | +0.158 |
| halo growth rate (+ = growing, /unit) | +0.0093 | +0.0117 | +0.0097 | +0.0097 | +0.0104 |
| halo pre-burst (t_b − 60 … − 4) | 0.88 | 0.87 | 6.4 | 5.0 | 4.7 |
| ruler: halo in comoving 1.5R–4R (areal measure) | 120, growth +0.035 | 92, +0.034 | = areal | = areal | = areal |

- **The stone is carried by the Hubble flow (Read).** D's deficit half-radius, areal/coordinate (A): t = 0 29.6/29.7,
  75 44.6/30.6, t_b 43.9/22.1, t_b + 125 107/36, t_b + 250 201/51. In the toy it stays 22–37 until the burst.
  The stone is not bound against expansion (gravity is weak: g ~ 1.6e-5).
- **Bookkeeping (Derived, arithmetic on Read values):** on main the D flux through 3R over the window (~9.8 × 244 ≈
  2.4e3) is under 10% of the loss of E3 (~2.8e4). Most of main's E3 loss is the stretched stone leaving the fixed
  physical 3R, not radiation. Toy: 1.1e4 of ~2.5e4.
- **Wall:** the D flux 40 units in from the wall never reaches 1% of the peak 1.5R flux (A, B). No reflections in the windows.
- **Kernel diagnostics (Read, worst over run):** exterior (r ≥ 5) |H| 3.5e-6 / 1.0e-6, |M| 1.0e-5 / 5.6e-6 (A / B).
  The origin cell (r = 0.1, inside the kernel's exclusion radius) grows to |H| 3.5e-2 (A) / 7.8e-3 (B), from 5e-10 at
  t = 0. It falls ×4.5 for dr ×0.75, so it is resolution error. Open: what it does to the centre readings. α(0) stays
  in 0.9936–1.0024.
- **Twin (A, S × (1 + 1e-10)), Read:** the RMS d difference over r < 3R goes 4.9e-9 (t 104) → 1.8e-8 (139) → 4.8e-6
  (244) → 4.5e-5 (418). That is ≈ 0.053/unit over 139–244 (toy #519: 0.058), slower after, not saturated by T.
  The snapshots are float32, so differences below ~1e-8 are not resolved. Its statistics equal A's (the twins have not
  parted by T). **Chaos persists on main (finite-time, one ruler).**

## Answer (0*)

- **Expansion does not drain the halo (Read, two grids).** The near-gap halo still forms after the burst, grows at the
  toy's rate (0.009–0.012 vs 0.010/unit) and drifts outward at the toy's speed (0.11–0.15 vs 0.09–0.16) in physical radius.
- **Expansion does change the averages (Read).** The burst is later by ~6–7 units and at a lower line (1.77 vs 1.96).
  The outgoing D flux is 4–8× smaller. Halo energy is 2.5–4× smaller. D_rms(0) is ~12% lower. E3 falls 12–14% before the
  burst and ~54% over the window (toy ~42%). The main driver is the stone being stretched with the flow (Derived from
  the bookkeeping above), not the halo being drained.
- **Grid dependence:** A vs B differ by up to 30% on the fluxes (toy: up to 30% across its rulers). The conclusions
  are differences of ×2.5–8, beyond this.

## Unpaid inputs

The start conversion (u at φ_today; d² = 2u − 1.5x), SP shape, R = 30, x(0) = 0.3, the coordinate profile on the
conformal slice, T = 420, dt = 0.2 dr, r_max = 8R, the windows (t_b + 6 … + 250, sub-windows, pre-burst), the
halo band 1.4–1.7 and shell 1.5R–4R (areal; comoving as ruler), W_HALO 32 / step 8, the on-shell energy
estimate ω²⟨D_band²⟩, the burst rule (#519's 1.5R flux peaks), the D-deficit size rule, plus #518/#519's own inputs.
Borrowed: λ, φ_today, β_D, the fluids (production, labelled there).

## Limits

One configuration. Spherical. Only two main grids, both coarser than the toy's finest. Main has gravity and
production's fluids; the toy has neither. The toy has a sponge; main has the reflecting wall. Times are the wall clock
(α(0) ≈ 1). The centre is exposed to the origin-cell residual.

Run data (npz, ~29 MB) are not committed; `run`/`toy` regenerate them deterministically, and `analyse` writes `summary.txt`.
