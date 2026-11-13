# Outer boundaries and return flow (HANDOFF to projects Nov 11 (2), Tasks 2 and 3)

**Research only. No production change.** Script: `research/sd_return_flow_v1.py`.
Set-up: PR #470's cold-start pocket. Flat space, no gravity, GEAR-67 quartic, S and D both evolved (full back-reaction),
R = 8, smooth edge (width 1), D seed 1e-10 inside, cold start (D relaxed with friction against the time-averaged S first).

**Unpaid inputs:** λ = 5000; D's seed 1e-10; cold start. **No φ / no vacuum** (no V0, no φ_today, no β_D). G_reset not
involved; G_N not involved (no gravity).

## Task 2: boundaries (map row 7), research settings

| Setting | What it does |
|---|---|
| absorb (sph) | sponge for r > 40, out to 70: outgoing energy removed and counted |
| reflect (sph) | no sponge, Neumann wall at r = 70: everything that leaves comes back |
| periodic (row) | planar slab pocket in a periodic box of length 40 (a row of identical pockets, 1D) |
| absorb (row) | the same slab with sponges at both ends (control for the row) |

Energy bookkeeping closes (box + absorbed − start) to ≤ 3e-4 of the start, except periodic x = 3.27 at dr = 0.04 (2e-3;
9e-4 at dr = 0.02). Energy leaving (absorbing, spherical, by t = 200): 0.57 of the start (x_in = 2.0), 0.69 (x_in = 3.27).

## Task 3: return flow. S's energy in the pocket (r < R + 3), relative to its start

| x_in | case | t = 40 | t = 80 | t = 130 | t = 160 | late (t = 350–400, mean) |
|---|---|---|---|---|---|---|
| 2.0 | sph absorb | 0.65 | 0.30 | 0.30 | 0.18 | – |
| 2.0 | sph reflect | 0.65 | 0.30 | 0.30 | 0.18 | – |
| 2.0 | row absorb | 0.87 | 0.10 | 0.015 | 0.03 | – |
| 2.0 | row periodic | 0.87 | 0.13 | 0.38 | 0.27 | 0.38 (dr 0.04) / 0.46 (dr 0.02) |
| 3.27 | sph absorb | 0.60 | 0.53 | 0.16 | 0.07 | – |
| 3.27 | sph reflect | 0.60 | 0.53 | 0.16 | 0.07 | – |
| 3.27 | row absorb | 0.36 | 0.03 | 0.04 | 0.06 | – |
| 3.27 | row periodic | 0.37 | 0.32 | 0.25 | 0.32 | 0.31 (dr 0.04) / 0.32 (dr 0.02) |

Periodic row, mean over each 50 time units, t = 0 to 400 (T = 400):

| x_in | dr | 0–50 | 50–100 | 100–150 | 150–200 | 200–250 | 250–300 | 300–350 | 350–400 |
|---|---|---|---|---|---|---|---|---|---|
| 2.0 | 0.04 | 1.03 | 0.33 | 0.27 | 0.33 | 0.33 | 0.26 | 0.36 | 0.38 |
| 2.0 | 0.02 | 1.03 | 0.37 | 0.27 | 0.34 | 0.47 | 0.44 | 0.40 | 0.46 |
| 3.27 | 0.04 | 0.67 | 0.29 | 0.31 | 0.28 | 0.28 | 0.31 | 0.29 | 0.31 |
| 3.27 | 0.02 | 0.68 | 0.28 | 0.31 | 0.30 | 0.29 | 0.30 | 0.31 | 0.32 |

D's energy by wavenumber band (fraction), late time:

| case | k < 1 | 1–3 | 3–10 | > 10 | top quarter of grid |
|---|---|---|---|---|---|
| sph reflect, x 2.0, t = 160 | 0.28 | 0.67 | 0.05 | 0.000 | 1e-7 (peak 2e-4 at t = 130) |
| sph reflect, x 3.27, t = 160 | 0.53 | 0.40 | 0.07 | 0.001 | 1e-4 |
| row periodic, x 2.0, t ≥ 350, dr 0.02 | 0.26 | 0.38 | 0.36 | 0.002 | 4e-19 |
| row periodic, x 3.27, t ≥ 350, dr 0.02 | 0.16 | 0.33 | 0.46 | 0.056 | 4e-19 |

## Readings (derived / read / open)

1. **Read:** a reflecting edge does not restore S. The spherical reflect run is identical to the absorbing run until the
   reflected wave returns (t ≈ 60–80) and drains to the same level (0.18 / 0.07 by t = 160). The returning energy stays in D.
2. **Read:** in the periodic row (1D, pockets 40 apart) S drains to 0.13–0.20 and then **settles into a fluctuating
   balance at ≈ 0.3 of its start (x_in = 3.27, resolution-independent to 0.01) and ≈ 0.3–0.46 (x_in = 2.0, the late level
   moves with resolution by ~0.08)**. Without return (row absorb) the same slab drains to 0.02–0.06.
3. **Read:** no cascade to the grid scale. The top quarter of the grid's wavenumbers holds ≤ 2e-4 of D's energy in every
   run; halving dr does not move energy there. D's energy does migrate to shorter wavelengths in the periodic row
   (3–10 band 0.36–0.46; k > 10 rises to 0.056 at x_in = 3.27 at dr = 0.02, from 0.026 at dr = 0.04), so the
   spectrum is still filling at k ~ 10 by t = 400. Whether that continues on longer runs is **open**.
4. **Open:** the 1D row is a strong return (nothing spreads in 1D). A 3D periodic lattice was not run.
5. max|D|/D* in the pocket stays 1–2 in every case: the pocket is filled with D (as in PR #470).
