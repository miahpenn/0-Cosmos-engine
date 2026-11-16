# Hand-floor removal, step 2: the engine change (for Jeremiah's review)

Step 1 (research, PRs #499-#502, `results/floor-swap-2026-11-15/`) showed the swap on production's own path through
research switches. This step moves it into `engine/`. Script: `research/hand_floor_step2_checks_v1.py`
(`fingerprint`, `lane`, `kernel N T`). Short runs only.

**Labels carried by every result here:**
- **Floor:** mechanism derived (D's vacuum at rest in the sea, felt through phi by GEAR-91); size set through lambda.
- **lambda = 1.56889e6 in production: a ruler**, computed from the constants: matched today so that the sea's floor
  equals the old hand floor, `KAPPA e^(-4 beta_D phi_today) / (2 lambda) = 0.10 V0`, through the borrowed "today"
  (V0, phi_today). C6's reset criterion gives lambda* ~ 10-1000 (conditional on C6). The gap of three to five orders of
  magnitude is **the relocated hierarchy**, recorded, not hidden.
- **beta_D = -0.04:** form derived; value not derived (the locked dark-matter coupling, borrowed). beta_S = 0 (switch).
- **All times are chart-dependent through phi_today = 0.179055.**
- **Held-step readings** with the coupling on: against D's phi-dependent rest, `engine.scalar_system.sd_rest_value`.
- **Unpaid inputs:** H/m_S (0.00567 today; the one unpaid number under the one-budget rule); V0 and phi_today
  (borrowed "today"); the floor fraction 0.10 (forced: the old hand floor's size, kept as the matching ruler); beta_D's
  value. Fixtures only: x = 2 of the rescaled pocket (by hand), the inherited amplitude 0.01, the D seed 1e-10, the
  cold start's T = 200, gamma = 0.3.

## Results

| # | check | result | status |
|---|---|---|---|
| 1 | old path kept as the reference is bit-identical | SHA-256 equal to main (2ba2da2): 25-step production evolution (N 60) and the COSMOS lane (`fingerprint_*.txt`); the tracefree-Ricci audit's report and profiles on the reference equal main's byte for byte | Read |
| 2 | new COSMOS lane: anchors | node / turnaround / crunch **2.9543 / 4.2604 / 5.9087** Hubble times from today (lane units, 1/m_S: 521.1 / 751.5 / 1042.3); ln(a/a_today) at turnaround 1.739. Reference lane 3.7526 / 5.4183 / 7.5883 (lane units 662.0 / 955.8 / 1338.6). Anchors 2.954 / 4.260 / 5.909 and 3.753 / 5.418 / 7.588 (`lane.txt`) | Read (floor form Derived) |
| 2 | new COSMOS lane vs the dynamic-D lane (`research/floor_swap_lane_v1.py`) | same events to 4 decimals; max relative difference over t <= 0.95 crunch: a 1.7e-7, H 6.4e-7, phi 1.3e-7, Pi 6.1e-7; scaled and unscaled identical; H today equal to 11 digits | Read |
| 3 | production kernel, background default, against the new lane to t = 200 (1.13 Hubble times), r_max 40, means over r < 20 | N 20: phi 2.4e-8, Pi 9.0e-8, H 3.3e-8 (worst relative); N 40: phi 6.3e-8, Pi 1.8e-7, H 5.2e-8. D / D_rest(phi) within 7.3e-5 (N 20) and 7.1e-5 (N 40): D rings about its moving rest (`kernel_N20.txt`, `kernel_N40.txt`) | Read |
| 4 | full suite | 1086 passed, 2 xfailed (main: 1072 passed, 2 xfailed; +14 new tests) | Read |

Rulers: the anchor times are given in Hubble times from today and in lane units (1/m_S); the kernel check in two
quantities (phi and H) at two resolutions. The step-1 expectation for check 3 was 1e-6 to 1e-5 (driver with a stale
lapse); with production's own initializer (lapse solved from alpha = 1) the agreement is ~1e-7.

**Also Read.** On the reference path the old default pocket (N 20, r_max 40, CFL 0.06) does not reach t = 50: D's
1e-10 seed grows ~e^t (lambda = 0: D is unbounded; max |D| 1.6e-9 at t = 3, 2.3e-3 at t = 18, 5e-2 at t = 24), the
lapse at the centre falls to 0.53 at t = 21 and 0.05 at t = 27, and the step fails between t = 25 and 50 (invalid
volume element). This is
the old production default, not a new defect; it is why the reference is kept for reproducibility only.

**Open.** (i) The COSMOS <-> local conversion (`research/lane_conversion_v1.py`) is pinned on the reference; on the
coupled path (D at rest in the local lane against the adiabatic floor in the COSMOS lane) it is checked here only
through check 3. (ii) The open-edge K-dot closure's undamped lapse mode (`results/floor-swap-2026-11-15/KERNEL_DRIFT.md`)
is not touched: task 2 (native zero-flux closure). (iii) Long runs on the new path (beyond t = 200) not made.
