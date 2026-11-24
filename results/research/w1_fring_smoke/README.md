# f_ring Stage 0 timing smoke: results

**Outcome: STOP.** Two stop rules hold on the dr 0.025 pair. Nothing is projected forward, and a new attempt needs Jeremiah's own go.

- **S0e (cost):** 2.70 s per step per arm at 2400 cells. The full run (52,650 steps) would take 79 serial arm-hours for the dr 0.025 pair, against the declared ceiling of 24. With the dr 0.05 pair added it is 99.
- **S0b (constraint growth):** three readings trip the limit.
  - Arm A's Hamiltonian residual.
  - The momentum residual in both arms.
  - Every one peaks at the origin cell.

Declarations: `PRERUN.md` (commit 9fd6ee8, merged with #583 before any step). Script: `research/w1_fring_smoke/run.py`. Full numbers are in `report.json`, `prechecks.json` and `arm_*.json`. Go: "Go and go" (03:54 UTC) and "Go" (04:02 UTC). The full run was never covered.

## Before any step: all PASS (`prechecks.json`)

- **Static-frame flux estimator** (areal R, H(t) = -<K>/3, f_R = f_r/R_r, G_st = f_R - v P_st), checked against the input-check integrand at step 0:
  - Correct estimator: max diff 5.1e-9 (D) and 5.1e-9 (S) of the window max. Tolerance 1e-3.
  - Fault, no cross term: misses by 0.36. Fault, H = 0: misses by 0.49. Both fail, as required.
- **Stop rules:** every fault case behaved as declared.
  - S0a detects a NaN cell.
  - S0b trips on 1.06 (H) and 1.08 (M) per-step growth from the floor, and passes 1.5 r0.
  - A 1e-10 checkerboard in K pushes M past its limit.
  - S0c, S0e and G_run trip on their faults and pass their clean cases.

## Stops, dr 0.025 pair (100 steps, t = 0.075)

| | r0 | step 10 | step 50 | step 100 | S0b limit | cell |
|---|---|---|---|---|---|---|
| A, H | 7.14e-12 | 5.8e-11 | 4.5e-10 | **1.25e-9** | 7.25e-10 | origin (R 0.013) |
| A, M | 4.0e-16 | 6.9e-12 | 5.3e-11 | **9.4e-11** | 4.0e-13 | origin |
| B, H | 2.91e-12 | 2.3e-12 | 1.05e-11 | 9.1e-12 | 2.95e-10 | origin; passes |
| B, M | 6.3e-16 | 2.6e-13 | 1.4e-12 | **5.8e-12** | 6.3e-13 | origin |

Each S0b limit is 1.54 r0 + 100 r_floor(H), or 1.54 r0 + 1000 r_floor(M). The floors were r0 itself in every case: eps x max|term| for H is about 7e-20, and eps max|K|/dr for M is 1.5e-16.

- **S0a:** all finite.
- **S0c:** max |S| 1.0000013 x step 0, passes.
- **G_run** (recorded): never near the guard. Max H is 1.1e-5 C_H and max M is 1.6e-6 C_M.
- **S0d** (recorded): arm B's non-uniformity reads exactly 0 at steps 1, 10 and 100.
- **H(t):** 5.66897e-3 to 5.66741e-3.

dr 0.05 pair (record only): S0b M trips in both arms. H passes in both, though arm A's H grows 60x and stays just under its limit.

Static-frame flux at R_f = 20, ring minus control (recorded): D is 5.28e-10 at step 0 and 2.89e-10 at step 100. Without the cross term it reads 4.64e-10 and 2.54e-10. These are instantaneous values, not cycle means.

## Observations (post hoc, not gates)

- ~~Growth is polynomial, not exponential.~~ Withdrawn by the reanalysis below: over the stored t range a power law and an exponential fit about equally well, so the t-law is undetermined. Every growing residual in arm A sits at the origin cell.
- **M starts at round-off, then jumps.** It begins at about 5e-16 and reaches 1e-13 to 1e-12 within the first step in both arms. The A4 floor sizes round-off, not the truncation that evolution switches on at step 1, so S0b M as written trips on a background-only arm. Arm B's H does not grow past its floor; arm A's does.
- **Cost:** 2.70 s per step at 2400 cells and 1.37 s at 1200, i.e. linear in N. With the two dr 0.025 arms in parallel, the full run is about 40 hours of wall-clock.

## Diagnostic reanalysis of the stored data (no new steps, nothing gated)

Asked for by the coordinator (Nov 24 04:11 UTC), ahead of Lumen²'s f_ring A5 S0b form. Script: `research/w1_fring_smoke/reanalyse.py`; numbers: `reanalysis.json`. The grids are compared at a common t, t in [0, 0.075]: step 2k on dr 0.025 against step k on dr 0.05.

**Limit:** each step stored only the max-norm H and M and the areal radius where each peaked. Cells 1 to 5, and pointwise ring-against-control values, were not stored, so they can't be read without new steps. The origin cell is read only where both grids peak there.

| arm, residual | t 0.0015 | 0.0075 | 0.015 | 0.0375 | 0.075 | where the max sits |
|---|---|---|---|---|---|---|
| A, H: dr 0.05 / dr 0.025 | 10.7 | 1.14 | 1.02 | 1.44 | **1.20** | origin on both grids from t 0.015 |
| A, M | 1.93 | 3.53 | 3.13 | 0.96 | 1.53 | mostly origin |
| B (control), H | 0.05 | 0.15 | 0.18 | 0.08 | 0.12 | edge or far cells (R 38 to 60) |
| B (control), M | 0.41 | 0.12 | 0.31 | 0.20 | 0.17 | mostly the edge (R near 60) |

- **Arm A, H at the origin:** the residual does not fall with refinement. The ratio is 1.0 to 1.4 from t 0.015 on, against about 16 for 4th-order truncation or 4 for a 2nd-order origin. By the replacement form's (iii) that reads as "no fall with refinement".
  - Fit over t in [0.0375, 0.075]: power law t^1.52 at dr 0.025 and t^1.17 at dr 0.05. An exponential fits about as well, at rates 28 and 22 per unit time.
  - A factor of 2 in t cannot separate the two laws, so the t-law is undetermined.
- **Arm A, M:** ratio 1 to 3.5, not a clean order.
- **Control (arm B):** its max-norm residuals are larger on the finer grid (ratio 0.05 to 0.4). They peak at the edge or at far cells, not the origin, so p cannot be measured on the control from the stored max-norms.
  - Control H at dr 0.025 runs 1.9e-12 to 1.05e-11. At dr 0.05 it runs 2.1e-13 to 1.1e-12.
- **Ring against control (max-norms, not pointwise):** at t 0.075 on dr 0.025, H is 1.25e-9 against 9.1e-12, a factor of 137, at different cells.

## Timing test (Jeremiah, "Go on the timing test", Nov 24 04:11 UTC; no gating)

Script: `research/w1_fring_smoke/timing.py`. Numbers: `timing_A025.json`, `timing_B025.json`. Each arm took one untimed warm-up step, then 10 timed steps at dr 0.025 on the corrected slice, one process at a time, single thread.

`kernel.step` always calls `kernel.diagnostics` on the accepted state, and that can't be switched off without an engine change. So the bare step is the step time minus the diagnostics time, with the diagnostics timed separately on the same state.

| per step (s, median) | arm A | arm B |
|---|---|---|
| step, including the kernel's diagnostics | 2.613 | 2.531 |
| the kernel's diagnostics alone | 0.389 (15%) | 0.385 |
| bare step | **2.224** | **2.146** |
| the smoke's extra readings (`measure`) | 0.116 | 0.114 |

**The diagnostics don't dominate: the step itself does.** One profiled step on arm B (cProfile; its overhead inflates the absolute times) puts the time here:
- **About 55-60%:** the radiation sector's primitive recovery, `valencia.recover_radiation`, called about 52,800 times per step as a per-cell Python loop.
- **About 30%:** the four CMC lapse solves.
- **About 17%:** the diagnostics.

Speeding up the radiation recovery would be an engine change, which is Jeremiah's call under the allowlist.

**Full-run cost for the dr 0.025 pair** (52,651 steps for 12 periods, 26,325 for 6). Wall-clock assumes the two arms run in parallel.

| | 12 periods | 6 periods |
|---|---|---|
| diagnostics every step (as now) | 38.3 + 37.1 = **75 arm-h** (about 38 h wall) | **38 arm-h** (about 19 h wall) |
| diagnostics only at sample times (every 34 steps; needs an engine change) | 32.7 + 31.6 = **64 arm-h** (about 33 h wall) | **32 arm-h** (about 16 h wall) |

Both are against the declared 24 arm-hour ceiling, so none of these four options fits.

## Plan estimate: Lumen²'s proposed origin diagnostic (NOT run; needs Jeremiah's go)

**Step 0 (t = 0):** pointwise H and M at cells 0 to 5 on dr 0.1, 0.05 and 0.025. This is a build plus one evaluation per arm: seconds for each.

**No dr 0.1 slice exists.** `native_static.npz`, Aster's `native.npz` and `members.json` hold only dr 0.05 and 0.025. Building one would take three things:
1. A fresh Newton solve on dr 0.1's own nodes, by the producer's own rule in `native.py`: `solver(0.1 sqrt(u), R)`, seeded from the same A1 member's even cubic, with the same gamma and the 1e-10 residual check. Aster's `native.npz` stores x only for its two grids.
2. The corrected-slice path of `research/w1_static_slice/build.py`: theta = omega_u H0 r^2/2 per solver node through the producer's `mapped()`, the theta = 0 positive control, and the S3 split for 'f.
3. The input checks before any step.

None of it is interpolated from the finer arrays (C9 G0.2). The build is minutes of compute.

**Open before the build:** dr 0.1 has no declared input-check bands. The edge tolerance, the 'd envelope and the 'f window and band are all declared per grid for 0.025 and 0.05 only, and would have to be declared first. At dr 0.1, S3's k dr is about 0.54 (about 11.7 cells per wavelength), so the 'f discretisation term is larger.

**Step 1:** both arms on all three grids to t 0.15, saving pointwise H and M at cells 0 to 5 every step plus the max and where it sits. Then dr 0.1 and dr 0.05 continue to one ring period, t 3.29.

Per arm, scaled linearly in N from the measured dr 0.025 numbers (dr 0.05 measured at 1.37 s per step in the smoke; dr 0.1 not measured):

| grid | steps to t 0.15 | steps to t 3.29 | hours per arm, diagnostics every step | at sample times only |
|---|---|---|---|---|
| dr 0.025 | 200 | (not continued) | 0.15 | 0.13 |
| dr 0.05 | 100 | 2,194 | 0.83 (to t 3.29) | 0.71 |
| dr 0.1 | 50 | 1,097 | 0.21 (to t 3.29) | 0.18 |

The total is about 2.4 arm-hours. With the six arm-runs spread over the 4 cores, wall-clock is about 1 hour, set by the dr 0.05 arms. Reading cells 0 to 5 every step uses a constraint evaluation like `measure`, already included above. "At sample times only" would still evaluate the cells-0-to-5 reading every step and skip only the kernel's full diagnostics, which needs an engine change.
