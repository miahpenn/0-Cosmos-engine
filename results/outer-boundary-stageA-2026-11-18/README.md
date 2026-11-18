# Outer boundary, Stage A: the absorbing sponge on S and D, gravity off (Nov 18 2026)

Script: `research/outer_boundary_stageA_v1.py` (reuses the toy's `Grid`, sponge profile and `run`; `stone_toy_v1.run` gained one optional key `rs`, default unchanged). Tables: `summary.txt`. Design and tolerances: project doc `outer-boundary-design`. Tolerances were fixed before the run; no sponge tuning (GAMMA_MAX=1, quadratic ramp, start r_s=52, r_max=88).

## Result (Measured, toy, gravity off)
| Test | Reading | Tolerance | Result |
|---|---|---|---|
| A1 pulse, D wave packet m=sqrt2 | returned energy fraction inside r<52, minus a no-sponge reference run | <1% | omega/m=1.05 (stress point): 7e-5; 1.2: 5e-7; 1.6, 2.0: <1e-7 |
| A1 rulers | sponge width 36 -> 72; dr 0.2 -> 0.1 (omega/m 1.2) | same | 1e-5 (1.05) / 7e-8; 5e-7 |
| A1 second ruler | sponge's take / initial energy | -> 1 | 0.990-0.994 (the rest is slow tail still inside r<52 at T, identical in the reference) |
| A2 boundary move, C1-type stone R12 x0.3, r_max 88 vs 132 | x_p, B-surrogate x_p/<r^2>, centre S^2 (t=50-100) | 1%, 1%, 2% | 1.3e-7, 7.7e-8, 7.2e-8 |
| A3 closure (box energy + sponge take - E0)/E0 | C1, dr 0.2 / 0.1 | <1% | -4.0e-3 / -9.2e-4 (second order) |
| A3 closure, pulse | dr 0.2 / 0.1, omega/m 1.2 | <1% | -4.8e-3 / -5.1e-3 |

All pass. The stress point (omega/m=1.05, lambda~14) returns 0.007%, far under the 1% line, so no closure-limit label is needed.

## How the returned fraction is measured (a correction made during the run)
The first version counted all energy left in r<52 at time T. That is dominated by slow components of the packet itself (raw 0.5-1%) and by the measure's own discretisation (it changed with dr), not by reflection. The reported number subtracts the same measure from the same pulse run in a domain wide enough (r_max = r0 + T) that nothing returns and with no sponge. A first reference domain was too small for omega/m=1.05 (the fast tail reflected off its wall and gave a spurious negative); widened to the light-speed bound. Raw values are in `summary.txt`.

## Limits (read before relying on this)
- A2 is weakly discriminating. At T=100 the C1 stone radiates only ~2e-4 of its energy into the sponge, so moving r_max changes almost nothing (1e-7). It shows the sponge does not disturb a quiet stone; it does not show the closure survives a strong outflow. The pulse test is the discriminating one.
- A3's pulse closure (about -0.5%) does not shrink with dr, unlike C1's. Within tolerance, cause not identified (candidates: energy bookkeeping of the damped leapfrog, the Neumann wall, the packet's bandwidth). Not chased.
- Linear D only, small amplitude, gravity off, one geometry. Nothing here says the sponge is clean with gravity on, with a nonlinear stone, or against the Newtonian exterior; that is Stage B (needs its own go).
