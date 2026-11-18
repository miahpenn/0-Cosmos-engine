# GEAR N-2: does the CLM threshold survive when the forcing is prescribed, not feedback? (Nov 18 2026)

Scripts: `research/clm_solver_v1.py` (CLM pseudospectral RK4, 2/3 dealiased), `research/clm_open_loop_v1.py` (run from repo root). Independent of the cosmos engine; no engine code touched.

## Why
GEAR N-1 forced with f_hat = gamma nu k^2 omega_hat in k = 10..20. That is proportional to omega: in the band it is just reduced viscosity nu(1-gamma), a feedback, not an external force. The 3D construction it compares itself to uses a prescribed forcing. Test: prescribe the forcing and see whether blowup (and a threshold) survives.

## Setup (Choices; N-1 did not state its initial condition)
omega0 = -sin x, N = 256, nu = 0.08, T = 35, dt = 5e-4, blowup = max|omega| > 1e3. Unforced nu = 0.08: bounded, peak 2.15, decays. Unforced nu = 0: blows up near t = 2.4 at this criterion (N-1: ~2.30).

## Result (Measured)
1. Closed-loop threshold with THIS initial condition: gamma_crit in [0.9547, 0.9558] (N-1: 0.977 with its own, unstated, initial condition). The 0.977 is not reproduced and depends on the initial data; the order of "a bit below 1" is.
2. Near the threshold the bounded run is violent: gamma = 0.95 peaks at max|omega| = 417 and then decays (so the 1e3 criterion sits only a factor 2.4 above a bounded run's peak).
3. Open loop (a): record f(x,t) from the bounded closed-loop run at gamma = 0.95 (replaying it at s = 1 reproduces that run, peak 417.44), then scale it by s. Blowup for s above 1.344-1.352. The closed-loop threshold is gamma/0.95 = 1.006. So the prescribed forcing must be about 35% stronger than the feedback forcing at the closed-loop threshold.
4. Open loop (b): a steady forcing, the time mean (t < 15) of the same record, scaled by A: blowup for A above 3.06-3.12. Not comparable in size to (a): the mean is dominated by the large transient and is held constant for all t.

## Reading
A threshold and finite-time blowup survive with a prescribed forcing, so "forcing converts a viscously stable system to blowup" does not rest on the feedback form. But the feedback forcing is more efficient: it follows the solution (large when omega is large), and the same trajectory-shaped forcing as a fixed function of time needs ~1.35x the amplitude. So N-1's gamma_crit ~ 0.977 is not a property of "forcing strength" alone; it is a property of the feedback loop plus the initial data. Open-loop forcing amplitude thresholds are the numbers a prescribed-force statement should use.

## Limits
- One initial condition, one band, one resolution (N = 256, T = 35). Not checked at N = 512.
- Blowup criterion 1e3 is close to bounded runs' peaks near threshold; the s and A brackets may shift with the criterion. Not tested.
- (a) forcing is trajectory-specific (recorded from one run) and linearly interpolated at RK stages; (b) is one arbitrary steady profile. Neither is the optimal forcing, so the open-loop thresholds are upper bounds on the minimal prescribed forcing needed.
- 1D scalar proxy only.

## Amendment 1: robustness check (same day), `research/clm_robustness_v1.py` (Measured)
| N | blowup criterion max\|omega\| > | closed-loop gamma_crit | replay s_crit |
|---|---|---|---|
| 256 | 1e3 | 0.9547-0.9558 | 1.344-1.352 |
| 256 | 1e4 | 0.9644-0.9662 | 1.794-1.800 (hit the search ceiling; may be higher) |
| 512 | 1e3 | 0.9531-0.9550 | 1.244-1.250 |

What this changes: the closed-loop threshold is stable in N at a fixed criterion (0.955) but moves 0.01 when the criterion goes 1e3 -> 1e4, and the open-loop ratio is NOT robust (1.25 to 1.8+). Near the threshold bounded runs peak at several hundred (417 at N=256, 457 at N=512, for gamma = 0.95) and then decay, so "max|omega| > X" counts transient amplification as blowup. The "35% stronger" in the main result is withdrawn as a number. What stands: with a prescribed forcing a threshold exists and the same trajectory-shaped forcing needs MORE amplitude than the feedback form (every row, at least 1.24x versus 1.006x). The size of the gap needs a proper blowup test (e.g. continue until omega_max follows ~1/(T-t) toward a finite-time singularity, or a much higher cut with the analytic CLM blowup-time fit). That also applies to GEAR N-1's 0.977, which was set by a cut of the same kind (its criterion is not stated).
