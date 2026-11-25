# #606 flat-vacuum centre growth: readback of stored saves

Go: Jeremiah, "Do the readback" (Lumen²'s channel, Nov 25 2026 21:47 UTC, relayed). Stored data only, no engine stepping.
Source: `results/research/control_tests_stage1/arms/1a_dr025_dt0_r8.{json,npz}` (#606). Script:
`research/control_tests/centre_readback.py`. Output: `readback.json`.

The verdict rule was fixed before any number was read (in the script's docstring). Team log not used.

## What #606 stored

Full H and M profiles only at its 8 saves, t = 0, 0.3, ..., 2.1. No per-step centre values. Steps 1, 10 and 100,
and t 0.5 and 1.0, were not stored. So the step-1 comparison with eps/dr² (3.55e-13) and the step-10 clause of
ROUNDING-LIKE cannot be read from stored data. #591 was not used (its centre carries the ring).

## Readings, dr .025 (cells 0, 1, 2; |H| then |M|)

| t | \|H\| | \|M\| |
|---|---|---|
| 0 | 6.3e-152 (all three) | 0, 0, 0 |
| 0.3 | 9.0e-11, 3.3e-11, 1.5e-11 | 4.6e-11, 4.9e-11, 8.8e-12 |
| 0.6 | 7.8e-10, 4.8e-10, 4.6e-10 | 2.0e-10, 9.2e-11, 6.9e-13 |
| 0.9 | 4.9e-9, 3.5e-10, 5.8e-10 | 4.7e-10, 5.0e-10, 2.0e-10 |
| 1.2 | 2.5e-9, 2.7e-10, 1.5e-9 | 1.2e-9, 2.0e-10, 4.8e-10 |
| 1.5 | 4.0e-8, 7.7e-9, 1.8e-9 | 1.1e-8, 6.1e-9, 9.4e-11 |
| 1.8 | 1.1e-8, 3.1e-9, 1.9e-9 | 1.0e-8, 5.3e-9, 9.8e-10 |
| 2.1 | 1.3e-7, 1.8e-8, 8.7e-10 | 1.4e-8, 5.3e-9, 3.3e-9 |

H at t 0 is the dust term (density 1e-150), not rounding.

## Envelope fits, running max of |cell 0| over t 0.3 to 2.1 (7 saves)

| | p ± SE | rho_pow | gamma ± SE | rho_exp | N_e | verdict |
|---|---|---|---|---|---|---|
| H0 | 3.64 ± 0.29 | 0.41 | 3.78 ± 0.44 | 0.59 | 6.8 | OPEN |
| M0 | 3.19 ± 0.38 | 0.53 | 3.38 ± 0.37 | 0.50 | 6.1 | OPEN |

Neither meets INSTABILITY (rho_exp is not below half of rho_pow). Neither can be ROUNDING-LIKE: the step-10 share
is not stored, and the exponent bound fails anyway (p + 2SE is 4.2 and 3.9, above 2).

Reading only, dr .05 (r_max 8): H0 p 3.6 ± 0.7, gamma 3.4 ± 1.1, rho 1.03 / 1.52; M0 p 2.5 ± 0.4, gamma 2.7 ± 0.3,
rho 0.56 / 0.42. Both OPEN by the same rule.

## What would settle it

The missing rows need the engine: about 100 steps of the 1a dr .025 arm (n 320, about 0.27 s/step, under a minute)
give steps 1, 10 and 100; a per-step envelope to t 2.1 is 5600 steps, about 25 minutes. Either needs Jeremiah's go.
