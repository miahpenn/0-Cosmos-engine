# Pumping under expansion: finite-k D mode, friction on S and D (Nov 18 2026)

Script `research/pumping_finite_k_v1.py` (run from repo root), table in `summary.txt`. Follows the k = 0 run (results/pumping-core-expansion-2026-11-18). Research only.

## Setup (fixed before running; the tolerance was set in the script header)
S: s'' = -3H s' - (4s + s^3/2), starting at x0 with friction. D: one linear comoving mode, d'' = -3H d' + d - 1.5 s^2 d - k0^2 e^(-2Ht) d (physical k falls as e^(-Ht)); no back-reaction on S. k0 = u * k_res(x0), k_res^2 = 5 - 0.375 x0, u scanned 0.7-1.4 (29 points); the reported gain is the maximum over u. G = rolling-max |d| at t_exit (x < 4/3) over d0. Paper estimate (chat, order of magnitude): G ~ exp(E), E = 0.64 Gamma0^2 / H with Gamma0 from this script's H = 0 scan (0.140, 0.172, 0.223 at x0 = 1.6, 2.0, 2.67). Tolerance: |ln G - E| <= max(1, 0.3 E).

## Result (Measured)
| x0 | H = 0.003 | 0.0057 | 0.01 | 0.02 | 0.03 | 0.05 | 0.11 |
|---|---|---|---|---|---|---|---|
| 1.6 | 14.6 | 4.0 | 2.3 | 1.6 | 1.3 | 1.1 | 0.98 |
| 2.0 | 255 | 24 | 6.4 | 2.5 | 1.9 | 1.3 | 1.0 |
| 2.67 | 2950 | 87 | 17 | 4.9 | 2.9 | 1.8 | 1.1 |

- The pump develops under expansion at production H (0.0057): gain 4, 24, 87 for x0 = 1.6, 2.0, 2.67, and at the pocket H (0.003): 15, 255, 2950. It is a transient gain of order 4-3000, not the e^5-e^23 of the rate comparison.
- G = 10 is reached below H ~ 0.004 (x0 1.6), 0.009 (x0 2.0), 0.014 (x0 2.67) (linear interpolation of ln G in H between the scanned rows; not bisected). Chat's paper range was 0.008-0.014.
- The paper estimate E agrees with the measured ln G within the tolerance in 20 of 21 rows; at x0 = 2.0 within 0.2 e-fold at every H >= 0.0057. The one failure is x0 = 1.6, H = 0.003: measured ln G = 2.68 vs E = 4.17 (difference 1.48 against a tolerance 1.25), measured lower than predicted. Reported as is.
- The best k0 sits above resonance (u* 1.05-1.25, growing with H), as chat's sweep-through-resonance picture says. At H >= 0.05, where G ~ 1, u* is not meaningful.

## Limits
- One linear mode, no back-reaction on S, constant H; Q2 (does the edge move with amplitude) needs a spatial run and is not answered.
- The k0 grid step is 0.025 in u; the maximum G is therefore a lower bound on the true maximum over k0.
- Gamma0 is this reading's own H = 0 value, so E is not independent of the code.
- S's decay is the only effect on x; the engine's real pocket also has density and a Newtonian background not represented here.
