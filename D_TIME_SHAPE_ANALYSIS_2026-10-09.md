# D-mode time-shape campaign: time-resolved analysis (Lumen), 2026-10-09

Scope: the four verified ledgers from run 37873115509 (commit 9f0afa0; summary c28888ea...,
ledger hashes D0 a1885657..., Dhalf db672802..., Dbase 9b0a701b..., Ddouble caebba2c...).
Grid: 3000 steps, dt = 0.0075, checkpoints 0.25 interval. r_max = 160 (NOT comparable
to the admitted r_max = 80 set).

## Rules (fixed before the onset numbers were examined in final form)
- Onset of a quantity q: the first sample at which |q(t) - q(0)| >= f * max_t |q(t) - q(0)|.
- Primary marker f = 0.25. Sensitivity f = 0.10 and 0.50 reported.
- D-active alignment: each lag = onset(q) - onset(D_active_at_source_probe), same case.
  Test: is the spread of lags across the three D cases smaller than the spread of coordinate-time onsets?
- A first-attempt rule (fraction of peak magnitude) was DISCARDED: quantities that start at or near
  their peak (cmc_kdot, tau_rate, the constant-ish K-sources) gave onsets at t = 0.007 as an artifact.
  This change is recorded here, not hidden.

## Controls for transients
- S_K_source_center and matter_K_source_center have identical early deviations (about 3.0e-3 by t = 1)
  in D0 and all D cases: an initialization transient, not D-driven. Their later onsets are unreliable
  (the transient falls below the marker threshold once the late D response raises the maximum).
  EXCLUDED from ordering.
- cmc_kdot: early deviation identical across cases (5.4e-7); total excursion very small; onset shifts
  with D (4.9 control, 4.4 / 4.3 / 4.1). Reported as weak, not used for ordering.

## Coordinate-time onsets, f = 0.25 (D-driven quantities)
| Quantity | Dhalf (5e-11) | Dbase (1e-10) | Ddouble (2e-10) | spread |
|---|---|---|---|---|
| D_center | 20.17 | 19.61 | 18.99 | 1.18 |
| D_at_source_probe | 19.98 | 19.22 | 18.67 | 1.31 |
| D_growth_rate_center | 20.54 | 19.91 | 19.30 | 1.24 |
| D_K_source_center | 20.33 | 19.60 | 18.86 | 1.46 |
| vacuum_K_source_center | 19.85 | 19.48 | 18.97 | 0.88 |
| alpha_r_max (lapse gradient) | 20.66 | 19.94 | 19.28 | 1.38 |
| tau_rate (central clock rate) | 20.66 | 20.03 | 19.40 | 1.27 |
| D_active_at_source_probe (reference) | 20.94 | 20.17 | 19.47 | - |

Larger D onsets earlier, monotonically, for every D-driven quantity.

## D-active alignment, f = 0.25: lag relative to D_active onset
| Quantity | Dhalf | Dbase | Ddouble | spread (lag) | spread (coord) |
|---|---|---|---|---|---|
| D_K_source_center | -0.62 | -0.57 | -0.61 | 0.045 | 1.46 |
| D_center | -0.77 | -0.56 | -0.48 | 0.29 | 1.18 |
| D_growth_rate_center | -0.41 | -0.26 | -0.17 | 0.23 | 1.24 |
| D_at_source_probe | -0.96 | -0.95 | -0.80 | 0.16 | 1.31 |
| alpha_r_max | -0.28 | -0.23 | -0.19 | 0.09 | 1.38 |
| tau_rate | -0.28 | -0.14 | -0.08 | 0.20 | 1.27 |
| vacuum_K_source_center | -1.10 | -0.69 | -0.50 | 0.59 | 0.88 |

Reading: for every D-driven quantity, the D-active alignment reduces the spread across amplitudes
(by roughly 4x to 30x, except vacuum_K at about 1.5x). Each D-driven onset also precedes the
D-active-probe onset by about 0.1 to 1.0 coordinate time units. So the D-active state is a lagging
indicator relative to these quantities at this marker.

## Sensitivity and ambiguity (reported as required)
- The sign and size of the lags depend on f. At f = 0.10, D_at_source_probe has a LARGER lag spread
  than its coordinate-time spread (0.73 vs 0.58); all other D-driven quantities keep the tighter
  alignment. At f = 0.50 the lags are larger (roughly -0.1 to -0.9) and the same ordering holds
  qualitatively.
- The fine ordering among D_center, D_K_source, D_growth, the lapse gradient and the clock rate is
  within about one time unit and is NOT resolved robustly across markers. The robust statements are:
  (a) D-driven onsets cluster within about one time unit of the D-active onset; (b) larger D onsets
  earlier; (c) D-active alignment tightens the cross-amplitude spread for most D-driven quantities.
- No one-way causality is inferred from these onsets (per the preregistration).

## Final-time values (verified against summary)
D0 22.343354; Dhalf 21.402742; Dbase 20.840071; Ddouble 20.219128. Monotone in amplitude.
Zero cycle events, handoffs, turnarounds, re-expansions, trapped roots in all four.

## Limits
- r_max = 160 only; not merged with the admitted r_max = 80 results.
- One resolution (N = 160). One probe location.
- Onset markers are analysis markers only; they never feed back into evolution.
