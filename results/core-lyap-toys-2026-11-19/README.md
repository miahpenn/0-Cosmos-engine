# #523 re-check with the toy's S (Nov 19 2026)

Scripts: `research/core2dof_lyap_toyS_v1.py` (#523's script, S mass term 4s -> 1s, nothing else), `research/core_lyap_toys_driver_v1.py` (process pool). Tables: `baseline.txt`, `claim.txt`, `oldrest.txt` (old S rerun where #523 has no published value), `rest.txt`. Design and the "moved" threshold were fixed in the coordination log (Nov 19) before running. Model, start, Benettin readings, T = 2000, 20 phases, rtol 1e-12 and 1e-13 as #523. Research only.

## Baseline (passed)
Old S, one phase, rtol 1e-12: n = 0.3: 0.092, 0.5: 0.049, 0.6: 0.016, 2/3: 0.0019. All inside #523's published per-n [min, max] (0.035-0.118, 0.029-0.068, 0.014-0.021, 0.0008-0.0032). The pipeline is unchanged.

## Result (Measured). "Moved" = |toy mean - old mean| > max(old sd over phases, old |1e-12 - 1e-13|)
| n | kappa | old lambda (1e-12 / 1e-13) | toy lambda (1e-12 / 1e-13) | moved |
|---|---|---|---|---|
| 0.1 | 0.92 | 0.072 | 0.077 / 0.078 | no |
| 0.2 | 0.84 | 0.077 | 0.088 / 0.084 | no |
| 0.3 | 0.742 | 0.080 / 0.074 | 0.085 / 0.078 | no |
| 0.4 | 0.63 | 0.064 | 0.076 / 0.084 | no (0.012 vs sd 0.014) |
| 0.5 | 0.500 | 0.049 / 0.052 | 0.074 / 0.072 | **yes, up 45%** |
| 0.55 | 0.42 | 0.048 | 0.028 / 0.028 | **yes, down 42%** |
| 0.6 | 0.316 | 0.017 / 0.018 | 0.038 / 0.040 | **yes, up 2.2x** |
| 0.62 | 0.26 | 0.0145 | 0.035 / 0.033 | **yes, up 2.4x** |
| 0.64 | 0.20 | 0.0066 | 0.0026 | yes (down to the control level) |
| 0.65 | 0.158 | 0.0033* | 0.0026 | not a published row |
| 2/3 | 0 | 0.0024 | 0.0027 | no |
| 0.7-1.3 | - | 0.0025-0.0028 | 0.0027-0.0029 | no (control level throughout) |
(*one phase, baseline run.)

- **lambda(n) changes shape in the middle, not at the ends.** n <= 0.4 and n >= 0.7 are unchanged within the stated spread. Between 0.5 and 0.64 the curve is no longer a smooth decline: 0.074 (0.5), 0.028 (0.55), 0.038 (0.6), 0.035 (0.62), 0.003 (0.64). A dip at 0.55 and a plateau at 0.6-0.62, then a sharp end. The old curve fell smoothly (0.049, 0.048, 0.017, 0.0145, 0.0066).
- **The edge of chaos did not move** (both between 0.62 and 0.64-0.65; the toy's drop is sharper). Chaos does not extend into the band's upper part (n up to 1.3 is at the control level).
- **The "lambda ~ 0.05" match moves from n = 0.5 to the toy's n ~ 0.55-0.6** (0.028 at 0.55, 0.038 at 0.6, 0.074 at 0.5; 0.05 is not hit by a scanned n). With the toy's S, n = 0.5 gives 0.07, not 0.05. lambda/kappa at 0.5 is 0.15 (old 0.10), at 0.6 it is 0.12 (old 0.055): "lambda/kappa ~ 0.055" is not preserved either.
- Phase spread is large (sd 0.004-0.024) and the dip at 0.55 sits on a curve with this noise; its exact depth is not resolved by 20 phases. Not scanned between 0.5 and 0.55 or 0.55 and 0.6.

## What it does not settle
- The effective n of a stone (chat's ITEM 3: n fluctuates across the band edge): unchanged question; the table only says where on the curve each n lands.
- Whether the toy's S equation is the right one for the core's purpose is chat's call (its review gave the three matches); this rerun only measures the consequence.
- The curve between 0.5 and 0.64 would need a finer n grid (0.52, 0.57, 0.58, 0.61, 0.63) to resolve the dip and the shoulder.
