# Toy-S rerun of the pumping and contraction rows (Nov 19 2026)

Why: chat's review (coordination log, Nov 18) found the pumping runs used #523's core S equation (S mass^2 = 4 at d = 0), while the stone toy has S mass^2 = 1 + 1.5 d^2 (4 only with D at the sea). The upper band edge x = 8/3 and the resonance k^2 = 2 - 0.375x need the toy's S. Scripts: `research/pumping_core_toyS_v1.py`, `pumping_finite_k_toyS_v1.py`, `contraction_crossing_toyS_v1.py` = the earlier scripts with the one term 4 s -> 1 s (nothing else changed). Tables: `summary.txt`, `contraction.txt`. The spatial back-reaction run (#536) is NOT rerun.

## Where the 4 came from (a read, no run)
#523's README: "Model (as specified; checked against V by differentiation): V = 2s^2 + s^4/8 + (3/4)s^2 d^2 - d^2/2 + d^4/8". Its 2s^2 is S's mass term with D frozen at the sea (1/2 + (3/4)(2) = 2) while D is also a dynamic field with the full (3/4)s^2 d^2 coupling; the sea's contribution is therefore counted in the S mass AND when d is at the sea. #523's lambda(n) used this potential, so its S frequency at d = 0 is 2, not 1. Whether lambda(n) shifts with the toy's S I did not test.

## Result (Measured), with chat's expectations
A. **Resonance and Gamma0 (toy S, H = 0, tiny D):** k_res^2 within 5% of 2 - 0.375x only at x = 1.6 (4%); at 2.0 it is 8% off, at 2.67 15% off (the measured resonance sits below the first-order formula; the pump is strong at large x, so higher-order shifts are plausible, not tested). Gamma0 = 0.24 / 0.28 / 0.35 at x0 = 1.6 / 2.0 / 2.67: inside 0.16-0.28 at the first two, above it at 2.67. So the expected 5% / window agreement holds in 1 of 3 rows on k_res^2 and 2 of 3 on Gamma0.
B. **Finite-k expansion gain (toy S):** at production H (0.0057) G = 678 (x0 2.0) and 1.2e5 (x0 2.67), against 24 and 87 with the core's S; at pocket H 2e5 and 3e9 (core 255, 2950). G = 10 now needs H of about 0.02 (x0 2.0). ln G is about 2x the core's, not 1.35x. The paper estimate E = 0.64 Gamma0^2 / H overpredicts: ln G 6.5 vs E 9.0 (x0 2.0), 11.7 vs 14.0 (2.67), 3.7 vs 5.1 at H = 0.01; the differences (2.5, 2.3, 1.4) are within max(1, 0.3E) in each row (2.7, 4.2, 1.5), but E is not independently accurate here.
C. **Contraction (toy S, d0 = 1e-8 unless stated):**
| Hdot | x_c 2.0 | x_c 2.6 | core (for comparison, x_c 2.0 / 2.6) |
|---|---|---|---|
| 1e-3 | fires at the crossing (t10 -18) | at crossing (t10 -17) | none / after |
| 1e-4 | way up (t10 -73) | way up (-78) | none / none |
| 1e-5 | way up (-204) | way up (-247) | way up |
| 4e-7 | way up (-1505) | way up (-2371) | way up |
All seeds (1e-8, 1e-4, 1e-2) give the same outcome in every row. So with the toy's S the pump is the trigger at every Hdot tested, 1e-3 included; the "needs a seed" bounce-rate result of the core run is gone. At Hdot 1e-5 and 4e-7 the fire again comes with x below 4/3 (about 1.1 and 0.7 at t10), now with the toy's S mass. The comparison E_b vs the measured crossing exponent could be read in two rows: 6.9 vs 7.5 (1e-3, 2.0) and 6.8 vs 10.3 (1e-3, 2.6; outside max(1, 0.3E)).

## What stands, what is superseded
- **Superseded:** the numbers of the core runs (pumping-core-expansion, pumping-finite-k, contraction-crossing tables, G = 24-87 at production H, H for G = 10, the 17-of-18 'no fire' count, E_b agreement) described the core's S, not the stone toy's. A note is added at the top of each of those READMEs.
- **Qualitatively survives:** the pump develops under expansion at production H; the k = 0 core has no pump; at machine-order Hdot (<= 1e-5) the pump fires before the bounce, for every seed and from below the band; large D saturates it.
- **Not tested:** the spatial back-reaction run with the toy's S (saturation scales shift with the gains); the 4e-7 rate itself; lambda(n) of #523 with the toy's S; whether the pumped mode stays trapped (omega_S at x = 2 is 1.32 < sea gap sqrt 2; not run).
