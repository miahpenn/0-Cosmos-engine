# Gradient-free 2-DOF core: does it give lambda ~ 0.055? (research only, Nov 18 2026)

Run request from chat's log entry (Jeremiah's go). Script: `research/core2dof_lyap_v1.py <n> <rtol>`. Readings: `summary.txt`.
Nothing in `engine/` or `vendor/` changed. No fitting. Forced = no spatial gradients, no radiation, no expansion.

Model (as specified; checked against V by differentiation): V = 2s² + s⁴/8 + (3/4)s²d² − d²/2 + d⁴/8, s'' = −(4s + s³/2 + 1.5 s d²), d'' = d − d³/2 − 1.5 s² d.
Start s = A (A² = 2n), d = 0, d' = 1e-3. Benettin: twin 1e-10 apart (fixed-seed direction), renormalised every 10, T = 2000, DOP853 at rtol 1e-12 and 1e-13 (atol 1e-14), 20 phases.

## Correction to the request's rulers (Derived)
The saddle rate κ = √(1 − 3n/2) is real only for n < 2/3. The requested n = 2/3 gives κ = 0 and n = 0.8 gives κ² = −0.2 (d = 0 is stable on average), so λ/κ is undefined there and the core should be regular. Run as negative controls; added n = 0.3 and 0.6 so three rulers lie below the threshold (n = 0.5 kept).

## Readings (Read; mean ± sd over 20 phases; two tolerances)
| n | κ | λ (rtol 1e-12 / 1e-13) | λ, 2nd half | λ/κ |
|---|---|---|---|---|
| 0.3 | 0.742 | 0.0795±0.023 / 0.0742±0.018 | 0.067 / 0.061 | 0.10–0.11 |
| 0.5 | 0.500 | 0.0494±0.011 / 0.0516±0.012 | 0.046 / 0.050 | 0.10 |
| 0.6 | 0.316 | 0.0169±0.002 / 0.0181±0.003 | 0.017 / 0.018 | 0.053–0.057 |
| 2/3 (control) | 0 | 0.0024±0.0006 | 0.0007 | undefined |
| 0.8 (control) | imaginary | 0.0026±0.0005 | 0.0007–0.0008 | undefined |

## Answer (0*)
- **At n = 0.5 the core alone gives λ ≈ 0.05 (0.049–0.052, phase spread ±0.011), close to main's 0.053 and the toy's 0.058** (Read). Weak chaos, as expected. Tolerance changes it by under 5% at n = 0.5; phase spread is the larger uncertainty.
- **But λ is not one number for the core**: it falls steeply with n (0.075 at 0.3, 0.05 at 0.5, 0.018 at 0.6, ≈ 0 at and above 2/3, where the controls are regular, as derived). So the match to 0.055 holds only if the stone's effective n is near 0.5; what n main and the toy actually carry is not read here (Open).
- λ/κ is 0.10 at n = 0.3 and 0.5, 0.055 at n = 0.6: the "λ/κ ~ 0.055" in the request is one point on a curve, not a constant.
- Reading per the request: ~0.05 → the core sets λ's size; the scaling in n is for chat to derive.

## Limits
Finite T (full-T mean exceeds the 2nd-half mean at n = 0.3, so there is a transient); 20 phases = s advanced by k·π/20 before measuring (a choice); one renormalisation interval; the 2/3 and 0.8 controls' small λ is finite-T residue.

## Pointer
Effective n: see chat's ITEM 3 entry and its close-out note (18 Nov 2026) in COORDINATION-LOG-0-star-work-chat-projects-Duet-read-first-each-session (Drive 1xj250uqmn2YnZXgH2M1LDME4eZUBsFz2TQEKy_LzlZk). Read there: the core's n fluctuates across the band edge during the chaotic phase (no single n; sampling coarse), so λ ≈ 0.055 is the core's λ(n) averaged over a fluctuating n. #520's own n has no stored core-x series: Open.
