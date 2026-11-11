# Floquet bands of D inside matter with the exact S oscillation (HANDOFF Nov 11, Task 1)

`research/sd_exact_floquet_v1.py` -> `output.txt`. Research only. Homogeneous S, linear D, no gravity.
**Unpaid input: λ (enters only through x = λA²). No V0, φ_today or β_D enters.** G_N plays no role; G_reset not involved.
Rates are per unit time (D's free rate is 1) and, where given, per radian of S's oscillation.

## Result: D is not held at large x

| x | S's frequency | uniform mode (k = 0) | peak rate | where |
|---|---|---|---|---|
| 1.0 | 1.171 | 0.472 | 0.472 | k = 0 |
| 4/3 | 1.222 | **0** | 0.204 | k² = 1.47 |
| 2.0 | 1.318 | **0** | 0.283 | k² = 1.17 |
| 2.45 | 1.379 | **0** | 0.331 | k² = 0.97 |
| 3.27 (bag surface) | 1.483 | 0.346 | 0.410 | k² = 0.60 |
| 10 | 2.150 | 0.690 | 0.690 | k = 0 |
| 100 | 6.077 | 0.804 | 0.804 | k = 0 |
| 10⁴ | 59.9 | 0.816 | 2.54 | k² = 0.808 x |
| 10⁶ | 599 | 0.8165 | 25.4 | k² = 0.808 x |

1. **The uniform mode is held only for 1.270 < x < 2.670** (exact S; the single-harmonic estimate gave 4/3 to 2.45).
   It is not held again at large x: PR #468's "held above x ≈ 60" was an artefact of the single-harmonic S and is
   withdrawn. Above the band the rate rises to √(2/3) = 0.8165 and stays there for every x (0.81525 at x = 10³,
   0.81648 at 10⁵).
2. **Inside the hold band, finite-wavelength modes still grow**: k² from about 0.5–1 up to 2, rates 0.20–0.33
   (PR #467 measured 0.10–0.26 in a confined geometry).
3. **Two instabilities survive to x → ∞:**
   - *Long wavelengths, k² < 2:* rate → √(2/3) = 0.8165 per unit time, independent of x. Per radian of S it falls as
     1.36/√x. Origin: in element variables x₁,₂ = (D ± S)/√2 the quartic is (λ/4)(x₁⁴ + x₂⁴), so at large amplitude the
     two elements are independent oscillators joined only by the bond −x₁x₂; matter is the antiphase state, and the bond
     pulls it toward in-phase (D) at a rate set by the bond alone.
   - *Short wavelengths, 0.75 x < k² < 0.866 x:* the Lamé equation D_uu + [κ² + 3cn²(u; ½)]D = 0, κ² = 2k²/x, u = √(x/2) t,
     with exactly one instability band 3/2 < κ² < √3 (found numerically: 1.5015 to 1.7310). Peak exponent 0.03598 per
     unit u at κ² = 1.614: rate 0.0254 √x per unit time, **0.0425 per radian of S (×1.31 per S period), for every large x.**
     This is the known λφ⁴ self-resonance (Greene, Kofman, Linde, Starobinsky 1997, g²/λ = 3).
4. **So there is no x at which matter holds D out completely.** For 1.27 < x < 2.67 only the uniform mode is held.

For x ~ 10¹⁵ (the handoff's estimate for tangible matter): long-wavelength D grows at 0.82 per unit time
(4×10⁻⁸ per S radian), and the short-wavelength band at 8×10⁵ per unit time (0.0425 per S radian) at k ≈ 2.8×10⁷.

## Limits

Homogeneous S (no pocket edge, no gradients of S); linear in D (says when D starts to grow, not where it ends);
no gravity, no expansion; classical. D's saturation and its back-reaction on S are not covered.
