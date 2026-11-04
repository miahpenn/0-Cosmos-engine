# Issue #4 — Full Bianchi predictor reuse comparison
Date: 2026-11-04
Status: diagnostic-only; one accepted-state class compared across three resolutions. No scientific admission.

## 1. Scope

This note records a reuse-only comparison of the source-derived full matter-coupled Bianchi predictor against the saved signed semi-discrete momentum tangent. It uses the archived accepted state at t=0.60 for N=60, 120, and 240. No workflow, observer, lapse solve, or evolution was run.

The predictor is
[
\dot M_{r,\mathrm{Bianchi}}=
\alpha K M_{r,\mathrm{tot}}
-\frac13(D_r\alpha)H_{\mathrm{tot}}
+\frac16\alpha D_rH_{\mathrm{tot}}
+D_j(\alpha W^j{}_r)
-8\pi\alpha\mathcal D_r^{\mathrm{tot}}.
]
The shift is zero on these accepted slices. The geometric constraints use metric-only Ricci; the evolved-connection tracefree offset is retained separately as W.

## 2. 0* parity audit of the connection offset

Keep the conformal factor \(\chi=X^2\) distinct from \(\chi_r=-D_1X/(2X)\). The connection difference \(C^r=\Lambda-\Lambda_{\mathrm{metric}}\) is odd at the origin, so its derivative uses the odd-parity stencil \(D_1^{(-)}\); derivatives of the even metric fields use \(D_1^{(+)}\):
[
\Delta_\Lambda R^r{}_r
=X^2\left[D_1^{(-)}C^r+\frac{D_1^{(+)}a}{2a}C^r\right],
\qquad
\Delta_\Lambda R=X^2D_1^{(-)}C^r.
]
Thus
[
W^r{}_r=X^2\left[\frac23D_1^{(-)}C^r+
\frac{D_1^{(+)}a}{2a}C^r\right].
]
Direct evaluation of the pinned Ricci brackets with evolved \(\Lambda\) and with \(\Lambda_{\mathrm{metric}}\) independently confirms these offsets:

| N | Scalar Ricci max gap | Radial Ricci max gap |
|---:|---:|---:|
| 60 | 5.52e-19 | 1.80e-19 |
| 120 | 6.99e-19 | 3.09e-19 |
| 240 | 7.15e-19 | 6.18e-19 |

This pins both the coefficient X² and the odd parity of C. The Bianchi connection divergence uses the even-parity derivative on F_C=alpha W^r_r and the pinned spherical tensor-divergence operator.

## 3. Operand validation

The accepted-state tangent artifacts include the independent total-stress divergence and the directional semi-discrete Mdot tangent. The total stress is independently reconstructed from six sectors: S, D, COSMOS phi, dark matter, baryons, and radiation.

The total-stress report's canonical projection cross-check passes for rho, pr, pt, and j. As a separate shared-state check, our reconstructed j_total agrees with the saved #423/#427 j_total operand to maximum absolute gaps:
- N=60: 9.76e-21
- N=120: 7.73e-21
- N=240: 8.47e-22

The six-sector raw divergence sum does not directly reproduce the total before time-density component accounting. The report explicitly accounts for the difference between the sum of per-sector finite-difference density rates and the total-projection finite-difference density rate. After that correction, the six-sector reconstruction is reported true, with unexplained remainder no larger than 5.70e-20 across the two epsilon branches and three resolutions. The uncorrected raw sum remains failed; it is not described as a direct pass.

## 4. Predictor versus measured tangent

Two independent discretizations of the total-stress divergence were retained:
1. pointwise covariant divergence;
2. finite-volume face-flux divergence.

Neither was selected based on residual size. Each was compared separately with the same measured directional tangent.

| N | Pointwise max, r<=10 | Pointwise annulus max, 5<=r<=10 | FV max, r<=10 | FV annulus max, 5<=r<=10 |
|---:|---:|---:|---:|---:|
| 60 | 1.8406e-5 | 2.0498e-6 | 9.6404e-6 | 1.8089e-6 |
| 120 | 1.5272e-6 | 1.6403e-7 | 1.2827e-6 | 1.3547e-7 |
| 240 | 3.9897e-7 | 1.0667e-8 | 4.0771e-8 | 8.9879e-9 |

The annular residuals exhibit approximately fourth-order refinement over this sequence:
- pointwise effective orders: 3.67, 3.90;
- finite-volume effective orders: 3.77, 3.88.

The centre cells remain anomalous and do not share a clean uniform order. Signed residuals at cells 0/1/2:

| N | Pointwise cells 0 / 1 / 2 | FV cells 0 / 1 / 2 |
|---:|---|---|
| 60 | +1.8406e-5 / +1.4017e-5 / +7.7777e-6 | -1.2950e-6 / +9.6404e-6 / +6.9830e-6 |
| 120 | +1.5272e-6 / +1.4480e-6 / +9.8044e-7 | -1.2827e-6 / +5.9546e-7 / +5.6046e-7 |
| 240 | +3.9897e-7 / +1.2524e-7 / +9.3559e-8 | +3.5862e-8 / +7.1513e-9 / +2.5872e-8 |

The directional tangent epsilon sweep (2e-5 versus 1e-5) is stable on r<=10. Maximum changes in the measured tangent are 1.38e-13, 1.93e-13, and 5.60e-13 for N=60, 120, 240; maximum changes in the pointwise predictor are 7.11e-15, 1.93e-14, and 4.36e-14.

## 5. Interpretation

This is the first full same-state comparison of the continuum Bianchi candidate—including normal trace pair, evolved-connection response, and independently measured total-stress divergence—against the signed semi-discrete tangent on the accepted-state class.

The annular comparison approaches fourth-order agreement. The centre-adjacent cells retain a distinct remainder consistent with the previously localized centre-row/stencil incompatibility. This does not prove that the centre row is the sole cause, and it does not establish full propagation closure across the other state classes.

The pointwise and finite-volume divergences remain separate comparators. The smaller finite-volume residual is not a reason to select it post hoc as the oracle.

No new observer or evolution is justified by this pass: all operands for this accepted-state comparison were already present in the archived artifacts. Next, apply the same frozen predictor and signed accounting to the remaining predeclared state classes, with identical conventions, parity, source projections, and operator ownership.

## 6. Provenance

Tangent/divergence artifacts are from workflow run 40236402303. Internal SHA256SUMS manifests were rechecked; all member hashes passed.

| N | Artifact ID | ZIP SHA-256 | NPZ SHA-256 |
|---:|---:|---|---|
| 60 | 12615136565 | 41d646090e05c08684806c5e39d57851f119cf0f38aa90f1d26ab5f72f87d1ce | 988ccd3c71ac6733bb47d013b32f8ceded9d132d2d21f7483c6ee4953993601b |
| 120 | 12615555594 | 9f5c2504f28290ff540c96f93e1fb569f5b960b7a2b3fc312dca9362c8e2eacf | 5751a16949b763ad5e96a4afac9db4fb1d92d55a0c6e7ddaabc4cc6b25652e0b |
| 240 | 12615550675 | e8a9cd26fbe3dbfb0bd62fc0606237c387a41efd14ade550df50da9842b38bea | 29a4bc1aeca718eff0dbc0de068809f8279aab94cb7023b004dea545589bf814 |

The tangent reports pin production source anchor 4aec88d0e3ff028e0a65a8ac36b30948625df792, vendor d6052d605673ce9d82cdc99b0df79855fa2ea215, Python 3.12.15, and NumPy 2.5.4. The total-stress helper is blob c10f9e8b9e654d5195fd55c493c3fe06ba5c7082; its report also records the historical source-reference commit and protected source blobs. CPU microarchitecture, BLAS core target, and NumPy SIMD dispatch are not captured, so no CPU-path-specific bitwise claim is made.

Research-only. No production changes, new workflow, observer, lapse solve, evolution, or admission-gate promotion.
