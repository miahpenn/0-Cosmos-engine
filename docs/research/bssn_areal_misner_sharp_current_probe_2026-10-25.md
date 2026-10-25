# Transformed fixed-areal-radius Misner–Sharp current probe — 2026-10-25

**Run:** [40160533041](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40160533041), all N=60/120/240 jobs passed. The probe used SHA-verified accepted \(t=0.60\) full-profile states and one copied kernel step at \(\Delta t\) and \(\Delta t/2\). Saved accepted states remained immutable; zero shift was checked explicitly; protected production/vendor/dependency files were unchanged.

This tests the coordinate map documented in [the areal-current map note](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/full-profile-propagation-candidate-ci-2026-10-24/docs/research/bssn_areal_current_coordinate_map_2026-10-25.md). It does not change production code.

## Method

From the pinned metric \(\gamma_{rr}=a/X^2\), \(\gamma_{\theta\theta}=br^2/X^2\), the probe computes \(R=r\sqrt b/X\), \(R_r\) with the pinned fourth-order odd-parity derivative, and \(R_t\) from the exact same-stage \(b_t,X_t\) RHS. It then forms
\[
\chi=\frac{X^2}{a}R_r^2-\frac{(R_t-\beta R_r)^2}{\alpha^2},\qquad
m_{\rm MS}=\frac R2(1-\chi)
\]
and transforms the mixed stress tensor to \(T_t{}^R\) with the full \((t,r)\to(t,R(t,r))\) Jacobian. The measured \(\Delta m_{\rm MS}/\Delta t\) compares interpolated start/end mass profiles at common fixed areal radii with the trapezoidal-in-time \(4\pi R^2 T_t{}^R\) RHS. No cell-index comparison is used.

## Results

The fixed-\(R\) identity residual decreases with spatial refinement, and halving \(\Delta t\) has little effect—especially at N=240—so the remaining error is not dominated by the finite time step.

| N | \(\Delta t\) | max \(|\dot m|_R-4\pi R^2T_t{}^R|\), \(R\le3\) | RMS residual, \(2\le R\le10\) |
|---:|---:|---:|---:|
| 60 | 0.040 | 2.745e-5 | 2.248e-5 |
| 60 | 0.020 | 2.752e-5 | 2.170e-5 |
| 120 | 0.020 | 9.508e-6 | 4.446e-6 |
| 120 | 0.010 | 9.593e-6 | 4.386e-6 |
| 240 | 0.010 | 5.330e-6 | 2.085e-6 |
| 240 | 0.005 | 5.296e-6 | 2.082e-6 |

At N=240 with \(\Delta t/2=0.005\), representative fixed-radius residuals are:

| Areal radius R | Mass rate | \(4\pi R^2T_t{}^R\) | Difference |
|---:|---:|---:|---:|
| 0.5 | +2.9572e-6 | +2.8312e-6 | +1.2599e-7 |
| 1 | +2.3174e-5 | +2.2942e-5 | +2.3190e-7 |
| 2 | +1.5251e-4 | +1.5236e-4 | +1.4962e-7 |
| 3 | +4.1187e-4 | +4.1208e-4 | -2.0685e-7 |
| 5 | +9.6631e-4 | +9.6668e-4 | -3.6677e-7 |
| 8 | +7.3039e-4 | +7.2967e-4 | +7.1499e-7 |
| 10 | +2.2505e-4 | +2.2539e-4 | -3.4146e-7 |

The relative mismatch is largest near the centre (about 4.5% at \(R=0.5\) for this N=240 sample), while it is around \(10^{-3}\) or smaller at most listed radii \(R\ge1\). The inner-cell residual remains visible and is why this is **not** marked as a closed conservation gate.

## Interpretation and next gate

This is the first direct test here that applies both pieces of the areal-coordinate transformation: the fixed-\(R\) time derivative and the mixed-tensor current \(T_t{}^R\). It is consistent with the archive's integrated-current target away from the centre and improves with resolution, but the inner residual still needs a source-aware regular-origin shell balance and an explicit truncation/error budget. The dt/2 comparison indicates that simply reducing the timestep is not the missing fix.

Next, split the inner shell balance into (i) geometric mass-rate contributions from \(R_t,R_r\), (ii) the transformed current contributions from \(E,j,p_r\), and (iii) the finite-volume shell integration/face reconstruction error. Keep this diagnostic independent of the pointwise-versus-expanded \(Q_r\) decomposition. Do not promote the source law or select a production operator from this result alone.

**Disposition:** transformed coordinate/current implementation exercised successfully; integrated-current identity remains open; full matter-coupled momentum propagation and total sector conservation remain open; no production changes; \`scientific_admission=false\`.
