# Transformed fixed-areal-radius Misner–Sharp current probe — 2026-10-25

**Corrected report:** [run 40160725682](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40160725682), all N=60/120/240 jobs passed. This supersedes the preliminary residual table from run 40160533041. The accepted \(t=0.60\) states remain SHA-verified and immutable; both runs use copied one-step probes at \(\Delta t\) and \(\Delta t/2\), with zero shift checked and protected production/vendor/dependency files unchanged.

## Method and correction

The pinned metric gives \(\gamma_{rr}=a/X^2\), \(\gamma_{\theta\theta}=br^2/X^2\), hence \(R=r\sqrt b/X\). The probe computes \(R_r\) with the pinned fourth-order odd-parity derivative, \(R_t\) from the same-stage \(b_t,X_t\) RHS, then
\[
\chi=\frac{X^2}{a}R_r^2-\frac{(R_t-\beta R_r)^2}{\alpha^2},\qquad
m_{\rm MS}=\frac R2(1-\chi).
\]
The mixed stress tensor is transformed to \(T_t{}^R\) using the full \((t,r)\to(t,R(t,r))\) Jacobian.

**Important interpolation correction:** the preliminary probe linearly interpolated \(m_{\rm MS}\) itself at each endpoint. That does not exactly preserve the nonlinear identity \(m=R(1-\chi)/2\) at a fixed interpolated \(R\). The corrected primary rate therefore interpolates \(\chi\) to common \(R\) first and computes
\[
\left.\dot m_{\rm MS}\right|_R=-\frac R2\frac{\Delta\chi|_R}{\Delta t}.
\]
The difference from the directly interpolated-mass rate is retained as a separate spatial interpolation diagnostic. The preliminary residual table should not be used as the primary result.

## Corrected results

| N | \(\Delta t\) | max \(|\dot m|_R-4\pi R^2T_t{}^R|\), \(R\le3\) | RMS residual, \(2\le R\le10\) | max mass-interpolation gap, \(R\le3\) |
|---:|---:|---:|---:|---:|
| 60 | 0.040 | 4.498e-5 | 3.424e-5 | 3.430e-5 |
| 60 | 0.020 | 4.517e-5 | 3.353e-5 | 3.424e-5 |
| 120 | 0.020 | 8.272e-6 | 7.275e-6 | 1.107e-5 |
| 120 | 0.010 | 8.375e-6 | 7.206e-6 | 1.106e-5 |
| 240 | 0.010 | 1.932e-6 | 2.948e-6 | 4.777e-6 |
| 240 | 0.005 | 1.953e-6 | 2.957e-6 | 4.777e-6 |

The corrected residual decreases under spatial refinement and changes very little when \(\Delta t\) is halved. The mass-interpolation gap is non-negligible, so it is kept separate instead of silently folded into the conservation residual. The corrected rate decomposition closes to \(<2\times10^{-13}\) absolute.

At N=240 with \(\Delta t/2=0.005\):

| Areal radius R | Corrected mass rate from \(\chi\) | \(4\pi R^2T_t{}^R\) | Difference |
|---:|---:|---:|---:|
| 0.5 | +2.84375e-6 | +2.83122e-6 | +1.25223e-8 |
| 1 | +2.18484e-5 | +2.29424e-5 | -1.09398e-6 |
| 2 | +1.50574e-4 | +1.52363e-4 | -1.78933e-6 |
| 3 | +4.10306e-4 | +4.12076e-4 | -1.77024e-6 |
| 5 | +9.67084e-4 | +9.66680e-4 | +4.04149e-7 |
| 8 | +7.32067e-4 | +7.29670e-4 | +2.39717e-6 |
| 10 | +2.26010e-4 | +2.25387e-4 | +6.22956e-7 |

The relative mismatch is largest near \(R=1\) (about 4.8% in this sample), falls to around 1.2% at \(R=2\), and is below about 0.4% at the listed radii \(R\ge3\). The \(R\le3\) residual and its refinement trend remain visible, so this is not a closed conservation gate.

## Interpretation and next gate

This probe now applies both parts of the areal-coordinate transformation: fixed-\(R\) differentiation and the transformed mixed-tensor current. The balance improves with spatial refinement, and the dt/2 comparison indicates timestep size is not the main remaining discrepancy. The explicit mass-interpolation gap was an additional numerical effect in the preliminary measurement and is now separately accounted for.

Next, split the corrected inner-shell residual into (i) spatial-gradient versus time/shift contributions to \(\chi\), (ii) transformed current contributions from \(E,j,p_r\), and (iii) finite-volume shell/face reconstruction. Keep this diagnostic independent of the pointwise-versus-expanded \(Q_r\) decomposition.

**Disposition:** coordinate/current transformation exercised; corrected fixed-\(R\) identity still open; full matter-coupled momentum propagation and total sector conservation remain open; no production changes; \`scientific_admission=false\`.
