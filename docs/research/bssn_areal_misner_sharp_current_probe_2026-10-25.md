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


## Component split of the corrected balance (N=240, \(\Delta t=0.005\))

The corrected rate splits exactly into the spatial-gradient and time/shift parts of \(\chi\), with decomposition closure below \(2\times10^{-13}\). The transformed current is split into energy-motion, radial-current, mixed-current and pressure-motion pieces.

| R | \(\dot m_{\rm spatial}\) | \(\dot m_{\rm time/shift}\) | \(4\pi R^2T_t{}^R\) | Residual |
|---:|---:|---:|---:|---:|
| 1 | +2.19798e-5 | -1.31432e-7 | +2.29424e-5 | -1.09398e-6 |
| 2 | +1.51945e-4 | -1.37149e-6 | +1.52363e-4 | -1.78933e-6 |
| 3 | +4.16339e-4 | -6.03329e-6 | +4.12076e-4 | -1.77024e-6 |
| 5 | +1.00621e-3 | -3.91301e-5 | +9.66680e-4 | +4.04149e-7 |
| 8 | +8.94602e-4 | -1.62535e-4 | +7.29670e-4 | +2.39717e-6 |
| 10 | +4.91893e-4 | -2.65883e-4 | +2.25387e-4 | +6.22956e-7 |

At \(R=1\), the transformed-current pieces are: energy-motion \(-4.10687\times10^{-6}\), radial current \(+2.59938\times10^{-5}\), pressure-motion \(+1.05449\times10^{-6}\), and mixed-current \(+9.81\times10^{-10}\). At \(R=5\), the corresponding pieces are \(-2.00704\times10^{-4}\), \(+1.10647\times10^{-3}\), \(+6.00394\times10^{-5}\), and \(+8.74\times10^{-7}\). Thus the total current is a cancellation of large signed pieces, just as the geometric mass rate is a cancellation between its spatial and time/shift contributions.

The residual is now localized more clearly: the coordinate/tensor transformation and component bookkeeping close internally, but the fixed-\(R\) identity still has \(O(10^{-6})\) residuals across \(R=1\ldots3\), while it is smaller at \(R=5\) and \(R=10\). This points to the inner spatial derivative / shell interpolation and regular-origin measure as the next numerical focus; it does not by itself imply a missing physical source term.


## Fixed-areal-radius shell balance

Using the same fixed-\(R\) endpoint values, I independently assembled each shell's mass-rate from the change in \(m(R_{\rm out})-m(R_{\rm in})\) and compared it with the net time-averaged boundary current \(4\pi R_{\rm out}^2T_t{}^R-4\pi R_{\rm in}^2T_t{}^R\). This is a telescoping shell balance, not yet the separate reconstructed-face FV observer.

For the smaller timestep at each resolution:

| Shell | N=60 residual | N=120 residual | N=240 residual | Relative residual at N=240 |
|---|---:|---:|---:|---:|
| R=1–2 | +1.113e-5 | -8.264e-6 | -6.965e-7 | 0.54% |
| R=2–3 | -3.344e-5 | +7.096e-6 | +2.459e-8 | 0.0095% |
| R=3–5 | +4.985e-5 | -7.978e-6 | +2.176e-6 | 0.39% |
| R=5–8 | +3.990e-5 | +1.647e-5 | +1.997e-6 | 0.84% |
| R=8–10 | -6.651e-5 | -5.770e-6 | -1.779e-6 | 0.35% |

The shell residual telescopes to the difference of the fixed-radius residuals with a closure gap below \(1.4\times10^{-19}\). All five shells improve substantially in absolute residual from N=60 to N=240; the R=2–3 shell is especially tight at N=240. The R=1–2 and R=5–8 shells remain the largest relative discrepancies. This supports the transformed integrated-current identity away from the centre, while keeping the inner-shell closure open.

The remaining follow-up is now narrower: compare this fixed-areal-boundary shell balance against a separately reconstructed face-flux FV calculation using the same \(R\)-measure. That will test face reconstruction and cell integration independently of the telescoping boundary arithmetic.

## Interpretation and next gate

This probe now applies both parts of the areal-coordinate transformation: fixed-\(R\) differentiation and the transformed mixed-tensor current. The balance improves with spatial refinement, and the dt/2 comparison indicates timestep size is not the main remaining discrepancy. The explicit mass-interpolation gap was an additional numerical effect in the preliminary measurement and is now separately accounted for.

Next, compare the fixed-areal-boundary shell balance with a separately reconstructed face-flux FV calculation in the same R-measure, and isolate the centre-cell reconstruction error. Keep this diagnostic independent of the pointwise-versus-expanded \(Q_r\) decomposition.

**Disposition:** coordinate/current transformation exercised; corrected fixed-\(R\) identity still open; full matter-coupled momentum propagation and total sector conservation remain open; no production changes; \`scientific_admission=false\`.
