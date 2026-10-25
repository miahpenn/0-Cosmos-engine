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


## Separate parity-aware FV face reconstruction

Run [40161253120](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40161253120) adds a separate face path: even metric/stress profiles are reconstructed as regular polynomials in \(r^2\), odd \(j,R_t\) profiles as \(r\,g(r^2)\), and \(R_r\) is differentiated from the reconstructed \(b,X\) face profiles. The transformed current is then formed from the independently reconstructed face primitives; it is not obtained by interpolating the already-combined cell-centred current.

For the smaller timestep, the face-shell residuals are:

| Shell | N=60 | N=120 | N=240 |
|---|---:|---:|---:|
| R=1–2 | -2.652e-5 | +5.569e-6 | +4.453e-7 |
| R=2–3 | +1.282e-6 | -9.001e-6 | -4.939e-7 |
| R=3–5 | +2.923e-5 | +8.682e-6 | -4.414e-6 |
| R=5–8 | +1.800e-5 | -1.061e-5 | -9.561e-7 |
| R=8–10 | +3.851e-6 | +9.991e-6 | +5.622e-6 |

The face-shell residual decreases overall with refinement, but it is not numerically interchangeable with the cell-centred fixed-R shell balance. At N=240, the largest face-shell residual is \(5.62\times10^{-6}\), compared with \(2.18\times10^{-6}\) for the cell-centred shell balance. The face-minus-cell residual difference reaches \(7.40\times10^{-6}\) in R=8–10.

The difference decomposes into the face reconstruction's shell-mass-rate change minus its net-current change. At N=240:
- R=3–5: face mass-rate shifts by \(-4.00\times10^{-6}\), while net current shifts by \(+2.59\times10^{-6}\), giving a \(-6.59\times10^{-6}\) residual shift.
- R=8–10: face mass-rate shifts by \(+6.62\times10^{-6}\), while net current shifts by \(-7.85\times10^{-7}\), giving a \(+7.40\times10^{-6}\) residual shift.

This isolates a meaningful numerical crosswalk: the independent face reconstruction changes the geometric mass-rate and flux differently. The residual does decrease with N, but this face path is not yet an accepted FV identity. The next focused check is to compare the face \(R_r\) and \(R_t\) reconstructed from metric profiles against the face reconstruction of the pinned cell-centred derivative/RHS fields, then identify which part accounts for the face-minus-cell shell shift.


## Derivative-vs-current decomposition crosswalk

Run [40161486243](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40161486243) passed all three resolutions after adding controlled alternatives: (i) metric-derived face \(R_r,R_t\), (ii) face \(R_r,R_t\) reconstructed from the pinned cell-centred derivative/RHS fields, (iii) face \(\chi\) reconstructed from cell-centred \(\chi\), and (iv) face current reconstructed from the already-combined cell-centred current. This distinguishes face geometry from current construction.

At \(N=240\), the metric-derived and cell-derivative-derived face shell residuals are nearly identical: maximum absolute residual \(5.62\times10^{-6}\) versus \(5.59\times10^{-6}\). Their maximum face-\(R_r\) difference for \(R\le10\) is \(3.19\times10^{-9}\), and maximum face-\(R_t\) difference is \(3.15\times10^{-11}\). The metric-\(\chi\)/cell-current, cell-\(\chi\)/primitive-current, and cell-\(\chi\)/cell-current alternatives also remain tightly clustered, with maximum residuals \(5.6213, 5.5914, 5.5910\times10^{-6}\).

At \(N=60\), the face geometry alternatives are visibly less consistent (maximum \(R_r\) gap \(6.06\times10^{-7}\), maximum \(R_t\) gap \(7.03\times10^{-9}\)), but those differences shrink strongly with refinement. This suggests the previously observed face-minus-cell shell discrepancy is not principally caused by the reconstructed \(R_r\) or \(R_t\) formulas. At high resolution, changing the current assembly path also barely changes the shell residual. The remaining residual is therefore more likely tied to the time-discrete fixed-\(R\) mass-rate/current identity itself, or to the accepted-state temporal sampling, than to face interpolation alone. This is a diagnostic localization, not a proof of the continuum conservation law.

Next isolate the temporal balance: compare a centred two-sided fixed-\(R\) mass derivative against the existing forward accepted-step difference, and separately compare instantaneous endpoint \(4\pi R^2T_t{}^R\) with the time-averaged current. Keep both checks on the same immutable accepted-state artifact and retain the existing \(dt\)/\(dt/2\) arms.
## Interpretation and next gate

This probe now applies both parts of the areal-coordinate transformation: fixed-\(R\) differentiation and the transformed mixed-tensor current. The balance improves with spatial refinement, and the dt/2 comparison indicates timestep size is not the main remaining discrepancy. The explicit mass-interpolation gap was an additional numerical effect in the preliminary measurement and is now separately accounted for.

Next, compare face (R_r) and (R_t) from the metric-profile reconstruction with face values reconstructed from the pinned cell-centred derivative/RHS fields, to localize the face-minus-cell shell shift. Keep this diagnostic independent of the pointwise-versus-expanded \(Q_r\) decomposition.

**Disposition:** coordinate/current transformation exercised; corrected fixed-\(R\) identity still open; full matter-coupled momentum propagation and total sector conservation remain open; no production changes; \`scientific_admission=false\`.


## Addendum — cubic fixed-areal-radius shell crosswalk (2026-10-25)

**Evidence:** shell component-budget workflow [40206277038](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40206277038), N=60/120/240, using the same SHA-verified accepted t=0.60 profiles, source/vendor pins, and copied `dt)/`dt/2) probes. The protected production/vendor/dependency diff was empty and accepted-state immutability passed.

The shell observer distinguishes a linearly reconstructed face path from a four-point cubic fixed-(R) query. The cubic face lane independently reconstructs face (chi) and transformed (F=4\pi R^2T_t{}^R), then evaluates shell mass-rate minus endpoint net current at fixed areal boundaries. It is an integrated Misner–Sharp balance diagnostic, not the fluid HLL operator and not the production scalar update.

At the smaller timestep, the maximum absolute shell residuals over (R=1ldots10) were:

| N | Cell-centred cubic query | Face-reconstructed cubic query | Face reconstruction, linear query |
|---:|---:|---:|---:|
| 60 | 4.543e-5 | 4.233e-5 | 2.923e-5 |
| 120 | 3.219e-6 | 2.638e-6 | 1.061e-5 |
| 240 | 2.078e-7 | 1.561e-7 | 5.622e-6 |

The cell-centred cubic maximum decreases at observed orders 3.82 and 3.95. The face-reconstructed cubic maximum decreases at orders **4.00 and 4.08**. By contrast, the face-reconstructed *linear*-query maximum does not display that convergence trend (orders about 1.46 and 0.92 over these three points).

At N=240, the per-shell cubic face residuals are (-4.03\times10^{-8}), (-4.79\times10^{-8}), (+1.03\times10^{-7}), (+1.55\times10^{-7}), and (-1.56\times10^{-7}) for (R=1\text{–}2,2\text{–}3,3\text{–}5,5\text{–}8,8\text{–}10), respectively. The maximum face-minus-cell discrepancy between the cubic-query residuals is (5.18\times10^{-8}) at N=240; the much larger linear-query discrepancy is therefore strongly interpolation-sensitive.

**Interpretation:** the earlier face-shell-versus-cell-centred mismatch is not sampling-independent. With a consistent cubic query on both mass and current, the independent integrated shell witness converges approximately fourth-order over N=60/120/240, in line with the target-map manufactured interpolation controls and the centered fixed-(R) mass/current trend. This is evidence for a finite-grid interpolation remainder in the tested annulus, not a universal error bound or a declaration that cubic interpolation is the production conservation convention.

This narrows the next scientific gate: do not continue tuning face geometry or treating the linear-sampled residual as a physical source. Return to the pinned full matter-coupled momentum-constraint propagation identity—specifically the normalized Hamiltonian-gradient, reference-metric connection response, and sector source/exchange mapping that remain unresolved at fixed coordinate radii. The integrated-current, full matter-coupled propagation, total stress-energy conservation, and scientific-admission flags remain false.
