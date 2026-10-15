# Manufactured control for the centre-Ricci derivative channel

**Status:** research-only numerical control. No production equation, vendor source, or campaign setting is changed. This is not scientific admission.

## Purpose

PR #254 found a large derivative-only response between the discrete derivative of a sampled metric-Lambda reconstruction and an analytic product/chain-rule expression on saved evolution stages. Before interpreting that result, test the same operator pair on a smooth, explicitly origin-regular manufactured metric.

## Manufactured profile and method

Use cell-centred radii on a fixed domain, \(r_{\max}=80\), and

\[
q(r)=1+0.1\left(1-e^{-(r/10)^2}\right),\qquad
a=q^{-4/3},\quad b=q^{2/3},\quad X=q^{1/3}.
\]

At the origin, the analytic continuation has \(q(0)=a(0)=b(0)=X(0)=1\); the metric fields are even and the reconstructed conformal connection is odd. No engine initialization or evolution is run.

For the pinned fourth-order cell-centred operators, compare

\[
\frac{2}{3}X^2\left[D_4(\Lambda_{\rm metric})-
\left.\partial_r\Lambda_{\rm metric}\right|_{\text{product/chain rule}}\right].
\]

The chain-rule operand uses the same \(D_1(a), D_1(b), D_2(a), D_2(b)\) construction as the PR #254 diagnostic. The \(2X^2/3\) factor is the pinned vendor TF-Ricci response coefficient for the \(D_4\Lambda\) operand. This is specifically a representation-difference control, not a claim that either discrete representation is the uniquely correct production choice.

## Result

| \(N\) | \(\Delta r\) | Signed cell-0 response | Maximum absolute response on fixed \(5\le r\le10\) |
|---:|---:|---:|---:|
| 40 | 2.0 | \(-1.148772\times10^{-5}\) | \(1.256354\times10^{-5}\) |
| 80 | 1.0 | \(-7.748350\times10^{-7}\) | \(8.240115\times10^{-7}\) |
| 160 | 0.5 | \(-4.950193\times10^{-8}\) | \(5.281102\times10^{-8}\) |
| 320 | 0.25 | \(-3.111749\times10^{-9}\) | \(3.328383\times10^{-9}\) |

The fixed-annulus maximum decreases by factors 15.25, 15.60, and 15.87 on successive refinements. That is consistent with fourth-order decay for this smooth manufactured profile. The first-cell magnitude also decreases by factors 14.83, 15.65, and 15.91; because cell 0 is at \(r=\Delta r/2\), these are origin-regularity diagnostics, not a same-coordinate pointwise convergence fit.

## Interpretation and stop condition

- The tested operator pair can agree rapidly under refinement on this smooth regular profile; a nonzero finite-grid representation difference alone is not evidence that either expression is incorrect.
- This control does **not** explain the substantially larger centre response in the saved production stages. The manufactured metric is not a substitute for the captured \(a,b,X,\Lambda\) profiles, and the fixed-annulus control cannot establish trajectory convergence.
- Next, if artifact provenance permits, compare the same profile-level identities on source-pinned saved stages and include a term-by-term signed closure. If the stored fields do not distinguish loss of regularity, steep gradients, and operator incompatibility, retain the blocker and do not edit production physics.
- No formal global convergence claim is made from this single manufactured family. The associated regression requires only that the tested response decreases by more than a factor of eight per halving, with the fixed-annulus ratios remaining in a broad stable range.

## Reproduction

The standalone regression is `tests/test_centre_ricci_chain_rule_manufactured.py`. It imports only the pinned derivative operators from `tools/analyze_centre_ricci_derivative_channels.py`, uses NumPy, and takes no artifact or production-engine inputs.
