# Manufactured momentum-constraint operator identity audit — 2026-10-12

## Result

- Workflow: [run 38288419513](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38288419513)
- Artifact: [manufactured-momentum-operator-identity](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38288419513/artifacts/11763654184)
- Artifact SHA-256: \`83c738421ca014b781d5a07d1218d3de116ddc43bea5bc4d762233bc1dfb2d77\`
- Production source imported read-only: \`6d8b3f2790644591326df57c55908863f994cb0c\`
- Grid: N=40, 60, 80, 120; \(r_{\max}=40\); fixed physical annulus \(3\le r\le6\)
- Evolution steps: none

The diagnostic samples smooth, regular, nontrivial metric and scalar profiles, computes the exact analytic radial derivatives independently, and compares them with the pinned production discrete operators. It records signed terms separately:
\[
M =
D_r A_a-\frac{2}{3}D_r K
+6A_a D_r\chi
+\left(A_a-A_b\right)\left(\frac{2}{r}+\frac{D_r b}{b}\right)
-8\pi j,
\quad A_b=-\frac{A_a}{2},\quad \chi_r=-\frac12\frac{D_rX}{X}.
\]
The scalar-current analytic reference uses \(j_S=-P_S S_r\), \(j_D=-P_D D_r\), and \(j_\phi=-\Pi\phi_r/(8\pi)\). Static fluid momentum is zero. The analytic connection \(\Lambda\) is built from the continuum metric derivatives so the discrete connection-constraint error can be measured independently.

## Refinement results

| N | \(\Delta r\) | \(M\) error RMS, \(3\le r\le6\) | \(C_\Lambda\) error RMS, \(3\le r\le6\) | signed \(M\) error, cell 0 | signed \(M\) error, cell 1 |
|---:|---:|---:|---:|---:|---:|
| 40 | 1.0000 | \(2.192\times10^{-6}\) | \(2.651\times10^{-6}\) | \(-6.138\times10^{-6}\) | \(-1.163\times10^{-5}\) |
| 60 | 0.6667 | \(5.117\times10^{-7}\) | \(6.234\times10^{-7}\) | \(-9.162\times10^{-7}\) | \(-2.222\times10^{-6}\) |
| 80 | 0.5000 | \(1.499\times10^{-7}\) | \(1.821\times10^{-7}\) | \(-2.277\times10^{-7}\) | \(-6.046\times10^{-7}\) |
| 120 | 0.3333 | \(3.009\times10^{-8}\) | \(3.658\times10^{-8}\) | \(-3.100\times10^{-8}\) | \(-8.804\times10^{-8}\) |

The descriptive fitted order from the fixed-annulus RMS errors is **3.927 for \(M\)** and **3.923 for \(C_\Lambda\)**. From N=40 to N=120, the corresponding RMS errors fall by factors of about 72.8 and 72.5. The centre-cell rows are at different physical radii as N changes, so they are reported as centre-adjacent witnesses, not as matched-radius order estimates.

The fixed-annulus RMS discretization errors in the largest-resolution case (N=120) split by signed term as follows:

| Signed term | RMS error |
|---|---:|
| \(D_r A_a\) | \(1.761\times10^{-8}\) |
| \(-\frac23 D_rK\) | \(8.806\times10^{-9}\) |
| \(6A_aD_r\chi\) | \(3.090\times10^{-10}\) |
| \((A_a-A_b)(2/r+D_rb/b)\) | \(2.056\times10^{-10}\) |
| \(-8\pi j\) | \(2.139\times10^{-8}\) |

Discrete signed-term reconstruction and the production momentum constraint agreed exactly at stored precision. The determinant residual was at most \(2.22\times10^{-16}\), and the static fluid-current check was zero.

## Interpretation and decision

This is a useful negative control: **the pinned production derivative/parity machinery and the momentum-constraint evaluator reproduce the analytic smooth-profile constraint terms with the expected near-fourth-order convergence on a fixed physical annulus.** It does not show a generic parity-sign defect or a broad inconsistency in the constraint evaluator.

It does **not** prove that the evolved fields preserve centre regularity, that the production RHS updates satisfy a discrete constraint-propagation identity, or that the full trajectory is scientifically admitted. The manufactured profiles are not claimed to solve the continuum Einstein–matter constraints, and no evolution was advanced.

Combined with the earlier one-step tangent and centre-series results, the remaining discriminator is more specific: test the evolved RHS compatibility of the centre-series regularity relations—especially the coupled \(A_a\), \(K\), metric, \(\Lambda\) and scalar-current tendencies—rather than repeat a generic parity or constraint-evaluator check. Keep production physics unchanged until that test isolates a particular discrete identity failure.
