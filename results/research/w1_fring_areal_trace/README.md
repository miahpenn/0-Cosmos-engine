# W1 accepted-slice curvature comparison

The comparison uses the same short runs as the [origin trace](../w1_fring_origin_trace/README.md), on the corrected ring (A) and background-only control (B). Every accepted slice through \(t=0.015\) is evaluated with production Ricci and with seven- and nine-point direct areal curvature. Each pair shares exactly one \(16\pi\rho\) and one kinematic array \(2K^2/3-3A_a^2/2\). The observer works from a copied metric and checks that its fields agree with the original after evaluation.

The evolved radial coordinate no longer equals the areal radius. The earlier zero-step formula assumes equality, so this comparison uses the actual \(R=r\sqrt{b}/X\) and \(B=X^2/a\) from each accepted metric. With \(d\ell=dr/\sqrt B\), the direct scalar is
\[
\mathcal R_{\rm metric}=\frac{2(1-R_\ell^2)}{R^2}-\frac{4R_{\ell\ell}}{R}.
\]
The implementation differentiates the small offsets \(R/r-1\) and \(B-1\), using an even extension across the origin and inward stencils at the edge. Both stencil widths recover zero for flat space in ordinary and non-areal coordinates, and \(6k\) for the constant-curvature profile \(R=r,\ B=1-kr^2\), within \(2\times10^{-14}\) on both grids.

The raw metric calculation is independent of the evolved connection. Production Ricci includes that connection. The trace also adds \(X^2D_r C\) to the direct scalar, where \(C=\Lambda-\Lambda_{\rm metric}\) and the derivative uses the production stencil. This term follows from the linear \(\Lambda'\) coefficient in the production scalar bracket. Comparing production with the connection-adjusted direct result measures the remaining difference between the two discrete curvature evaluations. Both direct forms are diagnostic readbacks.

## Ring centre at the common endpoint

All entries are signed values at cell 0 and \(t=0.015\). The nine-point result is used for the connection-adjusted column.

| Grid \(dr\) | \(\Delta t\) | Production \(H\) | Direct \(H\), 7 points | Direct \(H\), 9 points | \(X^2D_r C\) | Adjusted direct \(H\) |
|---|---:|---:|---:|---:|---:|---:|
| 0.05 | \(1.5\times10^{-3}\) | \(+1.4312\times10^{-10}\) | \(+8.4803\times10^{-11}\) | \(+8.5030\times10^{-11}\) | \(-3.5680\times10^{-11}\) | \(+4.9350\times10^{-11}\) |
| 0.025 | \(7.5\times10^{-4}\) | \(+1.4096\times10^{-10}\) | \(+8.7685\times10^{-11}\) | \(+8.9650\times10^{-11}\) | \(-3.1462\times10^{-11}\) | \(+5.8188\times10^{-11}\) |
| 0.025 | \(3.75\times10^{-4}\) | \(+2.7064\times10^{-11}\) | \(+1.4207\times10^{-11}\) | \(+1.4924\times10^{-11}\) | \(-6.6898\times10^{-12}\) | \(+8.2340\times10^{-12}\) |

The measured centre residual depends on both the accepted metric and the curvature evaluator. On the fine grid, halving the timestep reduces production \(H\) by \(5.21\times\) and nine-point metric-only \(H\) by \(6.01\times\). The direct metric-only residual retains \(7.473\times10^{-11}\) of the \(1.139\times10^{-10}\) production timestep difference, about 66%. The direct Ricci change itself is \(7.480\times10^{-11}\), while the fine-pair difference in \(16\pi\rho\) is \(7.582\times10^{-14}\) and the kinematic difference is \(6.058\times10^{-17}\).

The explicit connection term changes by \(-2.477\times10^{-11}\) from the half step to the nominal step. The production-minus-adjusted-direct gap changes by \(+6.394\times10^{-11}\). These two terms partially cancel, raising the production timestep difference by a net \(3.917\times10^{-11}\) relative to the metric-only comparison. At either timestep the explicit connection term reduces the direct residual at cell 0. It cannot account for the larger production value there.

At the same endpoint, the fine controls have cell-0 metric-only residuals smaller than \(7\times10^{-16}\). Away from the centre, the wider stencil can amplify edge noise. For the half-step background control, the nine-point maximum is \(7.06\times10^{-12}\) at the outer cell, while the seven-point maximum is \(3.00\times10^{-12}\) in the interior. The ring-centre conclusion uses signed, matched-cell values.

These measurements localize the difference to the evolved metric curvature and the remaining evaluator discrepancy. Two timestep sizes do not establish an order. The seven- and nine-point reconstructions use the same stored metric, so agreement between them is a stencil comparison. It does not provide an exact continuum value for the discrete slices. The Stage 0 stop still applies, including the separate runtime limit for the full-period measurement.

## Verification

All six arms reached the common endpoint with 140 accepted steps and no retries. Including initialization, 146 paired slice comparisons passed the signed residual check. The full pytest suite finished with 1,202 passed, two skipped and two expected failures.

## Reproduction

The code uses engine source `05c2a09aae5fc1224f110c055c0c5b4b3acf8b21` (PR #584), vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`, and corrected initial NPZ SHA-256 `339652a2a7e1379b19812893e43bafcf9ab0211422f2b8685c951107b863a295`. The engine, vendor and initializer metadata are checked against the smoke source. The original signed trace comes from `1377495f`. Every replayed row and first-step stage record is compared exactly against that trace. The JSON records the hashes of both diagnostic scripts and the prior trace for that arm. Runtime was Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, and one OpenBLAS 0.3.30 thread (SkylakeX).

From this research revision with the pinned vendor checked out, run:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
for arm in A B; do
    for dr in .05 .025; do
        python -m research.w1_fring_smoke.areal_trace "$dr" "$arm"
    done
    python -m research.w1_fring_smoke.areal_trace .025 "$arm" --half
done
```
