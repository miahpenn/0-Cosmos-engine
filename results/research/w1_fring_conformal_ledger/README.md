# W1 conformal metric determinant and tracefree ledger

The replay uses the same corrected ring, fine grid \(dr=0.025\), and endpoint \(t=0.015\) as the [X update ledger](../w1_fring_x_endpoint/README.md). The nominal, half and quarter trajectories take 20, 40 and 80 steps. The report compares the accepted readbacks with their archived counterparts and checks every predictor and corrector write for \(a,b\) bit for bit.

Let \(D=\ln a+2\ln b\) represent the conformal determinant and \(T=(\ln b-\ln a)/3\) the tracefree coordinate. With zero shift, \(\dot{\ln a}=-2\alpha A_a\) and \(\dot{\ln b}=\alpha A_a\), so \(\dot D=0\) and \(\dot T=\alpha A_a\). The inverse transformation is \(\ln a=D/3-2T\) and \(\ln b=D/3+T\).

For each corrector, the ledger uses four signed pieces in each coordinate. The tangent is the trapezoidal integral of the continuous logarithmic rate on the observed stage factors. The integration contribution is the difference between that tangent and the log change from the wide-precision, field-space trapezoidal update on the same stored stages. Increment arithmetic is the effect of rounding the computed binary64 increment. The write contribution is the effect of rounding the subsequent addition to the stored field. The four pieces telescope to the actual stored endpoint change from initialization. They use the stored predictor and primary-stage inputs, so the tangent includes any indirect effects of storage on those inputs.

The curvature attribution applies each difference in \(D,T\) to the finer endpoint \(a,b\) on a common frozen \(X\) profile, the arithmetic mean of both endpoint profiles. Both direct stencils are evaluated. The closure term measures the small difference between the sum of the individual curvature responses and the full frozen \(a,b\) response. The latter is also compared to the remainder from the earlier \(X\) ledger, which used a different frozen background for \(X\). The report combines the signed \(a,b\) and \(X\) write contributions as well as their sum of absolute component magnitudes.

## Determinant drift

All entries are maximum absolute changes from the initialized \(D\) across the grid at \(t=0.015\).

| Timestep | Measured drift | Corrector storage writes | Wide field-space integration |
|---|---:|---:|---:|
| Nominal | \(1.242\times10^{-15}\) | \(1.242\times10^{-15}\) | \(1.181\times10^{-21}\) |
| Half | \(2.126\times10^{-15}\) | \(2.126\times10^{-15}\) | \(1.477\times10^{-22}\) |
| Quarter | \(2.648\times10^{-15}\) | \(2.648\times10^{-15}\) | \(1.851\times10^{-23}\) |

The determinant tangent vanishes identically in the ledger. The field-space corrector produces a nonzero determinant defect on the frozen stage inputs, whose grid maximum falls by about eight per halving. The measured determinant drift grows with step count in this comparison and is dominated by the storage writes. The maximum storage-induced determinant drift over the first six cells stays near \(2.40\times10^{-17}\) for all three timesteps. Its differences have no resolved response at cell 0. Across the first six cells the largest nine-point determinant-write curvature difference is \(1.86\times10^{-14}\).

## Signed curvature at cell 0

These values use the nine-point direct stencil.

| Contribution | Nominal minus half | Half minus quarter |
|---|---:|---:|
| Determinant, all terms | Below readback resolution | Below readback resolution |
| Tracefree logarithmic rate | \(+6.5063\times10^{-11}\) | \(+1.8597\times10^{-11}\) |
| Tracefree storage writes | \(+3.5307\times10^{-12}\) | \(+1.5991\times10^{-12}\) |
| Closure | \(+4.75\times10^{-17}\) | \(+3.13\times10^{-19}\) |
| Total \(a,b\) curvature | \(+6.8594\times10^{-11}\) | \(+2.0196\times10^{-11}\) |

The integration and increment-arithmetic terms have no resolved curvature response at cell 0 for either coordinate. The trapezoidal tracefree rate accounts for about 95% and 92% of the two \(a,b\) curvature differences. The seven-point stencil reproduces the attribution. These shares identify the discrete rate contribution on recorded stages. They leave open the upstream cause of the difference in \(\alpha A_a\).

## Measured write contribution for A6

Adding the \(a,b\) and \(X\) write readbacks gives the following endpoint values. Cell 0 is the largest magnitude among the sampled centre cells in both comparisons. For these samples the absolute component sum has the same maximum as the signed sum.

| Direct stencil | Nominal minus half | Half minus quarter |
|---|---:|---:|
| Seven-point | \(+9.1168\times10^{-12}\) | \(+1.6572\times10^{-11}\) |
| Nine-point | \(+9.8005\times10^{-12}\) | \(+1.8516\times10^{-11}\) |

These are measured write contributions to signed direct-curvature differences at this grid, time and stencil. Their scope is narrower than a floor for the production Hamiltonian residual throughout an A6 run.

## Reproduction

All 140 steps completed without retries and reproduced the archived readbacks. All 37 targeted tests passed. The research commit leaves the engine and A6 settings intact.

    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
    output=$(mktemp -d)
    for divisor in 1 2 4; do
        python -m research.w1_fring_smoke.conformal_ledger "$output" --divisor "$divisor"
    done
    python -m research.w1_fring_smoke.conformal_ledger "$output"
