# W1 signed momentum budget through the PIRK stages

The diagnostic reads the unmodified production step on the pinned PR #584 engine and vendor d6052d6. It records the existing stage states from a research observer. It does not replace a derivative, modify a source, or add a projection.

The two matched grids are \(dr=0.025\) and \(dr=0.0125\) at \(\Delta t=0.000375\), using the same grid-aligned inputs as the prior spatial replay. The report traces steps 1, 20 and 40 through \(t=0.015\). It also restarts the two complete \(dr=0.05\) A6 snapshots at \(t=0.015\), advancing one step at \(\Delta t=0.0015\) and two at \(\Delta t=0.00075\), both to \(t=0.0165\). These earlier snapshots have different timestep histories and different starting residuals. The comparison reports changes from each run's own refreshed start.

The signed terms are \(D_r A_a\), \(3A_a/r\), \(-2D_rK/3\), the remaining metric couplings, \(-8\pi j_r^{\mathrm{scalar}}\), and \(-8\pi j_r^{\mathrm{fluid}}\). They are recorded before and after the boundary writes and at the relevant lapse and predictor stages. The report takes RMS only after summing signed changes at each sampled radius, so cancellation is preserved. The fine-grid band uses eight common coordinate shells spanning \(r=0.05\)–\(0.25\), with the fine cells averaged in pairs by volume. The \(dr=0.05\) companion uses its native cells in that radius interval.

The supplementary readback evaluates \(D_r(K_{\rm stage}-K_{\rm start})\) in 60 digits from the stored \(K\) profiles. This uses linearity of the fixed stencil to assess the numerically small change, while retaining the updates that were already stored. Comparing it to the difference of the two ordinary readbacks separates evaluation arithmetic from effects present in the saved \(K\) update. It does not identify which part of the CMC evolution or storage generated the spatial variation in \(K\).

## Matched grids at t = 0.015

The rows below are band RMS values of signed changes over the last sampled step on each grid, relative to its refreshed start. The no-\(K'\) change is a diagnostic sum of the other signed changes on the same states.

| Readback | \(dr=0.025\) | \(dr=0.0125\) |
|---|---:|---:|
| Explicit predictor \(\Delta M\) | \(3.1882\times10^{-9}\) | \(3.1943\times10^{-9}\) |
| Primary predictor \(\Delta M\) | \(1.3009\times10^{-12}\) | \(1.3036\times10^{-12}\) |
| Accepted change in \(D_rA_a+3A_a/r\) | \(3.1883\times10^{-9}\) | \(3.1942\times10^{-9}\) |
| Accepted change in total matter term | \(3.1881\times10^{-9}\) | \(3.1942\times10^{-9}\) |
| Accepted \(\Delta M\) without the \(K'\) change | \(1.4349\times10^{-13}\) | \(7.6506\times10^{-14}\) |
| Accepted change in \(-2D_rK/3\), 60-digit readback | \(1.3455\times10^{-13}\) | \(3.6458\times10^{-14}\) |
| Accepted \(\Delta M\), all terms | \(1.9772\times10^{-14}\) | \(5.2162\times10^{-14}\) |
| Total accepted band \(M\) | \(5.6363\times10^{-13}\) | \(8.1211\times10^{-13}\) |

The scalar momentum change accounts for nearly all the explicit predictor's temporary residual. The primary update brings the shear term into cancellation with it, and the corrector reduces the imbalance further. On the finer grid the accepted change *before* adding \(K'\) is smaller, but the \(K'\) change cancels less of it. The final signed update is therefore larger in band RMS. At the first native cell in the band, the coarse \(K'\) contribution is \(+4.637\times10^{-14}\) and the fine contribution is \(-2.825\times10^{-14}\). The report records the full shell profiles and preserves signed values. These first native cells sample slightly different radii. The comparison in the table uses common shells.

Reading the change in \(K\) before differencing it agrees with the ordinary two-readback subtraction to band RMS \(1.90\times10^{-17}\) and \(1.74\times10^{-17}\) in the final step on the reference and fine grids. The difference between the \(\Delta K'\) contributions is present in the stored stage fields at these steps. The source of that variation, including the spatial residual of the lapse solve and storage in the \(K\) update, remains unresolved.

The accepted, refreshed and boundary-bracketing momentum readbacks agree exactly at every sampled near-centre cell where only the lapse or the outer boundary is changed. The explicit predictor, primary predictor and corrector are where the local balance changes. The recorded steps sample the trajectory and do not telescope the complete forty-step residual history.

## A6 restart cross-check at dr = 0.05

The archived A1 and A2 runs repeat the same sequence. Starting from \(t=0.015\), their first explicit predictors have band RMS changes of \(1.265\times10^{-8}\) and \(6.324\times10^{-9}\), respectively. The corresponding fully accepted step changes are \(8.00\times10^{-13}\) and \(6.13\times10^{-14}\). The \(K'\) term contributes to the last cancellation in each case. Since the archived initial residuals differ, these are stage witnesses rather than a convergence-order test.

## Reproduction

The new observations match the archived cell-0 readbacks on both fine grids at \(t=0.0075\) and \(0.015\), and match all six A6 \(M\) samples at the start and end of each restarted step. All 83 research steps reached their targets without retries, and 37 targeted tests passed. The production sources and gate settings were preserved.

    git fetch origin refs/pull/598/head
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
    out=$(mktemp -d)
    for case in A1 A2 0.025 0.0125; do
        python -m research.w1_fring_smoke.stage_momentum "$case" "$out/$case.json"
    done
    python -m research.w1_fring_smoke.stage_momentum_report "$out"

The JSON files contain the source and input pins, script digest, 24 local-cell profiles at each sampled stage, signed differences, and budget closures. The fine input was created in 5b206037.
