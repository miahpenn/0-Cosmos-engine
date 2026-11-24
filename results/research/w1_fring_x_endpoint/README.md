# W1 X update ledger at t = 0.015

Lumen's prediction before the run was that the computed increments would dominate the centre difference at this endpoint, fall by roughly four to eight on successive timestep refinements, and leave the write contributions near \(10^{-12}\). The replay tests that prediction on the same corrected W1 ring and \(dr=0.025\) grid as the [first shared-time ledger](../w1_fring_x_update/README.md).

The three runs use nominal, half and quarter timesteps, with 20, 40 and 80 steps to \(t=0.015\). The nominal and half trajectories are compared against the full archived readbacks. The quarter trajectory is compared against its archived readbacks through \(t=0.0075\). The replay checks the \(X\) predictor and corrector writes bit for bit across the full grid at every step. The same native stage rates, binary64 increments and signed corrector storage errors are used in the telescoping account. The existing predictor write errors can influence later computed rates and are part of the actual trajectory.

Each signed \(X\) difference is read on a frozen common \(a,b\) background using both direct areal stencils. The grouped rate row includes all three factors and the increment arithmetic. The grouped write row includes all corrector storage errors. The residual after subtracting both rows from the \(X\)-only readback is retained as closure. Production \(H\) and direct \(H\) use the accepted states. Comparing them with the \(X\) ledger also exposes contributions from the other evolved metric fields and the evaluator.

## Accepted centre comparison

| Cell-0 Hamiltonian residual | Nominal | Half | Quarter |
|---|---:|---:|---:|
| Production \(H\) | \(+1.40959\times10^{-10}\) | \(+2.70636\times10^{-11}\) | \(-2.10732\times10^{-11}\) |
| Seven-point direct \(H\) | \(+8.76846\times10^{-11}\) | \(+1.42074\times10^{-11}\) | \(-2.07536\times10^{-11}\) |
| Nine-point direct \(H\) | \(+8.96496\times10^{-11}\) | \(+1.49238\times10^{-11}\) | \(-2.20199\times10^{-11}\) |

The ratios between successive signed timestep differences are 2.37 for production, 2.10 for the seven-point diagnostic, and 2.02 for the nine-point diagnostic. The quarter-step residual crosses zero in each readback. These are single-time ratios and do not establish a timestep order.

## Signed decomposition of the direct curvature difference

Nine-point results at cell 0 follow. The \(X\)-only and remaining-metric rows sum to the full metric difference. The remaining-metric row collects \(a,b\) and their small cross effects with \(X\) under the chosen frozen background.

| Contribution | Nominal minus half | Half minus quarter |
|---|---:|---:|
| Computed \(X\) increments | \(-6.2795\times10^{-14}\) | \(-1.4965\times10^{-13}\) |
| Accumulated \(X\) writes | \(+6.2698\times10^{-12}\) | \(+1.6917\times10^{-11}\) |
| \(X\)-only total | \(+6.2070\times10^{-12}\) | \(+1.6767\times10^{-11}\) |
| Remaining metric | \(+6.8594\times10^{-11}\) | \(+2.0196\times10^{-11}\) |
| Full direct curvature difference | \(+7.4801\times10^{-11}\) | \(+3.6963\times10^{-11}\) |

The forecast fails in the \(X\) ledger. The accumulated write term dominates \(X\) on both refinements, and the computed increment contribution grows in magnitude by about \(2.4\times\) on the second comparison. In the first comparison the earlier writes contribute \(+1.022\times10^{-11}\) and the final write partially offsets them by \(-3.948\times10^{-12}\). In the second comparison both write contributions have the same sign. The seven-point stencil gives the same attribution.

The \(X\)-only share of the full direct curvature difference is about 8% for nominal minus half and 45% for half minus quarter. The large nominal-to-half change in production \(H\) is still present in the direct metric comparison, where most of that change comes from \(a,b\). The origin of the \(a,b\) difference is outside this \(X\) account. Storing an \(X\) departure could address the identified storage contribution. This run does not establish its effect on the coupled trajectory or the rest of the metric difference.

## Verification and reproduction

All 140 short replay steps completed without retries. The 37 targeted tests passed. The diagnostic's endpoint time check allows the spacing expected from repeated binary64 additions. The prescribed step counts, timesteps and engine settings are the original ones. Stage 0 remains stopped.

With the pinned vendor checked out, run:

    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
    output=$(mktemp -d)
    for divisor in 1 2 4; do
        python -m research.w1_fring_smoke.x_endpoint "$output" --divisor "$divisor"
    done
    python -m research.w1_fring_smoke.x_endpoint "$output"

The resulting JSON includes the reference and script hashes, every sampled step contribution, and both stencil readbacks.
