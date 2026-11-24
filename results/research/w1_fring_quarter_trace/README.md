# W1 quarter-timestep comparison

The corrected W1 ring was replayed on the same fine grid, \(dr=0.025\), through \(t=0.0075\). The new timestep was \(1.875\times10^{-4}\), one quarter of the nominal \(7.5\times10^{-4}\). The run used the eight-node initializer, the open outer edge and KO off. The previously recorded nominal and half-step trajectories came from the [accepted-slice curvature comparison](../w1_fring_areal_trace/README.md). All three runs start from the same curvature and source readbacks. The signed differences below use cell 0 at the ten common positive times.

For a residual \(H\), define \(d_1=H_{\Delta t}-H_{\Delta t/2}\) and \(d_2=H_{\Delta t/2}-H_{\Delta t/4}\). A dominant second-order term gives \(d_1/d_2\) near four at a fixed time, with the same sign for \(d_1\) and \(d_2\). This comparison uses both direct areal stencils and the production Ricci readback.

## Matched endpoint

At \(t=0.0075\), the cell-0 residual changes direction between the two refinements in all three readbacks.

| Hamiltonian readback | Nominal | Half | Quarter | \(d_1/d_2\) |
|---|---:|---:|---:|---:|
| Production Ricci | \(+5.7594\times10^{-11}\) | \(+4.5659\times10^{-12}\) | \(+8.7744\times10^{-12}\) | −12.60 |
| Direct, seven-point | \(+2.8884\times10^{-11}\) | \(-1.3292\times10^{-12}\) | \(+1.2484\times10^{-12}\) | −11.72 |
| Direct, nine-point | \(+2.9005\times10^{-11}\) | \(-1.9655\times10^{-12}\) | \(+1.5997\times10^{-12}\) | −8.69 |

Both direct stencils show a sign reversal between \(d_1\) and \(d_2\) at five of ten matched times. Where the signs agree, the ratios are still dispersed. The production ratio is near four at three consecutive times, then becomes much larger when \(d_2\) approaches zero before reversing sign later. The data do not establish a second-order timestep contribution at this cell. The Stage 0 stop still applies.

## Reproduction

Run the following from this revision with the pinned vendor checked out.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.w1_fring_smoke.quarter_trace

The JSON records 40 accepted steps without retries, ten matched times, the original trace hashes and the shared source, vendor and input pins. The full pytest suite had 1,202 passes, two skips and two expected failures.
