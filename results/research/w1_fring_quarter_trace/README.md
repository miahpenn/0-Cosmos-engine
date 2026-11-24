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

## Predictor and corrector at the first shared time

Short replays bracketed the stages of the last substep approaching \(t=0.00075\) on each trajectory. The last substep starts at a different time on each trajectory, so this compares stage states targeting the same endpoint. The accepted endpoints reproduce the recorded readbacks exactly.

For the half-minus-quarter comparison, the nine-point direct scalar differs by \(-1.5592\times10^{-11}\) at cell 0 as soon as the explicit metric predictor has run. Frozen substitutions attribute \(+7.17185\times10^{-8}\) to the inverse radial metric and \(-7.17341\times10^{-8}\) to the areal-radius factor. Each contribution is averaged over the two orders of substitution, so they sum to the full signed difference. The near cancellation is already present at the predictor.

The corrector reduces the total difference to \(-3.7325\times10^{-12}\) before the boundary write. At that stage the inverse-radial and areal-radius contributions are \(+1.2587\times10^{-12}\) and \(-4.9911\times10^{-12}\), respectively. The nominal-minus-half difference reverses sign across the corrector, from \(-5.3733\times10^{-11}\) to \(+3.1466\times10^{-12}\). The seven-point readings show the same pattern. The outer writes leave the sampled centre curvature differences unchanged in these brackets.

This localizes the first shared-time signal to the explicit metric predictor and its correction. The comparison leaves open why the balance among the final accepted radial and angular metric components varies with timestep.

## Reproduction

Run the following from commit `92d31efa` with the pinned vendor checked out.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.w1_fring_smoke.quarter_trace

The JSON records 40 accepted steps without retries, ten matched times, the original trace hashes and the shared source, vendor and input pins. The full pytest suite had 1,202 passes, two skips and two expected failures.

For the stage trace, run from the revision containing the stage observer.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.w1_fring_smoke.quarter_stage_trace

The stage replays completed seven accepted steps without retries. All 35 targeted tests passed with the stage observer present.
