# W1 predictor and corrector metric balance

This diagnostic replays the last substep approaching the first shared time \(t=0.00075\) at \(dr=0.025\). The nominal, half and quarter trajectories are the same corrected ring, eight-node initializer, open edge and KO-off runs used in the [quarter-timestep comparison](../w1_fring_quarter_trace/README.md). The last substep starts at a different time on each trajectory. Its predictor and corrector target the same endpoint.

## Continuous rate and discrete writes

With zero shift, put \(F=\alpha A_a\) and \(G=\alpha K/3\). The native rates are \(\dot a/a=-2F\), \(\dot b/b=F\), and \(\dot X/X=G\). Therefore the inverse radial metric \(B=X^2/a\) and the areal-radius factor \(U=R/r=\sqrt b/X\) satisfy
\[
\dot{\ln B}=2F+2G,\qquad \dot{\ln U}=F/2-G .
\]
Uniform expansion enters these two rates with opposite signs. In the homogeneous flat case \(a=b=1\) and spatially uniform \(X\), the coordinated change preserves \(B(R')^2=1\) and \(R''=0\), so it produces no spatial curvature. On the ring, curvature reads the spatial variation of the updated factors as well as their local values.

The predictor writes each metric field \(y\) as \(y_0+\Delta t f_0\), using the refreshed starting lapse. The corrector writes it as \(y_0+(\Delta t/2)(f_0+f_1)\), where \(f_1\) is evaluated on the primary predictor with its solved lapse and updated \(A_a,K\). The corrector increment minus the predictor increment is \((\Delta t/2)(f_1-f_0)\). In all three replays, both formulas reproduce every stored \(a,b,X\) cell exactly in binary64 using the native rates recomputed on copies of the corresponding stages.

## Signed centre budget

The entries below are nominal-minus-half or half-minus-quarter differences in the direct scalar curvature at cell 0. The calculation uses the nine-point stencil and varies \(a,b,X\) on frozen metric copies. Each field's contribution is averaged across all six orders of substitution, and the contributions sum to the full difference. The seven-point readback shows the same pattern.

| Pair and stage | Share from \(a\) | Share from \(b\) | Share from \(X\) | Total |
|---|---:|---:|---:|---:|
| Nominal minus half, predictor | \(-2.345\times10^{-11}\) | \(-3.411\times10^{-11}\) | \(+3.833\times10^{-12}\) | \(-5.373\times10^{-11}\) |
| Nominal minus half, corrector | \(-5.075\times10^{-15}\) | \(-8.247\times10^{-15}\) | \(+3.160\times10^{-12}\) | \(+3.147\times10^{-12}\) |
| Half minus quarter, predictor | \(-6.439\times10^{-12}\) | \(-8.599\times10^{-12}\) | \(-5.541\times10^{-13}\) | \(-1.559\times10^{-11}\) |
| Half minus quarter, corrector | \(-5.414\times10^{-14}\) | \(-8.896\times10^{-14}\) | \(-3.589\times10^{-12}\) | \(-3.732\times10^{-12}\) |

The two geometric curvature contributions on the half-minus-quarter predictor are \(+7.17185\times10^{-8}\) for \(B\) and \(-7.17341\times10^{-8}\) for \(U\). These opposing contributions use the same metric. In the field-by-field budget, changing \(X\) moves both factors together. The difference that survives the predictor is dominated by \(a,b\). The corrector reduces their signed combined contribution from \(-1.504\times10^{-11}\) to \(-1.431\times10^{-13}\), leaving \(X\) responsible for about 96% of the corrected nine-point difference.

## Last-write rounding

The extended-precision comparison first re-evaluates curvature with the stored metric values promoted to a 64-bit significand. It then repeats the readback with only the last metric addition widened, keeping the binary64 starting fields and increments as operands. For half minus quarter at the corrector, the promoted stored difference is \(-3.73245\times10^{-12}\) and the pre-write difference is \(-4.70605\times10^{-12}\). Rounding the last write adds \(+9.73605\times10^{-13}\), partially cancelling the difference. This test preserves the earlier trajectory and the computed increments. It leaves open the effect of rounding in earlier updates and how the spatial \(X\alpha K\) rate develops along each trajectory.

All seven short replay steps completed without retries. The two stencil widths reproduce the earlier signed stage differences exactly in binary64. The full arrays are used for the bitwise write checks. The report records the first six cells. All 37 targeted tests passed. Stage 0 remains stopped.

From the research revision with the pinned vendor present, run:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.w1_fring_smoke.metric_balance
