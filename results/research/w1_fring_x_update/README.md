# W1 accumulated X update ledger

This read-only replay follows every predictor and corrector on the fine ring grid to the first shared time \(t=0.00075\). There are one, two and four accepted steps at the nominal, half and quarter timestep, respectively. The [previous last-write comparison](../w1_fring_metric_balance/README.md) established that the final addition partially cancels the corrected centre difference. This ledger measures the inherited write errors and the computed increments separately.

The observer checks every predictor and corrector write for \(X\) bit for bit using the native metric rates. Every accepted readback matches the archived trajectory exactly. The ledger subtracts the exact stored operands of the corrector addition from the value written by binary64. It telescopes these signed storage defects across the trajectory. Adding the computed increments and the storage defects to the common initial profile reproduces the final \(X\) profiles. The predictor write errors are recorded too. They can affect the later rate because the primary predictor is part of the actual trajectory.

The exact product \(X\alpha K/3\) of the stored stage factors is evaluated with a 64-bit significand. Changes from a common initial rate are assigned to \(X,\alpha,K\), with each interaction split equally among its factors. The ledger also retains the difference between the wide factor product and the binary64 trapezoidal increment. An identical initial-rate integral cancels between timestep runs at the same endpoint.

For each component of the endpoint \(X\) difference, the readback adds its entire radial profile to the finer trajectory's \(X\) on a frozen common \(a,b\) background. The common background is the arithmetic mean of the two endpoint \(a,b\) profiles. Curvature is read with both direct areal stencils. The sum includes a closure term for finite-precision and nonlinear cross effects. The complete trajectory rates are held fixed for this accounting. A different internal storage representation could change those rates and requires a separate evolved diagnostic.

## Signed curvature at cell 0

The nine-point entries are differences between timestep runs using the frozen common \(a,b\) background.

| Contribution | Nominal minus half | Half minus quarter |
|---|---:|---:|
| Earlier corrector writes | \(+2.417\times10^{-12}\) | \(-4.355\times10^{-12}\) |
| Last corrector write | \(+7.548\times10^{-13}\) | \(+7.658\times10^{-13}\) |
| Computed increments, all rate factors and increment arithmetic | \(-1.143\times10^{-14}\) | \(-1.536\times10^{-15}\) |
| Closure | \(-3.966\times10^{-16}\) | \(+1.609\times10^{-15}\) |
| Total \(X\)-only curvature difference | \(+3.160\times10^{-12}\) | \(-3.589\times10^{-12}\) |

For half minus quarter, the earlier corrector writes produce a larger negative difference than the endpoint value. The last \(X\) write offsets \(+7.66\times10^{-13}\) of it. This is the \(X\)-only result. The earlier diagnostic measured \(+9.74\times10^{-13}\) when widening the last writes to all three fields. The seven-point result has the same structure, with \(-4.035\times10^{-12}\) from earlier writes and \(+6.289\times10^{-13}\) from the last one.

Across the first six cells, the half-minus-quarter lapse and expansion-rate contributions to the \(X\) profile vary smoothly. The accumulated write difference ranges from \(-1.57\times10^{-16}\) to \(+1.77\times10^{-16}\) and changes sign across these cells. The curvature derivatives amplify this variation. The localization applies at the first shared time. It does not establish the dominant contribution at later times or determine what happens when metric storage is changed in the evolving kernel.

## Reproduction

All seven short replay steps completed without retries. The 37 targeted tests passed. The engine source, vendor, fixture and script hashes are recorded in the ledger. Stage 0 remains stopped.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.w1_fring_smoke.x_update_ledger
