# W1 live first step with stable Ricci second derivatives

The [stored-input precision witness](../w1_fring_first_rate_precision/README.md) identified the metric second derivatives in tracefree Ricci as the largest evaluated arithmetic difference in the first fine-grid corrector. This replay branches one live step of dt=0.000375 from each pinned aligned initializer. Every branch starts from the same geometry, scalar, fluid and clock arrays. The production branch reproduces the archived first-step stage ledger exactly. All three branches reach the same end time without retries.

The research observer redirects the selected even-parity second derivatives during the existing primary L2 Ricci calculation. The original five-point weights are applied to differences from the central value. The final two one-sided rows use the original evaluator. One arm applies this to X, and the other applies it to a, b and X. The lapse second derivative uses the original evaluation, and every intercepted primary L2 call verifies that its K output is identical on the same inputs. The production constraint readbacks use the original operators. The alternate rates are applied during live stage evolution, including the primary predictor and corrector.

## Accepted first step

Band statistics use the existing eight common coordinate shells from r=0.05 to 0.25, with the same volume-weighted quadratic fit used to define rough RMS. Counts refer to adjacent nonzero signs across shell means and across native cells.

| dr | Arm | Band M RMS | Rough M RMS | Shell / native reversals |
|---|---|---:|---:|---:|
| 0.025 | Production | 1.214e-14 | 3.811e-15 | 0 / 0 |
| 0.025 | Stable X | 1.233e-14 | 4.560e-15 | 1 / 1 |
| 0.025 | Stable a, b and X | 1.217e-14 | 3.872e-15 | 0 / 0 |
| 0.0125 | Production | 1.927e-14 | 1.522e-14 | 3 / 12 |
| 0.0125 | Stable X | 2.152e-14 | 1.746e-14 | 3 / 12 |
| 0.0125 | Stable a, b and X | 1.804e-14 | 1.225e-14 | 2 / 10 |

On the fine grid, stabilizing a, b and X together reduces the rough band norm by 19.5% and the complete band norm by 6.4%. Stabilizing X alone raises those norms by 14.8% and 11.7%. The combined arm changes the reference-grid rough norm by +1.6%. The fine-grid rough norm with the combined arm is still 3.16 times its reference-grid counterpart. The local scalar and fluid momentum terms are identical at the accepted endpoint on both grids, and all reported H statistics match.

The signed live corrector response agrees with the isolated frozen metric second-derivative correction to band RMS 3.32e-21 at dr=0.025 and 1.75e-20 at dr=0.0125 in the combined arm. The preceding full frozen arithmetic shadow also widened the Hessian and other Ricci terms, which accounts for testing a different combination of corrections in that report. The combined live arm retains two shell reversals and ten native-cell reversals on the fine grid. This first step establishes an arithmetic contribution to the rough residual and leaves its growth over multiple steps open.

## Reproduction

    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
    out=$(mktemp -d)
    for dr in 0.025 0.0125; do
        python -m research.w1_fring_smoke.live_ricci_precision "$dr" "$out/$dr.json"
    done

The JSON retains the signed momentum stage ledgers, constraint extrema, frozen-to-live comparison, checkpoint digest and source pins. All 59 targeted checks passed.
