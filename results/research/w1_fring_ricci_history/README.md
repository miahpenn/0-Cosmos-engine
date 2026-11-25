# W1 stable Ricci trajectory and timestep check

The [live first-step witness](../w1_fring_live_ricci_precision/README.md) reduced fine-grid rough momentum by stabilizing the evaluation of the existing stencil on a, b and X in the primary L2 Ricci calculation. This extension applies all three changes during each of 40 live steps at dt=0.000375 on dr=0.025 and 0.0125. Each replay starts from the same pinned state as the first-step witness. The first step reproduces its archived stages and endpoint digest. Comparison uses the 40 unmodified accepted M profiles from the [production momentum history](../w1_fring_momentum_history/README.md) at identical times.

The observer retains the signed momentum residual and changes in the shear, K-gradient, remaining metric coupling and matter terms at the explicit predictor, primary predictor before boundary, primary predictor after lapse, corrector before boundary, and accepted endpoint. Changes are measured relative to each step's lapse-refreshed start. Band and rough RMS use the same eight common shells in [0.05, 0.25) and the same volume-weighted quadratic fit as the preceding witnesses. The normal production constraint evaluator reads all endpoints.

## Forty-step result

| dr | Production band M at t=0.015 | Stable band M | Production rough M | Stable rough M |
|---|---:|---:|---:|---:|
| 0.025 | 5.636e-13 | 5.811e-13 | 3.409e-13 | 3.690e-13 |
| 0.0125 | 8.121e-13 | 1.044e-12 | 5.869e-13 | 7.396e-13 |

The fine-grid rough RMS improves on the first two steps and is larger on steps 3 through 40, peaking at 56.0% above production on step 33 and ending 26.0% above it. On the reference grid it is larger at every accepted step and ends 8.2% above production. At the fine endpoint the number of shell-mean reversals falls from six to four while native-cell reversals rise from nine to ten. The rough norm grows despite fewer shell reversals.

On step 40 at dr=0.0125, the accumulated production-to-stable difference in the signed band profile is 7.241e-13 after the explicit predictor and 7.396e-13 at acceptance. The differences between that step's own signed updates, starting from each branch's refreshed state, have RMS 1.596e-18 and 3.517e-14 at the same stages. By the last step most of the separation is inherited from the preceding trajectory. The largest difference during the first step appears at the primary/corrector evolution stages, and the primary-to-refreshed and corrector-to-accepted local M readbacks agree exactly throughout the stable replay.

The combined change isolates the contribution of metric derivative arithmetic. Its short-lived benefit and the larger late rough residual leave a production repair unsupported by this test. A fixed-dt control for the pairing of shear and K stage increments would address the timing question more directly.

## Matched-time timestep check

From the same initializer, both production and the stable Ricci arm take either one dt=0.000375 step or two dt=0.0001875 steps. All four arms end at t=0.000375. These are early-time controls. Changing the number of steps also changes the sampled stage operands and the sequence of binary64 operations.

| dr | Steps to same time | Production rough M RMS | Stable Ricci rough M RMS |
|---|---|---:|---:|
| 0.025 | One full step | 3.811e-15 | 3.872e-15 |
| 0.025 | Two half steps | 4.003e-15 | 3.391e-15 |
| 0.0125 | One full step | 1.522e-14 | 1.225e-14 |
| 0.0125 | Two half steps | 8.466e-15 | 1.099e-14 |

Halving dt reduces production fine-grid rough M by 44.4%, while increasing the reference-grid rough value by 5.0%. On the fine grid the stable Ricci change reduces rough RMS by 19.5% at the original timestep and increases it by 29.8% at the half timestep. Its signed profile effect has band RMS 7.794e-15 at the original timestep and 5.866e-15 at half timestep. The production signed profile difference due to timestep alone has band RMS 1.275e-14 on that grid. Production fine-grid reversal counts stay at 3 across shell means and 12 across native cells in both timestep arms. The fine-grid production timestep difference decomposes into band RMS values of 1.399e-14 in the shear increment, 9.686e-15 in the K-gradient increment and 3.183e-16 in the matter increment. Their signed sum with the remaining metric term reproduces the 1.275e-14 total difference. The local primary and corrector M profiles agree exactly across the bracketed boundary/lapse refreshes within each short arm. The timestep dependence enters through the step updates. Identifying a particular stage-pairing error would require a control that changes the pairing independently of dt.

## Reproduction

    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
    out=$(mktemp -d)
    for dr in 0.025 0.0125; do
        for mode in trajectory timing; do
            python -m research.w1_fring_smoke.ricci_history "$dr" "$mode" "$out/$dr-$mode.json"
        done
    done

Both grids, the source and vendor, the initial state, the first-step report and the production history are pinned in the JSON. All 59 targeted checks passed.
