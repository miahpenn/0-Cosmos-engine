# W1 momentum history and accumulated corrector writes

The [late audit](../w1_fring_early_lapse/README.md) found six sign reversals across the eight fine-grid band shells, already present before the last step. This replay traces the same signed profiles from initialization through all 40 accepted production steps on `dr=0.025` and `dr=0.0125`, at the common `dt=0.000375`, ending at `t=0.015`. It keeps the production lapse and the existing eight-node, open-edge settings. The band is the same fixed coordinate interval `[0.05, 0.25)`, with volume-weighted common shells.

For each corrector, the observer collects the four actual PIRK operands for `K` and `Aa`. It verifies each stored update exactly against the binary64 operation order, including the near-centre equality before and after the outer boundary write. The ledger separates the computed update, rounded increment arithmetic, and rounding in adding the increment to the stored field. Each piece is summed cell by cell in 60-digit arithmetic. Its total telescopes to the exact stored field difference.

The momentum budget applies the existing fourth-order even-parity derivative to each accumulated field contribution. The Aa term uses the metric coefficients at the current endpoint and accounts separately for the coefficient change acting on the initial Aa. Signed scalar and fluid momentum changes are taken from the same states. The budget reports the signed closure against the ordinary production readback, where differentiating the full background K introduces a small evaluation difference.

The write contribution is a direct attribution on the observed trajectory. All captured RHS values already contain the effects of earlier writes on the predictor, lapse, metric, and matter. The ledger also leaves those effects in any computed term. Subtracting the direct `K` and `Aa` write contribution from the readback is a diagnostic of its size on these states.

## Results

| Readback | dr=0.025 | dr=0.0125 |
|---|---:|---:|
| Step 1 accepted band M RMS | 1.214e-14 | 1.927e-14 |
| Step 1 reversals across shell means / native cells | 0 / 0 | 3 / 12 |
| Step 40 accepted band M RMS | 5.636e-13 | 8.121e-13 |
| Step 40 reversals across shell means / native cells | 1 / 1 | 6 / 9 |
| Step 40 accumulated K + Aa write response RMS | 1.543e-16 | 1.883e-16 |
| Step 40 write response / accepted band M RMS | 0.0274% | 0.0232% |
| Largest write response / accepted band M RMS over all 40 steps | 0.145% | 0.139% |

On the finer grid, the shell reversal count is five by step 16 and first reaches six at step 38. The finer grid already shows twelve native reversals at the first accepted step. Subtracting the accumulated direct write response leaves the sign of every native band cell and every shell mean at all 40 accepted steps on both grids. The direct write response also remains below 0.33% of the rough component on both grids, where rough is the remainder after removing the same volume-weighted quadratic shell fit used in the frozen derivative report.

The field increments telescope in extended precision. The reconstructed momentum budget agrees with the ordinary readback within 3.09e-17 on the reference grid and 6.17e-17 on the fine grid at their largest pointwise discrepancies. The independent late-witness last-write values are reproduced exactly. The measured residual starts early and persists independently of the direct K and Aa write attribution. The precision of the sampled RHS operands and feedback through the other stored fields remain open.

## Reproduction

    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
    out=$(mktemp -d)
    for dr in 0.025 0.0125; do
        python -m research.w1_fring_smoke.momentum_history "$dr" "$out/$dr.json"
    done

The JSON retains radius, provenance, 40 accepted samples on each grid, the signed cumulative term profiles, the per-step write response and shell signs. The production archive is checked at the four shared output times. Fifty-nine targeted tests passed.
