# W1 fixed-timestep corrector endpoint comparison

The [matched-time comparison](../w1_fring_ricci_history/README.md) found timestep sensitivity concentrated in the shear and K-gradient increments. This test keeps dt=0.000375 and isolates one feature of the PIRK corrector. Production evaluates the right-hand L2 rate on the final explicit metric copy and uses the primary predictor for the right-hand L3 rate. The research arms use the already computed primary L2 value for Aa, K, or both. The left-hand terms and right-hand L3 values follow the production calculation.

The first test runs one live step per arm from identical pinned initial data. The production endpoint digest reproduces the previous witness. The intervention verifies identical stage ledgers through the refreshed primary predictor and changes exactly one corrector evaluation. The late test works from the captured production operands of step 40 ending at t=0.015, reconstructs the two stored Aa/K writes exactly, and evaluates the changed momentum on copies of the same corrector. The late readback holds metric coefficients and matter momentum fixed. It measures a local response without changing the history leading into that step.

## Momentum response

Both grids use the existing eight common shells from r=0.05 to 0.25. Rough RMS is the residual after the existing volume-weighted quadratic shell fit. In the following table the step change is measured relative to the lapse-refreshed start. The paired arm changes both right-hand L2 values.

| Grid | Step | Production step change RMS | Paired step change RMS | Production rough M RMS | Paired rough M RMS |
|---|---|---:|---:|---:|---:|
| 0.025 | 1, live | 1.216e-14 | 5.083e-15 | 3.811e-15 | 2.668e-15 |
| 0.0125 | 1, live | 1.928e-14 | 1.135e-14 | 1.522e-14 | 1.087e-14 |
| 0.025 | 40, frozen | 1.977e-14 | 1.751e-14 | 3.409e-13 | 3.398e-13 |
| 0.0125 | 40, frozen | 5.216e-14 | 4.023e-14 | 5.869e-13 | 5.851e-13 |

The live fine-grid rough norm falls by 28.6%, while its shell reversals increase from three to five and its native reversals stay at twelve. The fine-grid late step change falls by 22.9%, while the rough norm of the total stored slice falls by 0.3%. Cell 0 moves from -1.679e-12 to -1.072e-12 on the fine first step and from +2.860e-11 to +2.922e-11 on the fine late slice. On the reference first step the cell-0 magnitude increases from 6.732e-14 to 2.292e-13. The band and centre therefore respond differently.

## Ricci contribution and comparison with halving the step

The operand copies share lapse, shift, evolved connection and primary Aa/K; their explicit metric components differ. Widening both L2 evaluations to 60 digits gives a fine-grid Ricci contribution to the pairing difference of 1.429e-14 on step 1 and 1.500e-14 on step 40. The respective residuals between binary64 and widened responses have RMS 4.583e-15 and 1.005e-14. The initial and late fine-grid K and shear-Hessian effects from changing operands are at or below 1.1e-20. The K-only stored correction is at or below 2.7e-17 across all four tests, consistent with sensitivity of the final addition at this scale. The substantial Ricci response survives higher arithmetic precision on the same stored inputs.

We compare the first-step pairing response against the ordinary change produced by taking two half steps to the same endpoint. Both differences use the same initializer and same production full-step reference. The normalized volume-weighted dot product across signed shell means is 0.970 at dr=0.025 and 0.906 at dr=0.0125. Equality would give +1 and zero residual difference. The differences between the response profiles have band RMS 4.691e-15 and 6.209e-15, respectively. The pairing therefore reproduces much of the direction of the early timestep response, with a substantial remaining difference. This is a single substitution in the archived scheme. Its multi-step response and temporal order remain open.

## Reproduction

    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
    out=$(mktemp -d)
    for dr in 0.025 0.0125; do
        python -m research.w1_fring_smoke.capture_lapse_stages "$dr" "$out/$dr.npz"
        python -m research.w1_fring_smoke.corrector_pairing "$dr" "$out/$dr.npz" "$out/$dr.json"
    done

The JSON retains the signed geometry/constraint responses, the widened Ricci/Hessian split, and the source, checkpoint and comparison hashes. All 59 targeted tests passed.
