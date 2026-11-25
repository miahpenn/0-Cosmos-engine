# W1 first-step rate precision

The [40-step momentum history](../w1_fring_momentum_history/README.md) places the first fine-grid alternation at the first accepted step. This comparison captures the first production step on both pinned aligned inputs at dt=0.000375. All production momentum stage readbacks reproduce the archived first-step rows exactly.

On copies of the start, explicit predictor and final L2 operand, the diagnostic evaluates the lapse Hessian and tracefree Ricci terms with 60-digit arithmetic and the existing fourth-order interior stencils. It holds the cell values, parity, evolved connection and grid fixed. The same L2 expression reproduces the production binary64 components bit-for-bit when fed the production derivative values. L3 is tested with the geometry and stress projections held at their captured values, and reconstruction of its production output is also checked bit-for-bit. The diagnostic therefore tests arithmetic after its operands have been stored. It leaves their prior discretization and storage errors in place.

The signed rate differences enter the same primary and corrector weights as production. Applying `-2 D_r / 3` to the K difference and the full linear shear readback to the Aa difference gives a local contribution to M. This is a frozen signed readback with the metric coefficients and matter momentum held fixed. It does not re-evolve the predictor or store a different corrector. The band uses eight common coordinate shells in [0.05, 0.25). The rough RMS is the residual after the existing volume-weighted quadratic shell fit.

## First corrector

All numerical values below are band RMS. Rate differences are expressed in momentum units after multiplying by the corresponding step weight and applying the momentum operator.

| Readback | dr=0.025 | dr=0.0125 |
|---|---:|---:|
| K lapse-Hessian correction | 7.435e-16 | 6.618e-15 |
| Aa lapse-Hessian correction | 1.015e-15 | 6.297e-15 |
| Signed joint lapse-Hessian correction | 5.349e-16 | 1.439e-15 |
| Joint Ricci metric second-derivative correction | 1.383e-15 | 7.794e-15 |
| X second-derivative contribution alone | 1.310e-15 | 9.565e-15 |
| All evaluated rate corrections, signed sum | 1.183e-15 | 7.558e-15 |
| Production M RMS | 1.214e-14 | 1.927e-14 |
| M RMS with all corrections | 1.195e-14 | 1.786e-14 |
| Production rough M RMS | 3.811e-15 | 1.522e-14 |
| Rough M RMS with all corrections | 3.631e-15 | 1.189e-14 |

The lapse-Hessian changes largely cancel between K and Aa. Adding only that correction leaves all fine-grid shell and native reversal counts at three and twelve, respectively. Including all evaluated rate terms reduces rough fine-grid M by 22%, and the sign reversals fall to zero across the shell means and ten across the native cells. The sign-count change is more pronounced than the reduction in rough amplitude. The largest individual correction comes from the second derivative of X in tracefree Ricci. The corrections from a, b and X must be summed with sign.

A binary64 control evaluates the same even-parity five-point second derivative using differences from the central value before combining weights. In the isolated metric second-derivative contribution to M, this control agrees with the widened derivative to within 1.1e-23 band RMS in both stages on both grids. An exactly representable quartic checks the even reflection and the interior stencil. This identifies a local, algebraically equivalent candidate for a live first-step replay.

## Reproduction

    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
    out=$(mktemp -d)
    for dr in 0.025 0.0125; do
        python -m research.w1_fring_smoke.first_step_rate_precision "$dr" "$out/$dr.json"
    done

The JSON includes the signed momentum responses for the primary predictor and corrector and pins the production source, vendor, input, observer script and reference data. All 59 targeted tests passed.
