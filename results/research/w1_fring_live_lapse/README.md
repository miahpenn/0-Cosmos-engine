# W1 one-step live lapse comparison

The [frozen corrector test](../w1_fring_matched_lapse/README.md) found that matching the CMC lapse stencil to the evolution operator removed most of the K-gradient contribution to the step. The resulting shear change canceled most of its effect on the complete momentum budget. This follow-up gives the predictor, metric, scalar and fluid updates one full PIRK step to respond.

The two W1 cases reproduce the original trajectory up to the last step from `t=0.014625` to `0.015`, with `dt=0.000375` on both grids. A1 restarts the archived `dr=0.05` state for one step from `t=0.015` to `0.0165`. For each case the production, refined three-point and matched fourth-order branches use the same checkpoint, checked by a digest over the geometry, scalar, fluid and clock arrays. The production readback is checked against the accepted archive. Each alternate branch starts a fresh step, reaches the same end time without a retry, and leaves the original checkpoint intact.

The alternate lapse is used in all four solves of the live step, including the final accepted refresh. Each solve recomputes the production proper-volume Kdot target on its own current branch state and uses the existing even centre reflection and unit outer lapse. The solver uses the departure from unit lapse with extended-precision residual refinement and returns a binary64 lapse, as in the frozen test. The three-point control uses that same solver and retains the production stencil. The outer boundary and all metric, scalar and fluid evolution code are shared across branches.

The near-centre band uses common coordinate shells from r=0.05 to 0.25 on the fine pair, with native cells over the same interval for A1. The stage ledger retains the six signed momentum terms, the independently evaluated total, their changes relative to each branch's refreshed start, and the signed differences between branches. The end ledger records the Hamiltonian and momentum constraints over the band, at cell 0 and as domain-wide maxima. The domain-wide raw field changes are in the JSON.

## Measured response

All momentum entries are band RMS values. The last two rows take the difference of signed profiles before computing RMS. The row immediately above them reports the change in the norm of the accepted profile as a percentage of production.

| Readback | `dr=0.025` | `dr=0.0125` | A1 `dr=0.05` |
|---|---:|---:|---:|
| Production step change in M | 1.977e-14 | 5.216e-14 | 7.996e-13 |
| Three-point control step change in M | 1.961e-14 | 5.156e-14 | 7.999e-13 |
| Matched step change in M | 2.061e-14 | 4.726e-14 | 7.999e-13 |
| Production accepted M | 5.636e-13 | 8.121e-13 | 8.749e-12 |
| Matched accepted M | 5.630e-13 | 8.109e-13 | 8.749e-12 |
| Relative change in accepted M norm | -0.118% | -0.148% | +0.0021% |
| Matched minus production, K-gradient term | 1.340e-13 | 3.698e-14 | 2.151e-12 |
| Matched minus production, total M | 3.008e-15 | 6.772e-15 | 2.197e-15 |

Each matched solve reduces the K-gradient contribution to the step, by factors of 81, 5.4 and about 2,900 in the order shown. The shear response offsets most of the difference in the signed budget. The accepted matter contribution to the between-branch momentum difference is 4.38e-18, 3.16e-18 and 7.89e-17 across the same grids. It is small relative to the net geometry response. The live and frozen matched momentum profiles differ in band RMS by 2.55e-15, 5.54e-15 and 1.22e-15, respectively.

At the finest grid, the matched branch reduces the RMS momentum change over the step by 9.4%, while the inherited residual dominates the accepted norm. The negative matched cell-0 momentum increment is smaller in magnitude and removes less of the positive residual already present at the start. The accepted cell-0 M therefore rises from 2.8602e-11 to 2.8752e-11 relative to production.

The finest-grid Hamiltonian readback changes more. In the matched arm, its band RMS falls from 2.970e-12 to 2.090e-12 and its signed cell-0 value moves from -1.679e-11 to +7.461e-13. The reference-grid band moves from 3.730e-12 to 3.336e-12, with the refined three-point control at 3.478e-12. These are late single-step Hamiltonian measurements. Persistence across checkpoints and separation of metric update from metric write rounding remain open.

## Reproduction

Run from this checkout with the pinned engine, vendor and fixtures.

    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
    out=$(mktemp -d)
    for case in A1 0.025 0.0125; do
        python -m research.w1_fring_smoke.live_lapse "$case" "$out/$case.json"
    done

All 59 targeted eight-node, W1 initializer and PIRK/CMC checks passed. This diagnostic samples one live step at each checkpoint. It measures immediate feedback and leaves long-time accumulation and gauge evolution open.
