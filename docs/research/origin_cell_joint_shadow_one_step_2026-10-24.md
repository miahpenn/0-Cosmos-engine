# Joint regular-origin shadow: one accepted production step

**Status:** diagnostic-only; no production authorization.  
**Workflow:** [run 40153240438](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40153240438)  
**PR:** [#416](https://github.com/miahpenn/0-Cosmos-engine/pull/416)  
**Artifact:** [12583486611](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40153240438/artifacts/12583486611)  
**Artifact ZIP SHA-256:** `d7130043ba7b69873544c747b717d30cf9ed45fc7ae18c52aa38f01852e9a652`

## Provenance and exact baseline check

- Source: `4aec88d0e3ff028e0a65a8ac36b30948625df792`
- Vendor: `d6052d605673ce9d82cdc99b0df79855fa2ea215`
- CI Python 3.12.15, NumPy 2.5.4, SciPy 1.18.1.
- The run's fail-closed checks verified all protected engine blobs, the vendor checkout and gitlink, and an empty production/vendor/dependency diff.
- Runtime configuration: (N=60,120,240), (r_{\max}=80), amplitude (0.01), width (7), (D_{\mathrm{amp}}=10^{-10}), radiation enabled, (Delta t=0.03\,\Delta r). One accepted step only.
- The baseline arm reproduces the archived first accepted (A_a) profile to roundoff: maximum full-grid absolute differences are (2.26\times10^{-17}), (3.61\times10^{-16}), and (1.69\times10^{-21}) at (N=60,120,240), respectively.

Original workflow report SHA-256: `40a2440274a0d1fea169c50eaf77cc6bb45d209528299ffc0b6abc74138a8424`.  
Original `profiles.npz` SHA-256: `244291ed596268d56225f2400d17261427d72e659aacf5f5e17a4365b1930e74`. These hashes are from the downloaded CI artifact's `SHA256SUMS.txt`.

## Four-arm one-step comparison

Residual maxima below use the full saved arrays and are reported over (r\leq3); in these runs, the same maxima are also the full-grid maxima. (C_\Lambda=\Lambda-\Lambda_{\mathrm{metric}}). (A_a/r_0^2) is the signed cell-zero value.

| N | Arm | (A_a/r_0^2) | (max|H|) | (max|M|) | (max|C_\Lambda|) |
|---:|---|---:|---:|---:|---:|
| 60 | Baseline | +9.7874e-5 | 3.9110e-3 | 1.5695e-4 | 6.3830e-6 |
| 60 | Metric shadow only | +2.9813e-6 | 1.0637e-5 | 6.3934e-6 | 2.6106e-7 |
| 60 | CMC-row shadow only | +9.4191e-5 | 3.9099e-3 | 1.6262e-4 | 6.6229e-6 |
| 60 | Joint metric + CMC row | -8.0528e-7 | 1.1681e-5 | 4.3440e-7 | 1.8000e-8 |
| 120 | Baseline | +1.9518e-4 | 4.0086e-3 | 1.5491e-4 | 3.1527e-6 |
| 120 | Metric shadow only | +7.2241e-6 | 1.6304e-6 | 6.0751e-6 | 1.2398e-7 |
| 120 | CMC-row shadow only | +1.8755e-4 | 4.0076e-3 | 1.6054e-4 | 3.2690e-6 |
| 120 | Joint metric + CMC row | -4.5879e-7 | 8.7683e-7 | 1.5207e-8 | 3.3217e-10 |
| 240 | Baseline | +3.9024e-4 | 4.0321e-3 | 1.5451e-4 | 1.5725e-6 |
| 240 | Metric shadow only | +1.5188e-5 | 1.0105e-6 | 6.0894e-6 | 6.2123e-8 |
| 240 | CMC-row shadow only | +3.7484e-4 | 4.0310e-3 | 1.6012e-4 | 1.6302e-6 |
| 240 | Joint metric + CMC row | -2.3364e-7 | 5.6846e-8 | 5.1304e-10 | 8.0313e-12 |

Cell-zero accepted (A_a) in the joint arm is (-3.5790\times10^{-7}), (-5.0977\times10^{-8}), and (-6.4899\times10^{-9}) at (N=60,120,240). Its magnitude falls rapidly as (\Delta r) is halved; meanwhile the baseline (|A_a|/r_0^2) doubles at each refinement. This is a strong one-step regularity signal, not a formal convergence/admission pass.

The single-factor controls discriminate the interaction: metric-only substantially reduces (H) but leaves a material momentum/connection residual; CMC-row-only leaves the large Hamiltonian and connection residual essentially unchanged. The joint arm improves all four requested centre measures.

## Initial-slice crosswalk used to select the combined shadow

This table is from the separate direct-Ricci four-arm initial-data shadow (same source/vendor pins). It is **not evolved**:

| N | Saved production (L_2) | Original metric + regular CMC row | Fine-integrated metric + pinned row | Fine metric + regular row |
|---:|---:|---:|---:|---:|
| 60 | +1.09282e-3 | -5.85526e-4 | +3.02673e-5 | -1.14673e-5 |
| 120 | +1.08826e-3 | -5.70900e-4 | +3.98145e-5 | -2.44855e-6 |
| 240 | +1.08751e-3 | -5.67554e-4 | +4.18149e-5 | -5.82833e-7 |

The fine integration was checked from 512 to 2048 substeps per coarse cell; the largest reported shift in the combined-shadow cell-zero (L_2) was approximately (2.03\times10^{-9}). The two factors interact: neither alone restores the expected centre source scaling.

## Interpretation and limits

1. This one-step run uses the untouched `V55ProductionKernel.step()`, including its production RHS, PIRK ordering, outer-boundary calls, and accepted-state construction.
2. Metric-shadow arms rebuild (a,b,X,\Lambda) from the continuous-source metric while holding (A_a,K), scalar fields and conserved matter variables fixed; matter projections are recalculated on the rebuilt metric. This is not a full metric-matter fixed-point initialization.
3. The CMC-row shadow changes only the central row in a copied lapse solve; it is a diagnostic sensitivity control, not a production gauge authorization.
4. One accepted step does not establish long-time stability, full matter-coupled momentum-constraint propagation, or total stress-energy conservation. Those admission gates remain open.

`scientific_admission=false`. Production and vendor sources remain unchanged. Do not merge this diagnostic as a production change.
