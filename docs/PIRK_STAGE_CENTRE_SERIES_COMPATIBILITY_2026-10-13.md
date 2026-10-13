# PIRK Stage Centre-Series Compatibility Audit — 13 October 2026

## Status and provenance

**Status:** diagnostic harness PASS; physics/admission result remains open.

- Exact action run: [38326054469](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38326054469)
- Audited diagnostic HEAD: `10c7101f67f43b4cebce5c9dc34cc8f582933565`
- Frozen production source: `6d8b3f2790644591326df57c55908863f994cb0c`
- `engine/production_kernel.py` source blob matched the frozen commit exactly.
- Vendor submodule: `d6052d605673ce9d82cdc99b0df79855fa2ea215`
- Artifact: [pirk-stage-centre-series-compatibility.zip](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38326054469/artifacts/11780133843)
- Artifact ZIP SHA-256: `16c6bc678b376828309f65bacbe28370c8985bc877b39a341448d5174f2e412c`
- `report.json` SHA-256: `39ed6e71f7c4c9fe65995e3b4183d6715df5472c47921152509bf3859d60e7fa`
- `profiles.npz` SHA-256: `2ad0806ecbd1ad0a0364fadbbc4ce4185b0a4e26d9c6a2d07e5ee7ee27c2b90e`
- `profiles.csv` SHA-256: `00a9579e3b17b3e282c75057efe5506add045e77af421951b151eb010789accf`

The diagnostic captured 35 finite field/RHS records across five resolutions. At each of N=40, 60, 80, 120, and 160 it applied **one ordinary production step only**, with (r_{\max}=80), (dt=0.03\,dr), radiation enabled, and the frozen kernel/vendor above. Input-state numerical digests were unchanged; accepted geometry passed the positivity assertion. No physics source, production kernel, coefficient, limiter, projection, or admission threshold was modified.

The report's `PIRK_STAGE_CENTRE_SERIES_AUDIT_PASS` means the pinned instrumentation and one-step capture completed with finite arrays and an internally consistent record. It does **not** mean the center series, constraints, conservation, convergence, or production admission passed.

## Primary finding: the central ẊAₐ/r² RHS witness grows under refinement

Let (r_0=dr/2=40/N) be the first cell centre. The initial data have (A_a(r_0)=0), but the initial total geometry RHS has a nonzero first-cell ẊAₐ. The table reports the signed RHS quotient ((\partial_t A_a)_{i=0}/r_0^2), the reconstructed un-divided RHS, and the related even-series derivative witness.

| N | (r_0) | ((\partial_t A_a)/r_0^2) at initial RHS | Reconstructed ((\partial_t A_a)_{i=0}) | (D_r(\partial_t K)/r) at (i=0) |
|---:|---:|---:|---:|---:|
| 40 | 1.0000 | +1.104394e-3 | +1.104394e-3 | +5.10704e-4 |
| 60 | 0.6667 | +2.456064e-3 | +1.091584e-3 | +1.19501e-3 |
| 80 | 0.5000 | +4.356510e-3 | +1.089128e-3 | +2.14877e-3 |
| 120 | 0.3333 | +9.791530e-3 | +1.087948e-3 | +4.86944e-3 |
| 160 | 0.2500 | +1.740206e-2 | +1.087629e-3 | +8.67665e-3 |

Across a fourfold refinement from N=40 to N=160, the divided ẊAₐ/r² witness grows by 15.76 times, while the raw first-cell ẊAₐ remains near (+1.09\times10^{-3}). In the captured split, initial (A_a)-L2 contributes about (+1.107\times10^{-3}/r_0^2) at N=40 and (+1.7405\times10^{-2}) at N=160; the L3 contribution remains small, from (-2.65\times10^{-6}) to (-1.76\times10^{-7}) in the same divided quotient.

The centre-rate witness is present at the initial RHS and changes only slightly in the explicit-predictor and primary-predictor RHS. It persists in the read-only RHS recomputed after the accepted step. This rules against the PIRK time advance suddenly creating the original (A_a) centre mismatch; the nonzero central RHS is already in the initial semi-discrete operator. The (D_r(\partial_t K)/r) witness also grows with N, agreeing with the independently recorded central K-rate series concern.

## Fixed physical-radius context

The CSV and JSON reports now interpolate at the requested physical radii rather than choosing a different nearest cell at each resolution. Initial ((\partial_t A_a)/r^2) at fixed (r=3,5,8) is:

| N | r=3 | r=5 | r=8 |
|---:|---:|---:|---:|
| 40 | -4.67715e-6 | -9.77680e-6 | -4.29725e-6 |
| 60 | -8.88033e-6 | -1.10177e-5 | -4.18582e-6 |
| 80 | -1.58916e-5 | -1.12436e-5 | -4.12433e-6 |
| 120 | -1.77996e-5 | -1.12656e-5 | -4.07824e-6 |
| 160 | -1.79169e-5 | -1.12720e-5 | -4.06197e-6 |

The fixed-radius values settle at finer resolutions while the first-cell divided quotient grows. That is consistent with a localized central non-regular layer rather than a uniformly growing rate across the whole radial domain. This is a descriptive one-step refinement comparison, not a formal convergence proof.

## Second stage-localized witness: ẊΛ/r

At the first cell the ẊΛ/r RHS is at roundoff on the initial slice. It is already a small positive value at the explicit predictor, then becomes negative at the primary-predictor input and remains at nearly the same value in the post-step read-only RHS.

| N | Initial slice ẊΛ/r | Explicit predictor ẊΛ/r | Primary predictor ẊΛ/r |
|---:|---:|---:|---:|
| 40 | ∼0 | +1.33859e-5 | -2.78868e-5 |
| 60 | ∼0 | +9.25284e-6 | -5.52506e-5 |
| 80 | ∼0 | +7.01719e-6 | -8.00481e-5 |
| 120 | ∼0 | +4.71292e-6 | -1.26914e-4 |
| 160 | ∼0 | +3.54368e-6 | -1.72407e-4 |

At N=160 on the primary-predictor input, the signed split is approximately (L2=-1.76052\times10^{-4}/r_0) and (L3=+3.64488\times10^{-6}/r_0), leaving a total of (-1.72407\times10^{-4}/r_0). This identifies a second **stage-sensitive** connection-rate term. The measurement does not by itself determine whether its geometrical L2 contribution is wrong, or which coupled term is responsible.

## Implementation closure and boundary distinction

The code-path weighted-increment reconstruction closes as follows. `Aa`'s full-domain maximum is a distinct outer-boundary correction, not a central closure loss.

| Field | Largest max-absolute reconstruction residual across tested N | Location / interpretation |
|---|---:|---|
| a | 7.40e-15 | Roundoff scale |
| b | 3.68e-15 | Roundoff scale |
| X | 3.68e-15 | Roundoff scale |
| K | 1.16e-16 | Roundoff scale |
| Lambda | 1.33e-17 | Roundoff scale |
| Aa | 3.48e-8 | Outer-most cell at N=40; centre cells i=0..5 are at or below 2.71e-20 |

The `Aa` residual excluding the outer two cells is at or below (2.71\times10^{-20}) for every resolution. Therefore this reconstruction supports that the integrator is applying the coded stage weighting correctly in the interior. It does not certify that the RHS itself is regular or constraint preserving.

## Interpretation and next discriminating check

This diagnostic reinforces the already documented finding that the central ẊAₐ/r² RHS is not regular under this frozen baseline source. The centre mismatch is present before a step, persists across the actual PIRK stages, and is not explained by a failure of the weighted-increment implementation. The L2 term dominates the first-cell ẊAₐ split, but the present hook records the combined L2 value; it does not itself separate all internal Ricci and lapse-Hessian contributions.

Next check should capture a read-only signed split of those individual geometric pieces at the **same actual RHS states** (initial slice, explicit predictor, primary predictor), preserving the original returned sum exactly. Cross-compare the baseline and already-proposed uniform-quadrature arm only if that same-source comparison is still uncovered by the coordination board. Do not fit a coefficient, project the state, edit the production RHS, or run a long campaign on this evidence alone.

## Explicit limits

- No Hamiltonian/momentum/connection-constraint admission, Misner–Sharp balance/conservation check, high-resolution radiation gate, formal convergence study, or long-horizon run was evaluated here.
- This audit uses the frozen baseline source commit listed above, not the separate draft scalar inverse-metric candidate PR #15.
- The result is a stage-localization diagnostic only. It supports neither a new physics law nor production readiness.
