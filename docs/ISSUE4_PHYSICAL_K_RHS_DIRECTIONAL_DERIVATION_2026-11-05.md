# Issue #4 physical CMC K-RHS directional derivation

**Date:** 2026-11-05  
**Branch:** `research/issue4-coupled-geometric-matter-control-2026-11-05`  
**Status:** analytic derivation; implementation and physical fixture not yet complete  
**Scope:** research-only, instantaneous, zero-step. No production RHS, production CMC solver, vendor RHS, or stepper calls.

## 1. Physical K equation

Define
[
\mathcal Q=\frac32 A_a^2+\frac13K^2+
4\pi(\rho+p_r+2p_t),\qquad \nu=\frac{X^2}{a},
]
and
[
c=-\frac{D_1a}{2a}+\frac{D_1b}{b}
-\frac{D_1X}{X}+\frac2r.
]
On the zero-shift CMC slice,
[
K_{\rm RHS}(\alpha)=-\Delta_g\alpha+\alpha\mathcal Q
=-\nu(D_2\alpha+cD_1\alpha)+\alpha\mathcal Q.
]
This is the trace of the pinned primary L2 plus L3 equation: L2 supplies (-\Delta_g\alpha), while L3 supplies (alpha(3A_a^2/2+K^2/3)) and the six-sector matter source (4\pi\alpha(\rho+p_r+2p_t)).

## 2. Stage-aware lapse contract: keep input and output distinct

The inspected CMC implementation is stage-aware. Its target routine projects the actual K RHS using the **incoming** `geometry.alpha`; the lapse solver then solves for a new **output** lapse. Thus write these as (\alpha_{\rm in}) and (\alpha_{\rm out}), respectively:
[
\lambda=P_w[K_{\rm RHS}(\alpha_{\rm in})],
\qquad
L(a,b,X,\mathcal Q)\alpha_{\rm out}=-\lambda,
]
with the fixed centre regularity row and outer normalization. Here (P_w[f]=\sum_iw_if_i/W).

The source contract was inspected at repository commit `0c40172f32b72a8149549adadf63d668a2f24f3e`, file `engine/cmc_gauge.py`; the same file is present on this research branch. The source confirms:
- `target_kdot` forms L2 plus L3 K RHS from the incoming geometry lapse and projects it over all proper-volume cells;
- `solve_cmc_lapse` then builds the lapse operator and solves for the output lapse;
- (w_i=r_i^2\sqrt{a_i}b_iX_i^{-3}\Delta r);
- (c) is assembled using the grid's fourth-order parity-aware first derivative;
- the interior lapse matrix uses the three-point second-order radial Laplacian row with the centred (cD_1) coefficient;
- the centre row is (\alpha_0-\alpha_1=0), and the outer row fixes (\alpha_{N-1}=1).

Therefore the three-point interior matrix in `tools/issue4_cmc_lapse_tangent.py` matches the inspected CMC solver's matrix structure; it is not a guessed fourth-order replacement. Its coefficient derivative must nevertheless use the same fourth-order parity-aware (D_1), and all row placements must match the source.

**Do not silently identify** (\alpha_{\rm in}) and (\alpha_{\rm out}). They are separate arguments/outputs of the stage-aware map. A self-consistent analytic fixed-point control may deliberately impose (\alpha_{\rm in}=\alpha_{\rm out}), but that is an additional control condition, not an automatic property of the stage-aware implementation. The direction must carry both (\dot\alpha_{\rm in}) and (\dot\alpha_{\rm out}) until that condition is explicitly established.

## 3. Directional derivative of the physical K RHS

For a fixed grid and fixed linear discrete derivative operators,
[
\frac{\dot\nu}{\nu}=2\frac{\dot X}{X}-\frac{\dot a}{a},
]
[
\dot c=-\frac12\left(\frac{D_1\dot a}{a}
-\frac{(D_1a)\dot a}{a^2}\right)
+\frac{D_1\dot b}{b}-\frac{(D_1b)\dot b}{b^2}
-\frac{D_1\dot X}{X}+\frac{(D_1X)\dot X}{X^2},
]
[
\dot{\mathcal Q}=3A_a\dot A_a+\frac23K\dot K
+4\pi(\dot\rho+\dot p_r+2\dot p_t).
]
The physical projected K-RHS direction is
[
\boxed{
\begin{aligned}
\dot K_{\rm RHS}^{\rm in}={}&
-\dot\nu(D_2\alpha_{\rm in}+cD_1\alpha_{\rm in})\\
&-\nu(D_2\dot\alpha_{\rm in}
+cD_1\dot\alpha_{\rm in}+\dot cD_1\alpha_{\rm in})\\
&+\dot\alpha_{\rm in}\,\mathcal Q
+\alpha_{\rm in}\dot{\mathcal Q}.
\end{aligned}}
]
Here (dot K) inside (dot{\mathcal Q}) is the physical base (K_{\rm RHS}^{\rm in}). The six-sector source directions and all other dotted geometry/primary fields must be derived independently from the pinned instantaneous equations. No free `k_rhs_dot` profile is permitted.

## 4. Differentiate the proper-volume projection

With
[
w_i=r_i^2\frac{\sqrt{a_i}b_i}{X_i^3}\Delta r,
\qquad W=\sum_iw_i,
\qquad
\frac{\dot w}{w}=\frac12\frac{\dot a}{a}
+\frac{\dot b}{b}-3\frac{\dot X}{X},
]
the target derivative is
[
\boxed{
\dot\lambda=\frac1W\sum_i
\left[\dot w_i((K_{\rm RHS}^{\rm in})_i-\lambda)
+w_i(\dot K_{\rm RHS}^{\rm in})_i\right].
}
]
The projection includes the complete proper-volume domain, including boundary cells, exactly as the source implementation does. It is differentiated from the physical K RHS built with (\alpha_{\rm in}), not from the solved output lapse unless the fixed-point identity is part of the control definition.

## 5. Differentiate the output CMC lapse solve

The discrete output equation is
[
(D_2+cD_1-\mathcal Q/\nu)\alpha_{\rm out}
=-\lambda/\nu
]
on interior rows, with the fixed centre and outer rows above. Let its assembled matrix and RHS be (A) and (b). Its directional system is
[
\boxed{A\dot\alpha_{\rm out}
=\dot b-\dot A\alpha_{\rm out}.}
]
The interior RHS derivative is
[
\dot b_i=-\frac{\dot\lambda}{\nu_i}
+\frac{\lambda\dot\nu_i}{\nu_i^2}.
]
The matrix derivative contains (\dot c), (\dot{\mathcal Q}), and (\dot\nu), with the exact three-point row placement. The centre and fixed outer-normalization rows have zero directional RHS.

This is the stage-aware coupling: the incoming lapse direction enters (dot K_{\rm RHS}^{\rm in}), then (dot\lambda) enters the output lapse tangent equation. If the analytic control additionally imposes (\alpha_{\rm in}=\alpha_{\rm out}) and (\dot\alpha_{\rm in}=\dot\alpha_{\rm out}), these equations become one coupled linear system in the common lapse tangent and (dot\lambda). Otherwise retain the two lapse directions and do not force a fixed point.

## 6. Implementation and validation gates

1. Reconstruct the incoming physical K RHS from the independent L2 trace, L3 trace, and six-sector projections. Cross-check its algebraic form against the formula in §1 without importing or calling production/vendor routines.
2. Reconstruct the proper-volume projection and its derivative from the incoming lapse, including (\dot w) and (\dot K_{\rm RHS}^{\rm in}).
3. Bind the output lapse matrix to the inspected source stencil and exact row convention. Retain base and differentiated row-wise residuals.
4. Solve for (\dot\alpha_{\rm out}) from the differentiated output system. Preserve the distinction between input and output lapse directions.
5. If and only if the chosen analytic base is a fixed point, impose and test the input/output lapse equality and solve the corresponding coupled block system.
6. Independently construct plus/minus states at both required epsilons and compare the incoming K-RHS projection, target rate, output lapse, and output lapse tangent. Check positivity, parity, determinant/trace-free identities, six-sector projections, and all CMC residuals.
7. Only after these pass, construct the new immutable coupled archive and run the existing zero-step two-channel observer.

No production code, accepted archive, existing comparison class, predictor, admission flag, or conservation flag is changed by this derivation.


## 7. Independent stage-aware implementation checkpoint (2026-11-05)

Implemented in `tools/issue4_cmc_lapse_tangent.py`:
- `physical_k_rhs_and_direction` independently transcribes the native warped-product lapse-Laplacian identity using (A=a/X^2), (R=r\sqrt b/X), and the parity-aware (D_1,D_2) operators. It differentiates (A), (R), the lapse Hessian, and the matter/primary potential analytically.
- `solve_stage_aware_cmc_tangent` constructs incoming (K_{\rm RHS}) and its direction internally, projects the incoming source over the full proper-volume measure, solves the distinct output lapse, and differentiates that output solve. Its API has no `k_rhs_dot` argument.
- Incoming `alpha_in_dot) and output `alpha_dot) remain distinct. The existing output operator keeps the source-matched three-point interior row and fixed centre/outer rows.
- No engine, production, or vendor imports are used.

The bounded workflow run for commit `a21363a3557f16006bacedb31691b89ac85856c2` completed the compile, all seven tests, and the static no-production/no-vendor import guard successfully:
[CI run #42092298042](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42092298042).

**Scope limitation:** these new plus/minus tests use a manufactured smooth control and combined ((\rho,p_r,p_t)) input profiles. They verify the independent directional algebra and stage-aware solve, not yet the actual pinned accepted-state six-sector arrays. The next binding gate is to feed the locked six-sector projection outputs and independently derived coupled geometry/primary directions into this API, compare the reconstructed incoming (K_{\rm RHS}) against the pinned L2+L3 component ledger, and run the required two-epsilon plus/minus checks on the accepted (t=0.60) state. Do not claim that physical-state binding or coupled-fixture admission has passed yet.


### Instantaneous six-sector and coupled-geometry direction transcription (2026-11-05)

A new no-import research transcription has been added to `tools/issue4_cmc_lapse_tangent.py`:
- `independent_scalar_and_geometry_rates` derives the scalar-field rates, three Valencia conservative RHS profiles, native CMC ((a,b,X)) rates, pinned vendor L2/L3 ((A_a,K)) rates, and the (Lambda) L2/L3 rates directly from their source equations;
- `_independent_fluid_rhs` transcribes cell-centred primitive recovery, minmod reconstruction, face-metric interpolation, HLL fluxes, radial geometric sources, and the DM-only covector exchange/rest-current source. Baryons and radiation receive no DM exchange source;
- `six_sector_projection_directions` analytically differentiates each scalar and fluid sector's own ((\rho,p_r,p_t,j)) projection, including the metric/volume response and the implicit radiation primitive relation.

The source pins used are diagnostic source `4aec88d0e3ff028e0a65a8ac36b30948625df792`, vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`, and the scalar/Valencia source files at the diagnostic pin. The vendor primary L2 curvature brackets, L2 lapse Hessian, L3 primary sources, Lambda L2 momentum term, and Lambda L3 matter term are transcribed explicitly. The CMC explicit metric rates are the beta=0 source equations (\dot a=-2\alpha aA_a), (\dot b=\alpha bA_a), (\dot X=\alpha XK/3). No production RHS, production kernel, adapter, or vendor module is imported or called.

The first regression now compares the analytic six-sector projection directions with centered plus/minus reconstruction at both (\epsilon=2\times10^{-5}) and (10^{-5}). CI run [#42097574717](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42097574717) completed successfully: 11 tests passed, including the two-epsilon six-sector projection-direction comparison; the accepted-state base crosswalk and static no-production/no-vendor import guard also passed.

**Accepted-state stage-chain result and remaining lapse-input boundary:** CI [#42097909692](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42097909692) passed the accepted-state independent rate reconstruction and both stage-aware output-lapse finite-difference checks; the follow-up scope-explicit rerun [#42098254533](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42098254533) also passed both jobs at (\epsilon=2\times10^{-5}) and (10^{-5}). The stage-0 output lapse tangent has max absolute value (4.0622815467394005\times10^{-2}); it is propagated as a nonzero incoming-lapse direction to stage 1. Stage-1 K-RHS tangent max absolute value is (4.8344873607288966\times10^{-3}), stage-1 output-lapse tangent max absolute value is (4.0608853426564964\times10^{-2}), and stage-1 CMC base residual max absolute value is (2.656566373132607\times10^{-16}).

The accepted archive `accepted-state-N60.zip` contains only `tangent_arrays.npz` and `report.json`; it has the accepted (alpha) profile but no preceding CMC-stage input lapse or its directional history. Accordingly, the first archived incoming-lapse direction was anchored at zero for the first refresh. The subsequent stage uses the computed, nonzero stage-0 output direction. This verifies the two-stage propagation after that boundary, but it is **not yet the complete no-frozen-lapse chain** requested. Do not promote this to the final coupled-tangent gate.

**Next source-fidelity requirement:** recover the immediately preceding CMC-stage input lapse (or a source-pinned stage record from which its direction can be reconstructed) for the accepted (t=0.60,N=60) state. Then derive its directional response from the same instantaneous equations and rerun the two-epsilon stage chain without the zero incoming-direction anchor. If the upstream stage record cannot be recovered, retain this as an explicit missing archive datum; do not set (dot\alpha_{\rm in}=0) as a physical approximation, identify it with (dot\alpha_{\rm out}), or impose a fixed point.
 
### Discrete-operator crosswalk and source-matched implementation (2026-11-05)

The first accepted-state crosswalk was run against the immutable `accepted-state-N60.zip` nested in workflow artifact `13442329240` from run `42084531386`. The archived base arrays are at (t=0.60), (N=60), with the locked vendor commit `d6052d605673ce9d82cdc99b0df79855fa2ea215`.

The pinned source chain is:
- `engine/production_kernel.py` at diagnostic source `4aec88d0e3ff028e0a65a8ac36b30948625df792`, blob `1720615628a98a8d7cfc9fcadcc984fe28c88be1`: native CMC explicit metric RHS and stage dispatch;
- `engine/v55_pirk_adapter.py` at the same source commit, blob `be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3`: `primary_l2_rhs` plus `primary_l3_with_matter`;
- vendor `vacuum_bssn.py` at the locked vendor commit, blob `099ce52ba9e6dcd38d50c2915248c76c4236f7f4`: `geometry_terms()["lap_alpha"]`, `primary_l2_rhs()["K"]`, and `primary_l3_rhs()["K"]`.

For the CMC slice, (\\beta=0), and the vendor has (mathrm{KO\_EPSILON}=0). The exact pinned physical K RHS therefore has the fieldwise form
[
K_{\rm RHS}^{\rm pin}=
-\frac{X^2}{a}\left[
D_2\alpha-D_1\alpha\left(
\frac{D_1a}{2a}-\frac{D_1b}{b}
+\frac{D_1X}{X}-\frac2r
\right)\right]
+\alpha\left[\frac32 A_a^2+\frac13K^2
+4\pi(\rho+p_r+2p_t)\right].
]
The independent warped-product transcription (A=a/X^2, R=r\sqrt b/X) was evaluated on the same archived geometry. Although continuum-equivalent, its discrete lapse-Laplacian differs from the pinned fieldwise stencil by max absolute (6.892277181380808\times10^{-11}) and RMS (1.0987951371436368\times10^{-11}), with the largest difference in centre cells 0–3. This directly confirms the product/quotient finite-difference non-commutation; it is not a physical residual or a missing matter term.

The independent implementation has now been changed to transcribe the exact pinned fieldwise (D_1/D_2) expression and its analytic directional derivative. It still imports no engine or vendor code, and it retains separate incoming/output lapse directions. A regression test asserts the stage-aware K RHS against the pinned L2+L3 algebraic form. CI run [#42094234083](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42094234083) completed successfully after compile, all tests, and the static no-production/no-vendor import guard.

A strict six-sector binding layer has also been added: `sum_six_sector_projections` and `solve_stage_aware_cmc_tangent_six_sector` require exactly `dark_matter`, `baryons`, `radiation`, `S`, `D`, and `phi`, with base `rho/pr/pt/j` and directional `rho_dot/pr_dot/pt_dot/j_dot` profiles for every sector. Missing or extra sectors fail closed. The synthetic wrapper tests completed successfully in [CI run #42094841283](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42094841283): 9 tests passed, compile passed, and the static no-production/no-vendor import guard passed.

**Accepted-state base binding now passed:** `project_six_sector_stress_from_arrays` independently reconstructs all six sectors directly from the immutable accepted archive's scalar fields/momenta and three Valencia conserved triples. The accepted-state CI job passes those six named profiles into `solve_stage_aware_cmc_tangent_six_sector`, checks the exact pinned K RHS, and solves the distinct output CMC lapse. [CI run #42095400227](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42095400227) passed both jobs, including the no-production/no-vendor import guard.

Accepted-state (N=60, t=0.60) results:
- Six sector names reconstructed and bound: dark matter, baryons, radiation, S, D, phi.
- Total projection maxima: (max|\rho|=5.489752198434982\times10^{-5}), (max|p_r|=1.4287188104970626\times10^{-5}), (max|p_t|=1.4331527744886261\times10^{-5}), (max|j|=4.222814213833781\times10^{-6}).
- Pinned K-RHS maximum absolute value: (6.261589268138539\times10^{-5}).
- Independently solved output lapse minimum: (0.9989137986988448).
- Base CMC row residual maximum: (2.949978586061497\times10^{-16}).

**Remaining directional boundary:** the accepted-state run intentionally supplies zero sector and geometry directions as a harness control. It validates the six-sector base projection, source assembly, and output CMC base solve; it does **not** validate the physical accepted-state tangent. The six sector directional profiles ((\dot\rho,\dot p_r,\dot p_t,\dot j)), together with the coupled geometry/primary directions and (dot\alpha_{\rm in}), still need to be independently derived from the pinned instantaneous equations and bound at both required epsilons. Do not treat the zero-direction accepted-state control as the physical coupled fixture or as a conservation-closure result.



### Superseding scope correction — incoming lapse direction is unavailable (2026-11-05)

The earlier accepted-state two-stage witness in CI runs [#42097909692](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42097909692) and [#42098254533](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42098254533) used a zero incoming-lapse direction as the first-stage anchor. Those runs remain valid demonstrations of the differentiated CMC algebra **conditional on that imposed input direction**, and they show that a nonzero stage-0 output direction propagates through the next refresh. They do **not** establish the physical accepted-state stage-0 tangent or the complete no-frozen-lapse chain. The earlier wording that described this as an anchored first direction is superseded: the archived direction is **unavailable**, not zero.

The accepted-state workflow is therefore narrowed. It now:
- reconstructs the accepted-state CMC output lapse as a base solve;
- independently reconstructs scalar, Valencia, metric, primary, and Lambda instantaneous directions from the pinned equations;
- analytically differentiates all six sector stress projections, including metric/volume dependence;
- checks each accepted-state sector direction against centered plus/minus projection reconstruction at both \\(\\epsilon=2\\times10^{-5}\\) and \\(10^{-5}\\);
- records the preceding incoming-lapse direction as unavailable and does not compute or certify an accepted-state stage-0 or stage-1 lapse tangent.

The stage-aware CMC tangent API remains independently tested on controls where the incoming lapse direction is supplied. The accepted-state physical CMC tangent remains blocked until an actual preceding-stage directional record exists, or until a separately labelled self-consistent fixed-point control is constructed and verified. No fixed point is inferred from the archive. No production equations, accepted arrays, source pins, or admission flags are changed.
