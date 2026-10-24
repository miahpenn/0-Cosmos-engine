# Conservation-ledger source crosswalk — 23 October 2026

**Status:** research-only derivation map; not a conservation proof or scientific admission.  
**Base inspected:** `main` at `4aec88d0e3ff028e0a65a8ac36b30948625df792`.  
**Coordination:** Issue [#5](https://github.com/miahpenn/0-Cosmos-engine/issues/5); density-normalization gate [#210](https://github.com/miahpenn/0-Cosmos-engine/issues/210); fixed-areal-radius lane [#35](https://github.com/miahpenn/0-Cosmos-engine/issues/35), [#79](https://github.com/miahpenn/0-Cosmos-engine/issues/79), [#89](https://github.com/miahpenn/0-Cosmos-engine/issues/89).

## Purpose and non-goals

Identify the existing equation/diagnostic sources from which a physically interpretable ledger can be derived, state the conditions under which each identity applies, and prevent a coordinate-fixed shell-tail proxy from being mistaken for conserved flux.

This change adds documentation only. It does not alter production equations, introduce a source law, rerun simulations, close Issue #4, or change scientific-admission flags.

## 1. Separate the three balances

These are distinct checks and must not be collapsed into one residual.

1. **Covariant total stress-energy balance:** establish the applicable identity for the *total* tensor, including any derived exchange/source law and all geometric terms.
2. **Finite-volume matter-stage budget:** reconcile stored conserved variables with signed face fluxes, metric sources, and explicit exchange sources for the actual discrete operator.
3. **Spherical Misner–Sharp/worldtube balance:** compare the geometric mass change along a specifically defined worldtube with the stress-energy flux/work expression appropriate to that same worldtube.

A passing result in one lane does not prove the other two. In particular, the existing total-`T` observer continues to report `independent_total_stress_energy_divergence_closed=false`; the Valencia stage budget explicitly disclaims being an independent proof of `div(T_total)=0`.

## 2. Existing source map

| Source | What it contributes | Limit to preserve |
|---|---|---|
| `docs/RESEARCH_SOURCE_CROSSWALK_2026-10-15.md` | Archive-to-implementation evidence and unresolved normalization/provenance gates. | It explicitly does not find an action-bearing primary source fixing the local S/D/φ stress-energy normalization and gravitational prefactor. |
| `tools/audit_total_stress_energy_divergence.py` (merged diagnostic PR #286) | Independently reconstructs total stress projections and evaluates a frozen-state covariant divergence using pinned production sources. | Its documented measured slice has `beta^r=0`; it is diagnostic-only and does not accept the physical exchange law. It is not by itself an evolving-state conservation proof. |
| `tools/audit_valencia_total_stage_budget.py` | Separates signed face-flux divergence, metric sources, and conditional exchange sources in the pinned matter RHS. | It expressly does not claim an independent proof of total covariant stress-energy closure. |
| `tools/audit_issue4_production_hll_operator_pairing.py` and the Issue #4 artifacts | Existing source-pinned operator/constraint evidence to reuse if the conservation ledger later feeds the momentum-propagation identity. | Issue #4 still lacks its complete convention-pinned matter-coupled propagation RHS; this ledger must not be presented as closing it. |
| Library archive note `GEAR_V5_5_CMC_BRIDGE_GATE_2026-10-05.md` and companion `GEAR_V5_5_CMC_BRIDGE_GATE_2026-10-05_v2.md` | Record an actual V5.5 stress projection `(rho, p_r, p_t, j_r)`, geometric Misner–Sharp mass, and a beta=0 spherical balance evaluated on the same state; report residuals at N=40/60/80. | These are valuable archived witnesses, not a substitute for rerunning/rehashing their artifacts against the current pinned source or resolving the fixed-coordinate versus fixed-areal-radius surface definition. |
| Issue #210 | Specifies the missing provenance/dimensional mapping between archive COSMOS densities and local `κ=8πG` Einstein equations. | Open. Do not infer the conversion from a convenient factor or wire in a reciprocal exchange source to make the ledger close. |


## 2A. Historical audit pin versus crosswalk-base source check

A read-only blob comparison was made against the exact crosswalk base `main` snapshot `4aec88d0e3ff028e0a65a8ac36b30948625df792`. The comparison uses the source blob pins declared by `tools/audit_total_stress_energy_divergence.py`; it is a provenance check, not a numerical replay.

| Pinned path | Historical audit blob | Blob at crosswalk base | Result |
|---|---|---|---|
| `engine/production_kernel.py` | `9685e76bb5a6f26d6bf383e240de24f37523f265` | `1720615628a98a8d7cfc9fcadcc984fe28c88be1` | changed |
| `engine/matter_system.py` | `e52a117690253ddad7b24d116145a65474a351d2` | `0969f4d996e4f8e0fc5b70a7dff147f024907302` | changed |
| `engine/matter_rhs.py` | `f41180efbcff2c72cf36477a765318ab539fac80` | `0153c232ee15b5a9d3dbac37ea858e4d2d730011` | changed |
| `engine/scalar_system.py` | `ffc119b1a6f2d495d972323901f76a85783b0a1b` | `e2c2dae9a016b4d6454cf74a9eae8b0167d5733d` | changed |
| `engine/v55_matter.py` | `0ee3c51fac7dc46c63feba2747a750542bcee3d0` | `0ee3c51fac7dc46c63feba2747a750542bcee3d0` | match |
| `engine/v55_pirk_adapter.py` | `eeee0492d80889d70e52b0946563eee5054c0f11` | `be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3` | changed |
| `engine/cmc_gauge.py` | `7b49b5751abebd4faf0130208032cd2080b71f64` | `7b49b5751abebd4faf0130208032cd2080b71f64` | match |
| `engine/v55_initial.py` | `7db2e1a10cf16fb40e4221b0df504fdf39aae611` | `7db2e1a10cf16fb40e4221b0df504fdf39aae611` | match |
| `engine/stress_energy.py` | `ac3e7518c62bb72e2218a61489fa486bb787809f` | `ac3e7518c62bb72e2218a61489fa486bb787809f` | match |
| `engine/valencia.py` | `03ed86398f424c03d9d6a3e646da9222a34fac5c` | `03ed86398f424c03d9d6a3e646da9222a34fac5c` | match |
| `engine/worldtube.py` | `aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6` | `aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6` | match |

**Interpretation:** five of eleven source blobs differ at the crosswalk base. This does not show that the newer sources are incorrect, nor that the old residuals are invalid for their original source. It does mean that the old audit's strict provenance guard will not accept the current tree as the historical experiment. Do not silently update the old pins or describe a run against changed sources as a reproduction.

**Initial semantic triage of those diffs:**

- `engine/matter_system.py` and `engine/matter_rhs.py`: the matter RHS now includes a variable-mass rest-current source `q_rest`; the momentum exchange projection also carries an explicit lapse factor. These directly affect conservative source bookkeeping.
- `engine/scalar_system.py`: the scalar–DM source contribution to coordinate-time evolution now includes `geometry.alpha`, changing the pinned RHS relative to the historical audit.
- `engine/production_kernel.py` and `engine/v55_pirk_adapter.py`: the production geometry stage now supplies a native CMC explicit metric RHS for `a/b/X`, while retaining other source blocks. This changes the geometry/time-derivative path used by a current-source witness.

These are equation/RHS or geometry-path changes, not merely report-format differences. The hash comparison alone cannot determine whether they improve or worsen physical closure; that requires a like-for-like term audit.

**Bounded next action:** preserve the old witness in its exact historical source/runtime environment and create a separately named current-source witness only after mapping these changed terms into the continuum balance. Give each witness its own full source/vendor/runtime pins. Keep the historical and current-source evidence records distinct.

## 2B. Equation-to-source reconciliation at the inspected source base

This is a read-only source reconciliation against the exact base recorded above, 4aec88d0e3ff028e0a65a8ac36b30948625df792. Links below pin that snapshot. No numerical audit was rerun for this addition. The archived scalar–dark-matter conditional projection is already derived in [conditional_scalar_dm_exchange_projection.md](conditional_scalar_dm_exchange_projection.md); reuse that derivation and its regression rather than create a competing version.

| Ledger term | Current implementation path at the inspected base | What is established in code | What is not established |
|---|---|---|---|
| Total covariant balance | [Frozen-state observer](../../tools/audit_total_stress_energy_divergence.py#L3-L15); [canonical stress adapter](../../engine/stress_energy.py#L56-L86); [total projection](../../engine/v55_matter.py#L138-L224) | The observer evaluates the mixed-tensor identity \(\nabla_\mu T^\mu{}_\nu=(\sqrt{-g})^{-1}\partial_\mu(\sqrt{-g}T^\mu{}_\nu)-\tfrac12T^{ab}\partial_\nu g_{ab}\). The adapter assembles the scalar and fluid projection, with an optional EM projection that must already be in consumer-normalized units. | The observer is frozen-state and documented for \(\beta^r_{\rm shift}=0\). Its historical source pins differ from five current source blobs; it is not a current-source replay or evolving-state proof. Its recorded closure/admission fields remain false. |
| Dark-matter exchange source | [Covector source](../../engine/valencia.py#L323-L331); [matter source projections](../../engine/matter_system.py#L500-L606); [conservative RHS](../../engine/matter_rhs.py#L13-L47); [scalar RHS](../../engine/scalar_system.py#L53-L92) | The implementation computes \(Q_t=\beta_{\rm DM}\rho_{\rm DM}\partial_t\phi\), \(Q_r=\beta_{\rm DM}\rho_{\rm DM}\partial_r\phi\), then projects energy and radial momentum with the coded lapse/shift factors. The matter rest-current source \(q_D\) is a separate term. The scalar momentum-variable RHS adds \(+\alpha\beta_{\rm DM}\rho_{\rm DM}\). These terms are documented and covered by the existing conditional projection regression, including nonzero shift and non-unit lapse. | This proves a conditional source-expression/implementation cross-check, not the action premise, the physical value of \(\beta_{\rm DM}\), sector-normalization provenance, or evolving total-\(T\) closure. Keep \(\beta_{\rm DM}\) (coupling) distinct from \(\beta^r_{\rm shift}\) (coordinate shift). |
| Finite-volume matter balance | [Face/source construction](../../engine/matter_system.py#L500-L606); [semi-discrete assembly](../../engine/matter_rhs.py#L13-L47) | The returned derivative is separated into signed face-flux differences, metric energy/momentum sources, explicit exchange energy/momentum sources, and the variable-mass rest-current source. The interior face flux uses the production reconstruction/HLL path; the boundary faces are separately assigned. | A cellwise RHS decomposition does not alone prove the integrated physical conservation law. The ledger still needs consistent cell measure, face orientation/area, boundary contributions and the cross-sector/geometric terms derived from one convention. |
| Stress projection and normalization | [Scalar projection](../../engine/scalar_system.py#L95-L112); [fluid projection](../../engine/matter_system.py#L381-L401); [total sector combination](../../engine/v55_matter.py#L138-L224) | In the scalar projection, S and D use their archive-scaled expressions while the COSMOS \(\phi\) projection is divided by \(\kappa=8\pi\). Conventional fluid projections are divided by \(\kappa\). The total projection sums scalar and fluid pieces, plus optional pre-normalized EM arrays. | This is an observable implementation asymmetry, not proof of a coding error. The action-bearing source and physical interpretation needed to authorize the S/D/\(\phi\) mapping remain open under Issue #210. Do not “correct” the factor by inspection or fit it to a residual. |
| Geometry time tangent | [Production stage wrapper](../../engine/production_kernel.py#L224-L250); [native CMC explicit RHS](../../engine/native_geometry/evolution.py#L305-L340); [PIRK adapter seam](../../engine/v55_pirk_adapter.py#L94-L135) | The production CMC path builds a matching native grid and supplies the native explicit RHS for \(a,b,X\), while the adapter retains the separate primary/Lambda and matter-source blocks. The native API is specifically CMC/zero-shift and does not return the gauge RHS. | A current-source \(d_t M_{\rm MS}\) calculation must differentiate the mass using this actual stage RHS and the same state/time level as the stress terms. It cannot silently reuse a historical geometry tangent or mix stage and accepted-state operands. |
| Worldtube mass/flux | [Worldtube observables](../../engine/worldtube.py#L23-L93) | The code calculates \(R=r\sqrt b/X\), \(\chi\), and \(m_{\rm MS}=R(1-\chi)/2\). The summarize helper chooses the cell nearest a requested **coordinate** radius; its stress-flux and work-term arguments default to zero and are not derived there. The current-residual helper differentiates supplied sampled masses and subtracts a supplied RHS. | The `worldtube.py` helper itself is not a fixed-areal-radius conservation instrument. However, the repository already has a more complete historical implementation in [audit_em_total_energy_budget.py](../../tools/audit_em_total_energy_budget.py#L169-L220) and its fixed-R rate snapshot (around lines 2247 onward): it interpolates at fixed areal radius, records the coordinate speed of the tube, includes mass-gradient/sweep terms, and reports both product-of-interpolated and composite-interpolated flux/work operators without choosing between them by residual size. That is prior work to reuse, not recreate. Crucially, this audit has a strict historical source-pin set and vendor pin, so its results cannot be presented as a witness of the newer five-file source base without a separate source reconciliation. The unresolved task is to map the existing fixed-R operands to the current coupled matter/scalar source terms and current geometry tangent, preserving both interpolation operators as diagnostics. |



### 2B. Existing fixed-\(R\) operand map: present formula versus current-source witness

The older tools/audit_em_total_energy_budget.py already supplies the needed *shape* of a fixed-areal-radius ledger. This is a code-path map, not a claim that the historical output is valid evidence for the current source base.

| Signed operand or operation | Existing audit expression | Current implementation correspondence | Status / guardrail |
|---|---|---|---|
| Areal geometry | \(R=r\sqrt b/X\), \(R_r=D_rR\), \(\chi=\gamma^{rr}R_r^2-(RK^\theta{}_\theta)^2\), \(m_{\rm MS}=R(1-\chi)/2\) | engine/worldtube.py computes \(R,\chi,m_{\rm MS}\); engine/native_geometry/evolution.py supplies the CMC explicit metric RHS used by engine/production_kernel.py and engine/v55_pirk_adapter.py. | Formula exists in current production code, but worldtube.py itself selects a coordinate-nearest cell and does not produce the fixed-\(R\) tangent. |
| Coordinate-time areal-radius tangent | \(\dot R|_r=R(\dot b/(2b)-\dot X/X)\); \(\dot R_r=D_r(\dot R)\) | Current CMC explicit RHS returns stage \((\dot a,\dot b,\dot X)\); derive \(\dot R\) from those same-stage operands and apply the repository's parity-aware radial derivative. | Mapping is available when the captured RHS and state are synchronized. The native CMC API requires \(\beta^r_{\rm shift}=0\); do not generalize this witness to a nonzero-shift stage. |
| Inverse-radial-metric tangent | \(\gamma^{rr}=X^2/a\), \(\dot\gamma^{rr}=\gamma^{rr}(2\dot X/X-\dot a/a)\) | native_cmc_explicit_metric_rhs produces \(\dot a,\dot X\) through production_kernel.py; the expression is reconstructed by the existing audit helper. | Algebraic reconstruction can be compared at one stage. A historical residual is not a current-source result until the complete source/vendor pins match. |
| Extrinsic-curvature tangent | \(K^\theta{}_\theta=K/3-A_a/2\), \(\dot K^\theta{}_\theta=\dot K/3-\dot A_a/2\) | engine/v55_pirk_adapter.py::geometry_stage_terms keeps primary_l2 and matter-coupled primary_l3 as separate blocks; the matter primary_l3 uses the normalized total \((\rho,p_r,p_t,j)\) from engine/v55_matter.py. | This is a particularly important stage boundary: sum the actual K and Aa RHS pieces at the same stage; do not substitute an accepted-state difference or only the vacuum RHS. |
| Coordinate-fixed mass tangent | \(\dot m|_r=\tfrac12\dot R(1-\chi)-\tfrac12R\dot\gamma^{rr}R_r^2-R\gamma^{rr}R_rD_r\dot R+R^2\dot R(K^\theta{}_\theta)^2+R^3K^\theta{}_\theta\dot K^\theta{}_\theta\) | Already expanded and retained term-by-term by _misner_sharp_chain_rule_rate_profiles in the historical audit. | This is a differentiated definition of mass, not by itself a conservation law. Keep all five signed terms and the decomposition residual visible. |
| Moving-surface correction | \(v_R\equiv (dr/dt)|_R=-\dot R/R_r\); \(\dot m|_R=\dot m|_r+v_R D_rm\) | Existing audit records tube_coordinate_speed_for_fixed_R, tube_mass_gradient, and tube_coordinate_sweep; it also keeps a piecewise-linear bracket/secant witness for the moving interpolant. | Requires \(R_r\ne0\). The interpolant/secant and fourth-order derivative are different discrete operators and must not be silently substituted for each other. |
| Signed stress flux and radial-pressure work | Factorized: \(F=-4\pi R_*^2\alpha_*\gamma_*^{rr}j_*R_{r,*}\), \(W=-4\pi R_*^2p_{r,*}\dot R_*\). The alternate operator interpolates the local profiles \(-4\pi R_i^2\alpha_i\gamma_i^{rr}j_iR_{r,i}\) and \(-4\pi R_i^2p_{r,i}\dot R_i\). | Current total projections are assembled by engine/v55_matter.py from scalar/fluid/optional normalized EM pieces; engine/stress_energy.py is the shared total-stress representation. | Preserve both interpolation orders already exposed by _fixed_R_composite_flux_work_profiles and _sector_fixed_R_rates. Neither is selected because it gives a smaller conservation residual. Normalization provenance under Issue #210 remains unresolved. |
| Matter finite-volume and exchange terms | Signed face differences plus metric source and explicit exchange source; the DM rest-current source is separate. | engine/matter_system.py::_species_rhs_components returns face flux, metric energy/momentum sources, q_e, q_s, q_rest; engine/matter_rhs.py::species_rhs assembles their signed semi-discrete derivatives. engine/valencia.py::dark_matter_covector_source supplies the conditional Q_t,Q_r. engine/scalar_system.py::scalar_rhs_arrays adds the lapse-weighted scalar source to dot(Pi). | These terms must appear in the matter/scalar sector ledger and any derivation of the evolving total-T identity. They are **not** extra terms to append directly to dot(m_MS)|_R by inspection: the mass-rate formula uses geometry and the total stress projection. Exchange cancellation remains conditional on the source law and consistent sector normalization. |
| Spatial mass-gradient identity / constraints | The audit compares \(D_rm_{\rm MS}\) with \(4\pi R^2(\rho D_rR+j\dot R/\alpha)\) and retains Hamiltonian/momentum-constraint contributions and an operator remainder. | Current stress projections and geometry constraint blocks exist, but the crosswalk has not established a same-stage, convention-pinned total identity on the current accepted state. | Use this as a named diagnostic decomposition, not as a source correction or a way to absorb the remaining residual into a fitted term. |

The corresponding fixed-\(R\) chain rule is therefore already encoded in the historical audit:

\[
\left.\frac{dm_{\rm MS}}{dt}\right|_R
=\left.\frac{\partial m_{\rm MS}}{\partial t}\right|_r
-\frac{\partial_t R}{\partial_r R}\,\partial_r m_{\rm MS}.
\]

The code keeps two distinct routes for the sweep: the sampled spatial derivative of \(m_{\rm MS}\), and the exact piecewise-linear secant implied by its interpolation bracket. Their difference is meaningful operator evidence, not permission to select whichever closes the flux ledger.

For the exchange algebra, reuse the pinned research note [conditional scalar–dark-matter exchange projection](./conditional_scalar_dm_exchange_projection.md) (blob `e793fc304908e118256adda832beba08e9876ea1`). It already records the conditional mass/action premise, \(Q_\nu\), the three matter source channels \(q_D,q_E,q_{S_r}\), and the reciprocal \(+\alpha\beta_{\rm DM}\rho_{\rm DM}\) scalar source. Its regression includes nonzero shift and non-unit lapse. This work should consume that derivation as a conditional source contract, not duplicate it or elevate it to proof of the action premise, physical coupling, or global evolving-state closure.

**Fail-closed provenance comparison for this fixed-\(R\) audit:** its own declared `SOURCE_PINS` were checked read-only against the exact crosswalk base `4aec88d0e3ff028e0a65a8ac36b30948625df792`. Eleven of the 17 pinned engine blobs match. Six differ:

| Pinned source | Audit-expected blob | Observed blob at crosswalk base |
|---|---|---|
| engine/matter_rhs.py | `f41180efbcff2c72cf36477a765318ab539fac80` | `0153c232ee15b5a9d3dbac37ea858e4d2d730011` |
| engine/matter_system.py | `05e23d2368120ea8a1a6aad3fa438625d5380750` | `0969f4d996e4f8e0fc5b70a7dff147f024907302` |
| engine/production_kernel.py | `9685e76bb5a6f26d6bf383e240de24f37523f265` | `1720615628a98a8d7cfc9fcadcc984fe28c88be1` |
| engine/scalar_system.py | `ffc119b1a6f2d495d972323901f76a85783b0a1b` | `e2c2dae9a016b4d6454cf74a9eae8b0167d5733d` |
| engine/v55_matter.py | `acc8f8ceae28d82511c9a5ad9d2ca670e09fb0f7` | `0ee3c51fac7dc46c63feba2747a750542bcee3d0` |
| engine/v55_pirk_adapter.py | `eeee0492d80889d70e52b0946563eee5054c0f11` | `be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3` |

The vendor gitlink still matches the audit's expected pin `d6052d605673ce9d82cdc99b0df79855fa2ea215`. Thus the known provenance break is localized to the six listed source files, not the vendor pin. This was a blob/pin comparison only: no simulation or numerical audit was run, and the six expected pins must not be overwritten to make the old audit accept a newer source state.

**What the crosswalk establishes now:** every major fixed-\(R\) geometric operand has an existing implementation or explicit reconstruction, and the fluid/scalar exchange components are separately recoverable from the current source path. **What it does not establish:** that the old audit's historical vendor/source pins, its stress normalization, the current captured geometry tangent, and the actual matter/scalar RHS all describe one identical stage and convention. That requires a separately named current-source witness after pin and operand reconciliation. Do not run or reinterpret the old numerical residual as that witness.

### 2C. Applying the map: one-step current-source witness

A separate research-only tool now applies the crosswalk without altering the historical audit:

- [\`tools/audit_current_source_fixed_r_ledger.py\`](../../tools/audit_current_source_fixed_r_ledger.py) has its own current-source blob manifest and vendor gitlink check. It records the actual checkout/runtime identity and runs exactly one \(N=24\), \(R_{\max}=24\), CFL \(=0.01\), \(t_{\rm final}=0.01\) step with the optional Maxwell carrier disabled.
- It reuses the existing fixed-\(R\) mass/flux/work and passive geometry-stage observation functions, while capturing the actual **start RHS stage**. It verifies exact array equality between the stage geometry and the saved geometry operands before computing the joint witness.
- The current observer adapts to the six returned matter components: face flux, geometric energy source, geometric momentum source, exchange energy, exchange momentum, and the separate rest-current \(q_D\). This matters because the historical observer's instrumentation destructures the earlier five-value interface. The historical audit and its old source pins remain unchanged.
- It reconstructs the live rest, energy, and momentum RHS channels from their named discrete terms, and checks the conditional scalar \(\Pi\) exchange using the same stage state against a no-exchange evaluation. It preserves the factorized and composite-interpolated fixed-\(R\) flux/work operators and reports both physical balance remainders **signed and unthresholded**.
- [The bounded witness contract test](../../tests/test_research_current_source_fixed_r_ledger.py) exercises that one-step path, checks the source pins and stage match, and asserts that neither scientific admission nor the acceptance oracle is promoted.

This is a real engineering application of the source map, but its matter RHS reconstructions are operator bookkeeping—not independent integrated conservation proofs. The scalar-source check verifies the conditional implementation—not the action premise. The fixed-\(R\) remainder is an observation, not a pass/fail target. This witness does not close Issue #210, the total-\(T\) gate, or Issue #4.

### 2D. First current-source artifact — observed, not admitted

The current-source witness was executed on PR #414's merge checkout \`e636deec2ba744ad40e38f599f8e1cf895e96d78\` with source/vendor pins verified. Exact-head CI for PR head \`b46912666d1847f453ae73cdf65c0c09f331f4f1\` passed the 952-test repository suite, the bounded witness contract, the one-step production smoke, and all three PR workflows. The signed JSON is retained as [artifact \`current-source-fixed-r-ledger-e636deec2ba744ad40e38f599f8e1cf895e96d78\`](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40038594192/artifacts/12537759261) (14-day retention at upload).

**Observed internal checks**
- Three species were captured (\`dark_matter\`, \`baryons\`, \`radiation\`): all nine rest/energy/momentum RHS reconstructions passed; each reconstructed profile had zero observed arithmetic residual in this run.
- The live scalar exchange difference matched the conditional \(+\alpha\beta_{\rm DM}\rho_{\rm DM}\) source within its operand-scaled roundoff bound: maximum absolute difference \`8.470329472543003e-22\`, bound \`6.481289684017647e-19\`.
- All seven saved stage-geometry arrays (\`a,b,X,alpha,beta,Aa,K\`) matched the RHS stage arrays exactly.
- This particular initial state has \`max |beta| = 0\`, radial Maxwell disabled, and one accepted step at \`t=0\` for the captured start RHS. It is **not** a shift-sensitive or active-EM witness.

**Signed fixed-areal-radius results at \`R=10\`**
- Factorized interpolation: flux \`−4.846324899808471e−19\`, radial-pressure work \`2.1020774473368874e−4\`, coordinate sweep \`−3.4869686256279903e−4\`, fixed-\`R\` mass rate \`−1.4094452669647614e−4\`, residual \`−2.455408867365369e−6\`.
- Composite interpolation: flux \`−4.827228557047549e−19\`, radial-pressure work \`2.0848819518718583e−4\`, composite fixed-\`R\` residual \`−7.358593208624318e−7\`.
- The two physical remainders differ by about \`1.72e−6\`. Neither is zeroed, thresholded, or selected as the answer. The report also retains the signed Hamiltonian residual at the tube (\`3.570668042581941e−7\`) and momentum residual (\`7.563650701049904e−19\`).

These results establish **stage-aligned operand traceability and local source assembly**, not total stress-energy conservation. In particular, the algebraic mass-rate decomposition closes at \`−2.51e−18\`, but that identity only checks the expansion of the chain rule; it does not make either physical flux/work remainder acceptable. The fixed-\`R\` continuum law and the stress normalization remain independent obligations.

### 2D.1 Latest exact-head current-source artifact

The later exact-head run is preserved at [workflow run 40039446398 / artifact 12539085382](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40039446398/artifacts/12539085382). PR head was `3b68e00d0617bec623bade5949ec0916b2a0a539`; the workflow checked out merge SHA `c6c5f448ea0c8ec1a932da023309fb19e12c9cc2`. All 25 production blob pins and the vendor gitlink matched. Runtime recorded Python 3.12.15 and NumPy 2.5.4; `OPENBLAS_CORETYPE` and `NPY_DISABLE_CPU_FEATURES` were unset. This pins the source witness, but does not claim cross-run bitwise identity across CPU-dispatch environments.

The one-step capture remains (N=24), (R_{\max}=24), CFL (=0.01), (t_{\rm final}=0.01), with the start RHS at (t=0), zero shift, and radial Maxwell disabled. All seven geometry arrays match their captured RHS-stage arrays exactly. All nine live matter RHS component reconstructions have zero observed arithmetic residual in this report. The conditional scalar exchange check has maximum absolute residual (8.470329472543003\times10^{-22}) against a roundoff bound (6.48128968401772\times10^{-19}).

At (R=10), the retained signed fixed-(R) residuals are:

- Factorized interpolation: (-2.455408867470726\times10^{-6}).
- Composite interpolation: (-7.358593209677621\times10^{-7}).
- Factorized versus composite radial-pressure work: (2.1020774473369088\times10^{-4}) versus (2.0848819518718792\times10^{-4}).
- Stress-flux values: (-4.846324899808519\times10^{-19}) versus (-4.827228557047598\times10^{-19}).

The residual difference is (1.7195495465\times10^{-6}), essentially the difference between the two recorded pressure-work operators; the flux difference is negligible on this scale. This identifies interpolation order as the immediate discrete sensitivity, **not** which operator is physically authorized. Neither result is thresholded or selected.

The same report records a separate mass-gradient witness at the active bracket: native discrete gradient (0.006134080955367126), local-composite RHS (0.006106681695669592), and their residual (2.7399259697534304\times10^{-5}). The Hamiltonian contribution is (9.013332209810083\times10^{-6}), the momentum contribution is (2.2018768110436523\times10^{-18}), and the remaining operator term is (1.8385927487722022\times10^{-5}). The chain-rule term expansion reconstructs to (-5.66\times10^{-18}); that arithmetic closure does not close the physical flux/work balance.

**Source mapping retained for the next derivation:** the live boundary observer consumes the canonical `engine/stress_energy.py` projection assembled by `engine/v55_matter.py`. It takes total `j` and `pr` from scalar plus fluid projections (EM is zero in this run); the archived bridge uses (T^r{}_t=-\alpha\gamma^{rr}j_r) and (T^r{}_r=p_r) under its stated zero-shift convention. The code defines scalar `j=-\Pi\,\partial_r f` per sector scale and fluid `j=hW^2\gamma_{rr}v_r/\kappa`; the radial stresses are assembled sector-wise. This is the implementation mapping, not yet independent authorization that every active S/D/\(\phi\) normalization and tensor convention matches the action-bearing archive source (Issue #210).

Limits remain unchanged: `scientific_admission=false`, `acceptance_oracle_satisfied=false`, no physical residual threshold, and no closure claim for total (T^{\mu\nu}) or Issue #4.

### 2E. Fixed-\(R\) continuum identity from the archived coordinate-worldtube balance

This is a **conditional continuum derivation** under the archived spherical CMC convention, \(G=c=1\), \(\beta^r_{\rm shift}=0\), coordinate time \(t\), a consistent total stress projection, and \(R_r\ne0\). It derives the moving-tube chain rule from the existing fixed-coordinate identity; it does not resolve the missing primary action/normalization provenance.

The in-repository Valencia covector balance is

\[
\partial_\mu(\sqrt{-g}\,T^\mu{}_\nu)
=\frac12\sqrt{-g}\,T^{ab}\partial_\nu g_{ab}
+\sqrt{-g}\,Q_\nu .
\]

For the total tensor, an ordinary zero-source conservation equation follows only if the sector exchange terms cancel under one authorized source law and compatible normalizations. The archived CMC bridge records, for its **fixed-coordinate-\(r\)** spherical worldtube,

\[
\left.\partial_t m_{\rm MS}\right|_r
=4\pi R^2\left(T^r{}_t\,R_r-T^r{}_r\,R_t\right)
=F_r+W_r ,
\]
where, using that archived projection/sign convention,
\[
F_r=-4\pi R^2\alpha\gamma^{rr}j_r R_r,\qquad
W_r=-4\pi R^2p_r R_t,\qquad p_r\equiv T^r{}_r .
\]

Here \(R_t=(\partial_t R)|_r\), \(R_r=\partial_rR\), and \(j_r\) is the signed radial momentum projection for which the archived relation is \(T^r{}_t=-\alpha\gamma^{rr}j_r\). \(p_r\) is the **radial stress projection** in the same tensor convention, not automatically the fluid rest-frame pressure. The index placement and signs are essential; if the source projection differs, the equation must be re-derived rather than repaired by choosing a sign that reduces the residual.

A fixed-areal-radius tube is the coordinate curve \(r=r_*(t)\) satisfying \(R(t,r_*(t))=R_*\). Differentiating that defining equation gives
\[
v_R\equiv\frac{dr_*}{dt}=-\frac{R_t}{R_r},
\quad (R_r\ne0).
\]
Applying the scalar chain rule to \(m_{\rm MS}\) along the curve, then substituting the archived fixed-coordinate balance, gives
\[
\left.\frac{dm_{\rm MS}}{dt}\right|_R
=\left.\partial_t m_{\rm MS}\right|_r
+v_R\,\partial_r m_{\rm MS}
=F_r+W_r+S_R ,
\]
with the signed moving-surface/sweep term
\[
S_R=v_R\,\partial_r m_{\rm MS}
=-\frac{R_t}{R_r}\,\partial_r m_{\rm MS}.
\]

This derivation fixes the **continuum term structure and signs**, conditional on the stated archived identity. It also explains why setting \(R_t=0\) inside \(W_r\) while leaving \(\partial_t m|_r\) on the left would mix two different worldtubes. \(F_r\), \(W_r\), and \(S_R\) must be evaluated at the same stage and surface; at \(R_r=0\), this coordinate parameterization fails and requires separate treatment.

The current witness already records all three terms, their operands, and both discrete interpolation orders. Its observed difference between factorized and composite physical remainder is almost entirely the radial-pressure-work interpolation difference, while the stress-flux values themselves are both about \(5\times10^{-19}\) for this zero-shift/EM-off start stage. That is evidence of discrete-operator sensitivity, **not** a criterion for selecting the composite operator. The remaining physical remainder is still signed and unthresholded. A next source-level review should now answer, term by term: (a) does the archive's continuum convention authorize the same \(j_r,p_r\) as the live total projection; (b) what spatial operator defines \(\partial_r m_{\rm MS}\) for the intended discrete fixed-\(R\) observable; and (c) which difference is truncation/operator remainder versus a genuinely missing continuum term? No source factor or operator is selected by fitting the residual.

### Structural pattern visible in the five changed blobs

The source changes cluster along two coupled paths rather than five unrelated locations:

- **Matter/scalar source path:** matter_system.py and matter_rhs.py now carry the variable-mass rest-current contribution plus lapse/shift-projected exchange terms; scalar_system.py carries the reciprocal lapse-weighted scalar source.
- **Geometry tangent path:** production_kernel.py and v55_pirk_adapter.py route the CMC \(a,b,X\) explicit metric RHS through the native geometry implementation while retaining the distinct primary/Lambda blocks.

That is a useful structural pattern for deciding what operands belong in the next ledger. It is **not** evidence that these changes cause the old residual or that either path is physically wrong. The exchange regression already covers its conditional local source algebra; the open task is to connect those terms, the normalized total projection and the geometry-derived mass tangent on the same accepted state and convention.

One caution: v55_matter.exchange_pair() returns opposite dm_* and phi_* entries and zero net_* arrays by construction. That helper is a representation of the assumed pairing, not an independently measured divergence and not a closure witness.

### Bounded next action from this map

The first current-source artifact is now captured and its operator checks pass. Do not repeat that run just to obtain a smaller residual. The next scientific task is a convention-pinned derivation of the fixed-areal-radius continuum mass/flux/work identity from the archived covariant balance, including the moving-tube orientation and metric/index signs. Then map each derived term to the already captured discrete operands and explicitly identify any continuum term that the present projection does not supply. Keep both interpolation operators as diagnostics; do not use their difference as a fit target.

No broad campaign rerun, historical-pin update, production-equation edit, runtime-default change, scientific-admission promotion, or Issue #4 status change is authorized by this step.


### 2F. Current-head operator discriminator and source/work identity — 24 October 2026

The exact-head [Lightweight validation run 40113174226](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40113174226) completed successfully and preserved [artifact 12567890444](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40113174226/artifacts/12567890444). Its report records PR source head \`5c7e57146cd13eb18b147f27938982a14f4f19e7\` and execution checkout \`7075f578af165dc1967e8ab17810face989314d1\`. All 25 production-file blob pins and the vendor gitlink match. The run is a single accepted \(N=24\) step; its captured RHS is at \(t=0\), \(R_*=10\), zero shift, and radial Maxwell disabled. \`scientific_admission=false\`, \`acceptance_oracle_satisfied=false\`, and no physical residual threshold is applied.

The nine live matter rest/energy/momentum reconstructions have zero observed arithmetic residual in this report. The conditional scalar exchange reconstruction has maximum absolute residual \(8.470329472543003\times10^{-22}\) against an operand-scaled bound \(6.48128968401772\times10^{-19}\). All seven stage arrays \((a,b,X,\alpha,\beta,A_a,K)\) match their captured RHS-stage arrays exactly. These are implementation and stage-alignment checks, not a proof of physical exchange or total-\(T\) conservation.

**Fixed-\(R\) mass-rate comparison, with the same piecewise-linear mass observable and its secant sweep used for both stress operators:**

| Operand at \(R_*=10\) | Factorized factors | Composite local product |
|---|---:|---:|
| Stress flux | \(-4.846324899808519\times10^{-19}\) | \(-4.827228557047598\times10^{-19}\) |
| Radial-pressure work | \(2.1020774473369088\times10^{-4}\) | \(2.0848819518718792\times10^{-4}\) |
| Piecewise-linear fixed-\(R\) sweep | \(-3.4668261280958277\times10^{-4}\) | \(-3.4668261280958277\times10^{-4}\) |
| Predicted fixed-\(R\) rate | \(-1.3647486807589237\times10^{-4}\) | \(-1.3819441762239534\times10^{-4}\) |
| Measured piecewise-linear fixed-\(R\) rate | \(-1.389302769433631\times10^{-4}\) | \(-1.389302769433631\times10^{-4}\) |
| Signed residual | \(-2.455408867470726\times10^{-6}\) | \(-7.358593209677621\times10^{-7}\) |

Because the measured mass tangent and secant sweep are held common in this table, the residual difference \(1.7195495465\times10^{-6}\) is the radial-pressure-work difference to the displayed precision; the stress-flux difference is negligible at this scale. This is a clean comparison of the two recorded stress/work interpolation operators, but it still does **not** authorize choosing one from the smaller physical remainder. The fourth-order sampled-\(D_r M\) sweep is also retained separately; it must not be mixed with the piecewise-linear tangent or secant sweep.

**The mass-gradient identity is even more operator-sensitive.** At the same target, the sampled discrete gradient is \(0.006134080955367126\). If the primitive factors are interpolated and then multiplied, the RHS is \(0.006140802495695930\), giving residual \(-6.7215403288046\times10^{-6}\). If the cell-local product \(4\pi R_i^2(\rho_i D_rR_i+j_i\dot R_i/\alpha_i)\) is formed first and that composite profile is then interpolated, the RHS is \(0.006106681695669592\), giving residual \(+2.7399259697534304\times10^{-5}\). The two RHS operators differ by \(3.412080002633804\times10^{-5}\), and the residual changes sign. The composite residual decomposes into a Hamiltonian contribution \(9.013332209810083\times10^{-6}\), a negligible momentum contribution \(2.2018768110436523\times10^{-18}\), and an unassigned operator remainder \(1.8385927487722022\times10^{-5}\).

The bracket cells independently show the same unresolved split: at \(R=9.5\), the mass-gradient residual is \(2.9585674663249253\times10^{-5}\), with Hamiltonian contribution \(7.475602277123505\times10^{-6}\) and operator remainder \(2.211007238612415\times10^{-5}\); at \(R=10.5\), the residual is \(2.5212844731819367\times10^{-5}\), with Hamiltonian contribution \(1.0551062142496657\times10^{-5}\) and operator remainder \(1.4661782589319905\times10^{-5}\). Do not confuse this cell-local composite residual with the distinct factor-interpolated target residual above. Neither remainder by itself establishes a missing continuum source.

**A local source/work identity does close algebraically.** The existing diagnostic reconstructs the matter-source contribution
\[
S=R^3K^\theta{}_\theta\,\dot K^\theta{}_\theta\big|_{\rm matter},\quad
W=-4\pi R^2p_r\dot R,\quad
D=\frac{4\pi}{3}R^2\rho\dot R
\]
and tests
\[
S-W+D=4\pi R^2(\dot R+\alpha R K^\theta{}_\theta)(\rho/3+p_r).
\]
At the target, the composite-profile kinematic closure is \(-3.469446951953608\times10^{-18}\), and the composite \(S-W+D\) residual is \(-4.81482486096809\times10^{-35}\) (its algebraic reconstruction residual is about \(5.86\times10^{-21}\)). Under factorized interpolation, the kinematic closure is \(-8.553779305628395\times10^{-7}\) and the identity's nonzero value is \(1.4121634228192285\times10^{-9}\), consistent with the independently computed kinematic right-hand side. This identifies a noncommuting interpolation effect in the source/work identity; it does not prove global conservation or pick the physical stress-flux operator.

**Next bounded action:** specify and test the discrete operator for the mass-gradient identity from the archived continuum derivation and the chosen fixed-\(R\) observable, using operator-matched manufactured profiles and refinement evidence. Keep both product orders visible until the discretization is justified independently of residual size. Then map the resulting signed identity to the already captured total-stress projection and geometry tangent. No production equation, source normalization, historical pin, or admission flag changes follow from this report.


### 2G. One-step current-source fixed-\(R\) refinement — 24 October 2026

The bounded three-resolution current-source report completed on exact PR source head \`e7339e89aa7241ff549f3bcdd2d9e01dd7b873f4\`. Workflow [40114701026](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40114701026) preserved the three-report artifact [12568460374](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40114701026/artifacts/12568460374). The synthetic execution checkout recorded in all reports is \`5781daf0c406b80e4c6a7c59578b9cae4c979821\`; all 25 current production-file blob pins and the vendor gitlink match in every report.

Each resolution performs exactly one accepted step and measures the start RHS at \(t=0\), with \(R_{\max}=24\), CFL \(=0.01\), \(R_*=10\), zero shift, and radial Maxwell disabled. The final time is \(0.01,0.005,0.0025\) for \(N=24,48,96\), respectively, so each step uses the corresponding \(dt=\mathrm{CFL}\,\Delta r\). All three runs report nine of nine matter RHS reconstructions within their operand-scaled roundoff checks, the conditional scalar-exchange reconstruction within its roundoff bound, and exact stage identity for \((a,b,X,\alpha,\beta,A_a,K)\). All reports retain \`scientific_admission=false\`, \`acceptance_oracle_satisfied=false\`, and \`physical_fixed_R_residual_threshold_applied=false\).

**Matched discrete mass-gradient identity at \(R_*=10\).** Define the measured composite residual as the sampled fourth-order mass gradient interpolated to the target minus the cell-local composite continuum RHS interpolated to that same target. Results are signed:

| \(N\) | \(\Delta r\) | Composite mass-gradient residual | Hamiltonian contribution | Momentum contribution | Remaining operator term |
|---:|---:|---:|---:|---:|---:|
| 24 | 1.00 | \(+2.73992597\times10^{-5}\) | \(+9.01333221\times10^{-6}\) | \(+2.20\times10^{-18}\) | \(+1.83859275\times10^{-5}\) |
| 48 | 0.50 | \(+6.46311726\times10^{-6}\) | \(+5.26292778\times10^{-6}\) | \(+6.59\times10^{-19}\) | \(+1.20018948\times10^{-6}\) |
| 96 | 0.25 | \(+1.59180115\times10^{-6}\) | \(+1.51607041\times10^{-6}\) | \(+4.25\times10^{-18}\) | \(+7.57307377\times10^{-8}\) |

The total composite mass-gradient residual decreases by factors \(4.239\) and \(4.060\) when the grid spacing halves, an approximately second-order trend over these three grids. The remaining operator term decreases by factors \(15.319\) and \(15.848\), consistent with the fourth-order spatial derivative used by the pinned production operator. Its share of the full residual declines from about \(67.1\%\) at \(N=24\), to \(18.6\%\) at \(N=48\), to \(4.8\%\) at \(N=96\). At the finest grid the recorded Hamiltonian contribution accounts for about \(95.2\%\) of the composite residual, while the momentum contribution is negligible at this stage.

This is the key discriminator so far: the unassigned discrete operator remainder displays the expected fourth-order trend, while the remaining full identity residual is increasingly dominated by the separately recorded Hamiltonian-constraint contribution. It does **not** establish that the physical Hamiltonian constraint is acceptable, nor does it by itself establish global stress-energy conservation. No term is fitted or dropped.

**Factorized versus composite target evaluation.** The alternative factor-interpolated mass-gradient residuals are \(-6.72154033\times10^{-6}\), \(-1.96033140\times10^{-6}\), and \(-5.07411507\times10^{-7}\) for \(N=24,48,96\). The factorized-minus-composite RHS gap is \(3.41208000\times10^{-5}\), \(8.42344866\times10^{-6}\), and \(2.09921266\times10^{-6}\): this gap decreases by factors \(4.051\) and \(4.013\), consistent with second-order interpolation. The sign difference at a given finite resolution therefore reflects noncommuting discrete evaluation order, not by itself evidence for or against a missing continuum term.

**The fixed-\(R\) mass-rate balance remains open.** The factorized signed residuals are \(-2.45540887\times10^{-6}\), \(-3.82257016\times10^{-7}\), and \(-8.03733469\times10^{-8}\). The composite signed residuals are \(-7.35859321\times10^{-7}\), \(+4.10208893\times10^{-8}\), and \(+2.50290952\times10^{-8}\). The factorized magnitudes decrease over the three grids; the composite residual changes sign between \(N=24\) and \(N=48\) and does not yet show a clean refinement ratio between the last two grids. Neither operator is authorized by residual size, and neither result is treated as a conservation pass.

The source/work identity evaluated from the same composite profiles continues to reconstruct algebraically near floating-point roundoff; the factorized kinematic mismatch decreases at approximately second order. These are internal identity/operator checks, not an integrated evolving-state total-\(T^{\mu\nu}\) closure.

**Disposition / next bounded step.** Keep this three-grid trend as the current-source operator evidence. Do not launch another fixed-R grid sweep just to improve a residual. The next substantive step is to reconcile the Hamiltonian contribution and its constraint convention against the archive-backed mass-gradient identity and the exact stress projection, while keeping the fixed-\(R\) mass-rate remainder signed. Issue #210 normalization provenance, full evolving total-stress-energy closure, and Issue #4 remain open; production equations and scientific admission remain unchanged.

## 3. Spherical mass identity: exact scope matters

The archived CMC bridge defines the areal radius `R(t,r)` and Misner–Sharp mass

```text
m(t,r) = R/2 [1 - g^{ab} (∂_a R)(∂_b R)]     (G = c = 1 convention)
```

For its **spherically symmetric, beta=0, fixed-coordinate-r worldtube**, the recorded balance is

```text
(∂_t m)|_r = 4π R² [T_t^r (∂_r R) - T_r^r (∂_t R)]
           = 4π R² [-α γ^rr j_r (∂_r R) - p_r (∂_t R)].
```

Use the second line only with the same stress-projection, index, metric, sign, and time conventions as the source derivation. The `-p_r ∂_t R` term is part of the balance; dropping it changes the observable. This formula is not a universal flat-space flux `4πr²J_r`, and `r` must not be silently identified with areal radius `R`.

### Fixed coordinate versus fixed areal radius

The repository's admission lane explicitly requires a **fixed-areal-radius** measurement surface. A fixed coordinate `r` generally has `∂ₜR ≠ 0`, so it is a different worldtube. To compare mass on a constant-`R` surface, first derive the moving coordinate location:

```text
dr/dt |_R = -(∂_t R)/(∂_r R),     provided `∂ᵣR ≠ 0`
(dm/dt)|_R = (∂_t m)|_r + (dr/dt)|_R (∂_r m).
```

This chain rule specifies the surface derivative; it does **not** by itself finish the physical flux law. The flux through the moving worldtube, its orientation, and all work/geometry terms must be transformed consistently from the same continuum identity. Do not set `∂t R=0` in the fixed-coordinate formula while retaining the fixed-coordinate derivative on the left. At points where `∂ᵣR = 0`, this coordinate conversion is not valid and needs a separate treatment.

### Conditional fixed-\(R\) identity available from the archive

The archived CMC bridge explicitly assumes \(\beta^r=0\) and states the fixed-coordinate balance
\[
\left.\partial_t M\right|_r=-4\pi R^2\left(\alpha\gamma^{rr}j_r\,\partial_rR+p_r\,\partial_tR\right).
\]
Combining that exact archived identity with the surface chain rule gives, only under those same conventions,
\[
\left.\frac{dM}{dt}\right|_R=
-4\pi R^2\left(\alpha\gamma^{rr}j_r\,\partial_rR+p_r\,\partial_tR\right)
-\frac{\partial_tR}{\partial_rR}\,\partial_rM.
\]
The last term is the moving-coordinate sweep, not a new matter source. It must use the mass gradient from the same state and a consistent interpolation at the target \(R\). The implementation already reports the flux, pressure-work, sweep, left-hand-side mass tangent, and their signed remainder under both factorized and composite interpolation.

**Convention gate:** this expression is not licensed for nonzero shift. The archived source fixes \(\beta^r=0\); for a nonzero-shift state, derive the balance directly from the correctly indexed \(T^r{}_t\) and \(T^r{}_r\), including shift terms, before evaluating a physical residual. The current-source witness therefore now fails closed if the captured stage has any nonzero shift. Its zero-shift artifact is a source-aligned diagnostic, not a general-shift conservation proof.

## 4. Derivation and implementation sequence

1. **Pin the source convention.** Record metric signature, units, index placement, `T^μ{}_ν` versus `T_μν`, definitions/sign of `j_r`, lapse/shift, and whether the derivative is coordinate-time or proper-time.
2. **Recover the archive equation, not just its implementation.** For each matter/scalar sector, locate the exact source equation and its provenance. Keep the unresolved S/D/φ action and density-normalization issue visible; a normalized archive model may be the most specific justified claim.
3. **Derive the continuum identities independently.** (a) total covariant stress-energy divergence and any conditional exchange law; (b) the matter conservative-variable equation; (c) Misner–Sharp mass/worldtube balance. Show which source terms cancel only after summing sectors and which represent external/non-conservative sources.
4. **Choose the physical measurement surface.** Specify constant areal radius, its coordinate trajectory, orientation, and time parameter before calculating a flux.
5. **Map continuum terms to discrete operators.** Reuse the existing Valencia flux/source budget and total-`T` observer. Identify face areas, cell measures, metric factors, and interpolation/quadrature. Avoid a second implementation of an already tested term unless it is explicitly an independent cross-check.
6. **Construct a signed ledger.** Report, separately, stored-quantity change, inner and outer face contributions, radial-pressure/geometric work, explicit matter-sector exchange, and the local and integrated closure residual. State sign conventions in the artifact itself.
7. **Use bounded controls before a production seam.** Manufactured solution with known source, zero-current control, orientation reversal, constant-areal-radius surface check, and refinement of local and integrated residuals. Vary the outer domain while keeping the physical surface fixed. No gain fitting, clipping, hidden source, or profile correction.

## 5. Admission gates / next work

- [ ] Recover exact archive file/revision/page or equation identifier for the action-bearing S/D/φ stress-energy convention, or explicitly narrow the claim to the normalized reference model (Issue #210).
- [ ] Reconcile the archived Misner–Sharp bridge source and artifacts against current pinned production blobs before treating its numerical residuals as current-main evidence.
- [ ] Derive and review the fixed-areal-radius worldtube expression from the same continuum identity; include moving-surface terms consistently.
- [ ] Demonstrate term-by-term correspondence between the continuum identity and existing stress projection / Valencia flux / metric-source / exchange code.
- [ ] Build a research-only ledger with signed local and integrated residuals, manufactured controls, orientation reversal, and refinement.
- [ ] Keep Issue #4, radiation near-cone behavior, and global total-stress-energy closure open until their own admission witnesses pass.

**Current conclusion:** the repository and archive already contain useful pieces of a real ledger, especially the source-pinned total-`T` observer, the signed Valencia stage budget, and the archived Misner–Sharp bridge. The correct next action is to reconcile their conventions and worldtube definitions—not to add a new feedback law or tune the recirculation probe.
