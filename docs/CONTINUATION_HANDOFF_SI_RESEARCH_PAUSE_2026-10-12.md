# 0-Cosmos Engine — SI Research Continuation Handoff

**Handoff date:** 12 October 2026 (UTC)  
**Purpose:** Pause cleanly; preserve current derivations, evidence, limitations, and the exact next restart gate.  
**Repository:** [miahpenn/0-Cosmos-engine](https://github.com/miahpenn/0-Cosmos-engine)  
**Production branch:** `main` at `2035acc1210db18da0cededea18750b7ed813260` at last verification.  
**Production disposition:** main unchanged; no diagnostic candidate merged; no full physical trajectory is admitted by this work.

## 1. Read this first — current conclusion

The strongest current numerical finding is a discrete CMC-gauge boundary-row defect on the cell-centred radial mesh. The grid stores values at \(r_i=(i+1/2)\Delta r\), so `alpha[0] = alpha[1]` equates values at two different positive radii; it is not the even-parity condition at the origin and it omits the CMC PDE at the first cell. Replacing that equality row with the first-cell CMC PDE causes a large improvement in the one-step momentum and connection constraint witnesses at N=60–160. A three-arm test further shows that most of that improvement comes from adding the missing PDE row, rather than from fourth-order versus second-order stencil order alone.

This is a **promising, well-supported numerical-layer correction candidate, not a full machine repair**. Important blockers remain:

- The centre-series witness `Aa_dot / r_0^2` is still finite and nonzero with refinement; regular-centre evolution is not admitted.
- At N=240 and N=320, both old and corrected CMC arms fail the radiation RHS admissibility check at cell 0 **before the first predictor CMC solve**. The row correction therefore does not resolve the high-resolution first-cell radiation failure.
- The new one-step Misner–Sharp trapezoidal residual decreases with resolution but is nonzero and has no preregistered acceptance tolerance. It is not a conservation pass, and its value is nearly unchanged by the CMC row switch.
- A separate prior late-time trace remains: near t=46.98, an outer cell crosses the radiation cone when the explicit predictor metric is applied to an almost saturated accepted state. That late-time failure is distinct from the high-resolution first-cell RHS failure.
- The lapse-weighted covariant scalar–DM source pair passes source-stage checks and the homogeneous-limit RHS regression, but a tested one-step source A/B changes global geometry/constraint observables only at tiny levels. This supports source accounting on the tested stages; it does not identify that source as the dominant global defect.

**Restart posture:** review the artifacts and preserve the branches. Do not start a long campaign, change a physical law, or merge the CMC row from the current evidence alone.

## 2. Primary durable report and repo state

The expanded technical derivation is recorded here:

- [Production-readiness derivation review — 12 Oct 2026](https://github.com/miahpenn/0-Cosmos-engine/blob/research/production-readiness-review-2026-10-12/docs/PRODUCTION_READINESS_DERIVATION_REVIEW_2026-10-12.md)
- Review branch: `research/production-readiness-review-2026-10-12`
- Latest report commit: `ceaeb06bae506102d33f4942b1bf05bbbae398d1`
- Report content blob SHA: `b8d5607accc6807421c3a0ae3d87dba6b5c076a6`

Main was last confirmed at `2035acc1210db18da0cededea18750b7ed813260`. New source and gauge edits are isolated branches. Do not mix branch-specific test counts and describe them as main CI.

The canonical GEAR archive and REV16 criterion remain authoritative. The later V5.5 implementation target is one spatial spacetime and one total stress-energy tensor. `CosmosState` is a homogeneous reference/control, not a second state to inject into the live spatial production kernel. The new homogeneous-limit test is an unadvanced matched RHS/constraint check, not an evolved homogeneous trajectory. The distinct open GEAR bridges (D-mode to ordinary matter and local two-lane source to a global cosmological Q) are not closed by the tests below. REV16 remains `0*_eff = FALSE` under the currently recorded scorecard; numerical progress does not change that criterion.

## 3. Key CMC centre-row evidence and provenance

### Frozen-state row audit

- Run: [38302913119](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38302913119)
- Branch: `diag/cmc-lapse-cellcenter-row-repair-2026-10-12`
- Commit: `9f7ebf1aab2076fe8925737591366b416427305e`
- Artifact ID: `11771098012`
- Artifact ZIP SHA-256: `2c64f2f5b8db074f40748e50f22cfa14da71c98398d52bb626282932a79f9d5c`
- N = 40/60/80/120/160/240/320; frozen profiles only; no time advance.

The parity-aware first-cell PDE row, interior fourth-order rows, outer interior closure, and outer normalization yield CMC and production K-dot target residuals near roundoff: maximum about `2.00e-14` at N=320. The first-eight-cell even-series fit residual falls from `2.554e-4` at N=40 to `2.046e-10` at N=320. Lapse remained finite and positive. The nonzero difference `alpha[0]-alpha[1]` is expected for an even profile sampled at different positive radii.

Limitation: even after correcting this lapse row, the frozen `Aa_dot/r_0^2` witness remains about `-1.43e-5` at N=80, `-5.59e-5` at N=160, and `-6.72e-5` at N=320. The row consistency result is not proof of regular-centre evolution or continuum correctness.

### Matched one-step A/B

- Run: [38303550121](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38303550121)
- Branch: `diag/cmc-lapse-cellcenter-one-step-ab-2026-10-12`
- Commit: `0de3fa09511c919941072390a24426d10681b1fd`
- Artifact ID: `11770494466`
- Artifact ZIP SHA-256: `9c4934ddafed682044aced6222d94caf7c8fcff2434ec00a72a423c4abccf19d`
- N = 40/60/80/120/160; r_max=80; dt = 0.03 dr; one normal production step per arm.

Initial non-lapse geometry, scalars and conservative matter were identical; only the CMC solver row/operator differed. The parity-aware candidate strongly improves M and C_Lambda at N=60–160, often by roughly 20–30 times in the fixed annulus, while H is mixed. Full-domain H at N=120/160 worsens; M and C_Lambda worsen at N=40. `Aa_dot/r_0^2` drops greatly in magnitude, but does not approach zero with refinement.

### Three-arm row-order separation

- Run: [38304709682](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304709682)
- Branch: `diag/cmc-lapse-three-arm-row-order-2026-10-12`
- Commit: `1fca3f2f03626346720c95f4d4c114b884d1be96`
- Artifact ID: `11771896263`
- Artifact ZIP SHA-256: `59918d53bc883fbc796bf79057666692a01c49a0a355bd26f89cbac22fbb08d8`
- N = 80/120/160; one-step A/B/C.

The arms were: A, historical equality row; B, second-order PDE row at cell 0; C, fourth-order parity-aware PDE row matched to the production operator. At N=80, A→B reduces admitted-inner momentum max by `3.70e-6`, while B↔C differs by only `2.64e-10` in annulus M, `1.27e-11` in annulus H and `8.30e-12` in annulus C_Lambda. The pattern holds at N=120 and N=160. This materially strengthens the conclusion that the principal one-step M/C improvement is due to replacing the invalid equality row with the missing first-cell PDE, rather than merely increasing stencil order.

B and C do have different centre `Aa_dot/r_0^2` values (difference about `3.96e-7` at N=80, `1.93e-7` at N=120, and `1.12e-7` at N=160); neither is centre-series admitted. Full H is also worse for B/C than A at N=120/160. Do not claim universal improvement.

### Regression and promotion status

The corrected row is tracked in open draft [PR #7](https://github.com/miahpenn/0-Cosmos-engine/pull/7), branch `fix/cmc-cellcenter-row-2026-10-12`, head `72483464c6008b4529985c8767b501b11f57fcc9`. Two branch runs passed 2 focused tests and the full 81-test suite: [38304045777](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304045777) and [38304170442](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304170442). The isolated CMC full-suite workflow [38303841261](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38303841261) passed 89 tests; it added a workflow only relative to the candidate source. Regression success is not physical admission. PR #7 remains draft and should not be merged from these results alone.

## 4. Combined radiation-cone and Misner–Sharp A/B — latest diagnostic

- Run: [38304860664](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304860664)
- Branch: `diag/cmc-lapse-rad-ms-one-step-ab-2026-10-12`
- Audit commit: `759c5ead2d1b591667ddbe8ae82dbdc622aee15d`
- Candidate source under test: `9f7ebf1aab2076fe8925737591366b416427305e`
- Vendor: `d6052d605673ce9d82cdc99b0df79855fa2ea215`
- Artifact ID: `11772240963`
- Artifact ZIP SHA-256: `b63218c17b91a85d329441d9252124f96b6b7095602009cd1bf5e0d24e29e045`
- `report.json` SHA-256: `0a938cfca543e838cf0415ecd37c450375c34f4c103faab9ba28f212b862d90f`
- [Local copy of the artifact](sandbox:/mnt/data/cmc-lapse-radiation-ms-one-step-ab.zip)

The fail-closed job verified no engine/vendor/dependency changes relative to the pinned candidate. Helper tests passed 2/2 in 0.27 s and the complete repository suite passed 91 tests in 5.41 s. The workflow ended successfully as a **capture**: both arms completed one step for N=40/60/80/120/160; N=240/320 deliberately report a captured radiation RHS failure. Do not call the run a seven-resolution physical-step pass.

### Radiation outcome

The diagnostic reconstructs the cone invariant from the actual conservative radiation state and metric: `C = U_E - sqrt(X^2/a)*abs(U_r)` and `|S|/E = abs(U_r)*sqrt(X^2/a)/U_E`.

At N=40–160, both arms complete one step with positive cone margin. This small one-step state has predictor max ratios near 0.0101 and is not representative of the later near-saturated outer-cell trajectory. At N=240/320, both arms fail inside the radiation RHS at cell 0 **before the first explicit predictor CMC solve**. The CMC row change therefore does not fix or cause that high-resolution event.

| N | Old-row RHS failure ratio | Candidate-row RHS failure ratio | Difference |
|---:|---:|---:|---:|
| 240 | 1.0202778 | 1.0201269 | -0.0001509 |
| 320 | 1.3602899 | 1.3601766 | -0.0001134 |

At N=240, the failing updated conservative state is approximately E=`8.81887e-9`, |S|=`8.99770e-9` in the old arm versus E=`8.81887e-9`, |S|=`8.99636e-9` in the candidate arm. At N=320 the ratios are 1.36029 and 1.36018. Both remain inadmissible. The failure record included individual RHS contributions; for N=240 the flux-energy tendency is zero, source-energy about `-1.4166e-12`, flux-momentum tendency `-1.2495e-9`, and source-momentum `+9.9955e-10` in the old arm. This is the right target for the next diagnostic: identify the actual term/stage producing the first-cell cone loss without clipping or changing the cone rule.

Keep this failure separate from the earlier late-time outer-cell failure. The earlier trace ([report](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/radiation-predictor-stage-ledger/docs/RADIATION_PREDICTOR_TRACE_2026-10-10.md)) captured t≈46.98, central proper time τ≈22.0161, outer cell i=79, ratio 1.000057. Flux transport and geometric sources improved the margin on the accepted metric, but the explicit predictor metric changed gamma^rr enough to cross the cone. That trace did not prove the metric law wrong; it localized the crossing stage.

### One-step Misner–Sharp balance outcome

The combined audit records a descriptive one-step trapezoidal residual:

`Delta M_MS - dt * 0.5 * [(flux_T + work_pR)_initial + (flux_T + work_pR)_final]`.

| N | dt | Old scaled integrated residual | Candidate scaled integrated residual | Candidate−old absolute residual difference |
|---:|---:|---:|---:|---:|
| 40 | 0.0600 | -5.9588e-2 | -5.9002e-2 | 8.11e-9 |
| 60 | 0.0400 | -1.9503e-2 | -1.9487e-2 | 1.31e-10 |
| 80 | 0.0300 | -6.7155e-3 | -6.7132e-3 | 1.34e-11 |
| 120 | 0.0200 | -9.8935e-4 | -9.8985e-4 | 2.21e-12 |
| 160 | 0.0150 | -3.8669e-4 | -3.8679e-4 | 3.20e-13 |

Absolute integrated residual falls from roughly `-8.48e-7` at N=40 to `-1.22e-9` at N=160. The normalized mismatch also falls across this sequence, but no tolerance was preregistered, the endpoint quadrature has only one step, and the residual is nearly the same for both CMC arms. Correct interpretation: this shows a resolution trend worth following and indicates the CMC row switch is not the main cause of the current balance residual. It does **not** pass conservation.

## 5. Covariant scalar–DM source and homogeneous-limit gate

The newest documented current-main candidate snapshot is [run 38303904514](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38303904514), commit `92bfab1e330c9754d9909f0b2cb066baef7e8944`. It passed 2 actual-stage source-pair tests and 82 full-suite tests. The new homogeneous limit test compares the spatial scalar, DM rest/energy, baryon/radiation evolution and BSSN trace equation against the archived homogeneous COSMOS RHS on a matched flat-FLRW slice, and checks unadvanced H/M constraints. It excludes Bianchi-I shear in agreement with the exact-spherical V5.5 scope. This is an unadvanced RHS/constraint regression, not an evolved cosmology.

A later one-step source A/B, [run 38304517568](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304517568) at candidate snapshot commit `cd835df3d87f30f92d29afcd9a4b50a54b2877cf`, passed 4 focused tests and 83 full-suite tests. Actual stage pair residuals were around 1e-19 for rest/energy and 1e-24 for radial current. In its N=40, rmax=16, dt=0.012 comparison, candidate-minus-baseline changes were tiny: H max about `-8.49e-11`, admitted-inner M about `-6.0e-16`, connection max about `-3.7e-17`, M_MS at R≈10 about `-2.27e-13`. Therefore source reciprocity closes at the tested stages, but its global one-step impact in that setup is not dominant. Do not infer full trajectory conservation or close D-to-matter/local-to-global Q from this gate.

An earlier comparison run [38304363097](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304363097) failed during collection because `tests/test_covariant_beta_source_one_step_ab.py` contained an unclosed parenthesis. Later diagnostic source A/B and suite runs succeeded; that earlier failure was a harness syntax error, not evidence against the physical source law.

## 6. H-difference attribution

The one-step CMC attribution run [38304326260](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38304326260), branch `diag/cmc-lapse-cellcenter-one-step-ab-2026-10-12`, commit `920114ca17d029f6389ab204f197cec0754027fa`, used a Shapley-style component decomposition to reconcile the H difference between the old and candidate rows. Closure errors were around `1e-21` or smaller, so the decomposition reproduces the endpoint H difference. At the full-domain peak, the matter contribution was much smaller than the geometry contributions (for example, at N=160, about `1e-11` versus X, K, Lambda, and a*b contributions at `1e-7` to `1e-6` scales). This supports attributing the H response primarily to the changed metric/curvature solution rather than to matter source updates; it does not make the H response uniformly better. Artifact ID `11771589087`, archive digest `sha256:31910170208d4eda43f95660393ebb74b45f262acbc3195bd1a1ef1814b37c07`.

## 7. Failures, controls, and what not to infer

- Initial quadrature corrections substantially improve the initial Hamiltonian residual but do not repair `Aa_dot = O(r^2)`.
- Discrete-consistent H/M/C_Lambda initial data still produces an M tangent under the unchanged RHS; the one-step quotient agrees closely with that tangent.
- Smooth manufactured Ricci/Aa controls converge regularly, disfavoring a generic vendor sign/parity error, but they do not validate the actual production initial profile.
- The CMC row candidate is a real equation/discretization issue, but it is not proven to be the only centre RHS defect.
- High-resolution first-cell radiation failure and long-run outer-cell predictor-metric failure are distinct failure families.
- A successful workflow can mean successful capture of an expected failure. Workflow success, test-suite success, completion at target time, and physics admission remain distinct statuses.
- The homogeneous COSMOS state remains a reference/control. Do not add a second live homogeneous engine to the single-spacetime step.
- No bounce, turnaround, completed cycle, or global closure is admitted by this batch of evidence.

## 8. Exact next work when resuming

1. **Start with the combined artifact**, not a new simulation: download and inspect `report.json` from artifact 11772240963. Verify its SHA-256 (`0a938cfca543e838cf0415ecd37c450375c34f4c103faab9ba28f212b862d90f`) before using its numbers. The local artifact copy is linked above.
2. **Continue the first-cell radiation derivation** on the pinned production source. Instrument pre-RHS U_E/U_r, accepted metric, the separate flux and geometric source increments, and the exact admissibility check ordering at cell 0 for N=240/320. Keep the CMC row identical across this diagnostic unless a separate source attribution requires a fully matched A/B. Do not clip, floor, damp, or relax E≥|S|.
3. **Keep the CMC row candidate isolated**. The three-arm run shows the principal M/C gain tracks replacing the invalid equality row with a PDE row; the fourth-order matched row is a sensible candidate, but H and centre-series behavior are not fully admitted. PR #7 remains draft; do not merge before the radiation and balance witnesses are joined.
4. **Close the Misner–Sharp budget methodologically**: carry a predeclared residual, scale, surface definition, and multi-step/time-refinement comparison, first on a bounded admitted state. Current one-step values trend downward but do not pass any threshold.
5. **Only if the radiation and local constraint witnesses improve together**, run a further bounded old-row versus candidate-row comparison with H/M/C_Lambda, centre Aa/K series, determinant/Lambda identities, full predictor radiation cone, and MS budget recorded on identical stage states. No long campaign until this gate is clean.
6. The source-pair and homogeneous-limit regressions can be retained as completed code-level gates. Revisit their physics relevance only using a saved production-relevant state with nonzero phi_r and matched stage snapshots; do not claim they explain the remaining primary defects.
7. Before any production action, freeze exact commit and vendor hashes, verify exact-commit CI, preserve artifacts and JSON checksums, and separately state code regression, physics-admission and campaign status.

## 9. Guardrails for the next agent

No fitted feedback coefficients, artificial damping, floors/clamps, centre projections, imposed bounce/reset, invented D-to-matter law, phenomenological global Q, or silent physics. Change one diagnostic layer at a time. Preserve the GEAR archive and REV16. When describing convergence, use fixed physical radii or masks and do not treat moving first-cell/first-admitted-cell values as fixed-radius sequences. Report negative results and workflow/harness failures; distinguish syntax/test-harness errors from physical failure captures.

**Pause-state:** No more physics/campaign work was started after the combined one-step radiation/MS capture. The research report and this handoff are the durable restart anchors. Main remains at `2035acc1210db18da0cededea18750b7ed813260`; CMC PR #7 is still draft/unmerged.