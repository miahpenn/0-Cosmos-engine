# Build Status — 2026-10-13



## Current production snapshot — 2026-10-13

**Production ref at this snapshot:** `main` commit `efd006e33e95a19132af5066df21321f473580fc` (documentation-only status sync). The exact production runtime-source snapshot validated below is `7776a9156c24303c47a57e171a2d5997f6676668`.

**Disposition:** software-contract corrections have been promoted and post-merge smoke validation is green. The engine is **not yet admitted by the full scientific campaign gate**.

### Software promotions present on `main`

- [PR #12](https://github.com/miahpenn/0-Cosmos-engine/pull/12) — lightweight unit/control/one-step validation workflow.
- [PR #13](https://github.com/miahpenn/0-Cosmos-engine/pull/13) — campaign configuration, floating-point time-progress, and checkpoint-name preflight.
- [PR #14](https://github.com/miahpenn/0-Cosmos-engine/pull/14) — shared SI-agent communications protocol.
- [PR #15](https://github.com/miahpenn/0-Cosmos-engine/pull/15) — scalar lapse-gradient contraction uses the declared physical inverse radial metric `gamma^rr = X^2/a`.
- [PR #10](https://github.com/miahpenn/0-Cosmos-engine/pull/10) — removes the inadmissibility guard from only the fictitious `dt=1` RHS-extraction intermediate, while retaining default validation for direct updates and validating actual production stages.
- [PR #11](https://github.com/miahpenn/0-Cosmos-engine/pull/11) — preserves exact-zero radiation rest-current through radiation-specific HLL/outer-face flux dispatch.
- [PR #20](https://github.com/miahpenn/0-Cosmos-engine/pull/20) — refreshes production status and preserves the scientific-admission boundaries.
- [PR #23](https://github.com/miahpenn/0-Cosmos-engine/pull/23) — integrates strict campaign admission, typed radiation retries, strict radiation primitive-recovery checks, and standards-compliant result/ledger JSON on current main.
- [PR #22](https://github.com/miahpenn/0-Cosmos-engine/pull/22) — extracts the semidiscrete matter RHS directly while retaining componentwise finite-step ordering, stage-aware radiation admissibility, and radiation's exact-zero rest-current channel.
- [PR #27](https://github.com/miahpenn/0-Cosmos-engine/pull/27) — makes fixed-radius diagnostics robust when a valid grid does not sample either side of r=5, with explicit coverage flags and an `r_max=4` regression.
- [PR #30](https://github.com/miahpenn/0-Cosmos-engine/pull/30) — recursively normalizes strict JSON and persists fail-closed campaign reports when the returned initial state or its diagnostics are invalid, without unsafe diagnostic re-entry.

These changes close discrete software and equation-consistency defects; they are not evidence of a stable long-horizon solution.

### Exact post-merge verification

The PR #22 exact updated-head gates all passed at `fad6ea8c5fef7abb0bcecb8e01fdaad9b08387d3`: the direct-RHS contract passed 4 focused + **204 repository tests**; strict JSON passed 204 tests; lightweight validation passed 204 tests plus control smoke and an ordinary N=40 production step; strict retry integration passed 204 tests plus the N=40, `t=0.03` campaign-admission smoke.

The exact-current-main lightweight validation workflow [run 38331006940](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331006940) checked out `main` commit `7776a9156c24303c47a57e171a2d5997f6676668` after PR #30 was merged and passed **208 repository tests**, the control-only smoke, and one ordinary N=40 production step to `t=0.03`. The three long numerical workflows were skipped by path filters, so this is software/reliability validation, not scientific admission.

PR #30's exact-head checks also passed at `8f35c1148db76c86bae38d5aba78b4d516f65a53`: strict-JSON [38330896951](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330896951), lightweight/control/production smoke [38330896913](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330896913), and strict retry/admission [38330896900](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330896900). Each workflow passed 208 repository tests; lightweight and integrated workflows also passed their bounded smoke. Duplicate [PR #31](https://github.com/miahpenn/0-Cosmos-engine/pull/31) was closed unmerged because it was based on the pre-PR-30 main tree and omitted the initial-diagnostic `RuntimeError` path/regression included in #30.

### Remaining blockers to production admission

1. **Strict campaign admission / retry integration — software contract now promoted:** [PR #23](https://github.com/miahpenn/0-Cosmos-engine/pull/23) is merged to `main` at `5b7284600d64eaa70772a9c1e47519872a807ac8`. Its exact updated-head gates [38330179302](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330179302), [38330179305](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330179305), and [38330179296](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330179296) each passed 200 tests; these include lightweight/control/one-step smoke, strict JSON, and strict retry/admission smoke. It preserves PR #13's campaign preflight and the PR #10/#11 radiation contracts. Do not merge the stale standalone [PR #3](https://github.com/miahpenn/0-Cosmos-engine/pull/3) or the diverged [PR #16](https://github.com/miahpenn/0-Cosmos-engine/pull/16) stack; PR #23 is the clean-main integration. This closes the software/reporting integration gap, not the known archived-radiation continuation or full scientific admission.
   The follow-up reporting failure path is also promoted in [PR #30](https://github.com/miahpenn/0-Cosmos-engine/pull/30), merged at `7776a9156c24303c47a57e171a2d5997f6676668`. It recursively normalizes nested NumPy/non-finite JSON leaves and persists a fail-closed ledger/summary without re-running diagnostics on an invalid returned initial state; initial-diagnostic and fallback final-diagnostic `RuntimeError` paths are recorded as numerical failures. Exact PR-head strict JSON, lightweight, and strict retry workflows each passed 208 tests; the post-merge current-main lightweight run [38331006940](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331006940) also passed 208 tests + control/one-step smoke. Duplicate [PR #31](https://github.com/miahpenn/0-Cosmos-engine/pull/31) was closed unmerged. This closes reporting reliability gaps only, not scientific admission.
2. **Archived outer-radiation continuation:** at the archived `N=80` state near `t=47.13`, one smaller step succeeds after a predictor rejection, but the next ordinary continuation approaches the finite-pressure inversion boundary and ultimately fails closed near `t=47.153`. The strict-interior validator/retry contract improves failure classification; it does not repair the trajectory. Do not relax the cone, clip/project conserved variables, or add an arbitrary timestep floor.
3. **Initial centre geometry / CMC series:** a four-arm factorial baseline / quadrature-only / CMC-row-only / combined comparison completed on current production *numerical code* snapshot `3e5aedab6c3632d446a8ee4814432d47e14e86fc` in [run 38330898524](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330898524), artifact [11783041214](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330898524/artifacts/11783041214). The candidate source is `c2379c80ce0d9dbb598b837253291f336d559cda`; vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`. Only `engine/v55_initial.py` and `engine/cmc_gauge.py` differ under `engine/`; all 12 arms (four variants × N=80/120/160) reached matched (t=0.12) in 4/6/8 steps. The source blobs are the same candidate implementation now held by draft PR #18 (candidate initial-builder blob `2e5c41c10d4bbb53df3d7299b898e97f3239c120`; CMC blob `5ac1b1ba5cc4b840e485b74841ceb2c2b0b6b69b)). At N=160, H/M/C_lambda maxima are baseline `2.947e-3 / 1.054e-3 / 1.393e-4`; quadrature-only `6.889e-5 / 4.332e-5 / 5.620e-6`; CMC-row-only `2.884e-3 / 1.091e-3 / 1.443e-4`; combined `6.657e-6 / 1.955e-6 / 2.641e-7`. Quadrature supplies most initial constraint improvement; CMC-row-only improves the lapse even-series and K-target residual but does not fix the primary constraints; the combined candidate performs best on these short-time maxima.

   This factorial used runtime code base `3e5aedab...`, before PR #30's *runner/report-only* fix. The production geometry/matter stepping layer is unchanged by PR #30; exact post-PR-30 main smoke passed on `7776a915...`. Treat the factorial as pinned short-horizon evidence, not a replay on the post-PR-30 reporting tree. The quadrature-only arm worsens the local-time Misner–Sharp current residual/integrated mismatch, and the combined candidate has not closed fixed-worldtube conservation. The growing N=160-to-refinement M/connection maxima, unresolved centre-series limit and archived radiation crossing prevent promotion. The earlier artifact-first attribution is at [Issue #5 comment 6135914328](https://github.com/miahpenn/0-Cosmos-engine/issues/5#issuecomment-6135914328).

   A source-pinned **current-runtime** follow-up isolates baseline versus quadrature-only initialization using current production blobs and vendor. The successful three-slice audit [run 38331732142](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331732142), artifact [11782838101](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331732142/artifacts/11782838101), pins runtime source `efd006e33e95a19132af5066df21321f473580fc`, baseline initializer blob `4cfbfe4be6ac3ce6ed492f4e435186997d94116a`, candidate initializer blob `2e5c41c10d4bbb53df3d7299b898e97f3239c120`, unchanged CMC blob `b1e6b3e95ff37702509e1b38ec6401d10e2482c3`, and vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`. All six arms (baseline/candidate × N=80/120/160) completed two steps, storing t=0, dt and 2dt. At coordinate r=10, `dM_MS/dt - (flux_T + work_pR)` is worse in magnitude for quadrature-only at essentially all sampled times; at N=120, for example, t=0 changes from `+4.70e-9` to `-2.04e-7`, and t=dt from `-5.93e-8` to `-2.88e-7`. The interval defects `ΔM_MS - ∫(flux_T + work_pR)dt` also increase in absolute value at all tested resolutions: N=80 baseline `-2.32e-8, -2.64e-8` vs candidate `-3.84e-8, -4.39e-8`; N=120 `-5.46e-10, -6.52e-10` vs `-4.92e-9, -5.44e-9`; N=160 `+7.25e-10, +7.77e-10` vs `-1.12e-9, -1.20e-9`. The recorded full-radial peak remains centre-dominated and changes little, while the recorded full-radial RMS is worse at 2dt. **These global peak/RMS values are provisional:** at capture time both audit paths used the now-identified even-parity derivative for the odd areal-radius field near the origin. Draft [PR #37](https://github.com/miahpenn/0-Cosmos-engine/pull/37) corrects both production diagnostic call sites to odd parity and passes flat-centre tests; replay the full-radial centre residual after that correction lands. Do not use the pre-correction centre peak/RMS as admission evidence. The coordinate-r=10 local witness is away from the centre and remains a separately recorded finding.

The follow-up cross-application attribution is archived at [run 38331997335](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331997335), artifact [11782619185](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331997335/artifacts/11782619185), latest diagnostic PR #35 head `0108e4b54649ce2d85b52753b273499cf85c6f28`. ZIP SHA-256 `e29acc0da9ac19fb47c46777ca1762bd3cd033d6c0c028b98e550c88bcca55a9`; report SHA-256 `0d962a33fade4be9a2dc5214be94791d04c1288352a88a03faba7b8d07cf1df4`. Applying both cumulative-integral rules to the captured baseline/candidate source profiles shows the fixed-r=10 actual integral delta at N=80/120/160 is `-1.1453929961e-4, -5.2961678860e-5, -3.0082893129e-5`; the rule effect on the candidate source is `-1.1456491049e-4, -5.2974980910e-5, -3.0090701324e-5`; the source-profile response under the baseline trapezoid is only `+2.5611e-8, +1.3302e-8, +7.8082e-9`. Decomposition closes to at most `3.55e-15`. This is exact cross-application on the captured final sources, not a fully iterated causal counterfactual. It attributes almost all integral shift to the rule itself, but does not improve the observed mass-current closure.

These diagnostics interpolate at **fixed coordinate r=10**: areal radius changes from 10 to about 10.00342/10.00228/10.00171 by 2dt at N=80/120/160. They are not a fixed-areal-radius R=10 tube; the two-step horizons are only t=0.06/0.04/0.03. The complementary direct-RHS one-step audit [run 38331565460](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331565460) finds secant-vs-mean-direct-rate discrepancy of `3.28e-8, 1.32e-8, 7.09e-9` at N=80/120/160, supporting the chain-rule witness; it likewise finds worse fixed-r=10 current closure for quadrature-only.

The longer matched-time fixed-r=10 attempt [run 38331709472](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331709472) produced interpolated t=0–0.12 reports but correctly ended **failure** in nearest-cell aggregation: it expects `M_MS` at the surface top level while the nearest-cell record nests fields under `values`. Its artifact [11783052605](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331709472/artifacts/11783052605) is diagnostic data, not a passing workflow; a strict rerun remains required. The audited scope is local semi-discrete versus current witnesses, not full physical admission or continuum convergence.

4. **Conservation and convergence:** full Misner–Sharp current/flux/work balance at a fixed physical worldtube, fixed-time spatial convergence of H/M/connection/determinant and centre-series witnesses, high-resolution radiation continuation beyond one step, and the long-horizon multi-resolution campaign remain open.
5. **Matter RHS extraction — software refactor now promoted:** [PR #22](https://github.com/miahpenn/0-Cosmos-engine/pull/22) merged at `a4183bdb1b10c9ae6953af2a064af8a25415d001` after reconciliation with PR #23's typed stage-aware radiation retry contract. Exact-head gates passed 204 repository tests in each of direct-RHS, strict-JSON, lightweight/control/production-smoke, and strict-retry integration workflows; the direct-RHS suite also passed 4 focused tests. The prior seven-resolution one-step A/B found only tiny radiation-momentum rounding differences on its pinned source snapshots. This validates the narrow software refactor at the tested level, not long-trajectory equivalence or scientific admission. Do not promote [PR #8](https://github.com/miahpenn/0-Cosmos-engine/pull/8) as-is because its generic radiation flux path predates PR #11; duplicate PR #21 was closed unmerged.

### Candidate-only diagnostics

- Draft [PR #18](https://github.com/miahpenn/0-Cosmos-engine/pull/18) combines origin/interior radial quadrature with the cell-centred CMC centre-row PDE. Its latest candidate head is `3a845ec1ec903785e8dd3230ecd231ca9221ed44`; the numerical implementation blobs match those used in the four-arm artifact above. Exact-head harness/strict-report workflows [38330648443](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330648443) and [38330652884](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330652884) passed. PR #18 remains stacked on the historical PR #16 integration branch and draft-only; the factorial A/B is source-pinned against main code snapshot `3e5aedab...`, not the merged PR head. No candidate source is promoted.
- Draft [PR #35](https://github.com/miahpenn/0-Cosmos-engine/pull/35) adds a source-pinned diagnostic harness and isolated initializer candidate; it does not edit the production kernel or CMC row. Its three-slice audit [run 38331732142](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331732142) and source-attribution rerun [38331997335](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331997335) are green as diagnostic jobs, but measure a coordinate-fixed r=10 tube, not fixed areal R=10. They attribute the initial integral delta almost wholly to the quadrature rule while finding worse quadrature-only local/current-integrated residuals. The full-radial peak/RMS from these runs uses the pre-PR-37 even-parity centre stencil and is provisional; rerun centre-localized comparisons with the corrected odd-parity derivative after PR #37 merges. Keep #35 diagnostic-only and #18 isolated; this evidence does not admit the equations.
- [PR #25](https://github.com/miahpenn/0-Cosmos-engine/pull/25) is the earlier candidate-branch matched baseline A/B on baseline `9f0bd16d...` and candidate `be6926d...`, with the same quadrature/CMC implementation blobs. It is useful historical attribution but superseded for current-main software-lineage questions by [run 38330898524](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330898524). Keep the candidate PR isolated; short-horizon improvements are not a formal convergence proof.
- The diagnostic harness hardening PR #19 was merged **only into the isolated candidate branch** at `3a845ec1...`, never into production `main`. It enforces required radiation witness presence/finite values and correct determinant residual reporting. Its diagnostic pass does not admit the underlying quadrature/CMC candidate.
- Draft [PR #9](https://github.com/miahpenn/0-Cosmos-engine/pull/9) remains an unadmitted radiation four-force counterfactual; [PR #2](https://github.com/miahpenn/0-Cosmos-engine/pull/2) retains the N=320 reconstruction audit's frozen-tolerance/host-variance and source-preregistration gaps; [PR #1](https://github.com/miahpenn/0-Cosmos-engine/pull/1) remains experimental strong-field/central-clock diagnostics. None is a production promotion candidate in its current state.
- The prior one-step sweep and four-arm factorial results are diagnostics, not a physical-law conclusion. No long campaign is authorized until the centre-series, archived outer radiation continuation, full H/M/connection constraints, and fixed-worldtube Misner–Sharp/convergence gates close.

---

## Architecture baseline from the 2026-10-11 snapshot

The canonical repository contains the unified production-facing state graph:

1. corrected archive COSMOS background equations;
2. local GEAR-03 S/D sector with GEAR-63/65 barrier;
3. scalar S/D/COSMOS evolution and stress projections;
4. conservative spherical dark-matter, baryon, and radiation matter sectors;
5. locked beta DM/scalar four-force exchange;
6. one total Einstein stress-energy assembly;
7. pinned reference-metric BSSN/PIRK numerical geometry;
8. synchronized geometry/scalar/matter evolution state;
9. dynamically derived proper time, H_eff, e-folds, trapping, and Misner-Sharp observables;
10. cycle/handoff ledgers without a synthetic lambda_cycle;
11. checkpoint and long-horizon multi-resolution campaign infrastructure.

## Explicit exclusions

- Bianchi-I shear remains outside the exact spherical local source.
- No D-to-matter identification law.
- No shell EOS/action.
- No phenomenological interface source.
- No fitted coefficient.
- No imposed clock conversion.
- No bounce/reset/branch flip/lapse floor/ejection threshold/physical stop.

## Historical validation snapshot (2026-10-11; superseded by the exact current-main check above)

An independent repository review reported the following checks on 2026-10-11:

- Full test suite: 79/79 passed, including tests in `tests/` and the vendor BSSN/PIRK operator suite.
- Production smoke campaign: `resolution=40`, `r_max=40`, `final_time=1.0`; initialization and a kernel step completed without NaNs or assertion failures.
- Default-path kernel initialization and stepping succeeded.

These results establish basic test and smoke-run functionality. They do **not** close the full production campaign gate.

Current classification:

**IMPLEMENTED / TESTED AT UNIT AND SMOKE-RUN LEVEL / NOT YET FULL-CAMPAIGN-VALIDATED**

## Earlier numerical-audit notes (2026-10-11; use Issues #4/#5 for current results)

- Draft PR #2 tracks the N=320 reconstruction audit and host/runtime-dependent outer Hamiltonian L2 discrepancies. Frozen tolerances and the missing preregistration source documents must be resolved transparently; do not relax tolerances solely to make runners pass.
- Draft PR #1 contains experimental strong-field representation and central-clock diagnostics. It is not production validation and must remain separate from `main`.
- The N=640 signed-pole audit through t=12 is running on `diag/early-signed-pole-splice-audit`; its numerical result is pending until its artifact is available and admitted.
- The bond-sweep diagnostic artifact completed successfully. Its coupling-dependent mode differences are diagnostic evidence, not proof that the bond explains the full evolving constraint discrepancy.

## Full scientific engineering gate

The remaining work is the complete regression and long-horizon multi-resolution campaign, followed by analysis of:

- Hamiltonian, momentum, connection, and determinant witnesses;
- invariant trapped-surface roots;
- Misner-Sharp current balance;
- proper-time and H_eff trajectories;
- dark-sector exchange closure;
- resolution convergence and host sensitivity;
- turnaround and re-expansion events;
- genuinely emergent multi-cycle behavior, if the solved equations produce it.

A numerical failure is classified as a numerical/gauge result unless invariant, resolution-independent evidence supports a physical interpretation. No production physics change is authorized by a passing unit suite or short smoke campaign alone.

## Numerical reproducibility note

Host/runtime-dependent reconstruction discrepancies have been reported in the N=320 audit. Preserve the preregistered tolerance and report failing diagnostics, runtime details, and source hashes. Do not silently loosen tolerances or omit failed fields.

Campaign implementation note: the radiation inversion bracket uses the adjacent representable value only when endpoint roundoff reverses the analytically positive residual; the physical EOS and bound remain unchanged.
