# Build Status — 2026-10-13

## Current engineering snapshot — 2026-10-13 17:20 UTC (radiation/constraint reconciliation merged; source baseline unchanged)

**Last non-status main commit before this status refresh:** `f7034cefaa4ca3eb7d3f4fa7eb7281af59cb5211` (merged PR #122; source-pinned radiation/constraint report and CSV evidence only). **Production source baseline:** `bfbe59c5dc844ce0ee365fdfcdc963dfd681d403` (last production-source promotion; PR #131; production engine/vendor/physics source is unchanged by #122/#136/#137/#139/#141). Current-source-tree software validation [run 38437180218](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38437180218) on the immediate parent `592d0fc06654ebd95be78b2b61d71fce23bc4696` passed **320 tests** and the bounded N=40 one-step smoke to `t=0.03`; #122 adds only evidence Markdown/CSVs on top of that tree, so the production/test tree is unchanged. PR #141 removed only `.chatgpt-write-test.md` and two byte-identical duplicate test files (their Git blobs matched retained canonical tests); PR #136 added only research-only EM preflight tests. Campaign/zero-star/worldtube jobs were skipped by the current `push.paths` filters; this is not a scientific-campaign result.

### Recent qualified promotions

- **[PR #124 — run the three stable CI workflows on every pull request](https://github.com/miahpenn/0-Cosmos-engine/pull/124), merge `c7a57482902df460343ba50b1c286f32514cfa7f`.** Removes only `pull_request.paths` filters from lightweight, strict JSON and integrated retry. The existing push filters and job bodies stayed unchanged. This fixes skipped PR-check triggering; it does not itself configure required branch-protection checks.
- **[PR #128 — reject overflowing JSON evidence numbers](https://github.com/miahpenn/0-Cosmos-engine/pull/128), merge `60ee4554f14f8f6ea9d23863729f6cde1300d77f`.** The common evidence parser rejects syntactically valid exponent tokens that convert to `inf`/`-inf`; six regressions cover unchecked fields in summary, ledger and retry evidence.
- **[PR #116 — fail closed on invalid CMC proper-volume weights](https://github.com/miahpenn/0-Cosmos-engine/pull/116), merge `5426d31dd9060141ba3300e2fd3a12aa86ed60f2`.** Validates grid/metric/weight inputs before CMC consumers and checks separate L2/L3 K-RHS shapes/finiteness before summing. No CMC equation or physical source changed.
- **[PR #127 — preserve the solved radiation velocity norm](https://github.com/miahpenn/0-Cosmos-engine/pull/127), merge `faf8b81f3b27666d3401b51fc04ab4b9b8c52c7c`.** Makes the radiation pressure residual energy-relative and keeps the inversion's solved `v²` rather than recomputing it from a rounded coordinate velocity. The strict `E > |S|` rule remains unchanged.
- **[PR #130 — propagate solved `v²` through flux and shifted tensor products](https://github.com/miahpenn/0-Cosmos-engine/pull/130), merge `86010d095a90c94d53eb39043b63fcffc386cb2d`.** Keeps `_valencia_flux`, `_radiation_flux`, `stress_tensor_contravariant` and `mixed_flux` consistent with the preserved norm. This is algebraic representation consistency, not a new physical law.
- **[PR #131 — pin the binary64 flux-rounding regression](https://github.com/miahpenn/0-Cosmos-engine/pull/131), merge `bfbe59c5dc844ce0ee365fdfcdc963dfd681d403`.** Test-only. Includes zero/nonzero-shift coverage and a density-rescaled fixture that proves the legacy rounded `S_r*v_r` result differs from the preserved-norm expression after final float64 rounding. The first fixture scale failed to discriminate; a corrected scale then passed all three applicable exact-head workflows (356 tests each).
- **[PR #141 — remove repository write probe and identical test duplicates](https://github.com/miahpenn/0-Cosmos-engine/pull/141), merge `2b96c5b2df5d72578520f99b8881a8bb075d3cc2`.** Deletes only the unused `.chatgpt-write-test.md` note and two test files whose Git blob SHA exactly matched retained canonical tests. No engine/source/unique assertion/physics/campaign change. Exact merged-tree lightweight [38436620041](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38436620041) passed 317 tests and N=40 smoke; the lower count is from duplicate cases no longer being executed twice.
- **[PR #136 — combined research-only EM preflight suite](https://github.com/miahpenn/0-Cosmos-engine/pull/136), merge `592d0fc06654ebd95be78b2b61d71fce23bc4696`.** Adds two research-only test files for E0 projection parity and symbolic radial-electric consumer gain; the three exact-head stable PR workflows passed (lightweight [38433869059](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38433869059), strict JSON [38433868976](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38433868976), integrated retry [38433868993](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38433868993)). This is not Maxwell evolution, a selected charge/current law, or normalization.
- **[PR #122 — radiation/constraint evidence reconciliation](https://github.com/miahpenn/0-Cosmos-engine/pull/122), merge `f7034cefaa4ca3eb7d3f4fa7eb7281af59cb5211`.** Five evidence/report/CSV files only. Exact-head lightweight [38434108105](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38434108105), strict JSON [38434107854](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38434107854), and integrated retry [38434107673](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38434107673) passed. Records 27/27 saved accepted-state tangents negative and downstream-only Cdot delta at roundoff scale; no new trajectory or admission result.

These promotions improve evidence integrity, numerical fail-closed behavior and consistency of near-cone primitive consumers. They do **not** establish a convergent or scientifically admitted trajectory.

### Exact-head software evidence

- PR #116 refreshed candidate head `2775712cfc24fd99a70db042d669da155f7e3383`: three workflows, **340 tests each** — lightweight [38428811891](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38428811891), strict JSON [38428811916](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38428811916), integrated retry [38428811886](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38428811886).
- PR #127 exact tested head `5d2bda38b7172d3dee22ee5cd37eee7855f2a5ba`: four workflows, **353 tests each** — lightweight [38430157937](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38430157937), strict JSON [38430157952](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38430157952), integrated retry [38430158033](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38430158033), direct RHS [38430157929](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38430157929).
- PR #130 exact tested head `ece35e1f1066b7d3e1c4733cf34b4fb83bc28a35`: four workflows, **354 tests each** — lightweight [38430967360](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38430967360), strict JSON [38430967292](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38430967292), integrated retry [38430967422](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38430967422), direct RHS [38430967404](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38430967404).
- PR #131 exact tested head `5c5c7cba9253bac072b6d05aec9f98a5771e1c23`: three applicable workflows, **356 tests each** — lightweight [38431806877](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38431806877), strict JSON [38431806868](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38431806868), integrated retry [38431806965](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38431806965). It changed only a test file; the dedicated direct-RHS workflow did not trigger, but its test is included in the full repository suite.

### Current radiation and constraint evidence

- The checksummed #61 N=80 replay and source-pinned offline joins remain the baseline. On outer cell 79, accepted state + accepted metric gives `|S|/E = 0.998857696555`; putting that same conserved state on the first rejected predictor metric gives `1.006073578641`, while putting predictor conserved variables on the accepted metric gives `0.993226051837`; the combined actual predictor is `1.000401250166`. The immediate crossing is metric-dominated in this factorization, but that does not establish why the accumulated geometry reaches this state.
- The saved-state energy-face-flux identity is `F_E = alpha*gamma_rr_inv*U_r - beta*U_E`. Before PR #127, the archived outer flux was `2.8454937981219525e-4` versus the conserved-state identity `3.082618281298782e-4`, because a one-ULP `v²` round-trip changed `W²` near the cone. PRs #127/#130 correct the representation through the flux, stress and projection consumers without relaxing admissibility.
- The post-#127/#130 current-consumer frozen-state audit is now complete on the checksummed #61 replay: exact-key joins cover 27 saved accepted-state RHS/geometry stages, and the cone tangent remains negative in 27/27. It ranges from `Cdot_total = -1.2010510e-5` to `-1.3546275e-5` as the outer ratio approaches `0.9999999999999998`; a paired calculation changing only the downstream quadratic products alters the tangent by at most `1.63e-19`, so that consumer rewrite alone does not reverse the tendency. See [PR #122](https://github.com/miahpenn/0-Cosmos-engine/pull/122) and the [27-row exact-join CSV](https://github.com/miahpenn/0-Cosmos-engine/blob/435712e9d09bcdcfdf9c1fdcb38c715c91f722f1/diagnostics/radiation_accepted_state_tangent_current_consumers_2026-10-13.csv). This is a same-saved-state offline formula recomputation, not a replay using the pre-#127 primitive-recovery path and not next-step admissibility evidence. No new long trajectory has been run; scientific admission remains NOT ADMITTED.
- The completed-step implicit `Aa/K` arithmetic ledger was reconstructed offline for all 26 successful replay steps: recomputed `L2_final` matched saved stages to about `5.88e-15` for `Aa` and `4.11e-15` for `K`; the implemented update reconstructed saved interior `Aa` and all accepted `K` profiles exactly at binary64 precision. This closes a capture/reconstruction gap, not physical correctness, boundary accuracy or convergence.

### Scientific admission status

**NOT ADMITTED for production scientific campaigns.** Independent gates still open:

1. Radiation continuation and multi-resolution convergence after the near-cone consistency corrections. The consumer-complete frozen-state tangent was recalculated on the current production-source consumers and remains negative in 27/27 saved accepted-state joins; the downstream-consumer-only delta is at roundoff scale. This removes a stale diagnostic task, not the underlying gate: no next-step admissibility or long-horizon recovery has been established.
2. Centre regularity and H/M/connection/determinant constraint propagation, including the initialized high-resolution strict failure at N=120/160. Do not treat the short CMC/radiation diagnostics as closure.
3. Fixed-areal-radius Misner–Sharp mass/current/flux/work balance and long-horizon spatial/temporal convergence. Existing `r=10` coordinate-fixed results are not a fixed-areal-radius `R=10` witness.
4. Scalar–DM exchange: PR #120's source signs/index/lapse/shift projection are internally consistent conditional on the archived variable-mass law, but the origin of `m(phi)=m0 exp(-beta phi)` in 0* microphysics and the exact value/sign of `beta=-0.04` remain unpaid premises. Do not merge this production source candidate on algebra/tests alone.
5. Electromagnetism remains research-only. PR #136 is merged as a test-only preflight: E0 projection parity plus symbolic radial-electric consumer gain. It selects no Maxwell evolution, charge/current law or physical normalization and is not scientific admission. PR #101 remains a design/gate document.

### Governance and coordination

- PR #124 corrected PR-trigger coverage, but branch protection is still an owner/admin task. The connector has not verified required checks/reviews are enforced, and the blocked-merge negative test remains outstanding.
- The live open-PR queue reports **25 open PRs** at the latest query; the older #115 count of 29 is stale. Confirm each candidate's actual current base/head and exact-head checks before promotion.
- Continue using [Issue #5 — cross-agent coordination](https://github.com/miahpenn/0-Cosmos-engine/issues/5) and [Issue #4 — momentum-constraint propagation](https://github.com/miahpenn/0-Cosmos-engine/issues/4). Software tests, observer identities and frozen-state counterfactuals must remain explicitly separate from scientific admission.

---

## Historical production snapshot — 2026-10-13 (superseded by the current snapshot above)

**Production ref at this snapshot:** `main` commit `2454f4094f09f94343a3ba890b69d51bbb2defc4`, including the merged diagnostic odd-parity correction in PR #37. Exact current-main software validation: [run 38332385030](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38332385030).

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
- [PR #37](https://github.com/miahpenn/0-Cosmos-engine/pull/37) — corrects the regular-centre odd-parity derivative used by worldtube and production invariant diagnostics, with flat-centre regressions; no evolution equations or physical sources change.

These changes close discrete software and equation-consistency defects; they are not evidence of a stable long-horizon solution.

### Exact post-merge verification

The PR #22 exact updated-head gates all passed at `fad6ea8c5fef7abb0bcecb8e01fdaad9b08387d3`: the direct-RHS contract passed 4 focused + **204 repository tests**; strict JSON passed 204 tests; lightweight validation passed 204 tests plus control smoke and an ordinary N=40 production step; strict retry integration passed 204 tests plus the N=40, `t=0.03` campaign-admission smoke.

The post-PR-30 historical validation [run 38331006940](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331006940) checked out `main` commit `7776a9156c24303c47a57e171a2d5997f6676668` and passed **208 repository tests**, the control-only smoke, and one ordinary N=40 production step to `t=0.03`.

The exact-current-main lightweight workflow [run 38332385030](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38332385030) checked out `main` commit `2454f4094f09f94343a3ba890b69d51bbb2defc4` after PR #37 was merged. It passed **210 repository tests**, the control-only smoke, and one ordinary N=40 production step to `t=0.03`. The three long numerical workflows were skipped by path filters; this validates software/reliability at the tested scale, not the full scientific campaign.

PR #30's exact-head checks also passed at `8f35c1148db76c86bae38d5aba78b4d516f65a53`: strict-JSON [38330896951](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330896951), lightweight/control/production smoke [38330896913](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330896913), and strict retry/admission [38330896900](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330896900). Each workflow passed 208 repository tests; lightweight and integrated workflows also passed their bounded smoke. Duplicate [PR #31](https://github.com/miahpenn/0-Cosmos-engine/pull/31) was closed unmerged because it was based on the pre-PR-30 main tree and omitted the initial-diagnostic `RuntimeError` path/regression included in #30.

### Remaining blockers to production admission

1. **Strict campaign admission / retry integration — software contract now promoted:** [PR #23](https://github.com/miahpenn/0-Cosmos-engine/pull/23) is merged to `main` at `5b7284600d64eaa70772a9c1e47519872a807ac8`. Its exact updated-head gates [38330179302](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330179302), [38330179305](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330179305), and [38330179296](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330179296) each passed 200 tests; these include lightweight/control/one-step smoke, strict JSON, and strict retry/admission smoke. It preserves PR #13's campaign preflight and the PR #10/#11 radiation contracts. Do not merge the stale standalone [PR #3](https://github.com/miahpenn/0-Cosmos-engine/pull/3) or the diverged [PR #16](https://github.com/miahpenn/0-Cosmos-engine/pull/16) stack; PR #23 is the clean-main integration. This closes the software/reporting integration gap, not the known archived-radiation continuation or full scientific admission.
   The follow-up reporting failure path is also promoted in [PR #30](https://github.com/miahpenn/0-Cosmos-engine/pull/30), merged at `7776a9156c24303c47a57e171a2d5997f6676668`. It recursively normalizes nested NumPy/non-finite JSON leaves and persists a fail-closed ledger/summary without re-running diagnostics on an invalid returned initial state; initial-diagnostic and fallback final-diagnostic `RuntimeError` paths are recorded as numerical failures. Exact PR-head strict JSON, lightweight, and strict retry workflows each passed 208 tests; the post-merge current-main lightweight run [38331006940](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331006940) also passed 208 tests + control/one-step smoke. Duplicate [PR #31](https://github.com/miahpenn/0-Cosmos-engine/pull/31) was closed unmerged. This closes reporting reliability gaps only, not scientific admission.
2. **Archived outer-radiation continuation:** at the archived `N=80` state near `t=47.13`, one smaller step succeeds after a predictor rejection, but the next ordinary continuation approaches the finite-pressure inversion boundary and ultimately fails closed near `t=47.153`. The strict-interior validator/retry contract improves failure classification; it does not repair the trajectory. Do not relax the cone, clip/project conserved variables, or add an arbitrary timestep floor.
3. **Initial centre geometry / CMC series:** a four-arm factorial baseline / quadrature-only / CMC-row-only / combined comparison completed on current production *numerical code* snapshot `3e5aedab6c3632d446a8ee4814432d47e14e86fc` in [run 38330898524](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330898524), artifact [11783041214](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330898524/artifacts/11783041214). The candidate source is `c2379c80ce0d9dbb598b837253291f336d559cda`; vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`. Only `engine/v55_initial.py` and `engine/cmc_gauge.py` differ under `engine/`; all 12 arms (four variants × N=80/120/160) reached matched (t=0.12) in 4/6/8 steps. The source blobs are the same candidate implementation now held by draft PR #18 (candidate initial-builder blob `2e5c41c10d4bbb53df3d7299b898e97f3239c120`; CMC blob `5ac1b1ba5cc4b840e485b74841ceb2c2b0b6b69b)). At N=160, H/M/C_lambda maxima are baseline `2.947e-3 / 1.054e-3 / 1.393e-4`; quadrature-only `6.889e-5 / 4.332e-5 / 5.620e-6`; CMC-row-only `2.884e-3 / 1.091e-3 / 1.443e-4`; combined `6.657e-6 / 1.955e-6 / 2.641e-7`. Quadrature supplies most initial constraint improvement; CMC-row-only improves the lapse even-series and K-target residual but does not fix the primary constraints; the combined candidate performs best on these short-time maxima.

   This factorial used runtime code base `3e5aedab...`, before PR #30's *runner/report-only* fix. The production geometry/matter stepping layer is unchanged by PR #30; exact post-PR-30 main smoke passed on `7776a915...`. Treat the factorial as pinned short-horizon evidence, not a replay on the post-PR-30 reporting tree. The quadrature-only arm worsens the local-time Misner–Sharp current residual/integrated mismatch, and the combined candidate has not closed fixed-worldtube conservation. The growing N=160-to-refinement M/connection maxima, unresolved centre-series limit and archived radiation crossing prevent promotion. The earlier artifact-first attribution is at [Issue #5 comment 6135914328](https://github.com/miahpenn/0-Cosmos-engine/issues/5#issuecomment-6135914328).

   A source-pinned **current-runtime** follow-up isolates baseline versus quadrature-only initialization using current production blobs and vendor. The successful three-slice audit [run 38331732142](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331732142), artifact [11782838101](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331732142/artifacts/11782838101), pins runtime source `efd006e33e95a19132af5066df21321f473580fc`, baseline initializer blob `4cfbfe4be6ac3ce6ed492f4e435186997d94116a`, candidate initializer blob `2e5c41c10d4bbb53df3d7299b898e97f3239c120`, unchanged CMC blob `b1e6b3e95ff37702509e1b38ec6401d10e2482c3`, and vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`. All six arms (baseline/candidate × N=80/120/160) completed two steps, storing t=0, dt and 2dt. At coordinate r=10, `dM_MS/dt - (flux_T + work_pR)` is worse in magnitude for quadrature-only at essentially all sampled times; at N=120, for example, t=0 changes from `+4.70e-9` to `-2.04e-7`, and t=dt from `-5.93e-8` to `-2.88e-7`. The interval defects `ΔM_MS - ∫(flux_T + work_pR)dt` also increase in absolute value at all tested resolutions: N=80 baseline `-2.32e-8, -2.64e-8` vs candidate `-3.84e-8, -4.39e-8`; N=120 `-5.46e-10, -6.52e-10` vs `-4.92e-9, -5.44e-9`; N=160 `+7.25e-10, +7.77e-10` vs `-1.12e-9, -1.20e-9`. The original full-radial max/RMS used the pre-PR-37 even-parity derivative and is not a valid centre-localized balance witness. A same-trajectory parity-counterfactual audit [run 38332132436](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38332132436), artifact [11783331787](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38332132436/artifacts/11783331787), computes the corrected odd-parity residual without changing or rerunning the trajectories. It drops the t=0 full-domain peak from about `9.48e-4 / 6.32e-4 / 4.74e-4` (N=80/120/160 old-even) to `1.08e-6 / 2.40e-7 / 1.57e-7`. At 2dt, corrected odd-parity peak and RMS favor quadrature-only at all three N: max residual baseline/candidate is `2.863e-5 / 2.004e-5`, `1.272e-5 / 9.251e-6`, `7.156e-6 / 5.285e-6`; RMS is `4.123e-6 / 3.212e-6`, `1.494e-6 / 1.196e-6`, `7.285e-7 / 5.872e-7`. At t=0, candidate peak is larger at all N, so preserve time dependence. The fixed-coordinate r=10 residual is unchanged by the centre-stencil correction and remains worse in the candidate. PR #37 is now merged at `2454f4094f09f94343a3ba890b69d51bbb2defc4`; it routes production diagnostic call sites through the same odd-parity helper. Use the same-trajectory odd-parity measurements as the corrected diagnostic values and do not cite the old-even global values.

The follow-up cross-application attribution is archived at [run 38331997335](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331997335), artifact [11782619185](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331997335/artifacts/11782619185), latest diagnostic PR #35 head `0108e4b54649ce2d85b52753b273499cf85c6f28`. ZIP SHA-256 `e29acc0da9ac19fb47c46777ca1762bd3cd033d6c0c028b98e550c88bcca55a9`; report SHA-256 `0d962a33fade4be9a2dc5214be94791d04c1288352a88a03faba7b8d07cf1df4`. Applying both cumulative-integral rules to the captured baseline/candidate source profiles shows the fixed-r=10 actual integral delta at N=80/120/160 is `-1.1453929961e-4, -5.2961678860e-5, -3.0082893129e-5`; the rule effect on the candidate source is `-1.1456491049e-4, -5.2974980910e-5, -3.0090701324e-5`; the source-profile response under the baseline trapezoid is only `+2.5611e-8, +1.3302e-8, +7.8082e-9`. Decomposition closes to at most `3.55e-15`. This is exact cross-application on the captured final sources, not a fully iterated causal counterfactual. It attributes almost all integral shift to the rule itself, but does not improve the observed mass-current closure.

These diagnostics interpolate at **fixed coordinate r=10**: areal radius changes from 10 to about 10.00342/10.00228/10.00171 by 2dt at N=80/120/160. They are not a fixed-areal-radius R=10 tube; the two-step horizons are only t=0.06/0.04/0.03. The complementary direct-RHS one-step audit [run 38331565460](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331565460) finds secant-vs-mean-direct-rate discrepancy of `3.28e-8, 1.32e-8, 7.09e-9` at N=80/120/160, supporting the chain-rule witness; it likewise finds worse fixed-r=10 current closure for quadrature-only.

The longer matched-time fixed-r=10 attempt [run 38331709472](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331709472) produced interpolated t=0–0.12 reports but correctly ended **failure** in nearest-cell aggregation: it expects `M_MS` at the surface top level while the nearest-cell record nests fields under `values`. Its artifact [11783052605](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331709472/artifacts/11783052605) is diagnostic data, not a passing workflow; a strict rerun remains required. The audited scope is local semi-discrete versus current witnesses, not full physical admission or continuum convergence.

4. **Conservation and convergence:** full Misner–Sharp current/flux/work balance at a fixed physical worldtube, fixed-time spatial convergence of H/M/connection/determinant and centre-series witnesses, high-resolution radiation continuation beyond one step, and the long-horizon multi-resolution campaign remain open.
5. **Matter RHS extraction — software refactor now promoted:** [PR #22](https://github.com/miahpenn/0-Cosmos-engine/pull/22) merged at `a4183bdb1b10c9ae6953af2a064af8a25415d001` after reconciliation with PR #23's typed stage-aware radiation retry contract. Exact-head gates passed 204 repository tests in each of direct-RHS, strict-JSON, lightweight/control/production-smoke, and strict-retry integration workflows; the direct-RHS suite also passed 4 focused tests. The prior seven-resolution one-step A/B found only tiny radiation-momentum rounding differences on its pinned source snapshots. This validates the narrow software refactor at the tested level, not long-trajectory equivalence or scientific admission. Do not promote [PR #8](https://github.com/miahpenn/0-Cosmos-engine/pull/8) as-is because its generic radiation flux path predates PR #11; duplicate PR #21 was closed unmerged.

### Candidate-only diagnostics

- Draft [PR #18](https://github.com/miahpenn/0-Cosmos-engine/pull/18) combines origin/interior radial quadrature with the cell-centred CMC centre-row PDE. Its latest candidate head is `3a845ec1ec903785e8dd3230ecd231ca9221ed44`; the numerical implementation blobs match those used in the four-arm artifact above. Exact-head harness/strict-report workflows [38330648443](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330648443) and [38330652884](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330652884) passed. PR #18 remains stacked on the historical PR #16 integration branch and draft-only; the factorial A/B is source-pinned against main code snapshot `3e5aedab...`, not the merged PR head. No candidate source is promoted.
- Draft [PR #35](https://github.com/miahpenn/0-Cosmos-engine/pull/35) adds a source-pinned diagnostic harness and isolated initializer candidate; it does not edit the production kernel or CMC row. Its three-slice audit [run 38331732142](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331732142) and source-attribution rerun [38331997335](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331997335) are green as diagnostic jobs, but measure a coordinate-fixed r=10 tube, not fixed areal R=10. They attribute the initial integral delta almost wholly to the quadrature rule while finding worse quadrature-only local/current-integrated residuals. The same-trajectory counterfactual [run 38332132436](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38332132436) computes global residuals with the corrected odd-parity derivative: its t=2dt full-domain peak/RMS favor quadrature-only at all N, while its t=0 peak is higher and fixed-coordinate r=10 current closure remains worse. The pre-correction even-parity global values are invalid near the centre; use the odd-parity numbers and keep the spatial/time-local qualifications explicit. PR #37 is merged at `2454f4094f09f94343a3ba890b69d51bbb2defc4` and applies that helper to the production diagnostic code paths. Keep #35 diagnostic-only and #18 isolated; this evidence does not admit the equations.
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
