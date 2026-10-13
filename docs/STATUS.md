# Build Status — 2026-10-13



## Current production snapshot — 2026-10-13

**Production ref at this snapshot:** `main` commit `1db7d554ae43e447e8b4859b68ff8d42e224eb61`.

**Disposition:** software-contract corrections have been promoted and post-merge smoke validation is green. The engine is **not yet admitted by the full scientific campaign gate**.

### Software promotions present on `main`

- [PR #12](https://github.com/miahpenn/0-Cosmos-engine/pull/12) — lightweight unit/control/one-step validation workflow.
- [PR #13](https://github.com/miahpenn/0-Cosmos-engine/pull/13) — campaign configuration, floating-point time-progress, and checkpoint-name preflight.
- [PR #14](https://github.com/miahpenn/0-Cosmos-engine/pull/14) — shared SI-agent communications protocol.
- [PR #15](https://github.com/miahpenn/0-Cosmos-engine/pull/15) — scalar lapse-gradient contraction uses the declared physical inverse radial metric `gamma^rr = X^2/a`.
- [PR #10](https://github.com/miahpenn/0-Cosmos-engine/pull/10) — removes the inadmissibility guard from only the fictitious `dt=1` RHS-extraction intermediate, while retaining default validation for direct updates and validating actual production stages.
- [PR #11](https://github.com/miahpenn/0-Cosmos-engine/pull/11) — preserves exact-zero radiation rest-current through radiation-specific HLL/outer-face flux dispatch.

These changes close discrete software and equation-consistency defects; they are not evidence of a stable long-horizon solution.

### Exact post-merge verification

The lightweight validation workflow [run 38329612360](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38329612360) checked out the exact `main` SHA above and passed:

- **126 repository tests**;
- the control-only smoke;
- one ordinary production step at `N=40`, `r_max=40`, `dt=0.03`, reaching `t=0.03` with finite diagnostics and positive geometry.

This is a current-main integration check. The workflow's three long numerical campaigns were not run; it does not close scientific admission.

### Remaining blockers to production admission

1. **Strict campaign admission / retry integration:** draft [PR #3](https://github.com/miahpenn/0-Cosmos-engine/pull/3) must not be merged standalone; its older runner can emit non-standard `NaN` for undefined short-history Misner–Sharp summaries. Draft [PR #16](https://github.com/miahpenn/0-Cosmos-engine/pull/16) has the strict-JSON and stage-aware retry changes, but is stacked on PR #3 and still needs reconciliation against the current `main` tree (including PRs #10/#11 and the checkpoint preflight) with exact merged-source CI. Do not bypass its admission checks just to obtain a completed status.
2. **Archived outer-radiation continuation:** at the archived `N=80` state near `t=47.13`, one smaller step succeeds after a predictor rejection, but the next ordinary continuation approaches the finite-pressure inversion boundary and ultimately fails closed near `t=47.153`. The strict-interior validator/retry contract improves failure classification; it does not repair the trajectory. Do not relax the cone, clip/project conserved variables, or add an arbitrary timestep floor.
3. **Initial centre geometry / CMC series:** the combined initial-quadrature + CMC centre-row candidate in draft [PR #18](https://github.com/miahpenn/0-Cosmos-engine/pull/18) passes focused/full software tests and bounded candidate-only checks. The CMC-row-only control on the default initializer did not fix the dominant centre H/M/connection witnesses. A matched-time baseline-versus-combined-candidate comparison is the required next discriminator; candidate-only trajectories are not a comparative admission result.
4. **Conservation and convergence:** full Misner–Sharp current/flux/work balance at a fixed physical worldtube, fixed-time spatial convergence of H/M/connection/determinant and centre-series witnesses, high-resolution radiation continuation beyond one step, and the long-horizon multi-resolution campaign remain open.
5. **Competing RHS extraction design:** the matched finite-`dt` A/B is complete in [run 38329555957](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38329555957): seven resolutions, identical initialized states/vendor/common-source pins, all 21 one-step species updates completed; DM/baryon fields and radiation energy were bit-identical, with only a tiny radiation-momentum difference. Do not promote [PR #8](https://github.com/miahpenn/0-Cosmos-engine/pull/8) as-is: its stale generic HLL/outer-face path predates PR #11 and would regress the exact-zero radiation rest-current invariant. A clean current-main API port is being prepared to retain PR #10's finite-time arithmetic, its actual-stage guards, and PR #11's radiation-specific dispatch; that port still needs focused/full-suite validation.

### Candidate-only diagnostics

- Draft [PR #18](https://github.com/miahpenn/0-Cosmos-engine/pull/18) was refreshed at candidate head `be6926d079a4f8c9c919554fedb2a2a92182e987`. Exact-head workflow [38329679022](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38329679022) passed 13 focused tests, 165 repository tests, and the bounded N=40/`t=0.03` smoke. The corrected candidate-only matched-time trajectory [38329573898](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38329573898) reached `t=0.12` at N=80/120/160 (4/6/8 steps); H falls with refinement, but M/connection are not uniformly monotone and fixed-worldtube time-integrated Misner–Sharp balance remains open.
- Draft [PR #19](https://github.com/miahpenn/0-Cosmos-engine/pull/19) hardens the multiresolution diagnostic's determinant and radiation witnesses. Its exact-head harness tests and five-resolution one-step run passed on its pinned earlier candidate base, but the candidate initializer has since advanced; **do not use that older artifact as validation of the latest initializer blob**.
- The latest default-initializer baseline-versus-combined-candidate matched-time A/B is assigned on the coordination board and remains the next discriminator. Keep PR #18/#19 isolated from `main` until the exact-source A/B, centre-series/constraint attribution, archived radiation continuation and Misner–Sharp balance are independently reconciled.
- No bounce, turnaround, re-expansion, completed cycle, or multi-cycle physical claim is admitted by the current software/smoke evidence.

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
