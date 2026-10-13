# Build Status — 2026-10-13



## Current production snapshot — 2026-10-13

**Production ref at this snapshot:** `main` commit `7776a9156c24303c47a57e171a2d5997f6676668`.

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
- [PR #30](https://github.com/miahpenn/0-Cosmos-engine/pull/30) — preserves fail-closed per-resolution ledgers and campaign summaries when initialization/initial diagnostics fail, and recursively normalizes non-finite NumPy values before strict JSON serialization.

These changes close discrete software and equation-consistency defects; they are not evidence of a stable long-horizon solution.

### Exact post-merge verification

The PR #22 exact updated-head gates all passed at `fad6ea8c5fef7abb0bcecb8e01fdaad9b08387d3`: the direct-RHS contract passed 4 focused + **204 repository tests**; strict JSON passed 204 tests; lightweight validation passed 204 tests plus control smoke and an ordinary N=40 production step; strict retry integration passed 204 tests plus the N=40, `t=0.03` campaign-admission smoke.

The post-PR-30 lightweight workflow [run 38331006940](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38331006940) checked out exact production `main` commit `7776a9156c24303c47a57e171a2d5997f6676668` and passed **208 repository tests**, the control-only smoke, and one ordinary N=40 production step to `t=0.03`. PR #30’s exact candidate head `8f35c1148db76c86bae38d5aba78b4d516f65a53` also passed all three [strict JSON](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330896951), [lightweight](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330896913), and [strict retry integration](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330896900) workflows with 208 tests in each; control and bounded N=40 smokes passed. The three long numerical workflows on current main were **skipped by path filters**, not executed. Software CI/smoke is not scientific admission.

### Remaining blockers to production admission

1. **Strict campaign admission / retry integration — software contract now promoted:** [PR #23](https://github.com/miahpenn/0-Cosmos-engine/pull/23) is merged to `main` at `5b7284600d64eaa70772a9c1e47519872a807ac8`. Its exact updated-head gates [38330179302](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330179302), [38330179305](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330179305), and [38330179296](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330179296) each passed 200 tests; these include lightweight/control/one-step smoke, strict JSON, and strict retry/admission smoke. It preserves PR #13's campaign preflight and the PR #10/#11 radiation contracts. Do not merge the stale standalone [PR #3](https://github.com/miahpenn/0-Cosmos-engine/pull/3) or the diverged [PR #16](https://github.com/miahpenn/0-Cosmos-engine/pull/16) stack; PR #23 is the clean-main integration. This closes the core admission/retry integration gap. PR #30’s follow-up reporting hardening is merged at `7776a9156c24303c47a57e171a2d5997f6676668` after its strict-JSON, lightweight, and strict-retry exact-head checks passed 208 tests each. These software changes do not repair the known archived-radiation continuation or close full scientific admission.
2. **Archived outer-radiation continuation:** at the archived `N=80` state near `t=47.13`, one smaller step succeeds after a predictor rejection, but the next ordinary continuation approaches the finite-pressure inversion boundary and ultimately fails closed near `t=47.153`. The strict-interior validator/retry contract improves failure classification; it does not repair the trajectory. Do not relax the cone, clip/project conserved variables, or add an arbitrary timestep floor.
3. **Initial centre geometry / CMC series:** the previously required matched-time baseline-versus-combined-candidate comparison has completed at [run 38329810735](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38329810735), artifact [11782600093](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38329810735/artifacts/11782600093). It compared baseline source `9f0bd16d403ada22a4f21413048c93b37c9924a7` against the combined-candidate source `0c58922db287c458dabad8582891afb2fd4d46d0`, vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`, N=80/120/160, `r_max=80`, CFL 0.03, and common coordinate time `t=0.12` (4/6/8 steps). Both arms completed at all three resolutions. Full-grid maxima changed baseline → combined candidate as follows: N=80, H `3.692e-3 → 3.591e-5`, M `5.997e-4 → 1.651e-6`, C_lambda `7.437e-5 → 2.014e-7`; N=120, H `3.387e-3 → 1.047e-5`, M `8.512e-4 → 1.784e-6`, C_lambda `1.083e-4 → 2.318e-7`; N=160, H `2.947e-3 → 6.657e-6`, M `1.054e-3 → 1.955e-6`, C_lambda `1.393e-4 → 2.641e-7`. The centre `Aa_rhs/r0²`, lapse even-series fit and CMC K-target mismatch also improved strongly at these short-time endpoints. Report/ZIP SHA-256 are `0fc22a110dd1981aa299f1014f36d88835dc499260c4623a81d2aa5647d16ca8` and `fd286435b167765b2071921f653d9efe8782ec624b2484f89097374b5f50a2e3` respectively.
   
   Interpretation remains bounded: this supports the **combined** initializer-quadrature + CMC-row change on the matched short trajectory, but does not isolate each component, prove continuum convergence, close fixed-worldtube time-integrated Misner–Sharp balance, or establish archived radiation continuation. The A/B production source predates the latest PR #18 head `be6926d079a4f8c9c919554fedb2a2a92182e987` by four commits; the only numerical-source change between them is extraction of the same origin half-cell weighted dot product into a guarded helper (same weights/arithmetic expression), with an added exact even-quartic regression. Treat the reported trajectory as matched evidence for that formula, while retaining exact latest-head validation and all scientific gates.
4. **Conservation and convergence:** full Misner–Sharp current/flux/work balance at a fixed physical worldtube, fixed-time spatial convergence of H/M/connection/determinant and centre-series witnesses, high-resolution radiation continuation beyond one step, and the long-horizon multi-resolution campaign remain open.
5. **Matter RHS extraction — software refactor now promoted:** [PR #22](https://github.com/miahpenn/0-Cosmos-engine/pull/22) merged at `a4183bdb1b10c9ae6953af2a064af8a25415d001` after reconciliation with PR #23's typed stage-aware radiation retry contract. Exact-head gates passed 204 repository tests in each of direct-RHS, strict-JSON, lightweight/control/production-smoke, and strict-retry integration workflows; the direct-RHS suite also passed 4 focused tests. The prior seven-resolution one-step A/B found only tiny radiation-momentum rounding differences on its pinned source snapshots. This validates the narrow software refactor at the tested level, not long-trajectory equivalence or scientific admission. Do not promote [PR #8](https://github.com/miahpenn/0-Cosmos-engine/pull/8) as-is because its generic radiation flux path predates PR #11; duplicate PR #21 was closed unmerged.

### Latest source-pinned quadrature/CMC factorial (diagnostic only)

The latest four-arm replay [run 38330898524](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330898524), artifact [11783041214](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38330898524/artifacts/11783041214), completed all 12 arms (baseline, quadrature-only, CMC-row-only, combined × N=80/120/160) to the common coordinate time `t=0.12`. Exact pins: baseline code `3e5aedab6c3632d446a8ee4814432d47e14e86fc`; candidate `c2379c80ce0d9dbb598b837253291f336d559cda` (only `engine/v55_initial.py` and `engine/cmc_gauge.py` differ under `engine/`); vendor `d6052d605673ce9d82cdc99b0df79855fa2ea215`; diagnostic head `672bc5d45296ebc432e6338768b1ea2432d21c61`. ZIP SHA-256 `d9a085fed0540906d60d39d1bd074537f8de400204fc2c273a26ce3f42d0b244`, report SHA-256 `f8251a6932849a6a9bcdc806df2f146235f1dfa079df21a2072ad4d395e51c9b`. This physics-source baseline predates PR #30's report-only hardening, but includes the PR #22/#27 numerical integration; PR #30 changed runner/reporting tests only.

At N=160, peak H/M/C_lambda maxima for baseline / quadrature-only / CMC-only / combined were `2.947e-3 / 1.054e-3 / 1.393e-4`, `6.889e-5 / 4.332e-5 / 5.620e-6`, `2.884e-3 / 1.091e-3 / 1.443e-4`, and `6.657e-6 / 1.955e-6 / 2.641e-7`, respectively. The local production-defined Misner–Sharp current residual max/RMS was baseline `1.057e-7 / 4.600e-8`, quadrature-only `1.397e-7 / 1.108e-7`, CMC-only `1.058e-7 / 4.601e-8`, and combined `1.397e-7 / 1.108e-7`. Thus the short-time H/M/C improvement is associated with quadrature (combined is best), while the CMC row mainly improves the lapse/K-target witness; the mass-current residual increase is also associated with quadrature. The fixed-radius sample called `r=10` used nearest cell centres 9.5 (N=80), 10.333… (N=120), and 9.75 (N=160), not a common interpolated physical tube.

No quadrature polynomial-moment defect has been demonstrated. However, the increase in local and integrated mass-current mismatch means the initializer change is **not admitted for production**: current fixed-worldtube conservation, direct semi-discrete `dM_MS/dt`, and fixed-time spatial convergence remain open. The current artifact-first Misner–Sharp lane is investigating existing records; no new evolution run is authorized by this short-time factorial alone.

### Candidate-only diagnostics

- Draft [PR #18](https://github.com/miahpenn/0-Cosmos-engine/pull/18) currently points to candidate branch head `3a845ec1ec903785e8dd3230ecd231ca9221ed44`. The earlier exact-head workflow [38329679022](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38329679022) passed 13 focused tests, 165 repository tests, and the bounded N=40/`t=0.03` smoke at candidate head `be6926d079a4f8c9c919554fedb2a2a92182e987`; do not treat that run as exact-head validation of the later branch tip. The corrected candidate-only matched-time trajectory [38329573898](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38329573898) reached `t=0.12` at N=80/120/160 (4/6/8 steps); H falls with refinement, but M/connection are not uniformly monotone and fixed-worldtube time-integrated Misner–Sharp balance remains open.
- Draft [PR #19](https://github.com/miahpenn/0-Cosmos-engine/pull/19) hardens the multiresolution diagnostic's determinant and radiation witnesses. Following rebase to the latest candidate initializer, exact-head workflow [38329976316](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38329976316) passed the five-resolution one-step harness contract. This remains candidate-only validation and does not promote either PR #18 or #19 to production.
- The matched default-initializer baseline-versus-combined-candidate A/B is complete at [run 38329810735](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38329810735), with its report/artifact and source/vendor pins recorded above. It supports the combined quadrature + CMC-row candidate's short-horizon improvement, but does not isolate the two components, prove continuum convergence, close fixed-worldtube Misner–Sharp balance, or repair archived radiation continuation. Keep PR #18/#19 isolated from `main` until these independent admission gates close.
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
