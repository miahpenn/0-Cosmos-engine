# Matter RHS Rate-Extraction Correction — SI Handoff

**Date:** 12 October 2026 (UTC)  
**Status:** Isolated engineering correction; bounded diagnostic re-run passed. Not a production admission or merge approval.  
**Branch:** `fix/matter-rhs-generator-2026-10-12`  
**Audit workflow:** [Run 38307235876](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38307235876)  
**Audit commit:** `149273559da17d5966ae196a38a10348a8263818`  
**Engine correction source commit:** `d7f964657df2cb8aabe4bd30c04eaff9e65d635d`  
**Geometry/PIRK vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`

## Executive result

The archived CMC-row/radiation/Misner–Sharp A/B audit previously captured an apparent first-cell radiation-cone violation at N=240 and N=320, before the first predictor CMC solve. Code review located a false-failure path: the matter RHS wrapper implemented a rate query by advancing the conservative state with an artificial `dt=1.0`, then subtracting the input state. The lower-level `evolve_species` routine checked radiation admissibility on that artificial finite-time result. Thus the captured failure was not the production integrator's requested update (dt = 0.01 at N=240; dt = 0.0075 at N=320).

An isolated correction now calculates the conservative rate directly from the existing fluxes and source terms. Finite-time radiation checks remain enabled on real updates; primitive recovery still validates incoming radiation states. The original matched audit has been re-run. Both historical and candidate CMC-row arms completed one real step at all seven resolutions (N=40, 60, 80, 120, 160, 240, 320); all recorded radiation states were inside the cone at the four captured stages.

This confirms the earlier N=240/320 *one-step RHS-wrapper* failure was an artefact of the artificial unit step. It does **not** close the separate long-run outer-cell radiation predictor failure, centre-regularity admission, Misner–Sharp conservation, constraint propagation, or long-horizon campaign gates.

## Root cause and data path

The previous `engine/matter_rhs.py::species_rhs` called `evolve_species(..., 1.0)`, then returned `advanced - state`. In `engine/matter_system.py::evolve_species`, the argument controls a finite-volume state update. After that update the radiation branch verifies $E \geq 0$ and $|S| \leq E$ and raises if the updated state is outside the cone.

A unit-step update is not mathematically identical in floating-point arithmetic to extracting a differential rate: it both advances the state by a physical-sized interval unrelated to the integrator's `dt` and subjects that artificial state to a finite-step admissibility test. At high resolution, the tiny first-cell cone margin was crossed by that synthetic update even though the actual small-`dt` step retained a positive margin.

This distinction is independently visible in the archived data. Before the correction, N=240 reported `|S|/E=1.0202778` in the wrapper-created state and N=320 reported `1.3602899`; these were values after the unit-time surrogate. In the corrected matched run, the first explicit predictor's maximum ratio is about (0.01014) at those same resolutions and the minimum cone margins remain positive.

## Change made (no physical-law change)

1. `engine/matter_system.py` now has `_species_rhs_components(...)`, which reconstructs the input primitives, computes the existing HLL face fluxes and Valencia/DM source terms, and returns the conservative differential rates plus separate diagnostic flux/source arrays. It does not choose a timestep or test a fictitious updated state.
2. `evolve_species(..., dt)` now applies `dt) to those rate components and evaluates the radiation cone on the resulting requested finite-time state. This keeps validation at the actual state update.
3. `engine/matter_rhs.py::species_rhs` now calls the rate-only component helper. It no longer calls `evolve_species` with an arbitrary unit step.
4. A regression test, `test_species_rhs_does_not_emulate_a_unit_time_update`, forbids the old rate-wrapper call path.

The change does not alter the stress-energy law, the radiation EOS, the dark-matter exchange law, the CFL, timestep selection, boundary condition, CMC target, physical parameters, floors, clamps, or damping. The existing pinned geometry/PIRK vendor is unchanged.

## Validation and provenance

Workflow [38307235876](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38307235876) completed successfully on the isolated branch.

- Python compilation: passed.
- Focused audit-helper tests: **3 passed** in 0.28 s (the pre-existing two plus the new rate-wrapper regression).
- Full repository suite: **92 passed** in 7.14 s.
- Matched CMC A/B diagnostic: **7/7 resolutions completed both arms; zero partial captures and zero setup failures**.
- Source guard: engine/vendor/dependencies were pinned and unchanged after source commit `d7f9646...`; vendor submodule pin verified by the workflow.

Artifact provenance:
- Artifact ID: `11773485002`, name `cmc-lapse-radiation-ms-one-step-ab`.
- Downloaded archive SHA-256: `787bdc721c3345d19bde8290b73f0329775c266a5965119beed53fcdd23dcd63`.
- Internal `report.json` SHA-256: `e9e2dbd55ebe4cc7895f3fa6c7338e8221bc4e6f95a4f1a1b5e74dbe3bfb6e98`.
- Machine-readable `report.json` records `schema=cmc_lapse_radiation_ms_one_step_ab_v1`, `status=completed_capture`, `diagnostic_only=true`, audit commit `1492735...`, and candidate engine source `d7f9646...`.

### Actual-step radiation-cone witnesses

The following table gives the maximum full-domain $|S|/E$ recorded at the first explicit predictor (the maximum among the captured pre-solve stages) and the minimum cone margin over those stages. The cone margin is reported in the audit's densitized conservative variables $C=U_E-\sqrt{X^2/a}\,|U_r|$. All old/candidate states, including the completed step, passed the audit's positivity and cone checks.

| N | Actual dt | Max (|S|/E), old row | Max (|S|/E), candidate | Min cone margin, old | Min cone margin, candidate |
|---:|---:|---:|---:|---:|---:|
| 40 | 0.0600 | 1.017536e-2 | 1.012768e-2 | 8.778783e-9 | 8.779200e-9 |
| 60 | 0.0400 | 1.016172e-2 | 1.013901e-2 | 3.901587e-9 | 3.901675e-9 |
| 80 | 0.0300 | 1.015426e-2 | 1.014115e-2 | 2.194688e-9 | 2.194717e-9 |
| 120 | 0.0200 | 1.014812e-2 | 1.014217e-2 | 9.754500e-10 | 9.754558e-10 |
| 160 | 0.0150 | 1.014581e-2 | 1.014244e-2 | 5.487021e-10 | 5.487039e-10 |
| 240 | 0.0100 | 1.014412e-2 | 1.014261e-2 | 2.438732e-10 | 2.438736e-10 |
| 320 | 0.0075 | 1.014352e-2 | 1.014267e-2 | 1.371804e-10 | 1.371805e-10 |

At N=240 and N=320 both arms completed the actual update; neither generated an exception, and no recorded cell exceeded (|S|/E=1). The candidate/old cone differences are small. This is a one-step diagnostic, not proof of radiation admissibility over a full physical trajectory.

### Misner–Sharp balance remains open

The script evaluates the one-step residual
[
$\Delta M_{MS} - \Delta t\,\frac{(F_T+W_{pR})_0+(F_T+W_{pR})_1}{2}$.
]
No preregistered acceptance tolerance exists for this witness, so the following normalized residuals are descriptive—not pass/fail results.

| N | Scaled residual, old row | Scaled residual, candidate |
|---:|---:|---:|
| 40 | -5.958767e-2 | -5.900167e-2 |
| 60 | -1.950282e-2 | -1.948668e-2 |
| 80 | -6.715467e-3 | -6.713182e-3 |
| 120 | -9.893546e-4 | -9.898456e-4 |
| 160 | -3.866915e-4 | -3.867917e-4 |
| 240 | -5.647967e-5 | -5.647726e-5 |
| 320 | -1.684872e-5 | -1.686494e-5 |

The magnitude trends down with resolution in this one-step comparison, while the CMC row switch makes a very small difference. Without preregistered tolerances and a multi-step interval/convergence study, this is not a conservation pass.

## Open gates and safe interpretation

- The parity-aware CMC first-cell PDE-row candidate remains isolated. Earlier matched one-step evidence showed strong M and (C_\Lambda) improvement for N=60–160, but mixed H response. At N=240/320 the candidate centre (\dot A_a/r_0^2) remains about (-6.38\times10^{-5}) and (-6.68\times10^{-5}), respectively—not a demonstrated (\dot A_a=O(r^2)) centre-series limit. Do not promote the CMC row based on this run.
- The full-machine review still classifies the machine as implemented / unit- and smoke-tested / not full-campaign validated. This run advances the code-rate/admissibility diagnostic only.
- The distinct long-run radiation predictor crossing near (t\approx46.98), (\tau\approx22.0161), remains a separate issue and must be revisited with its archived predictor-stage ledger. The current one-step result does not erase it.
- Keep the canonical GEAR archive and REV16 status intact. Do not invent a global Q, D-to-matter law, clock conversion, bounce, reset, or radiation limiter.

## Recommended handoff

1. Have an independent SI agent review the diff for algebraic equivalence of the extracted conservative rate and the original flux/source RHS, including dark-matter reciprocal source terms and floating-point ordering.
2. Keep this result and correction isolated; do not merge the CMC candidate or any related diagnostic branch to `main`.
3. Use the newly successful N=240/320 actual-step captures to replace the erroneous wrapper-failure interpretation in subsequent reports.
4. Return to the separate long-run radiation predictor failure and centre-series/constraint-tangent gates, using existing artifacts before requesting any new expensive trajectory.
5. Only after the relevant invariant witnesses close together should the exact candidate commit be frozen for a preregistered multi-resolution campaign.

**Conclusion:** a real rate extraction must return the semi-discrete RHS, not simulate a one-unit-time update. That implementation error created the earlier high-resolution one-step cone report. The corrected code passes the expanded regression suite and the seven-resolution *one-step* diagnostic, while full physical and numerical admission remains open.
