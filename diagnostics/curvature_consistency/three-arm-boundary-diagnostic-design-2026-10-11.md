# Three-arm boundary diagnostic design and admission gates — 2026-10-11

## Decision

Hold all new full trajectories. The admitted evidence supports a domain-dependent CMC lapse as a candidate contributor to the R=80/R=160 divergence, but the current boundary override would also change matter boundary closure inputs. No causal attribution is admitted yet. This plan is diagnostic-only; production equations, gauge, initial data, and physical interpretation remain unchanged.

## Evidence already admitted

- R=80 and R=160 evolved trajectories begin with shared physical initial arrays equal at stored precision, but the initial centre lapse differs by about 5.669e-4.
- After one step, the centre physical-state differences remain small (largest reported among listed fields: X about 1.205e-8; PS about 2.383e-8), which identifies a candidate early seed but does not prove causality.
- At t=12, independent-trajectory lapse difference reaches 0.116139 at r=39.75.
- In the frozen-state crop audit, matching the full-domain CMC operator reproduces the crop solution to near numerical precision; this does not explain the larger independently evolved mismatch.
- The R=160 reference lapse at the shared R=80 outer-cell centre r=79.75 is domain-normalized by the R=160 outer boundary. It is about 0.9896 at t=12, versus the R=80 CMC Dirichlet value 1.0. Treat this as a domain-relative reference, not an absolute lapse.
- Linear interpolation holdout over 0.5 time units produced max absolute error 1.48936e-4 and RMS 5.62025e-5 (24 points). Scaling the interval from 0.5 to 0.25 suggests a max error around 3.7e-5 only under a smooth quadratic-error assumption; that is an estimate, not a measurement or stage-level bound. The 49-slice record cannot establish accurate boundary values at internal stages dt=0.00375.
- The one-step CI workflow passed: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38153562804.

## Required diagnostic arms

All arms must begin from the same pinned baseline commit, same initial physical state, timestep/CFL, solver tolerances, and output checkpoints. Each arm is a counterfactual diagnostic, never a production physics change.

1. **Baseline control**: native R=80 CMC normalization and native matter boundary closure. Reproduce the existing admitted trajectory first if a new code path could affect baseline execution.
2. **Joint intervention**: at every CMC solve in every integrator stage, use the time-dependent R=160 lapse reference at R=80's outer cell centre r=79.75 as the CMC Dirichlet value. Let the resulting boundary metric feed through to the native radiation/dust closure. This measures only the combined intervention (CMC boundary value + changed closure inputs).
3. **CMC-only separated intervention**: use the R=160 reference in the CMC solve, but retain baseline R=80 outer-face metric inputs for matter boundary characteristic speeds and flux closure. In particular, audit reads of geom.alpha[-1], geom.beta[-1], and gamma_rr at the outer face. Freeze only the specifically identified boundary closure inputs; do not alter interior matter equations or source terms. Label this a split diagnostic counterfactual, not a physically self-consistent production evolution.
4. **Radiation-closure-only intervention**: retain native R=80 CMC normalization, but apply the explicitly defined reference/baseline-separated boundary radiation closure inputs. Do not quietly change dust closure too; if dust also needs a corresponding intervention, report it as a separate factor or expand the design transparently. This arm should isolate the already-open radiation boundary concern rather than silently bundle radiation and dust.

The three requested intervention arms are joint, CMC-only, and radiation-only, but an unchanged baseline control is mandatory for comparison. If inspection shows that dust closure also changes under the joint intervention, record it as a separate confound and do not describe a radiation-only arm as a complete decomposition of all boundary effects.

## Dense reference history

Preferred path: replay the pinned R=160 reference trajectory once with the same admitted settings and record the solved lapse at index 159 (r=79.75) at every accepted step and at every CMC/integrator stage where the boundary value is consumed. Include timestamps, stage identifiers, values, accepted/rejected status, and a content hash. Do not interpolate across rejected stages or assume accepted-step output equals all internal stage values.

Before replay, verify the R=160 replay reproduces the admitted reference trajectory at existing saved checkpoints within declared tolerances. If exact stage-level capture is not supported by the integrator, add a diagnostic recorder only; demonstrate it does not alter solver state or step acceptance. Preserve the baseline path unchanged.

## Additional reference-domain check

The R=160 value at r=79.75 is itself domain-dependent. Before treating it as a sufficiently stable reference, compare the same-state solved CMC lapse at r=79.75 using a larger-domain solve, holding the evolved state over the shared region fixed and documenting each outer normalization. This is a frozen-state sensitivity check, not another full trajectory. Report the difference as reference-domain uncertainty; do not call the R=160 profile domain-independent.

## Mandatory pre-trajectory validation

1. **Initializer audit**: show that R=80 and R=160 initialization algorithms differ only by domain size/grid length for the shared region; compare arrays, lapse solve inputs, boundary conditions, and solver settings.
2. **Call-site audit**: enumerate every CMC lapse solve and every consumer of outer alpha/beta/gamma in matter fluxes/characteristic speeds. Prove the selected override reaches every intended stage and no unintended call site.
3. **Baseline identity test**: all diagnostic flags off must reproduce the pinned baseline outputs/hashes or declared floating-point tolerances exactly.
4. **CMC override unit test**: with a fixed state, verify the requested boundary value is imposed at every solve; verify interior operator rows and source terms are unchanged; verify the native path remains byte-for-byte or numerically identical when disabled.
5. **Closure separation tests**: for CMC-only, prove closure inputs equal the recorded baseline boundary inputs while the CMC solve uses the reference. For radiation-only, prove CMC normalization remains native and only the declared radiation boundary inputs change. Assert no interior cells or matter source terms are modified.
6. **Reference-data tests**: monotonic timestamps; correct index 159/r=79.75 mapping; no duplicate/missing stages; matching reference commit/configuration; SHA-256 recorded; fail closed on missing timestamps or incompatible grids.
7. **Short diagnostic smoke test**: tiny grid/short duration only, to validate flags, recorded inputs, and state hashes. No expensive full trajectory until all checks pass.
8. **Admission and comparison**: report baseline vs each arm at identical checkpoints: alpha(r), a, b, D, PD, scalar fields, radiation/dust boundary fluxes, constraints/residuals, accepted/rejected steps, and tau accumulator (without relabelling it proper time). Compare both absolute differences and spatial locations; no front/bounce/turnaround inference from thresholds alone.

## Interpretation rules

- If joint and CMC-only suppress late divergence while radiation-only does not, that supports CMC normalization as the dominant boundary lever, subject to the split-counterfactual limitation.
- If only joint suppresses divergence, the effect is interaction-dependent; do not attribute it to CMC alone.
- If radiation-only suppresses divergence, the radiation closure is implicated; distinguish the known prior radiation boundary concern from a new defect.
- If none suppresses divergence, boundary treatment is not sufficient to explain it; investigate interior evolution and initialization history next.
- Early-seed claims require the dense history and are not supported by the current coarse interpolation record.
- None of these outcomes alone establishes a physical front, bounce, turnaround, or production-gauge defect.

## Prior reported radiation issue

An earlier diagnostic was reported to show a 19.7% outer-flux mismatch. Revalidate its artifact, exact definition, sign, denominator, and configuration before using the percentage as admitted evidence in this experiment. Do not assume that number directly predicts the result of the new intervention.

## Current status

- Full R=80/R=160 domain comparison: admitted.
- Frozen-state crop audit: admitted.
- One-step domain comparison CI: passed.
- Interpolation holdout: admitted as a coarse empirical check only.
- Dense stage-level reference history: not yet produced.
- Larger-domain frozen-state reference check: not yet run.
- Three-arm code paths and unit tests: not yet implemented/validated.
- Full boundary-matched trajectory: **not launched**.
- No production physics, gauge, or matter equation changes.
