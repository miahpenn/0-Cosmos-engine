# Independent Confirmation — Radiation RHS Unit-Step Guard

**Date:** 2026-10-12 UTC  
**Status:** Diagnostic-only; no production merge recommendation.  
**Branch:** diag/radiation-unit-step-extraction-audit-2026-10-12  
**Audit commit:** 9dfd88065dfe0219feae9b0e55e9b54f34294976  
**Workflow:** https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38307129998  
**Artifact:** https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38307129998/artifacts/11772998975  
**Artifact ZIP SHA-256:** 2dd08350f0dc0ccee19478b727f2c93d78dc2312d5c2cc9c85696a822637d42d  
**Internal report.json SHA-256:** 09943b1b18a0ec66878b72f329d743d73167760c087f6eb77daa4505b0cb0a2f  
**Frozen production source:** 6d8b3f2790644591326df57c55908863f994cb0c  
**CMC candidate:** 9f7ebf1aab2076fe8925737591366b416427305e  
**Pinned vendor:** d6052d605673ce9d82cdc99b0df79855fa2ea215

## Finding

The prior N=240/320 radiation failure is reproduced by the production species_rhs wrapper but not by the same conservative matter update at the production CFL increment. The pinned wrapper called evolve_species(dt=1.0) to estimate a differential rate, then subtracted the input state. The radiation-cone guard inside evolve_species was evaluated after this artificial finite update. The guard therefore rejected a unit-time surrogate that the integrator never intended to take.

The unchanged admissibility inequality is E >= |S|. This audit neither clips nor projects the state and does not relax the cone check.

## Method

The diagnostic reused the existing uniform-fourth-order initial-profile builder and parity-aware CMC candidate, pinned matter update and source equations, and the production timestep rule dt=0.03*dr. At each frozen state it compared (1) the production species_rhs wrapper, (2) evolve_species at actual production dt, and (3) evolve_species at dt/10. It also compared rate quotients from the last two increments. No coupled geometry/matter trajectory was advanced. The fail-closed source/vendor check passed; the complete repository regression suite passed 91 tests in 10.67 s.

## Results

The earlier combined audit reports its minimum cone margin in densitized conserved variables, C = U_E - sqrt(X^2/a)|U_r|. This report uses the undensitized physical margin E - |S|. Those margins are not directly comparable as magnitudes; multiplying the physical margin by local sqrt(gamma) recovers the earlier values. Dimensionless |S|/E ratios and exception values match.

| N | Production dt | Initial max |S|/E | species_rhs | Direct unit update | Actual-dt update | Actual-dt max |S|/E |
|---:|---:|---:|---|---|---|---|---:|
| 160 | 0.0150 | 0 | returned | returned | admitted | 0.0101433191 |
| 240 | 0.0100 | 0 | ValueError | ValueError | admitted | 0.0101431965 |
| 320 | 0.0075 | 0 | ValueError | ValueError | admitted | 0.0101431053 |

The unit-increment exceptions reproduce the earlier combined artifact's candidate-arm ratios:
- N=240: |S|/E = 1.0201269083.
- N=320: |S|/E = 1.3601765554.

After the production-sized increment, all cells remain admissible. The minimum undensitized physical cone margins are 8.77939356e-9 at N=240 and 8.77952058e-9 at N=320; the corresponding densitized margins recorded in the original audit are 2.46386633e-10 and 1.38592128e-10. They agree after accounting for the volume factor.

### RHS quotient agreement

| N | Energy rate max relative difference, dt vs dt/10 | Momentum rate max relative difference |
|---:|---:|---:|
| 160 | 7.06e-12 | 4.08e-16 |
| 240 | 1.06e-11 | 5.24e-16 |
| 320 | 1.41e-11 | 4.08e-16 |

These are equivalence witnesses on a frozen input, not a physical-law proof.

## Cross-agent reconciliation

The existing SI-agent branch fix/matter-rhs-generator-2026-10-12 already implements the corresponding isolated correction in commit d7f964657df2cb8aabe4bd30c04eaff9e65d635d. Its matched A/B workflow, run 38307235876, passed on that branch; the agent's report is:

https://github.com/miahpenn/0-Cosmos-engine/blob/fix/matter-rhs-generator-2026-10-12/docs/MATTER_RHS_RATE_EXTRACTION_CORRECTION_2026-10-12.md

The correction computes a rate-only _species_rhs_components() result; evolve_species(dt) continues to check the actual requested finite-time update, and primitive recovery continues to check incoming radiation states. The corrected matched audit completed both CMC-row arms at N=40/60/80/120/160/240/320 with all four recorded radiation stages admitted. The expanded full suite passed 92 tests on the isolated branch.

This audit is independent corroboration, not a competing engine patch.

## Limits and disposition

This resolves the specific high-resolution first-step wrapper failure. It does not resolve the separate long-run outer-cell radiation predictor crossing near t≈46.98, proper time τ≈22.0161; the regular-centre condition A_a_dot = O(r^2); the nonzero momentum-constraint tangent from a discrete-consistent initial state; or the Misner-Sharp balance tolerance/convergence gate. The CMC row candidate also remains unmerged because its matched constraint response is mixed for H and does not pass centre regularity.

Do not infer a physical-source failure from the superseded unit-increment exception. Do not relax the radiation cone, introduce floors/clamps, change a physical source law, or start a long horizon before the remaining gates and their existing stage artifacts are reconciled.

**Disposition:** Independent confirmation complete. Use the existing isolated rate-extraction correction for further review. Keep main, the canonical GEAR archive, and draft PR #7 unchanged.
