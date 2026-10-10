# GEAR-to-kernel state ownership contract — draft

**Status:** source audit / draft only; not a physics or architecture approval  
**Date:** 2026-10-10  
**Repository:** miahpenn/0-Cosmos-engine  
**Inspected source commit:** d407ff2ff3197f435346d05c4112116f16725698  
**Authoritative reference:** the Library copy of CYCLIC-COSMOS ENGINE-.docx and the archived GEAR/REV16 records. The archive itself is not modified by this draft.

## 1. Purpose and non-goals

This document distinguishes what the checked-in production kernel actively evolves from what exists only as a separate homogeneous helper or legacy coupled-state container. It is intended to stop an unproved field identification, source term, or clock conversion from being inserted merely to make the machine appear fully coupled.

This is not a claim that the source is fully validated. The V5.5 kernel class explicitly labels itself as an unvalidated production-facing implementation and sets its validation state to implemented_not_campaign_validated.

No production equations, acceptance criteria, cycle rules, radiation boundaries, or archive files are changed by this document.

## 2. Source evidence at the inspected commit

- [Production state and V5.5 step path](https://github.com/miahpenn/0-Cosmos-engine/blob/d407ff2ff3197f435346d05c4112116f16725698/engine/production_kernel.py)
- [Archive-locked homogeneous COSMOS lane](https://github.com/miahpenn/0-Cosmos-engine/blob/d407ff2ff3197f435346d05c4112116f16725698/engine/cosmos.py)
- [Legacy coupled-state container](https://github.com/miahpenn/0-Cosmos-engine/blob/d407ff2ff3197f435346d05c4112116f16725698/engine/coupled.py)
- [Cycle and handoff ledger](https://github.com/miahpenn/0-Cosmos-engine/blob/d407ff2ff3197f435346d05c4112116f16725698/engine/handoff.py)
- [Production campaign / output writer](https://github.com/miahpenn/0-Cosmos-engine/blob/d407ff2ff3197f435346d05c4112116f16725698/engine/run_production.py)
- [Cycle coordinate decision](https://github.com/miahpenn/0-Cosmos-engine/blob/d407ff2ff3197f435346d05c4112116f16725698/docs/CYCLE_COORDINATES.md)
- [Post-turnaround failure adjudication](https://github.com/miahpenn/0-Cosmos-engine/blob/d407ff2ff3197f435346d05c4112116f16725698/docs/POST_TURNAROUND_FAILURE_2026-10-07.md)

The authoritative archive names local S/D modes and a COSMOS sector with scalar, dark matter, baryons, radiation, shear, expansion/contraction, and separate local/cosmic clocks. It deliberately leaves several derivations open. The existence of code fields or helper functions does not close those gates.

## 3. State ownership as implemented

| State owner | Variables / role | Source-level status | Limitation |
|---|---|---|---|
| ProductionState | Radial grid and geometry; scalar-field bundle; matter bundle; coordinate time t; accumulated central proper time tau; e-fold diagnostic; cycle/history/handoff ledgers | Actively advanced by V55ProductionKernel.step() | This is the live spatial evolution path. It does not contain a CosmosState field. |
| Spatial geometry and scalar/matter bundles | Shared radial spacetime; scalar and species RHS, metric/lapse, total stress-energy, radiation fluxes and existing outer-light-boundary call | Active mutual coupling inside the production kernel | This establishes local spatial coupling, not a completed global cyclic-cosmos model or a D-to-baryonic-matter conversion law. |
| CosmosState | Homogeneous a, H, phi, pi_phi, rho_dm, rho_b, rho_r, rho_shear | Separate homogeneous state defined in engine/cosmos.py | It is not advanced from inside the current ProductionState.step() call path. |
| CoupledState | A CosmosState, a LocalState, and an InterfaceState | Legacy/top-level container exists in engine/coupled.py | A container and pack() serializer do not establish an active, conserved coupled evolution. |
| CycleLedger / HandoffPoint | H_eff crossing records, t/tau, local phases and geometric/flux observables | Event and diagnostic bookkeeping on the production path | These records are observations, not additional field equations or a physical bounce rule. |

## 4. The COSMOS bridge: implemented helper, not a closed reciprocal source

The homogeneous lane implements the archive-locked equations, including:

- da/dt = H a;
- dH/dt = -0.5 (rho + pressure), with H evolved in the standalone homogeneous RHS;
- scalar evolution with the beta-coupled dark-matter source;
- dilution terms proportional to 3H for dust-like components, 4H for radiation, and 6H for the effective shear density.

In engine/cosmos.py, cosmos_rhs() explicitly rejects nonzero exchange_Q. The code does not silently treat an arbitrary phenomenological Q as a physical source.

The geometry-driven helper consumes already-solved H-driver values and advances the non-H homogeneous state. It is one-way: no homogeneous stress-energy or derived Q is returned to the radial production geometry. The source review found no call to geometry_driven_step() in the current V55ProductionKernel.step() or run_production.py path, and ProductionState has no homogeneous CosmosState member.

Therefore these three propositions must remain distinct:

1. A geometry-driven COSMOS helper exists.
2. The live production kernel derives and records geometric H_eff observables.
3. A fully reciprocal, energy-conserving production-to-COSMOS-to-production coupling has been derived and is evolved.

The source supports the first two; it does not establish the third.

## 5. Archive gates that remain open

| Gate | Required evidence before promotion | Prohibited shortcut |
|---|---|---|
| Local proper time versus cosmic time | A parameter-free clock relation or an explicit single-spacetime time definition with convergence evidence | Arbitrary time rescaling or equating tau and t by convenience |
| D-mode energy to localized ordinary matter | A derived conversion/identification law with a ledger that closes | An invented D-to-matter coefficient or classifying D flux as matter creation |
| Local source/current to global exchange Q | Derived sequence from physical work/power to event energy/rate and then Q, plus independent conservation and target checks | Setting Q from instantaneous force or fitting a normalization |
| Local phi versus homogeneous phi | A field-ownership/matching rule, including which state is authoritative | Evolving duplicate copies of the same field without a matching equation |
| Shear and global anisotropic curvature | A clear owner and matching/curvature equation appropriate to the chosen model | Adding homogeneous shear pressure to the local spherical stress tensor without derivation |
| Cycle and turnaround interpretation | Event rows and neighboring states showing the same trajectory, signs, coordinate/proper times, expansion measure and constraints | Calling one H_eff crossing a completed cycle or bounce |

The archive and prior GEAR work explicitly distinguish a diagnostic bridge from a promoted physical source. Those labels must be preserved.

## 6. Architecture decision still required

Two representations are conceptually possible; this draft does not silently select one.

### Candidate A — one solved spatial spacetime

Local S/D and matter fields evolve on one metric and contribute to the shared stress-energy. Cosmological expansion and its cycle history must emerge from the solved solution or from a defined asymptotic/averaged observable. The homogeneous CosmosState code would need a precisely declared role—such as a reference-limit regression, a projection of the spatial solution, or a non-production comparison. It must not be integrated as a second physical copy of a field without an ownership rule.

### Candidate B — homogeneous background plus local region

The homogeneous COSMOS and local spatial states remain distinct. Promotion requires explicit matching/boundary data, a derived local-to-global clock rule, and conserved source/current equations in both directions. The existing geometry-driven helper alone is not enough because it does not close the reciprocal source or exchange_Q.

The decision must come from the archive's governing equations and the intended physical representation, not from whichever code path is easiest to connect.

## 7. Event language and the two-turnaround observation

The active CycleLedger records a turnaround when consecutive H_eff samples satisfy H_before > 0 and H_after <= 0, and a re-expansion crossing when H_before < 0 and H_after >= 0. A second positive-to-nonpositive crossing in the same continuous sample sequence entails an intervening recorded negative-to-nonnegative crossing, unless the rows originate from different trajectories or reset/duplicate data.

That still does not, by itself, prove a scale-factor turnaround, a physical contraction interval, or a completed cycle. Admission of the user's remembered two-turnaround observation requires the exact originating run/artifact and the raw event rows, neighboring checkpoints, same-run configuration, coordinate t and named-worldline tau, H/H_eff signs, and the actual expansion measure. The previously admitted D-mode ledgers have zero turnaround/re-expansion events; the earlier True-CMC post-turnaround adjudication records one crossing and no re-expansion before failure. Neither has been identified as the user's two-event source.

## 8. Gate order

1. Complete and admit the currently running read-only P-minus boundary audit. Separate diagnostic success from physical trajectory success or failure.
2. Decide the state-ownership architecture from the archive before wiring or promoting the homogeneous lane.
3. For every coupling edge, name its exact state owner, equation, inputs, output/recipient, time level, clock, conservation identity, constraints, and validation evidence.
4. Add only deterministic unit tests and diagnostics for that contract before any new long integration.
5. Promote a coupling only when its source equation and conservation/constraint tests are derived and pass across resolution. Keep the exact source commit, logs, artifacts, and SHA-256 digests for every numerical claim.

## 9. Standing constraints

- No artificial radiation floors, conservative-variable clamps, or fitted damping.
- No manufactured bounce, branch reset, or hand-forced sign change.
- No fitted or arbitrary local-to-global time scaling.
- No invented matter-identification or D-to-matter law.
- No arbitrary exchange_Q; no silent changes to the REV16 criterion.
- No edits to the authoritative archive and no merge to main as part of this draft.

**Draft conclusion:** the local spatial kernel contains genuine matter–geometry coupling, but the inspected source does not demonstrate that the archive's homogeneous COSMOS lane is part of the active production evolution or that a reciprocal physical Q/clock bridge is closed. The radiation P-minus audit and the archive-to-kernel state-ownership decision are separate gates and must be resolved separately.
