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

## 4. The homogeneous COSMOS lane is a control, not a second live state

The homogeneous implementation in `engine/cosmos.py` preserves the archive equations, including the signed-H homogeneous RHS, the scalar/DM beta exchange, and the dilution laws for dust, radiation and effective shear. It rejects unapproved nonzero `exchange_Q`.

The later `GEAR_V5_5_STRONGFIELD_CLOSURE_AUDIT_2026-10-05.md`, Sections 10–11, selects a V5.5 target with **one spatial metric and one total stress-energy tensor**. The intended live cosmological scalar is `ProductionState.scalars.phi(r,t)` and shares that geometry with the local S/D fields and radial matter species. The standalone `CosmosState` is not expected to be a second state evolved in the same production step.

The existing `geometry_driven_step()` consumes a supplied `H_eff` and advances the homogeneous helper in one direction. It is useful only where an explicitly driven homogeneous comparison is intended; it is not the unifying mechanism for the selected one-spacetime target. The absent production call is therefore not by itself a defect. The required test is that, with local perturbations/spatial gradients suppressed and parameters, units and initial data matched, the live radial scalar/matter evolution reproduces the homogeneous COSMOS equations to the order expected of the numerical method.

The source supports: (1) a separate homogeneous reference/control; (2) a live spatial evolution with local `phi(r,t)` and conservative matter; and (3) a one-spacetime V5.5 design target. It does **not** yet establish full campaign validation, homogeneous-limit agreement, all strong-field continuation, global shear ownership, D-to-matter conversion, or a completed 0* gate.

## 5. Archive gates that remain open

| Gate | Required evidence before promotion | Prohibited shortcut |
|---|---|---|
| Local proper time versus cosmic time | A parameter-free clock relation or an explicit single-spacetime time definition with convergence evidence | Arbitrary time rescaling or equating tau and t by convenience |
| D-mode energy to localized ordinary matter | A derived conversion/identification law with a ledger that closes | An invented D-to-matter coefficient or classifying D flux as matter creation |
| Local source/current to global exchange Q | Derived sequence from physical work/power to event energy/rate and then Q, plus independent conservation and target checks | Setting Q from instantaneous force or fitting a normalization |
| Spatial phi versus homogeneous phi | Homogeneous-limit regression with matched units, parameters and initial data | Injecting CosmosState as a second live copy instead of testing the spatial equations against the control |
| Shear and global anisotropic curvature | A clear owner and matching/curvature equation appropriate to the chosen model | Adding homogeneous shear pressure to the local spherical stress tensor without derivation |
| Cycle and turnaround interpretation | Event rows and neighboring states showing the same trajectory, signs, coordinate/proper times, expansion measure and constraints | Calling one H_eff crossing a completed cycle or bounce |

The archive and prior GEAR work explicitly distinguish a diagnostic bridge from a promoted physical source. Those labels must be preserved.

## 6. Architecture decision from the later V5.5 closure audit

The original 2 October harness remains the canonical two-lane reference and must not be silently rewritten. For the later V5.5 unified production target, the 5 October Strong-Field and Covariant-Closure Audit, Sections 10–11, already selects **Candidate A: one solved spatial spacetime**:

- local S/D, the COSMOS scalar, dark matter, baryons and radiation share one radial geometry and one total stress-energy accounting system;
- the standalone `CosmosState` is a control for homogeneous-limit regression, not a second live scalar state;
- the exact spherical branch does not emulate the archive's Bianchi-I shear as an isotropic fluid;
- the covariant DM exchange law is given in Sections 8–9 of that audit; its source terms still require implementation-level conservation checks.

Candidate B (homogeneous background plus local region with reciprocal interface) describes a different architecture and would need its own derived interface contract. It must not be created by calling the current one-way geometry-driven helper or by adding `CosmosState` to the production state.

The next decision is therefore not whether to wire two state objects. It is whether the active V5.5 spatial equations pass the locked equation/source accounting, constraints, and homogeneous-limit tests needed for admission.

## 7. Event language and the two-turnaround observation

The active CycleLedger records a turnaround when consecutive H_eff samples satisfy H_before > 0 and H_after <= 0, and a re-expansion crossing when H_before < 0 and H_after >= 0. A second positive-to-nonpositive crossing in the same continuous sample sequence entails an intervening recorded negative-to-nonnegative crossing, unless the rows originate from different trajectories or reset/duplicate data.

That still does not, by itself, prove a scale-factor turnaround, a physical contraction interval, or a completed cycle. Admission of the user's remembered two-turnaround observation requires the exact originating run/artifact and the raw event rows, neighboring checkpoints, same-run configuration, coordinate t and named-worldline tau, H/H_eff signs, and the actual expansion measure. The previously admitted D-mode ledgers have zero turnaround/re-expansion events; the earlier True-CMC post-turnaround adjudication records one crossing and no re-expansion before failure. Neither has been identified as the user's two-event source.

## 8. Gate order

1. Complete and admit the currently running read-only P-minus boundary audit. Separate diagnostic success from physical trajectory success or failure.
2. Preserve the original independent-lane archive; follow the later V5.5 single-spacetime target; keep the homogeneous CosmosState as a control and verify the spatial homogeneous limit rather than wiring a duplicate state.
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

**Draft conclusion:** V5.5's target architecture is a single radial spacetime, while the original 2 October two-lane engine remains the preserved archive reference. The standalone homogeneous COSMOS implementation should be used as a control, and its equations must be recovered in a homogeneous-limit regression of the spatial kernel. A missing call to the helper does not itself establish a broken interface. The source-accounting witness found a lapse-weighting mismatch in the spatial beta-exchange implementation; a minimal correction is being tested on an isolated physics branch. The radiation P-minus failure remains a separate diagnostic stream.
