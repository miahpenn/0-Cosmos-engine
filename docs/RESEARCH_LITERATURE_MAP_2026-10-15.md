# COSMOS scientific literature map — 15 October 2026

**Disposition: research aid only. Production scientific admission remains NOT ADMITTED.**

This is a focused, source-linked reading map for the existing COSMOS Engine research lanes, not a systematic review and not a claim that any paper has identified a defect in this code. It records (1) the result actually established by the cited work, (2) the method that may transfer, (3) what does *not* transfer without a source-level derivation, and (4) one falsifiable next discriminator in the current task lane.

## Provenance and coordination

- The latest Issue [#5 cross-agent handoff](https://github.com/miahpenn/0-Cosmos-engine/issues/5) recorded main as 8350f23fc919b5f2129d20f4671769667a86ccb7 for the literature-synthesis claim. Verify the actual base SHA again when reviewing this document; this is a claim-time anchor, not a permanent statement of current main.
- The checked-in [docs/STATUS.md](STATUS.md) and [Issue #135 task board](https://github.com/miahpenn/0-Cosmos-engine/issues/135) contain older snapshots. Do not reuse their embedded SHA or PR counts as live state without checking the newest evidence.
- Operating protocol: [docs/RESEARCH_WORKFLOW.md](RESEARCH_WORKFLOW.md); ownership and cross-lane corrections: [Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5); canonical task inventory: [Issue #135](https://github.com/miahpenn/0-Cosmos-engine/issues/135).
- No engine/vendor/production-physics changes, tests, or campaigns were performed to prepare this literature map.

## Evidence map

| Open question | Literature-supported point | Transferable method | Boundary / what it does not establish | Existing lane and next falsifiable gate |
|---|---|---|---|---|
| Centre regularity and momentum/constraint propagation | Spherical BSSN evolutions can be formulated with careful coordinate-origin treatment; continuum BSSN constraint propagation is a distinct system from the evolution equations. | Compare exact production and diagnostic discrete operators, parity rules, and signed contributions on the same frozen state. Test a smooth regular manufactured profile before evolving production data. | Neither a successful spherical code in a paper nor continuum well-posedness names the defective term in this implementation. Numerical constraint growth can originate in an operator/parity mismatch or in an unbalanced source; it should not be presumed to be one or the other. | [Issue #4](https://github.com/miahpenn/0-Cosmos-engine/issues/4); existing current-head Ricci formula test [PR #231](https://github.com/miahpenn/0-Cosmos-engine/pull/231) is separate. Next: signed term-by-term M-tangent budget with exact closure; report cells i=0,1 separately from a fixed physical annulus. |
| COSMOS/local density convention and scalar–DM transfer | Coupled-field literature specifies a model/action and then defines how component energy-momentum transfers cancel in the *combined* system. Different actions/coupling currents are different models. | Trace action, sign, units, normalization, and component stress-energy definitions back to the archive; derive a reciprocal transfer identity from that source. Add separate homogeneous-limit/Friedmann checks for the two code lanes. | A paper about coupled quintessence or interacting scalar fields is not evidence that COSMOS uses its interaction law, coupling coefficient, density normalization, or action. No conversion factor may be inferred just because it makes two equations numerically agree. | [Issue #210](https://github.com/miahpenn/0-Cosmos-engine/issues/210). Next: establish primary-action and dimensional provenance first; then document the mapping or explicitly keep the archive as a separately normalized reference model. |
| Radiation cone / continuation admissibility | Moment closures have realizability constraints: numerical states must continue to correspond to a physically possible underlying radiation distribution. A discretization can violate a continuum invariant unless its actual update preserves it. | Derive the cone inequality from the code's own variable normalization and closure; evaluate signed cone margin and its tangent on source-pinned saved states at the exact consumer. Include a synthetic negative control and all relevant stages. | The cited M1 paper uses a particular planar model and RKDG scheme. Its limiter, CFL restriction, or variables cannot be copied into this code without a formulation match. A cone margin affected at roundoff scale is not trajectory recovery; do not relax the cone or add clipping to hide a failing state. | Radiation lane tracked through [Issue #135](https://github.com/miahpenn/0-Cosmos-engine/issues/135) and [Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5). Next: source-matched frozen-state RHS/cone tangent; no long continuation until the current formula/tangent preflight yields a discriminating result. |
| Fixed-areal-radius Misner–Sharp balance | In spherical symmetry, the unified first law separates changes in quasi-local gravitational energy into energy-supply and work terms. | Evaluate all signed balance terms on one source-pinned sequence at a genuinely fixed areal radius R, mapping to the coordinate radius separately at each state. Independently recompute cell/face flux pairing and show closure. | Coordinate-fixed r is not generally fixed areal radius R in an evolving metric. A small residual or a mass-gradient identity alone is not proof that flux, work, and metric terms are all represented once. | Coordinate #35 / #79 / #89 dependencies are recorded on [Issue #135](https://github.com/miahpenn/0-Cosmos-engine/issues/135). Next: one same-state fixed-R witness including the required flux and work terms, reconciled with the centre/constraint checks. |
| Conservative matter and coupled source accounting | Flux-conservative relativistic-hydrodynamics formulations make the conserved variables, characteristic fluxes, and source terms explicit. | Test primitive-to-conserved recovery, signed face fluxes, and geometry/matter source contractions as separate terms before summing a global ledger. Keep every contribution tagged by stage and counted once. | The literature is a numerical-method precedent, not the project's equation of state, carrier/current selection, or a missing physical coupling law. Exact conservation of a local exchange pair does not imply global energy conservation in evolving spacetime. | Existing matter/EM tasks and research-only gates are cross-linked from [Issue #135](https://github.com/miahpenn/0-Cosmos-engine/issues/135). Preserve the opt-in/research-only boundary until the source law and independent conservation review are complete. |

## Selected peer-reviewed references

### 1. Origin regularity in spherical BSSN

**Pedro J. Montero and Isabel Cordero-Carrión (2012).** “BSSN equations in spherical coordinates without regularization: Vacuum and nonvacuum spherically symmetric spacetimes.” *Physical Review D* **85**, 124037. [https://doi.org/10.1103/PhysRevD.85.124037](https://doi.org/10.1103/PhysRevD.85.124037)

The paper describes a spherical BSSN and flux-conservative relativistic-hydrodynamics code, using a second-order partially implicit Runge–Kutta method, with regular-spacetime tests that remain stable without a separate origin-regularization algorithm. **Useful transfer:** the origin treatment and benchmark discipline are relevant to centre-cell diagnostics. **Limit:** a different formulation and implementation do not diagnose this repository's stencil, parity, or stress-source balance.

### 2. BSSN hyperbolicity and constraint-system reasoning

**Olivier Sarbach, Gioel Calabrese, Jorge Pullin, and Manuel Tiglio (2002).** “Hyperbolicity of the Baumgarte-Shapiro-Shibata-Nakamura system of Einstein evolution equations.” *Physical Review D* **66**, 064002. [https://doi.org/10.1103/PhysRevD.66.064002](https://doi.org/10.1103/PhysRevD.66.064002)

The authors analyze the relationship between BSSN and strongly/symmetric hyperbolic first-order systems under stated gauge conditions. **Useful transfer:** use continuum constraint propagation to frame a discrete audit, then independently measure the discrete defect. **Limit:** hyperbolicity is conditional on formulation/gauge assumptions and does not prove the discretized coupled RHS preserves the implementation's constraints.

### 3. Conservative relativistic hydrodynamics

**Francesc Banyuls, José A. Font, José M. Ibáñez, José M. Martí, and Juan A. Miralles (1997).** “Numerical 3+1 General Relativistic Hydrodynamics: A Local Characteristic Approach.” *The Astrophysical Journal* **476**, 221–231. [https://doi.org/10.1086/303604](https://doi.org/10.1086/303604)

This work expresses the relativistic-hydrodynamics equations in conservation form and develops a local characteristic approach. **Useful transfer:** audit conserved/primitive definitions and interface fluxes separately from geometric source terms. **Limit:** its fluid formulation does not establish the project's COSMOS density convention, interaction law, or production stress-energy choices.

### 4. Quasi-local mass, energy supply, and work in spherical symmetry

**Sean A. Hayward (1998).** “Unified first law of black-hole dynamics and relativistic thermodynamics.” *Classical and Quantum Gravity* **15**, 3147–3162. [https://doi.org/10.1088/0264-9381/15/10/017](https://doi.org/10.1088/0264-9381/15/10/017) · [arXiv:gr-qc/9710089](https://arxiv.org/abs/gr-qc/9710089)

This derives a spherical-GR relation for the gradient of active gravitational energy with distinct energy-supply and work terms. **Useful transfer:** use the split to design an independently recomputed fixed-areal-radius worldtube ledger. **Limit:** the continuum identity does not by itself determine the discrete face quadrature, moving-coordinate mapping, sign convention, or exact source terms in the engine.

### 5. Realizability-preserving radiation moments

**Edgar Olbrant, Cory D. Hauck, and Martin Frank (2012).** “A realizability-preserving discontinuous Galerkin method for the M1 model of radiative transfer.” *Journal of Computational Physics* **231**(17), 5612–5639. [https://doi.org/10.1016/j.jcp.2012.03.002](https://doi.org/10.1016/j.jcp.2012.03.002)

The authors study realizability constraints on M1 radiation moments and provide a numerical method with a proof of preservation under specified conditions. **Useful transfer:** distinguish the continuum admissible set from what the fully discrete update actually preserves. **Limit:** their geometry, closure, discretization, and CFL assumptions differ; use the paper to frame a test, not to import a limiter or relax a project's admission threshold.

### 6. Model-specific scalar-field interaction

**Luca Amendola (2000).** “Coupled quintessence.” *Physical Review D* **62**, 043511. [https://doi.org/10.1103/PhysRevD.62.043511](https://doi.org/10.1103/PhysRevD.62.043511)

This analyzes a particular cosmological scalar-field/matter model with specified exponential potential and linear coupling. **Useful transfer:** coupled-field claims need an explicit model, rather than a coupling inferred from numerical residuals. **Limit:** its action and coupling assumptions cannot be assigned to COSMOS by analogy.

**Gabriela Garcia-Arroyo, L. Arturo Ureña-López, and J. Alberto Vázquez (2024).** “Interacting scalar fields: Dark matter and early dark energy.” *Physical Review D* **110**, 023529. [https://doi.org/10.1103/PhysRevD.110.023529](https://doi.org/10.1103/PhysRevD.110.023529)

The paper explores specified scalar-field models coupled through their conservation equations and discusses background and perturbation-level consequences. **Useful transfer:** make the transfer current and total conservation statement explicit at equation level. **Limit:** this is not a first-principles derivation of the project's GEAR/COSMOS action or a basis for inserting a reciprocal source absent from the archive.

## Recommended sequence (least expensive first)

1. **Momentum/centre:** continue the active source-only preflight on [Issue #4](https://github.com/miahpenn/0-Cosmos-engine/issues/4). Reuse pinned states/artifacts, compare RHS and diagnostic operators term by term, report signed terms and algebraic closure, and isolate centre cells i=0,1 from a fixed annulus. No production equation edit until an exact discrete identity fails and a minimal profile distinguishes competing explanations.
2. **Density/action provenance:** resolve [Issue #210](https://github.com/miahpenn/0-Cosmos-engine/issues/210) from archived equations/source documents. State dimensions, units, normalization, action, and sign conventions in writing before any conversion or exchange implementation.
3. **Radiation:** derive the exact implemented cone formula and source-pinned tangent from current main. Keep saved-state arithmetic separate from full trajectory continuation.
4. **Mass balance:** coordinate existing #35/#79/#89 work into one fixed-areal-radius test with independently computed mass, flux, and work terms. Do not launch a new expensive campaign until this local identity is discriminating.
5. **Independent review:** for each result require a second reviewer, exact source/head SHA, actual workflow run links (mark skipped/missing honestly), reproducible artifact/hash, what the result establishes, what it does not establish, and one next gate.

## Owner question requiring primary-source confirmation

**Which exact archived GEAR/COSMOS document or version defines the action and dimensional normalization for the archived COSMOS density and its relationship to the local κ = 8πG Einstein-source convention?** The repository issue identifies the discrepancy, but the generic journal literature above cannot decide this project-specific mapping. Until the primary source is identified and reconciled, keep the mapping open and do not add a coupling or conversion based only on numerical fit.

## Decision boundary

This map proposes tests and source questions only. It does not promote diagnostic PRs, alter physics, clear governance checks, or change the repository's **NOT ADMITTED** scientific status. Any material update should be linked back to the specific task and then summarized on Issue #5; do not start a duplicate tracker or copy transient CI events into a second status log.
