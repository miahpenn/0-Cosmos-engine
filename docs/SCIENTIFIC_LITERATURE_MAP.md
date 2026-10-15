# Scientific literature map: open COSMOS Engine admission questions

**Purpose.** A focused literature map for the existing scientific-admission blockers. This is a source guide and research handoff, **not** a new task board, a physics decision, or evidence that a campaign is admitted. Keep ownership and queue state on [Issue #135](https://github.com/miahpenn/0-Cosmos-engine/issues/135), cross-agent findings on [Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5), and technical acceptance criteria in the existing task issue/PR.

**Repository snapshot used for orientation:** live `main` was `8350f23fc919b5f2129d20f4671769667a86ccb7` (commit timestamp 2026-10-15 03:49:40 UTC). This is a literature-map anchor only, not a verified production-source baseline. At capture, `docs/STATUS.md` identified `7860d11...` and Issue #135 identified `6aa7a9d...`; both are older than this main tip. Treat those as historical snapshots until the canonical status document is refreshed. Verify live refs before launching work.

## How to use the papers

Use journal papers to establish equations and numerical methodology, not to infer which unrecorded physical model COSMOS intends. For each open question, the evidence chain should be:

1. **Model provenance:** original action/equations, sign convention, dimensions, units, and archive/source revision.
2. **Discrete identity:** derive the continuous identity under those conventions, then derive the exact semi-discrete version implemented by the production RHS and diagnostic.
3. **Discriminator:** source-pinned analytic/manufactured profile or frozen-state test, including signed terms, closure and a negative control.
4. **Admission:** only after the prior gates, independent review plus resolution and full-horizon/campaign evidence. A passing unit suite, a manufactured oracle, or one-step smoke does not suffice.

## Existing open questions mapped to literature

### 1. Central regularity and momentum/constraint propagation — [Issue #4](https://github.com/miahpenn/0-Cosmos-engine/issues/4), Ricci diagnostics [PR #220](https://github.com/miahpenn/0-Cosmos-engine/pull/220) and [PR #231](https://github.com/miahpenn/0-Cosmos-engine/pull/231)

**Relevant journal evidence.** Baumgarte et al. (2013) discuss BSSN evolution in spherical-polar coordinates, coordinate singularities, and stability/accuracy/convergence tests. Montero & Cordero-Carrión (2012) study the spherical BSSN/PIRK formulation. Calabrese et al. (2002) show that apparently reasonable discretizations of weakly hyperbolic formulations can be nonconvergent; extra numerical dissipation does not cure the underlying issue.

**What transfers to this code.** Separate the centre-adjacent cells from a fixed physical annulus; test parity and the exact derivative operators used by the constraint and RHS; compare *signed* contributions and the measured constraint tangent; use a smooth regular manufactured profile with a deliberate parity-breaking negative control. The existing Issue #4 acceptance criteria already specify this bounded test; this map does not claim to have run it.

**Limits.** The cited formulations and gauges are not identical to the current engine's complete coupled RHS. PR #231's nonzero trace-free-Ricci oracle is formula/test coverage, not proof of production-centre correctness. A manufactured profile verifies its declared oracle and sampled operators, not trajectory stability or continuum convergence.

### 2. Radiation admissibility and outer continuation — use the current radiation lane cross-linked from [Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5) and [Issue #135](https://github.com/miahpenn/0-Cosmos-engine/issues/135)

**Relevant journal evidence.** Liebendörfer, Mezzacappa & Thielemann (2001) derive a conservative general-relativistic radiation-hydrodynamics formulation in spherical symmetry. Park (2006) derives the covariant energy/momentum exchange between matter and radiation. Xu, Ding & Wu (2022) prove physical-constraint preservation for a particular special-relativistic hydrodynamics scheme; their theorem depends on that scheme's fluxes, limiter and assumptions.

**What transfers to this code.** Keep matter and radiation exchange terms equal-and-opposite in the *same declared frame and stage*, and distinguish the conserved total stress-energy statement from admissibility of the radiation moments or continuation geometry. The next useful gate is to identify the first saved state and first discrete term responsible for a negative cone margin, on an exact source/input pin, then check whether the same signed margin is predicted by the local discrete tangent. Record source cone definition, consumer, stage, cell, RHS terms and time-stepper treatment separately.

**Limits.** The relativistic-fluid positivity/sub-luminality constraints in a published PCP scheme are not the engine's radiation-cone condition. Do not import a limiter, relax a cone threshold, add a floor/clamp or launch a long continuation merely because the paper uses one. A consumer-only change at roundoff scale is not trajectory repair.

### 3. Fixed-areal-radius mass balance — existing worldtube/Misner–Sharp work in [Issue #135](https://github.com/miahpenn/0-Cosmos-engine/issues/135) and [Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5)

**Relevant journal evidence.** Misner & Sharp (1964) derive the dynamical spherical-fluid mass relation and pressure-work term; their stated fluid model excludes heat flow and radiation. Hayward (1998) formulates a spherical-symmetry “unified first law” separating energy-supply and work terms. Together these are useful equation-level references, but differ in assumptions and notation.

**What transfers to this code.** Re-derive the balance for the actual metric variables, lapse/shift convention, complete stress-energy tensor and selected worldtube. Hold *areal radius* (R) fixed when the target is (R=10); do not replace it with coordinate radius (r=10). List separately the boundary flux, pressure/work, radiation exchange, charged-matter/Lorentz terms where enabled, and metric/geometry terms; check orientation and signs from the actual implementation. Independently recompute the fixed-(R) observable from raw operands.

**Limits.** Do not transplant the ideal-fluid Misner–Sharp formula unchanged into a system with radiation, scalar sectors, electromagnetic stress, or a moving coordinate grid. Roundoff closure of a local face/cell pairing is evidence for that local operator only, not global conservation in a dynamical spacetime.

### 4. COSMOS/local density normalization and units — [Issue #210](https://github.com/miahpenn/0-Cosmos-engine/issues/210)

**Repository premise to verify, not assume.** Issue #210 records that the archival homogeneous lane uses (H^2=\rho/3) and \,\(\dot H=-(\rho+p)/2\), while the local production lane declares \,\(\kappa=8\pi G\); it also records hardcoded archival densities and an apparently unused `CosmosParams.G`. These expressions can reflect different normalizations, but the intended mapping cannot be selected from numerical convenience.

**Relevant journal evidence.** Reviews of modified gravity and explicit matter–geometry coupling show that conservation equations depend on the action and on how matter couples to geometry. Wu et al. (2018) derive non-conservation in a specific Palatini (f(R,T)) theory. That is a model-specific result, not a universal permission to use non-conserved matter or to rescale \,\(\rho\).

**Next gate.** Trace each archived input to the primary archive/paper/action, state whether \,\(\rho\) is physical energy density or a rescaled variable, audit dimensions of (G), (H), \,\(\rho\), and (p), then derive the conversion to the local Einstein-source convention. Add the two separate homogeneous/Friedmann regressions specified in Issue #210 using an initial state not reconstructed from the hardcoded present-day Hubble value. Do not add a conversion or couple `exchange_Q` before source provenance and independent review.

### 5. Scalar–dark-matter exchange and full anisotropic-shear completeness — existing GEAR/COSMOS archive research

**Relevant journal evidence.** Garcia-Arroyo, Ureña-López & Vázquez (2024) and Liu, Tsujikawa & Ichiki (2024) analyze interacting scalar/dark-matter scenarios with energy and, in some models, momentum transfer. They demonstrate that exchange structure and phenomenology are model-dependent; they do not establish the coupling, sign, or species assignment in this repository. Modified-gravity examples likewise cannot substitute for the archived GEAR action.

**Next gate.** Find the exact primary action/field equations used by the archive; vary them (or cite the already-derived equations), derive each continuity/exchange equation and its dimensions, and verify reciprocal source signs from those equations. Separately check whether a spherical/isotropic specialization discards shear terms present in the archived Bianchi-I lane. The repository explicitly warns that the exact spherical branch alone does not represent full Bianchi-I shear. No coupling coefficient, fitted law, or “missing RHS call” should be inferred from a mismatch.

### 6. Electromagnetic field–matter–geometry exchange — existing opt-in research lane; see [EM face-work note](RESEARCH_EM_ADJOINT_FACE_WORK_V1.md), [stress/geometry note](RESEARCH_EM_STRESS_GEOMETRY_FEEDBACK.md), and [Issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5)

**Methodological connection.** The same total-stress-energy accounting principle applies: sources transferred between EM field and charged matter must cancel in a common convention, while gravitational mass/geometry terms require their own derivation. The radiation-hydrodynamics papers above are analogies for bookkeeping discipline, not EM derivations.

**Next gate.** Continue the already-claimed bounded face-work derivation rather than opening a parallel calculation. Keep the charge carrier/current law, face flux, cell source, Lorentz power, metric-weight tangent and fixed-(R) mass terms distinct, signed and counted once. First establish the chosen action/Maxwell conventions and a closed discrete identity; only then consider a bounded counterfactual. Existing research-only probes do not admit a physical carrier law or full global energy conservation.

## Peer-reviewed reference shelf

1. **Misner, C. W. & Sharp, D. H. (1964).** “Relativistic Equations for Adiabatic, Spherically Symmetric Gravitational Collapse.” *Physical Review* **136**, B571. [DOI: 10.1103/PhysRev.136.B571](https://doi.org/10.1103/PhysRev.136.B571). Dynamical spherical-fluid mass and pressure-work relation; ideal-fluid assumptions exclude heat flow/radiation.
2. **Hayward, S. A. (1998).** “Unified first law of black-hole dynamics and relativistic thermodynamics.” *Classical and Quantum Gravity* **15**. [DOI: 10.1088/0264-9381/15/10/017](https://doi.org/10.1088/0264-9381/15/10/017). Spherical active gravitational energy; energy-supply/work separation.
3. **Liebendörfer, M., Mezzacappa, A. & Thielemann, F.-K. (2001).** “Conservative general relativistic radiation hydrodynamics in spherical symmetry and comoving coordinates.” *Physical Review D* **63**, 104003. [DOI: 10.1103/PhysRevD.63.104003](https://doi.org/10.1103/PhysRevD.63.104003). Conservative radiation/matter formulation under a specific frame choice.
4. **Park, M.-G. (2006).** “Equations of general relativistic radiation hydrodynamics from a tensor formalism.” *Monthly Notices of the Royal Astronomical Society* **367**, 1739–1745. [DOI: 10.1111/j.1365-2966.2006.10081.x](https://doi.org/10.1111/j.1365-2966.2006.10081.x). Covariant matter–radiation four-force and equal/opposite exchange.
5. **Montero, P. J. & Cordero-Carrión, I. (2012).** “BSSN equations in spherical coordinates without regularization of the origin.” *Physical Review D* **85**, 124037. [DOI: 10.1103/PhysRevD.85.124037](https://doi.org/10.1103/PhysRevD.85.124037). Origin handling in a particular spherical BSSN/PIRK formulation.
6. **Baumgarte, T. W., Montero, P. J., Cordero-Carrión, I. & Müller, E. (2013).** “Numerical relativity in spherical polar coordinates: Evolution calculations with the BSSN formulation.” *Physical Review D* **87**, 044026. [DOI: 10.1103/PhysRevD.87.044026](https://doi.org/10.1103/PhysRevD.87.044026). Regularity, stability, accuracy and convergence tests in spherical-polar numerical relativity.
7. **Calabrese, G., Pullin, J., Sarbach, O. & Tiglio, M. (2002).** “Convergence and stability in numerical relativity.” *Physical Review D* **66**, 041501(R). [DOI: 10.1103/PhysRevD.66.041501](https://doi.org/10.1103/PhysRevD.66.041501). Warns that standard-looking refinement tests and added dissipation may not expose weak-hyperbolicity-induced nonconvergence.
8. **Xu, L., Ding, S. & Wu, K. (2022).** “A physical-constraint-preserving finite volume WENO method for special relativistic hydrodynamics on unstructured meshes.” *Journal of Computational Physics* **466**, 111398. [DOI: 10.1016/j.jcp.2022.111398](https://doi.org/10.1016/j.jcp.2022.111398). Rigorous positivity/subluminality preservation for the paper's own special-relativistic scheme; not the COSMOS radiation cone.
9. **Wu, J. et al. (2018).** “Palatini formulation of (f(R,T)) gravity theory, and its cosmological implications.” *The European Physical Journal C* **78**, 430. [DOI: 10.1140/epjc/s10052-018-5923-9](https://doi.org/10.1140/epjc/s10052-018-5923-9). Example of action-dependent matter non-conservation in one modified-gravity theory.
10. **Garcia-Arroyo, G., Ureña-López, L. A. & Vázquez, J. A. (2024).** “Interacting scalar fields: Dark matter and early dark energy.” *Physical Review D* **110**, 023529. [DOI: 10.1103/PhysRevD.110.023529](https://doi.org/10.1103/PhysRevD.110.023529). Specific interacting-scalar construction; useful comparison, not source selection for this engine.
11. **Liu, X., Tsujikawa, S. & Ichiki, K. (2024).** “Observational constraints on interactions between dark energy and dark matter with momentum and energy transfers.” *Physical Review D* **109**, 043533. [DOI: 10.1103/PhysRevD.109.043533](https://doi.org/10.1103/PhysRevD.109.043533). Illustrates that conclusions are conditional on the chosen interaction model and observable/data assumptions.

## Owner question that blocks a physics decision

**Can you point us to the authoritative primary source for the COSMOS archived homogeneous equations and the S/D-to-matter exchange premise** (paper, derivation, notebook, or archive file with its version)? Specifically, we need the action or defining equations and the meaning/units of \,\(\rho\), (G), and (H). Without that provenance, the literature can guide consistency checks but cannot determine the correct density conversion, exchange sign, or whether the archived COSMOS lane is intended as a separate normalized reference model.

## Safety of interpretation

This document contains no engine-code change, new physical law, fitted coefficient, tolerance relaxation, campaign run, merge, or scientific-admission decision. Reconcile against the latest main and existing task artifacts before using it operationally; keep admission **NOT ADMITTED** until the canonical tests and independent review gates are satisfied.
