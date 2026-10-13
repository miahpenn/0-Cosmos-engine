# Electromagnetism literature review and model-selection gates

**Date:** 2026-10-13  
**Source baseline reviewed:** `main=6aa7a9d3107582606b129f6b5ed3d04ea5f27bf3`  
**Status:** research/design note; not an approved production model  
**Production effect:** none

## Executive disposition

The live V5.5 runtime has no Maxwell state, charge/current evolution, reciprocal Lorentz-force source, or EM stress-energy projection. Current-main PRs #136, #144 and #145 supply reference/preflight tests only; they are not an EM implementation.

Published theory provides a consistent way to build a conventional spherical Einstein-Maxwell extension. It does **not** decide what, if anything, carries electric charge in the 0-Cosmos/GEAR model. The archive's unresolved bridges do not establish a charge-generation law. Therefore this note does not assign charge to radiation, baryons, dark matter, or the existing real scalars S, D, and phi.

Two different goals must remain separate:

1. **Numerical/architecture extension:** a standard spherical Maxwell sector with a deliberately specified charged source, used to test and then couple the field equations. It is an added model sector, not a derivation from GEAR.
2. **Theory derivation:** demonstrate that the existing 0* dynamics imply a current/charge law. No literature reference alone completes that project-specific derivation.

Recommended implementation default *if the model owner approves a conventional extension*: introduce an explicit, separately named charged matter sector with its own declared number/charge current and source law. Keep all existing species neutral unless the model explicitly says otherwise. Do not make the new-sector choice silently, and do not enable it in production until its action/EOS/charge parameters and initial/boundary data are approved.

## Literature reviewed

### 1. Maxwell equations in curved spacetime; Eulerian field definitions

Jai-chan Hwang and Hyerim Noh, “Maxwell equations in curved spacetime,” *The European Physical Journal C* **83**, 969 (2023). DOI: https://doi.org/10.1140/epjc/s10052-023-12135-0 . Open-access article: https://link.springer.com/article/10.1140/epjc/s10052-023-12135-0

Relevance: the electric and magnetic fields depend on the observer. For this engine's 3+1/ADM architecture, use a clearly defined normal/Eulerian observer and derive charge/current projections using the same observer. Do not treat coordinate components of F as if they were the physical orthonormal electric field. This source is useful for conventions; it does not select the engine's charged species.

### 2. Spherical Einstein-Maxwell-scalar evolution

“Characteristic initial value problems for the Einstein-Maxwell-scalar field equations in spherical symmetry,” *Physical Review D* **111**, 124015 (2025), published 11 June 2025. DOI: https://doi.org/10.1103/f5fc-zpld . Journal page: https://journals.aps.org/prd/abstract/10.1103/f5fc-zpld

Relevance: establishes a published spherical Einstein-Maxwell-scalar system and studies its initial/boundary formulation. Its affine-null formulation is not a drop-in replacement for this engine's CMC/PIRK 3+1 formulation; use it as an existence/constraint/boundary reference, not copy-and-paste source code.

### 3. Charged scalar fields and Noether current

Flavio Rossetti, “Strong Cosmic Censorship for the Spherically Symmetric Einstein–Maxwell-Charged-Klein–Gordon System with Positive Lambda: Stability of the Cauchy Horizon and H1 Extensions,” *Annales Henri Poincaré* **26**, 675–753 (2025; published online 7 June 2024). DOI: https://doi.org/10.1007/s00023-024-01454-x . Open access: https://link.springer.com/article/10.1007/s00023-024-01454-x

Relevance: its charged scalar is a complex Klein-Gordon field coupled with a U(1) gauge-covariant derivative; the current is derived from that field/action. This is materially different from the current engine's independent real S, D and phi equations. A charged-scalar route means a model change to a complex/gauge-coupled field, not merely adding qA terms to the current real scalar RHS.

### 4. Conservative relativistic fluid equations (Valencia)

José A. Font, “Numerical Hydrodynamics and Magnetohydrodynamics in General Relativity,” *Living Reviews in Relativity* **11**, article 7 (2008). DOI: https://doi.org/10.12942/lrr-2008-7 . Open review: https://link.springer.com/article/10.12942/lrr-2008-7

Relevance: documents the Eulerian Valencia conserved variables and conservative finite-volume approach. The production code evolves densitized quantities; a covariant four-force must be lowered/projected and converted to these exact active variables before entering an RHS. A mixed-component helper or a continuum sign argument is not a substitute for an active-state test.

### 5. Independent charged-fluid conservation

Bekenstein and Oron, “New conservation laws in general-relativistic magnetohydrodynamics,” *Physical Review D* **18**, 1809 (1978). DOI: https://doi.org/10.1103/PhysRevD.18.1809 . Journal page: https://journals.aps.org/prd/abstract/10.1103/PhysRevD.18.1809

Relevance: illustrates that relativistic charged-fluid dynamics couples Maxwell equations, fluid equations, energy-momentum conservation and particle/baryon number conservation. It reinforces that a current law and matter law must be specified together; it does not supply a 0*-specific charge assignment.

## Conditional reference equations for a conventional spherical extension

The following are a consistent *candidate convention set*, not a claim about GEAR. Use Gaussian units with c=1, signature (-,+,+,+), F=dA, and the Eulerian unit normal n^mu. Use one convention consistently throughout Maxwell, stress divergence, force projection and tests.

- Maxwell: ∇_mu F^(nu mu) = 4 pi J^nu (equivalently ∇_mu F^(mu nu) = -4 pi J^nu, with this index/sign convention).
- Charge continuity: ∇_mu J^mu = 0.
- EM stress-energy: T_EM^(mu nu) = (1/4 pi)[F^(mu alpha) F^nu_alpha - (1/4) g^(mu nu) F_alpha beta F^alpha beta].
- Reciprocal exchange: ∇_mu T_EM^(mu nu) = -F^(nu lambda) J_lambda; the charged matter must satisfy ∇_mu T_matter^(mu nu) = +F^(nu lambda) J_lambda.

For a purely radial electric field with no magnetic field or Poynting flux, the orthonormal Eulerian projections are:

- rho_EM = Ehat_r^2/(8 pi)
- p_r,EM = -rho_EM
- p_t,EM = +rho_EM
- j_r,EM = 0

That anisotropic stress must not be routed through the radiation-fluid equation of state.

### Spherical charge state

For the code's radial grid, let sqrt(gamma) denote the spatial metric volume density including its spherical r^2 factor, and let R be areal radius. Define enclosed charge Q(r) by the radial Gauss constraint:

- dQ/dr = 4 pi sqrt(gamma) rho_e
- Ehat_r = Q/R^2

For Eulerian charge density rho_e and Eulerian radial current j_e^r, charge continuity in the active ADM variables is:

- d_t(sqrt(gamma) rho_e) + d_r[sqrt(gamma)(alpha j_e^r - beta^r rho_e)] = 0

Thus a symmetry-reduced dynamical field can be represented by the enclosed-charge profile Q(t,r) updated by the same stage's radial charge flux, with Ehat_r reconstructed from Q and areal radius. This is not a propagating transverse-wave solver; exact spherical symmetry excludes transverse electromagnetic radiation degrees of freedom. The center requires Q=O(R^3) for finite central charge density, hence Ehat_r=O(R). Boundary charge flux must be explicit.

The discretization must update charge and evaluate Q/E from a metric and current belonging to the same actual beginning/predictor/corrector/accepted stage. A post-step Gauss projection or clipping is not a substitute for compatible charge continuity and Gauss evolution.

### Active matter-source map (conditional)

Under the convention above, for a pure radial field E_r = sqrt(gamma_rr) Ehat_r and Eulerian charge/current split J^t = rho_e/alpha, J^r = j_e^r - beta^r rho_e/alpha, the lowered force covector Q_nu = F_nu lambda J^lambda is:

- Q_r = rho_e E_r
- Q_t = -alpha E_r j_e^r + beta^r rho_e E_r

For the active Eulerian energy/momentum variables used in the existing Valencia RHS, the conditional matter increments are:

- q_E = +alpha sqrt(gamma) E_r j_e^r
- q_Sr = +alpha sqrt(gamma) rho_e E_r

The electromagnetic field must carry the exact opposite exchange at the same stage. Re-derive all expressions if the convention, current definition, or state representation changes. In particular, the current-main matter momentum source conversion issue already tracked by draft PR #143 is a separate scalar-DM candidate: reconcile its review before layering an EM source, but do not conflate its source law with Lorentz force.

### Normalization

If the owner selects the standard Gaussian Maxwell kinetic action -1/(16 pi) integral sqrt(-g) F_mu nu F^mu nu d^4x with G=c=1, then the above field stress normalization follows. The physical-to-code field map must still be derived from the declared code length/time units and the approved charged-matter action/current. Existing scalar/fluid projections have sector-specific normalization paths. Do not divide an already consumer-normalized EM energy by kappa again, and do not fit a field-scale coefficient to cosmological output.

## Architecture fit and staged implementation gates

Current-main source paths reviewed:
- Owned step state and stage order: `engine/production_kernel.py::ProductionState`, `V55ProductionKernel._rhs`, and `step`.
- Existing fluid variables/source conversion: `engine/matter_system.py`, `engine/valencia.py`.
- Geometry/gauge projection: `engine/v55_matter.py::total_matter_projection`, `engine/v55_pirk_adapter.py::primary_l3_with_matter` and `lambda_l3_with_matter`, and `engine/cmc_gauge.py::solve_cmc_lapse`.
- Independent diagnostics/boundary projection: `engine/stress_energy.py::assemble_total_stress_energy` and `engine/production_kernel.py::_apply_outer_light_boundary`.
- Capability contract: `engine/production_contract.py`; draft PR #149 now marks `radial_maxwell=False` and adds an explicit opt-in fail-closed guard.

Implement only in an isolated, feature-gated candidate, in this order:

1. **E0 — Zero-field parity:** exactly zero Maxwell/charge state must preserve every baseline RHS and projection within a declared arithmetic tolerance.
2. **E1 — Charge/Gauss:** manufactured compact smooth charge fixture; analytic flat-space enclosed-charge/Coulomb oracle; curved frozen-metric Gauss residual, central regularity, explicit boundary charge flux, and spatial convergence. This is a manufactured source, not a GEAR-origin claim.
3. **E2 — Canonical stress projection:** calculate EM rho/pr/pt/j in one shared helper. At actual start/predictor/corrector/accepted states, require parity between `total_matter_projection` and `assemble_total_stress_energy`, including EM exactly once. Check the radial stress identities independently.
4. **E3 — Charge law and reciprocal exchange:** only after an owner-approved carrier/action/current law exists, test charge continuity and the signed matter↔EM energy and momentum sums at every actual integrator stage, including non-unit lapse and nonzero shift. A synthetic source may test this plumbing; it cannot count as physical charge selection.
5. **E4 — Geometry/gauge/boundary:** use the same EM projections in CMC target/source, K/Aa L3, connection/Lambda, constraint diagnostics, initialization (only for approved charged initial data), and boundary/admissibility paths. Fail closed if a required consumer is missing.
6. **E5 — Convergence and campaign admission:** source-off/on signed tests, spatial/time refinement, H/M/connection/determinant and fixed-areal-radius Misner–Sharp balance. No long campaign before the short discrete identities pass.

Do not add the physical source, a charge-to-mass ratio, coupling parameter, nonzero initializer, or field scale as an implementation convenience. Do not label a static Gauss fixture as a dynamically coupled trajectory.

## Decisions still requiring model-owner authorization

1. Is the aim a **conventional added spherical Einstein-Maxwell sector** or a derivation of electromagnetism from the existing GEAR dynamics? These are distinct deliverables.
2. If a conventional extension is allowed, should charge belong to a new explicitly specified charged-fluid species, or to a modified complex scalar with U(1) gauge symmetry? Neither follows from the present real scalar/fluid state.
3. Which action, current sign, charge normalization/units, equation of state/number law, and initial/boundary data are selected? These are related parts of one model, not independent knobs to fit.
4. Is the intended scope exactly spherical radial electric/Coulomb coupling, or is a symmetry-relaxed Maxwell-wave architecture also intended?

Until those are answered, proceed with conventional reference tests, consumer parity, explicit capability gating and documentation only; keep production physics and scientific admission unchanged.
