# GEAR Blueprint → Production Dynamic-Interface Audit

**Date:** 2026-10-12  
**Repository:** \`miahpenn/0-Cosmos-engine\`  
**Audited production snapshot:** \`main\` at \`2035acc1210db18da0cededea18750b7ed813260\`  
**Report branch:** \`research/gear-dynamic-interface-audit-2026-10-12\`  
**Scope:** Source and blueprint audit only. No production physics/code was changed and no simulation was launched for this audit.

## Executive conclusion

The proposed “dynamic interface” hypothesis is worth pursuing, but the GEAR history requires a distinction between two architectures.

1. **The early GEAR archive had two lanes and open bridges.** It explicitly left local proper time ↔ cosmic time, D-mode energy → localized matter, and local source → global cosmological source unresolved. Those are not established laws to copy into production.
2. **The later V5.5 blueprint deliberately replaced the old separately prescribed local/COSMOS interface with one global spherical spacetime and one total stress-energy tensor.** The transition region is produced by the fields and geometry, rather than by an invented shell or source term.
3. **The current \`main\` call graph largely implements that unified-state architecture.** It advances one geometry, one scalar bundle carrying local \`S\`, inverted \`D\`, and spatial COSMOS \`phi\`, plus conservative dark matter, baryons and radiation. These all enter the shared geometry/stress-energy evolution.
4. Therefore, **the fact that \`cosmos_rhs()\` is not called by each production step is not, by itself, a missing-call bug.** \`engine/cosmos.py\` is an independent homogeneous control lane; its \`exchange_Q\` argument is explicitly ignored. Directly wiring that lane into V5.5 without deriving field ownership and matching would risk duplicate descriptions of \`phi\` and conflict with the single-metric design.
5. The strongest near-term hypothesis is narrower and testable: **the synchronized discrete RHS, CMC lapse/gauge, centre stencils, geometry stages and conservative matter updates may not preserve the identities required for one coupled physical system.** This is not yet proof that the continuum model lacks a physical equation or that any one equation is uniquely faulty.
6. Global completeness remains open: homogeneous Bianchi-I shear/anisotropic three-curvature is not part of the exact spherical production geometry. If the model instead requires two distinct spacetime regions, a physical moving-worldtube matching/averaging treatment remains a separate architecture choice.

## 1. Blueprint chronology

### 1.1 GEAR full cyclic-cosmos archive (2026-10-02)

The archived harness separates:

- **Local lane:** GEAR-03 stable \`S\`, unstable/ejected \`D\`, GEAR-63/65 geometric barrier, local energy and transport.
- **COSMOS lane:** scalar, dark matter, baryons, radiation, shear, expansion/contraction and cycle diagnostics.
- **Explicitly open gates:** local proper time ↔ cosmic time; D-mode energy → localized matter source; local source → global cosmological source; full-GR oscillaton endpoint; first-principles bounce.

The archive prohibits arbitrary D-to-matter coefficients and arbitrary local-to-cosmic clock scaling. It says the cosmological lane was independently evolved rather than artificially forced to consume D energy, and the canonical baseline must not be silently rewritten.

Source: [GEAR full cyclic-cosmos archive snapshot in Library](https://github.com/miahpenn/0-Cosmos-engine/tree/main).

### 1.2 GEAR-66: force-to-Q is an accounting chain, not a production source law

GEAR-66 specifies:

\`local force → dissipated power → event energy → event/ejection rate → cosmological exchange Q\`

It rejects equating instantaneous force with \`Q\` by definition. The derived \`Q_mech(a)\` is supposed to match an independently obtained cosmological \`Q(a)\` without a fitted comparison-stage normalization. The event/node population and production/ejection history, complete force law and normalization match remained unpaid. “Closed / evidentiary” means the test/accounting chain was made explicit; it does not mean the full mechanism was shown to generate cosmological \`Q\`.

**Decision:** GEAR-66 does not authorize inserting a guessed \`Q\` into the V5.5 production kernel.

### 1.3 V5.4: covariant worldtube for a genuinely two-region model

The V5.4 covariant audit establishes that if local and cosmological regions remain distinct, a physical interface must be a dynamically determined timelike worldtube. It needs complete energy **and normal-momentum** transfer, constructed from the local \`T_mu_nu\`, worldtube tangent/normal and surface motion. Moving-surface power is not generally the same as a fixed-coordinate flux. The fixed outer numerical boundary is not, by itself, a physical interface.

The previous V5.4 Landau-tube bidirectional diagnostic was useful but resolution-sensitive and did not satisfy the required two-sided Israel matching. It remained diagnostic-only. The old \`Q=-dot(m_R)/V_c\` insertion was not a covariant boundary flux by itself; a donor/recipient transfer must be included consistently in the total budget.

Sources:
- [GEAR V5.4 covariant interface audit](https://github.com/miahpenn/0-Cosmos-engine/tree/main)
- [GEAR V5.4 current coupling log](https://github.com/miahpenn/0-Cosmos-engine/tree/main)

### 1.4 V5.5: the later single-spacetime architecture governs production

The later V5.5 master checkpoint and bidirectional production gate lock:

- one global spherical metric;
- one total stress-energy tensor;
- local \`S\`, inverted \`D\`, COSMOS scalar \`phi\`, beta-coupled DM, baryonic dust and radiation;
- a finite transition region generated by the coupled fields;
- no separate shell EOS, reservoir, phenomenological source, fitted coefficient or imposed bounce.

The intended interface is a **derived geometric/current ledger in the unified spacetime**, not a guessed second source term. This does not close all validation gates: long post-turnaround constraint preservation, full radiation admissibility, convergence and the final multi-cycle campaign remain open.

Sources: GEAR V5.5 master checkpoint and GEAR V5.5 bidirectional production gate, archived in the user's Library.

## 2. Current \`main\` call-graph observations

Audited against commit \`2035acc1210db18da0cededea18750b7ed813260\`.

| Layer | Source evidence | Finding |
|---|---|---|
| Evolved state | [\`engine/production_kernel.py\`](https://github.com/miahpenn/0-Cosmos-engine/blob/main/engine/production_kernel.py), blob \`cc1214ab899e43919bfbd197113b0696eab2958d\` | \`ProductionState\` has one geometry, one \`ScalarFields\`, one \`V55MatterState\`, coordinate/proper-time bookkeeping and ledgers; it does not carry a separate homogeneous \`CosmosState\`. |
| Scalar evolution | [\`engine/scalar_system.py\`](https://github.com/miahpenn/0-Cosmos-engine/blob/main/engine/scalar_system.py) | Local \`S\`, inverted \`D\`, and spatial COSMOS \`phi\` are evolved on the same geometry. The archive-corrected scalar potential is retained. |
| DM–phi exchange | [\`engine/matter_system.py\`](https://github.com/miahpenn/0-Cosmos-engine/blob/main/engine/matter_system.py), \`BETA_DM=-0.04\` in the production kernel | Scalar and dark-matter RHS both include the beta-coupled interaction. The source-level reciprocal channel exists; same-stage discrete exchange closure still needs explicit measurement. |
| Stress-energy | [\`engine/stress_energy.py\`](https://github.com/miahpenn/0-Cosmos-engine/blob/main/engine/stress_energy.py) | Scalar and conservative-fluid projections are assembled into one \`TotalStressEnergy\`; no independent interface reservoir is assembled here. |
| Geometry RHS | [\`engine/v55_pirk_adapter.py\`](https://github.com/miahpenn/0-Cosmos-engine/blob/main/engine/v55_pirk_adapter.py), blob \`f671b838d9aefd8926702d9c2235209f52cd6815\` | Total matter/stress projections contribute to geometry RHS components including \`Aa\`, \`K\` and connection evolution. |
| Synchronized step | \`V55ProductionKernel.step()\` in \`production_kernel.py\` | Geometry, scalar and matter are staged together; the CMC lapse is re-solved on stage states. The runner advances with \`kernel.step(state, dt)\`. |
| Global observables | [\`engine/cosmology_observables.py\`](https://github.com/miahpenn/0-Cosmos-engine/blob/main/engine/cosmology_observables.py) | \`H_eff=-<K>/3\`, e-folds and outer scalar/stress values are derived from the unified solution. These are diagnostics, not a second homogeneous evolution lane. |
| Homogeneous control | [\`engine/cosmos.py\`](https://github.com/miahpenn/0-Cosmos-engine/blob/main/engine/cosmos.py), blob \`01d36844ff1126305cc0b419b330791b469863cc\` | The independent homogeneous control has \`CosmosState\`; \`cosmos_rhs(exchange_Q)\` explicitly ignores \`exchange_Q\`. Production does not call it, consistently with the V5.5 architecture. |
| Interface helper | \`engine/interface.py\`, \`engine/coupled.py\` | Holds interface state and donor/recipient ledger bookkeeping, but is not a closed production receiving equation. Existence of these helpers does not demonstrate a production global-Q closure. |
| Validation state | \`production_kernel.py\` | Declares \`implemented_not_campaign_validated\`. This matches the current admission status. |

### Consequence

The present machine is not simply several independent toy engines running next to each other. Its spatial scalar/matter/geometry interactions are present. However, source-code presence alone does not establish that each coupled RHS is discretely compatible, that exchange cancels at every true integration stage, or that the global cosmological completion is physically complete.

## 3. Match the hypothesis to current numerical evidence

### A. Discrete RHS / centre / CMC compatibility

The October 12 continuation records a nonzero momentum-constraint tendency from an initially discrete constraint-consistent state; the actual one-step change follows the frozen RHS tangent. Improved radial quadrature greatly reduces the baseline first-cell Ricci contribution to \`Aa_dot\`, but a smaller trace-free lapse-Hessian contribution and a centre-to-neighbour \`K_dot\` jump remain on the uniform-quadrature profile. Signed component reconstruction closes to roundoff—this validates the decomposition of coded terms, not the intended continuum regularity identity.

The frozen CMC lapse centre-row audit shows good operator closure on the archived initial slice. The CMC target-rank probe indicates near-redundancy on its frozen slice. Neither proves the lapse equation is wrong; the unresolved issue is compatibility of the evolving gauge response, centre-series constraints and staged curvature RHS.

### B. Radiation cone margin through real predictor stages

Keep distinct the initial one-step failures at \`N=240/320\` in the first cell and the separate longer-path outer-domain/predictor failure. The latter moves outward when the domain is expanded, indicating an outer-domain/global-CMC/geometry-stage interaction under investigation—not a demonstrated radiation-law defect.

For each actual predictor stage, track the conserved radiation variables, metric coefficient in the cone condition, and margin \`E-|S|\` (or its correctly densitized equivalent). Close the budget to roundoff before attributing loss to flux, source, geometry or boundary treatment.

**Do not** clip conserved radiation variables, relax \`|S|/E ≤ 1\`, introduce a floor, or invent a new source to force admissibility.

### C. Global geometry scope

The exact spherical production branch excludes homogeneous Bianchi-I shear and anisotropic three-curvature as dynamical degrees of freedom. This is a known global-completeness boundary in the GEAR/REV16 record. It must be kept separate from the current local centre and radiation-stage failures.

### D. Distinct-regions alternative

If the intended model is still one global solved spacetime, the V5.5 field-generated transition region and total stress-energy are the selected architecture. If the intended model instead has a separately solved local region joined to an FLRW exterior, the V5.4 moving-worldtube/Israel-matching or a controlled averaging closure is required. These are mutually distinct model-ownership choices; silently mixing a second homogeneous state into V5.5 is not a derived interface.

## 4. Competing hypotheses and verdict

| Hypothesis | Verdict |
|---|---|
| Production forgot to call \`cosmos_rhs()\` | **Not supported as a wiring bug.** The homogeneous control is separate by design; direct insertion risks duplicate \`phi\` ownership and would still not supply reciprocal stress-energy closure. |
| A local-to-global \`Q\) source law is missing | **Open only if the two-region model is restored.** The old archive and GEAR-66 do not give a production-ready fitted-free source coefficient/population closure. |
| Coupled semi-discrete RHS/gauge/centre updates violate constraint tangency | **Strongest testable near-term hypothesis.** Evidence supports investigating it, but no single continuum equation has yet been proved wrong. |
| A global degree of freedom is outside the spherical model | **Known completeness limitation.** It does not by itself explain the current local predictor failures. |
| The system needs a physical moving worldtube | **Conditional on a two-region architecture.** It should not be inserted into the present single-spacetime kernel without resolving the blueprint choice. |

## 5. Recommended bounded diagnostic

Do not launch another long campaign or change production physics yet. The repository already has active centre/gauge, RHS-tangent, constraint and radiation-stage investigations. Reconcile their current artifacts and avoid duplicating active agent work.

Run a frozen-state **reciprocity and constraint-tangent audit**, with no time advance:

1. On admitted baseline and uniform-quadrature states, perturb spatial COSMOS \`phi\` and DM density independently. On the same frozen metric, recompute scalar RHS, DM RHS, scalar/fluid stress projections and geometry source blocks. Verify that each response required by the existing equations is present, appears at the same declared stage, and is not duplicated.
2. Evaluate the beta-coupled DM/phi donor-recipient exchange pair at identical predictor states and conventions. Report the continuum cancellation and the actual discrete stagewise balance separately.
3. Compare the full-RHS derivative of the discrete Hamiltonian, momentum and connection witnesses with the corresponding constraint-propagation identities, using the same centre parity, stencil, metric and boundary treatment on both sides.
4. For radiation, capture conserved variables and cone margin before and after each genuine predictor metric/flux/source stage; the stage budget must reconstruct the observed margin change.
5. Only after a specific signed source-to-RHS dependency or discrete identity is shown missing/inconsistent, test one isolated diagnostic-layer counterfactual. Require improvement in all implicated constraints and radiation admissibility, not just one residual.

This discriminates a missing implemented response from a response that exists but fails discrete compatibility. It tests the dynamic-interface idea without inventing new physics.

## 6. Guardrails and disposition

- Keep production \`main\` and canonical GEAR baselines unchanged.
- Do not add \`exchange_Q\`, \`Q_friction\`, shell pressure, a fitted clock scale, D-to-matter coefficient, arbitrary damping, floor/clamp, or forced bounce.
- Keep centre-series, radiation-cone, global-shear and moving-worldtube questions distinct.
- Do not interpret workflow success, a short completed run or passed unit tests as full scientific admission.
- Current full production admission remains open: constraints, centre regularity, radiation admissibility, Misner–Sharp balance and required convergence gates have not all passed.

**Bottom line:** the Library does identify unresolved interface questions, but the current V5.5 machine is not simply missing a call between two otherwise complete engines. The best immediate test is whether the *existing* coupled RHS, gauge and geometry stages preserve the identities needed for one physical system. That is the cleanest test of the “dynamic interface” hypothesis without forcing the answer.
