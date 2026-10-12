# GEAR Blueprint → Production Dynamic-Interface Audit

**Date:** 2026-10-12  
**Repository:** `miahpenn/0-Cosmos-engine`  
**Audited production snapshot:** `main` at `2035acc1210db18da0cededea18750b7ed813260`  
**Report branch:** `research/gear-dynamic-interface-audit-2026-10-12`  
**Scope:** Source and blueprint audit only. No production physics/code was changed and no simulation was launched for this audit.

## Executive conclusion

The proposed “dynamic interface” hypothesis is worth pursuing, but the GEAR history requires a distinction between two architectures.

1. **The early GEAR archive had two lanes and open bridges.** It explicitly left local proper time ↔ cosmic time, D-mode energy → localized matter, and local source → global cosmological source unresolved. Those are not established laws to copy into production.
2. **The later V5.5 blueprint deliberately replaced the old separately prescribed local/COSMOS interface with one global spherical spacetime and one total stress-energy tensor.** The transition region is produced by the fields and geometry, rather than by an invented shell or source term.
3. **The current `main` call graph largely implements that unified-state architecture.** It advances one geometry, one scalar bundle carrying local `S`, inverted `D`, and spatial COSMOS `phi`, plus conservative dark matter, baryons and radiation. These all enter the shared geometry/stress-energy evolution.
4. Therefore, **the fact that `cosmos_rhs()` is not called by each production step is not, by itself, a missing-call bug.** `engine/cosmos.py` is an independent homogeneous control lane; its `exchange_Q` argument is explicitly ignored. Directly wiring that lane into V5.5 without deriving field ownership and matching would risk duplicate descriptions of `phi` and conflict with the single-metric design.
5. The strongest near-term hypothesis is narrower and testable: **the synchronized discrete RHS, CMC lapse/gauge, centre stencils, geometry stages and conservative matter updates may not preserve the identities required for one coupled physical system.** This is not yet proof that the continuum model lacks a physical equation or that any one equation is uniquely faulty.
6. Global completeness remains open: homogeneous Bianchi-I shear/anisotropic three-curvature is not part of the exact spherical production geometry. If the model instead requires two distinct spacetime regions, a physical moving-worldtube matching/averaging treatment remains a separate architecture choice.

## 1. Blueprint chronology

### 1.1 GEAR full cyclic-cosmos archive (2026-10-02)

The archived harness separates:

- **Local lane:** GEAR-03 stable `S`, unstable/ejected `D`, GEAR-63/65 geometric barrier, local energy and transport.
- **COSMOS lane:** scalar, dark matter, baryons, radiation, shear, expansion/contraction and cycle diagnostics.
- **Explicitly open gates:** local proper time ↔ cosmic time; D-mode energy → localized matter source; local source → global cosmological source; full-GR oscillaton endpoint; first-principles bounce.

The archive prohibits arbitrary D-to-matter coefficients and arbitrary local-to-cosmic clock scaling. It says the cosmological lane was independently evolved rather than artificially forced to consume D energy, and the canonical baseline must not be silently rewritten.

Source: [GEAR full cyclic-cosmos archive snapshot in Library](https://github.com/miahpenn/0-Cosmos-engine/tree/main).

### 1.2 GEAR-66: force-to-Q is an accounting chain, not a production source law

GEAR-66 specifies:

`local force → dissipated power → event energy → event/ejection rate → cosmological exchange Q`

It rejects equating instantaneous force with `Q` by definition. The derived `Q_mech(a)` is supposed to match an independently obtained cosmological `Q(a)` without a fitted comparison-stage normalization. The event/node population and production/ejection history, complete force law and normalization match remained unpaid. “Closed / evidentiary” means the test/accounting chain was made explicit; it does not mean the full mechanism was shown to generate cosmological `Q`.

**Decision:** GEAR-66 does not authorize inserting a guessed `Q` into the V5.5 production kernel.

### 1.3 V5.4: covariant worldtube for a genuinely two-region model

The V5.4 covariant audit establishes that if local and cosmological regions remain distinct, a physical interface must be a dynamically determined timelike worldtube. It needs complete energy **and normal-momentum** transfer, constructed from the local `T_mu_nu`, worldtube tangent/normal and surface motion. Moving-surface power is not generally the same as a fixed-coordinate flux. The fixed outer numerical boundary is not, by itself, a physical interface.

The previous V5.4 Landau-tube bidirectional diagnostic was useful but resolution-sensitive and did not satisfy the required two-sided Israel matching. It remained diagnostic-only. The old `Q=-dot(m_R)/V_c` insertion was not a covariant boundary flux by itself; a donor/recipient transfer must be included consistently in the total budget.

Sources:
- [GEAR V5.4 covariant interface audit](https://github.com/miahpenn/0-Cosmos-engine/tree/main)
- [GEAR V5.4 current coupling log](https://github.com/miahpenn/0-Cosmos-engine/tree/main)

### 1.4 V5.5: the later single-spacetime architecture governs production

The later V5.5 master checkpoint and bidirectional production gate lock:

- one global spherical metric;
- one total stress-energy tensor;
- local `S`, inverted `D`, COSMOS scalar `phi`, beta-coupled DM, baryonic dust and radiation;
- a finite transition region generated by the coupled fields;
- no separate shell EOS, reservoir, phenomenological source, fitted coefficient or imposed bounce.

The intended interface is a **derived geometric/current ledger in the unified spacetime**, not a guessed second source term. This does not close all validation gates: long post-turnaround constraint preservation, full radiation admissibility, convergence and the final multi-cycle campaign remain open.

Sources: GEAR V5.5 master checkpoint and GEAR V5.5 bidirectional production gate, archived in the user's Library.

## 2. Current `main` call-graph observations

Audited against commit `2035acc1210db18da0cededea18750b7ed813260`.

| Layer | Source evidence | Finding |
|---|---|---|
| Evolved state | [`engine/production_kernel.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/main/engine/production_kernel.py), blob `cc1214ab899e43919bfbd197113b0696eab2958d` | `ProductionState` has one geometry, one `ScalarFields`, one `V55MatterState`, coordinate/proper-time bookkeeping and ledgers; it does not carry a separate homogeneous `CosmosState`. |
| Scalar evolution | [`engine/scalar_system.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/main/engine/scalar_system.py) | Local `S`, inverted `D`, and spatial COSMOS `phi` are evolved on the same geometry. The archive-corrected scalar potential is retained. |
| DM–phi exchange | [`engine/matter_system.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/main/engine/matter_system.py), `BETA_DM=-0.04` in the production kernel | Scalar and dark-matter RHS both include the beta-coupled interaction. The source-level reciprocal channel exists; same-stage discrete exchange closure still needs explicit measurement. |
| Stress-energy | [`engine/stress_energy.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/main/engine/stress_energy.py) | Scalar and conservative-fluid projections are assembled into one `TotalStressEnergy`; no independent interface reservoir is assembled here. |
| Geometry RHS | [`engine/v55_pirk_adapter.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/main/engine/v55_pirk_adapter.py), blob `f671b838d9aefd8926702d9c2235209f52cd6815` | Total matter/stress projections contribute to geometry RHS components including `Aa`, `K` and connection evolution. |
| Synchronized step | `V55ProductionKernel.step()` in `production_kernel.py` | Geometry, scalar and matter are staged together; the CMC lapse is re-solved on stage states. The runner advances with `kernel.step(state, dt)`. |
| Global observables | [`engine/cosmology_observables.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/main/engine/cosmology_observables.py) | `H_eff=-<K>/3`, e-folds and outer scalar/stress values are derived from the unified solution. These are diagnostics, not a second homogeneous evolution lane. |
| Homogeneous control | [`engine/cosmos.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/main/engine/cosmos.py), blob `01d36844ff1126305cc0b419b330791b469863cc` | The independent homogeneous control has `CosmosState`; `cosmos_rhs(exchange_Q)` explicitly ignores `exchange_Q`. Production does not call it, consistently with the V5.5 architecture. |
| Interface helper | `engine/interface.py`, `engine/coupled.py` | Holds interface state and donor/recipient ledger bookkeeping, but is not a closed production receiving equation. Existence of these helpers does not demonstrate a production global-Q closure. |
| Validation state | `production_kernel.py` | Declares `implemented_not_campaign_validated`. This matches the current admission status. |

### Consequence

The present machine is not simply several independent toy engines running next to each other. Its spatial scalar/matter/geometry interactions are present. However, source-code presence alone does not establish that each coupled RHS is discretely compatible, that exchange cancels at every true integration stage, or that the global cosmological completion is physically complete.

## 3. Match the hypothesis to current numerical evidence

### A. Discrete RHS / centre / CMC compatibility

The October 12 continuation records a nonzero momentum-constraint tendency from an initially discrete constraint-consistent state; the actual one-step change follows the frozen RHS tangent. Improved radial quadrature greatly reduces the baseline first-cell Ricci contribution to `Aa_dot`, but a smaller trace-free lapse-Hessian contribution and a centre-to-neighbour `K_dot` jump remain on the uniform-quadrature profile. Signed component reconstruction closes to roundoff—this validates the decomposition of coded terms, not the intended continuum regularity identity.

The frozen CMC lapse centre-row audit shows good operator closure on the archived initial slice. The CMC target-rank probe indicates near-redundancy on its frozen slice. Neither proves the lapse equation is wrong; the unresolved issue is compatibility of the evolving gauge response, centre-series constraints and staged curvature RHS.

### B. Radiation cone margin through real predictor stages

Keep distinct the initial one-step failures at `N=240/320` in the first cell and the separate longer-path outer-domain/predictor failure. The latter moves outward when the domain is expanded, indicating an outer-domain/global-CMC/geometry-stage interaction under investigation—not a demonstrated radiation-law defect.

For each actual predictor stage, track the conserved radiation variables, metric coefficient in the cone condition, and margin `E-|S|` (or its correctly densitized equivalent). Close the budget to roundoff before attributing loss to flux, source, geometry or boundary treatment.

**Do not** clip conserved radiation variables, relax `|S|/E ≤ 1`, introduce a floor, or invent a new source to force admissibility.

### C. Global geometry scope

The exact spherical production branch excludes homogeneous Bianchi-I shear and anisotropic three-curvature as dynamical degrees of freedom. This is a known global-completeness boundary in the GEAR/REV16 record. It must be kept separate from the current local centre and radiation-stage failures.

### D. Distinct-regions alternative

If the intended model is still one global solved spacetime, the V5.5 field-generated transition region and total stress-energy are the selected architecture. If the intended model instead has a separately solved local region joined to an FLRW exterior, the V5.4 moving-worldtube/Israel-matching or a controlled averaging closure is required. These are mutually distinct model-ownership choices; silently mixing a second homogeneous state into V5.5 is not a derived interface.

## 4. Competing hypotheses and verdict

| Hypothesis | Verdict |
|---|---|
| Production forgot to call `cosmos_rhs()` | **Not supported as a wiring bug.** The homogeneous control is separate by design; direct insertion risks duplicate `phi` ownership and would still not supply reciprocal stress-energy closure. |
| A local-to-global `Q` source law is missing | **Open only if the two-region model is restored.** The old archive and GEAR-66 do not give a production-ready fitted-free source coefficient/population closure. |
| Coupled semi-discrete RHS/gauge/centre updates violate constraint tangency | **Strongest testable near-term hypothesis.** Evidence supports investigating it, but no single continuum equation has yet been proved wrong. |
| A global degree of freedom is outside the spherical model | **Known completeness limitation.** It does not by itself explain the current local predictor failures. |
| The system needs a physical moving worldtube | **Conditional on a two-region architecture.** It should not be inserted into the present single-spacetime kernel without resolving the blueprint choice. |

## 5. Recommended bounded diagnostic

Do not launch another long campaign or change production physics yet. The repository already has active centre/gauge, RHS-tangent, constraint and radiation-stage investigations. Reconcile their current artifacts and avoid duplicating active agent work.

Run a frozen-state **reciprocity and constraint-tangent audit**, with no time advance:

1. On admitted baseline and uniform-quadrature states, perturb spatial COSMOS `phi` and DM density independently. On the same frozen metric, recompute scalar RHS, DM RHS, scalar/fluid stress projections and geometry source blocks. Verify that each response required by the existing equations is present, appears at the same declared stage, and is not duplicated.
2. Evaluate the beta-coupled DM/phi donor-recipient exchange pair at identical predictor states and conventions. Report the continuum cancellation and the actual discrete stagewise balance separately.
3. Compare the full-RHS derivative of the discrete Hamiltonian, momentum and connection witnesses with the corresponding constraint-propagation identities, using the same centre parity, stencil, metric and boundary treatment on both sides.
4. For radiation, capture conserved variables and cone margin before and after each genuine predictor metric/flux/source stage; the stage budget must reconstruct the observed margin change.
5. Only after a specific signed source-to-RHS dependency or discrete identity is shown missing/inconsistent, test one isolated diagnostic-layer counterfactual. Require improvement in all implicated constraints and radiation admissibility, not just one residual.

This discriminates a missing implemented response from a response that exists but fails discrete compatibility. It tests the dynamic-interface idea without inventing new physics.

## 6. Guardrails and disposition

- Keep production `main` and canonical GEAR baselines unchanged.
- Do not add `exchange_Q`, `Q_friction`, shell pressure, a fitted clock scale, D-to-matter coefficient, arbitrary damping, floor/clamp, or forced bounce.
- Keep centre-series, radiation-cone, global-shear and moving-worldtube questions distinct.
- Do not interpret workflow success, a short completed run or passed unit tests as full scientific admission.
- Current full production admission remains open: constraints, centre regularity, radiation admissibility, Misner–Sharp balance and required convergence gates have not all passed.

**Bottom line:** the Library does identify unresolved interface questions, but the current V5.5 machine is not simply missing a call between two otherwise complete engines. The best immediate test is whether the *existing* coupled RHS, gauge and geometry stages preserve the identities needed for one physical system. That is the cleanest test of the “dynamic interface” hypothesis without forcing the answer.


## 7. Cross-agent source reconciliation — verified equation mismatch, not the dominant defect

**Read-only source snapshot:** main at 2035acc1210db18da0cededea18750b7ed813260.  
**Coordination thread:** [GitHub issue #5](https://github.com/miahpenn/0-Cosmos-engine/issues/5).  
**Library save:** /0* Cosmos Engine/GEAR_DYNAMIC_INTERFACE_BLUEPRINT_AUDIT_2026-10-12.md.

### 7.1 A real scalar equation mismatch is confirmed

The spatial metric declared by engine/valencia.py and used by the matter kernel is:

\[
\gamma_{rr}=a/X^2,\qquad \gamma^{rr}=X^2/a.
\]

In engine/scalar_system.py, the scalar momentum RHS used X**-2 * alphap * fp for the lapse-gradient contraction, despite the same function already defining invr = X**2/a and using it in the scalar Laplacian. The covariant contraction is \(\gamma^{rr}\partial_r\alpha\,\partial_r f\), so invr * alphap * fp is the internally consistent expression.

An isolated candidate changed only this factor and added a regression test:
- [Candidate assessment](https://github.com/miahpenn/0-Cosmos-engine/blob/fix/scalar-inverse-metric-lapse-gradient-strict-2026-10-12/docs/SCALAR_INVERSE_METRIC_CANDIDATE_ASSESSMENT_2026-10-12.md)
- [Matched test/campaign run 38282790322](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38282790322), successful; artifact 11762490670, reported SHA-256 a3ede0fd7bb3b74851ccfadb244b5765a000b5e45b5f688131acc20218d5f2e2.

The matched N=40/60/80 histories all reached coordinate time 12, but the trajectory changes were tiny: absolute delta H_eff only about 2.6–3.9e-11, and proper-time endpoint differences about 1.1–1.4e-6. Centre-dominated constraints and Misner–Sharp residuals remained substantially open. Thus this is a demonstrated equation-correctness issue worth keeping isolated, **not** a demonstrated repair of the principal centre/constraint failure. It is not merged into production main.

### 7.2 The DM–COSMOS-scalar source coupling already exists

Source tracing shows a reciprocal beta-coupled channel in the spatial production path:
- scalar_system.scalar_rhs_arrays() adds + beta_dm * rho_dm to the COSMOS-scalar momentum RHS.
- production_kernel._rhs() supplies the same-stage spatial phi derivatives and dark-matter density to matter_system.species_rhs().
- valencia.dark_matter_covector_source() constructs Q_nu = beta rho_DM partial_nu phi, and the dark-matter conservative update includes the corresponding energy and momentum source.
- The scalar, dark matter, baryons and radiation stress projections are assembled in stress_energy.py and feed the shared geometry RHS through v55_pirk_adapter.py.

This refutes the broad hypothesis that the production kernel has no COSMOS-scalar/DM coupling at all. It does **not** prove same-stage discrete donor/recipient cancellation or full conservation through every PIRK stage. That remains an auditable numerical identity. Nor does this coupling identify inverted D with ordinary matter or solve local-to-global source closure.

### 7.3 The strongest remaining numerical lead is still centre/RHS compatibility

The independent manufactured regular-centre Ricci/Aa controls pass, while the production-initial centre profile does not satisfy the expected bounded Aa_dot/r^2 scaling. The constrained-initial-state momentum tangent also remains nonzero and closely matches the ordinary one-step change. Better initial Hamiltonian closure therefore does not restore full constraint propagation.

The radiation failures remain a separate gate: high-resolution first-cell failures occur in the one-step campaign, while the longer R80 failure is associated with a near-null outer-cell margin and moves outward under larger-domain control. Do not merge these into one presumed cause or weaken the |S|/E <= 1 condition.

### 7.4 Updated next action — reuse active artifacts, then test stagewise reciprocity

Issue #5 is the cross-agent coordination thread. The latest visible active tracks already cover centre-series/RHS tangent, CMC operator/target closure, bond sweep, scalar inverse-metric comparison and radiation predictor/P-minus tracing. No agent response to issue #5 has yet been verified; don't claim a human/agent acknowledgment until it is present.

The next useful bounded comparison is therefore not a new long campaign:

1. On identical admitted frozen states, calculate the beta-coupled scalar and DM source pieces at each actual predictor stage using identical metric, rho_DM, partial_t phi, and partial_r phi. Show the continuum donor/recipient cancellation and the discretized/stagewise budget separately.
2. On those same states, compare the full discrete constraint tangents with the production geometry/matter RHS and the already archived signed centre-series decomposition. Keep the centre row, first formally admitted cell, and fixed physical radii separate.
3. Reconcile the radiation predictor artifact before starting a further radiation run; verify whether its exact state/stage is already captured in the current artifact set.
4. Make no code change until the source-to-RHS dependency or signed identity that fails has been pinpointed. The scalar inverse-metric candidate remains a separate equation-correction branch, not evidence that centre regularity or physical interface closure is solved.

**Revised bottom line:** one metric-contraction inconsistency is confirmed and the local DM–COSMOS scalar exchange channel is already wired. The more consequential open question is whether the existing coupled source, gauge and geometry updates preserve the required continuum identities at the discrete stages. The early GEAR global-source bridge and anisotropic global completeness remain different, still-open physical-model questions.


## 8. Scalar–DM exchange audit reconciled across parallel branches

### 8.1 Governing identity

The archived V5.5 covariant exchange law is

\[
\nabla_\mu T_{DM}^{\mu\nu}=+\beta\rho_{DM}\nabla^\nu\phi,\qquad
\nabla_\mu T_{\phi}^{\mu\nu}=-\beta\rho_{DM}\nabla^\nu\phi.
\]

The spatial scalar kinematics are \(\partial_t\phi=\alpha\Pi+\beta^r\partial_r\phi\). The conservative DM source in the production variable conventions therefore has the source-only energy increment

\[
\Delta E^{DM}_{cons}=-\sqrt{\gamma}\,\alpha\,\beta\rho_{DM}\Pi,
\]

and the radial momentum increment

\[
\Delta S^{DM}_{cons}=\sqrt{\gamma}\,\alpha\,\beta\rho_{DM}\,\partial_r\phi.
\]

The reciprocal scalar source in the coordinate-time \(\Pi\) RHS is \(+\alpha\beta\rho_{DM}\). Removing either lapse factor makes the corresponding pointwise source pair fail when \(\alpha\ne1\).

### 8.2 What the older October 10 branch did—and did not—close

The existing [GEAR lapse accounting audit](https://github.com/miahpenn/0-Cosmos-engine/blob/physics/spatial-beta-lapse-source-correction/docs/SPATIAL_BETA_EXCHANGE_LAPSE_ACCOUNTING_AUDIT_2026-10-10.md) derived and corrected the scalar lapse factor. Its test passed 111 tests on [run 38060315441](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38060315441), but checked source-pair energy cancellation only. A later cross-agent review found the dark-matter radial source still used sqrt(gamma) Q_r rather than alpha sqrt(gamma) Q_r. The older test was therefore a useful, but incomplete, source-pair gate.

### 8.3 The concurrent October 12 audit covers both components

The newer isolated branch [diag/covariant-beta-lapse-source-2026-10-12](https://github.com/miahpenn/0-Cosmos-engine/tree/diag/covariant-beta-lapse-source-2026-10-12) changes both the scalar Pi source and the DM radial source to carry the lapse, and its [audit report](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/covariant-beta-lapse-source-2026-10-12/docs/COVARIANT_BETA_LAPSE_SOURCE_AUDIT_2026-10-12.md) contains a synthetic, frozen-state test of pointwise energy and radial-current cancellation. It did not alter main or run a trajectory.

### 8.4 Current-main source-pair gate: passed

To validate the same source-only candidate against the precise October 12 production snapshot, a second isolated branch was created from main at SHA 2035acc1210db18da0cededea18750b7ed813260:

- Candidate branch: [diag/current-main-beta-source-gate-2026-10-12](https://github.com/miahpenn/0-Cosmos-engine/tree/diag/current-main-beta-source-gate-2026-10-12)
- Candidate scalar source: [engine/scalar_system.py](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/current-main-beta-source-gate-2026-10-12/engine/scalar_system.py), using alpha * beta_dm * rho_dm
- Candidate DM radial source: [engine/matter_system.py](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/current-main-beta-source-gate-2026-10-12/engine/matter_system.py), using alpha * sqrt(gamma) * Q_r
- Current-main regression: [test_spatial_beta_exchange_current_main.py](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/current-main-beta-source-gate-2026-10-12/tests/test_spatial_beta_exchange_current_main.py)
- Passing workflow: [38302168862](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38302168862), head commit 88e4d2e809c8c4b30ed4e37dd33b84342502e1be

The regression initializes through the current production kernel and creates an unadvanced diagnostic state with non-unit lapse and nonzero shift. It toggles the beta coupling on/off on the same state, then checks the source increments and their pointwise energy and radial-current cancellation. The observed latest run results are:

- Direct current-main source-pair test: **1 passed in 0.34 s**.
- Full repository suite on this candidate: **80 passed in 7.04 s**.
- Pinned numerical vendor: d6052d605673ce9d82cdc99b0df79855fa2ea215.
- No time advance or trajectory campaign was part of this workflow.

This is a useful result: the exact current-main candidate source pair closes under the synthetic, unadvanced test. It is not evidence that the candidate has been evolved successfully, that the full constraint set is repaired, or that the radiation-admissibility failure is cured. Nothing was merged into main.

## 9. Coordination outcome and updated next gate

Coordination thread: [Issue #5 — GEAR dynamic-interface hypothesis and current diagnostics](https://github.com/miahpenn/0-Cosmos-engine/issues/5).

The parallel CMC agent reported a bounded frozen uniform-quadrature CMC comparison, run [38301808650](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38301808650). Its workflow **failed closed at the source-provenance guard** because the isolated branch’s checked-out source did not equal the declared frozen production SHA 6d8b3f2790644591326df57c55908863f994cb0c; the scientific comparison steps were skipped. This is a provenance refusal, not a CMC or physics result. Reconcile source/vendor provenance with that agent before re-running it; do not interpret the workflow’s failure as a physical failure.

### Recommended bounded next experiment

The source-pair regression is now complete for both energy and radial-current components on one synthetic state. The next discriminating step is a stage-local source-only comparison on already admitted, frozen production states at the real beginning, predictor and corrector levels. Reuse captured states/artifacts where possible; do not launch a long trajectory.

Record for each stage and cell: alpha, rho_DM, Pi, radial phi gradient, scalar source increments, DM conservative energy and radial-momentum source increments, signed cancellation residuals, and how those increments enter the discrete stress/current and constraint budgets. Distinguish continuum/source algebra from the actual staged update. Then:

1. If the source pair remains pointwise balanced at every actual stage, preserve that result and continue the centre/RHS compatibility and radiation investigations separately.
2. If it fails only after stage mixing/time-level mismatch, isolate that staging mechanism before making an equation change.
3. Only after the frozen-stage budget is clear should a bounded one-step A/B be considered; require evaluation of the full implicated H/M/C, centre-series, radiation cone and Misner–Sharp witnesses.

Do not merge either candidate branch, wire the independent homogeneous CosmosState into production, or add a new phenomenological Q law. The current-main beta-source candidate changes only the two derived lapse factors and the test; the separate scalar inverse-metric candidate and the dominant centre/RHS defect remain separate work.

**Updated bottom line:** two concrete spatial scalar–DM source omissions have now been identified—the scalar lapse weight and the DM radial-source lapse weight—and both are corrected on the isolated current-main diagnostic branch. The source-pair and full unit-test gates pass there, but the full evolution impact remains unmeasured. The strongest open numerical question is still whether the complete staged geometry, CMC gauge, scalar and conservative-matter updates preserve all required constraints and admissibility conditions. The original GEAR local-to-global-Q, D-to-matter, clock-ownership and global anisotropic-geometry questions remain separate physical-model gates.
