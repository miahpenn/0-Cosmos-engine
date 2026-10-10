# GEAR-to-kernel architecture and integration gates

**Status:** source-grounded architecture gate; not a physics approval, not a coupling implementation, and not permission to merge to `main`  
**Date:** 2026-10-10  
**Research branch:** `research/gear-kernel-ownership-contract`  
**Production source snapshot reviewed:** `0star-central-clock`, commit `409c4b9e4374cc05ffea0086efef613044028c2f`  
**Companion draft:** [GEAR_TO_KERNEL_STATE_OWNERSHIP_CONTRACT_DRAFT_2026-10-10.md](GEAR_TO_KERNEL_STATE_OWNERSHIP_CONTRACT_DRAFT_2026-10-10.md)

## 1. Decision for the current work

The source evidence requires us to keep three configurations distinct:

1. **Canonical GEAR archive reference — two independently evolved lanes.** The preserved `CYCLIC-COSMOS ENGINE-.docx` runner calls the local sector, runs homogeneous cosmic expansion and contraction, and compares their outputs. The lanes use their natural clocks. The archive explicitly leaves unsupported inter-lane bridges open. This is the reference baseline and must not be silently rewritten.
2. **V5.5 spatial production kernel — one radially resolved evolution.** Its `ProductionState` owns spatial geometry, spatial scalar fields, conservative matter species, coordinate time, accumulated central proper time, and diagnostic ledgers. Geometry, scalars, and species interact within that solve. Its state does not contain a homogeneous `CosmosState`.
3. **A future matched background-plus-local model — not yet derived or implemented.** It would need explicit matching, clock, source, and conservation equations. It is not created merely by putting both state objects in a container or calling the current geometry-driven helper.

This classification settles the near-term operating mode, not the ultimate physical architecture: **preserve the archive reference; diagnose the V5.5 spatial solve on its own terms; do not wire the homogeneous lane into production until the governing design supplies the missing ownership/matching equations.** It also avoids pretending the present spatial solver alone has solved every global GEAR degree of freedom.

## 2. Provenance and scope

### Authoritative archive

The 2 October master archive names `CYCLIC-COSMOS ENGINE-.docx` as the canonical baseline and says not to rewrite/reconnect it without a new version. Its archive snapshot describes the local GEAR-03 S/D modes and GEAR-63/65 barrier on one lane, and the homogeneous scalar/DM/baryon/radiation/shear expansion–contraction sector on the other. It explicitly states that the cosmological lane is independently evolved rather than forced to consume D-mode energy.

The archive leaves these bridges open:

- local proper time ↔ cosmic time;
- D-mode energy → localized ordinary matter;
- local source/current → global cosmological source (Q);
- full-GR oscillaton endpoint;
- first-principles bounce.

REV16 remains the authority for the 0* criterion and its status. A kernel code test, a trajectory that reaches its numerical target, or an interface test does not change REV16's open global gates or the result (0^*_{\mathrm{eff}}=\mathrm{FALSE}).

### V5.5 source

The source files fetched from `0star-central-clock` include:

- `engine/production_kernel.py` — blob `d9a4740592f2521f7dc04d66a718a854c70a8cde`;
- `engine/cosmos.py` — blob `465f38653ae864455987f3a82583b651c3f07faa`;
- `engine/coupled.py` — blob `a584c4058637801e2765b15470fb592087f607bf`;
- `engine/run_production.py` — blob `6bfcfb73ddb4d1299e88c80c36d4f6fde8bbf9f9`.

The kernel itself labels its validation state `implemented_not_campaign_validated`. Passing component tests must not be promoted to a validated full-machine claim.

### Existing radiation diagnostic — separate evidence stream

Workflow run [38028306184](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38028306184) finished with its diagnostic job successful and produced artifact `0star-radiation-proper-time-domain-control`. The downloaded ZIP SHA-256 was `f4f06a54bf66f0e5eed597efe0b6b0e605e31375edd7303bcb4edd4cdb97fac5`; the embedded `trace.json` SHA-256 matches the artifact's `SHA256SUMS.txt) and summary record: `f418a023c44dcd51802696e6468e3e09691d01ea44ecbb5a9463881be6b14e76`. The trace recorded 956 stage snapshots, 152 predictor budgets, 88 completed-state budgets (zero budget errors), 153 geometry-predictor budgets (zero errors), 174 boundary snapshots (zero errors), and zero instrumentation errors.

Its captured terminal state is a radiation-admissibility timestep underflow at coordinate time (t=61.046503364247606), central proper time (\tau=21.9130694389674), cell (i=119), (r=119.5). This is a successfully captured diagnostic failure, **not** a successful physical trajectory and not evidence that the homogeneous COSMOS lane is coupled. Its conclusions remain scoped to the specific diagnostic source/configuration in that artifact.

## 3. State ownership as currently implemented

| State / quantity | Current owner and role | Status and boundary |
|---|---|---|
| Radial geometry ((a,b,X,A_a,K,\Lambda,\alpha,\beta,ldots)) | `ProductionState.geometry` | Evolved in the V5.5 spatial solve. |
| Local S/D and spatial φ fields | `ProductionState.scalars` (`ScalarFields`) | Radial fields in the spatial evolution. The S/D mode is not, by that fact alone, ordinary baryonic matter. |
| Dark matter, baryons, radiation | `ProductionState.matter` (`V55MatterState`) | Conservative species evolved on the spatial metric. The existing spatial φ–DM beta coupling is an equation-level interaction; it is not a D-to-matter conversion law. |
| Coordinate time (t) and central proper time (\tau) | `ProductionState.t`, `ProductionState.tau` | Both are recorded, but their presence does not define an arbitrary conversion between clocks. |
| Homogeneous (a,H,φ,π_φ,\rho_{DM},\rho_b,\rho_r,\rho_{shear}) | `CosmosState` in `engine/cosmos.py` | Independent homogeneous state. It is not a member of the live V5.5 `ProductionState` and is not advanced in its `step()` path. |
| `CoupledState` | `engine/coupled.py` container of Cosmos, Local, Interface states | Construction/serialization is not proof of a conserved, active coupled evolution. |
| (H_{eff}), Misner–Sharp flux/work, cycle and handoff events | Production diagnostics / ledgers | Observables and records. An observable does not become a receiving evolution equation simply because another helper can read it. |
| Homogeneous shear | `CosmosState.rho_shear), with its archived (a^{-6}) dilution/effective pressure | It is explicitly not inserted as a homogeneous shear-pressure term into the exact spherical local stress tensor. Global anisotropic curvature remains a separate global-geometry issue. |

The local spatial φ field and homogeneous `CosmosState.phi` have no approved identity/matching rule. They may not be blindly equated, copied, or evolved as duplicate physical copies. The architecture decision must state which is authoritative or provide a derived background/perturbation/matching relation.

## 4. What the existing homogeneous helper does — and does not do

The homogeneous equations in `engine/cosmos.py` include:

- (ẏa=Ha);
- signed (H) evolution in the standalone homogeneous RHS;
- scalar/DM beta exchange;
- dilution terms for baryons, radiation, and effective shear.

`cosmos_rhs()` rejects nonzero `exchange_Q`. `cosmos_rhs_driven()` and `geometry_driven_step()` advance the homogeneous non-(H) state using externally supplied, solved (H_{eff}) values. The helper does not return homogeneous stress-energy to the spatial geometry, provide a D-to-matter law, or close local-source-to-global-(Q) accounting. The reviewed production `step()`/runner call graph does not invoke `geometry_driven_step()`.

Therefore:

- a geometric (H_{eff}) handoff is a one-way driver, not reciprocal source coupling;
- an active local spatial solve is not automatically the archive's two-lane reference;
- a combined state container is not an integrated physical machine;
- adding a call to the helper would be insufficient and is not authorized by this contract.

## 5. Coupling-edge contract

No cross-lane edge is promoted by this document. Every proposed edge must be written down with all of the following fields before implementation:

| Required field | Question it must answer |
|---|---|
| Source state and owner | Which unique state owns the source variable? |
| Governing equation | Which archived/derived equation defines the source term? |
| Inputs and units | Which quantities enter, in what units and frame? |
| Recipient state and owner | Which equation consumes the returned quantity? |
| Clock / time level | Is it evaluated at accepted, predictor, or corrector state? Is time (t), proper time (\tau), or another derived parameter used? |
| Conservation identity | What source/recipient balance must cancel or remain in the total ledger? |
| Constraint set | Which Hamiltonian, momentum, constitutive, cone, or geometric constraints must remain satisfied? |
| Evidence and falsifier | What trace/test would demonstrate the edge is active, and what result would disprove the proposed closure? |
| Baseline/control | Which exact source commit, artifact and control establish the comparison? |

### Required open edges

1. **Local proper time ↔ cosmic time.** Derive a parameter-free relation or define a single-spacetime time coordinate and its proper-time observables. No fitted scale or convenient identification (\tau=t).
2. **D-mode energy → localized ordinary matter.** Derive the identification/conversion law. D amplitude, D flux, or generated-energy proxy cannot simply be relabelled as baryon density.
3. **Local work/current → global (Q).** GEAR-66 supplies an evidentiary accounting route: local force → dissipated power → event energy → event/ejection rate → (Q). Its reduced event-energy result is not proof of a first-principles force/carrier law, node abundance, ejection history, or a normalization-free (Q(a)) match. Do not set (Q=F_{local}), add an arbitrary coefficient, or guess a recipient.
4. **Spatial φ ↔ homogeneous φ.** Establish ownership and matching, including the relation of a radial solution to a homogeneous observable, if such a relation is intended.
5. **Homogeneous shear / global anisotropic curvature ↔ spatial geometry.** Derive the representation. Do not inject the homogeneous (a^{-6}) effective shear term into the spherical stress tensor by convenience.
6. **Cycle/turnaround interpretation.** Require event rows and neighboring states from the same continuous run, with signs, (t), (\tau), named expansion measure and constraints. (H_{eff}) sign-crossings alone do not establish a scale-factor turnaround or completed physical cycle.

## 6. Numerical and conservation gates, in order

### Gate 0 — provenance and immutability

- Preserve the canonical GEAR document and archive files byte-for-byte.
- Pin branch, commit, configuration, source hashes and artifact hashes for every comparison.
- Keep the radiation audit and spatial-kernel campaigns as separate evidence streams.
- No long run is launched while the equation/ownership gate is open.

### Gate 1 — state ownership and accidental-coupling guard

Add deterministic tests that identify the active state owners and call graph. They must ensure that the archive reference remains independently evolved and that the V5.5 production path does not accidentally acquire a second copy of homogeneous φ, (a), or (H). These are scope guards, not a claim that non-coupling is the final physics.

### Gate 2 — equation-level source accounting

For the archive-locked homogeneous scalar/DM pair, the equations imply the exact exchange identity
[
\dot\rho_\phi+3H(\rho_\phi+p_\phi)=+\beta\rho_{DM}\pi_\phi,\qquad
\dot\rho_{DM}+3H\rho_{DM}=-\beta\rho_{DM}\pi_\phi.
]
Their sum must cancel the internal beta-exchange source. Test the symbolic/source identity and a numerical residual at the order justified by the integrator. This validates the pair inside the homogeneous lane; it does **not** connect that lane to D-mode energy or to the spatial production state.

Separately, verify the spatial implementation's source signs and conservative-variable accounting from its own RHS at matched metric/time levels. Do not infer a spatial conservation result merely from the homogeneous equations or identical parameter names.

### Gate 3 — local geometry and stage admissibility

- Verify exact algebraic/center regularity identities independently from evolution residuals.
- Record the accepted state and exact predictor/corrector RHS at the same stage/time level.
- Check Hamiltonian and momentum constraints and the existing Misner–Sharp flux/work balance with convergence-based tolerances.
- For radiation, reconstruct (C=U_E-\sqrt{\gamma^{rr}}|U_r|) on the exact stage metric; record flux, source, and metric contributions and a closure residual. An inadmissible stage is observed and captured, never clipped or repaired by a diagnostic.

### Gate 4 — derive a receiving equation before adding any interface code

For every cross-lane edge, derive both source and recipient sides and specify an independent conservation identity. In particular, the GEAR-66 (Q(a)) comparison must not have a fitted final normalization. Unknown carrier identity, event population/rate, or clock matching remains an explicit blocker rather than a guessed parameter.

### Gate 5 — deterministic integration test on a research branch

Only after Gates 1–4 are documented and individually pass:

- exercise one small, reproducible step sequence;
- prove that every declared edge is actually invoked and that source and recipient states respond at the intended time level;
- reconcile the transfer in the total ledger with no double count;
- verify the appropriate constraints and the separate homogeneous/spatial field ownership rule;
- compare against the exact independent-lane and spatial baselines.

A pair of trajectories plotted together, a populated `CoupledState`, or a one-way (H_{eff}) driver is not an integration pass.

### Gate 6 — resolution and full-machine admission

First demonstrate the predicted convergence of source accounting and constraints across predeclared resolutions. Tolerances must follow the method's truncation order and algebraic precision; do not pick a tolerance after seeing a desired answer. Only then preregister and run one full controlled campaign from the approved commit, with complete checkpoints, ledgers, hashes and fail-closed admission. N=320 and other expensive campaigns remain blocked until their own existing gates are satisfied.

## 7. Prohibited shortcuts

- No arbitrary/fitted local-to-cosmic clock conversion.
- No guessed (Q), no instantaneous-force-as-(Q) assumption, and no fitted final normalization.
- No invented D-to-matter or radiation-carrier law.
- No silent identification of spatial φ and homogeneous φ.
- No inserting homogeneous shear into the spherical local stress-energy without derivation.
- No clipping, radiation floors, artificial damping, branch resets, manufactured bounce, or hand-forced sign changes.
- No edits to the canonical archive; no merge of this documentation branch to `main`; no change to REV16 or existing campaign criteria.

## 8. Initial deterministic gates — executed 2026-10-10

The first gate set is now implemented on the isolated `research/gear-kernel-integration-gates` branch. No production physics file or archive baseline was changed.

- New `tests/test_architecture_gates.py`: checks homogeneous scalar–DM continuity-source cancellation at the initialized archive state, confirms unpromoted nonzero (Q) fails closed, and uses an AST/call-graph scope guard to ensure the live V5.5 production step/runner do not silently call the homogeneous COSMOS RHS/helper or introduce `CosmosState` as a production call.
- Existing spatial checks were included in the same suite: `tests/test_center_regularity.py` tests initialized/one-step center regularity and determinant identities; `tests/test_outer_constraint_audit.py` checks the finite full-grid constraint profile, Hamiltonian reconstruction/decomposition and regional coverage; `tests/test_v55_reference.py` checks finite reference initial constraints; the remaining source/geometry tests also ran.
- GitHub Actions run [38059261015](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38059261015) passed **110 tests in 32.14 s** at commit `c0e096d68b14a929a0df37676102ced82a720db9`. The normal test job succeeded; both optional long-campaign jobs were skipped.
- This is a source-scope / deterministic component-gate pass, not a full-machine coupling pass, convergence proof, or campaign admission.

## 9. Source-accounting review — source fingerprint complete; physical law still open

The source-only active-RHS witness is recorded in [SPATIAL_BETA_EXCHANGE_LAPSE_ACCOUNTING_AUDIT_2026-10-10.md](SPATIAL_BETA_EXCHANGE_LAPSE_ACCOUNTING_AUDIT_2026-10-10.md) and implemented in [tests/test_spatial_beta_exchange_audit.py](https://github.com/miahpenn/0-Cosmos-engine/blob/research/gear-kernel-integration-gates/tests/test_spatial_beta_exchange_audit.py). It confirms that the current scalar RHS source adds +beta*rho_DM without alpha, while the live DM conservative energy source reduces to -sqrt(gamma)*alpha*beta*rho_DM*Pi. On a synthetic, non-evolved source probe with non-unit lapse and nonzero shift, the source-pair residual matches sqrt(gamma)*(1-alpha)*beta*rho_DM*Pi.

The normal suite passed **111 tests in 30.22 s** at [run 38059772092](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38059772092); the optional long-campaign jobs were skipped. This confirms test repeatability and the current source fingerprint, **not** physical exchange conservation.

The first October 5 V5.5 inhomogeneous gate omitted the beta scalar-DM interaction pending a local momentum equation. The later same-date `GEAR_V5_5_STRONGFIELD_CLOSURE_AUDIT_2026-10-05.md`, Sections 8–9, locks the covariant source signs for signature (-,+,+,+) and specifies `Q_DM^nu = + beta rho_DM nabla^nu phi`, with the opposite scalar source. So the archive does supply the covariant interaction law; the remaining requirement is to derive its exact 3+1 coordinate-time implementation (lapse, shift, Pi convention, sign, and density normalization). The current source witness shows the present code's algebra, but is not itself a physical conservation pass. Do not patch the alpha factor by analogy alone. First write and test the 3+1 projection of the locked covariant equation; if that confirms the missing lapse, make the minimal correction on a separate physics branch, then rerun the deterministic gates.

Keep this separate from the still-open cross-lane local-source-to-global-Q bridge. No homogeneous helper wiring, other production-equation edit, or long trajectory is authorized until the source law and its paired balance are resolved.

