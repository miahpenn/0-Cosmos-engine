# Native component migration contract

**Status:** design baseline; no runtime behavior changed.  
**Baseline:** `main` at `6f63303d90d3af61867b0d29f4eaf7273725154d` (17 October 2026).  
**Purpose:** replace borrowed engine-specific runtime machinery with purpose-built `engine/` components in small, reversible increments while preserving the validated research archive as formula authority.

## Policy

- Replace borrowed **domain/runtime implementation** where a native equivalent can preserve the established equations, state ownership, stage ordering, grid conventions, boundary behavior, and diagnostics.
- Keep generic scientific infrastructure such as NumPy, SciPy, SymPy, mpmath, and pytest unless a separate, evidence-based need justifies changing it. Reimplementing general-purpose numerical libraries would increase risk without advancing engine ownership.
- Treat the existing validated research/archive as the source of truth for established formulas. Do not repeat scientific review of unchanged derivations. Apply focused implementation checks for transcription, sign/index, parity, interface, regression, provenance, and state-safety defects.
- Keep production equations, the research archive, vendor sources, and runtime behavior unchanged until a specific migration PR passes applicable software checks and its stated integration gate.
- Keep `scientific_admission=false` until the existing full-system admission criteria are met. Native ownership or green CI is not scientific admission.

## Borrowed-boundary inventory

| Boundary | Current evidence / native work | Disposition |
|---|---|---|
| Reference-metric radial grid and numerical operators | The adapter loads `r0_operators.SphericalCellGrid`; `engine/native_geometry.py` provides `NativeSphericalCellGrid` and native geometry operators. PRs #337/#339 exercise matched-grid and regional stage comparisons; #340/#341 cover frozen CMC slices and primary-L2 refinement. | **Replace incrementally.** Require matched cell centers, derivative/parity and origin behavior, fixed-areal-radius comparisons, and unchanged-state witnesses before any production caller switches. |
| Geometry RHS blocks (explicit, primary L2/L3, Lambda L2/L3) | Native RHS functions exist in `engine/native_geometry.py`. PR #340 separates geometry-only terms from primary-L3 stress and Lambda-L3 current increments. | **Native implementation exists; production substitution is not yet authorized by this document.** Reuse existing source-contract/refinement evidence; add only the missing stage/gauge contract. Keep matter increments separately owned by the local adapter. |
| PIRK2 ordering and predictor/corrector machinery | `native_pirk2_step` is available and merged in #329; #331 locks stage-order/regression controls; #339 captures actual stage snapshots. | **Keep as the native moving-puncture implementation/reference.** It is not a drop-in replacement for the production CMC stage sequence. Build and test a CMC-specific stage contract before considering production integration. |
| CMC lapse solve and gauge ownership | `engine/cmc_gauge.py` and `V55ProductionKernel` own the CMC lapse solve. Production re-solves lapse on stage-local slices and resets `beta/B); native moving-puncture evolution handles gauge variables differently. | **Keep local and explicit.** Any native integration must preserve each documented CMC solve/reset point rather than importing the vendor's moving-puncture gauge update. |
| Outer incoming-light/characteristic boundary | `V55ProductionKernel` currently imports `apply_light_constraint_boundary` from the pinned vendor module. | **Borrowed runtime boundary still present.** Do not remove or rewrite it in the first increment. Before porting, pin and document its exact source, formulas, indexing, parity/outer-face conventions, and negative controls; then create a native implementation with matched-output tests. |
| Constraint helper exposed by the adapter | `geometry_constraints` delegates to vendor `vacuum.constraints`. | **Inventory as a separate port candidate.** Establish which production callers require it and compare each returned constraint component against the native definitions before switching any caller. |
| Matter transport, scalar RHS, stress-energy and conservation ledgers | README code map places these implementations in local `engine/` modules; the adapter's matter source increments are assembled locally. | **Keep native/local ownership.** Do not route matter or scalar source ownership through a geometry vendor; preserve the existing source splits and archived normalization. |
| Generic Python numerical/test packages | `requirements.txt` lists NumPy, SciPy, SymPy, mpmath, and pytest. | **Retain by default.** These are general infrastructure dependencies, not borrowed engine equations or runtime-stage logic. Audit actual import/use before considering any dependency reduction. |

This inventory is a starting boundary map, not a claim that every vendor call has been eliminated or that a native component is already production-equivalent.

## Production CMC stage contract

The following sequence is read from `V55ProductionKernel.step`. The preflight must observe/document this contract without calling `.step`, advancing an accepted state, or adding a runtime selector.

| Stage | Input / ownership | CMC and gauge behavior | Candidate native substitution boundary |
|---|---|---|---|
| Accepted slice | Caller-owned accepted `state`; production first copies geometry to `g0`. Scalar/matter data are read as inputs. | The accepted lapse is read for the left-endpoint proper-time quadrature. | No writes to accepted state; hash accepted and copied inputs before/after the preflight. |
| Stage 0 | `g0`, `s0`, `m0`, and stage-local EM state. | Re-solve CMC lapse on `g0`; reset `beta/B=0`; validate radiation; evaluate local scalar/matter RHS and geometry blocks. | Same-input geometry RHS terms only. Keep primary-L3 stress and Lambda-L3 current increments separate from geometry-only terms. |
| Explicit predictor | Copy `g0`; advance only `a/b/X`; advance scalar/matter fields with stage-0 RHS. | Reset `beta/B`; re-solve CMC lapse on the predictor slice before predictor RHS evaluation. | Compare native explicit terms only where field and gauge semantics match; do not import moving-puncture gauge evolution. |
| Primary predictor | Apply L2 trapezoidal and stage-0 L3 updates to `Aa/K/Lambda`. | Apply the existing outer boundary after Lambda is current; re-solve CMC lapse on the conditioned primary predictor; reset `beta/B`; validate. | Native L2/L3 outputs are candidates only with same-stage inputs, same source ownership, and matching boundary treatment. |
| Stage 1 | Primary-predictor geometry plus stage-local scalar/matter/EM state. | Re-evaluate the complete split on this refreshed CMC slice. | Reuse #339's actual stage-observation approach; do not repeat its moving-puncture N=48/96/192 shadow matrix. |
| Final explicit state | Form `a/b/X` by trapezoidal pairing of stage-0 and stage-1 explicit terms; form scalar/matter trapezoidal state. | Reset `beta/B`; carry the stage-1 CMC lapse pending final implicit/boundary work. | Preserve state and source ownership; no native full-step substitution in this preflight. |
| Final implicit/boundary/CMC finish | Evaluate final L2 blocks on the final explicit state; assemble final `Aa/K/Lambda` from the existing split; apply outer boundary. | Re-solve CMC lapse on the completed final slice; reset `beta/B`. | Porting the boundary or final implicit stage is a later, separately tested increment. |

## First bounded increment: observation-only stage contract

1. Reuse #339's stage-hook pattern and #340's source decomposition; do not duplicate their existing resolution matrices or matter-source contract.
2. Record stage names, frozen inputs, CMC lapse solve/reset points, boundary placement, and source ownership. Compare only identical input slices and quantities with matching gauge conventions.
3. Report primary-L3 stress and Lambda-L3 current increments independently from geometry-only RHS outputs.
4. Fingerprint accepted and copied input states before/after. Do not invoke production `V55ProductionKernel.step` in this preflight.
5. Add no centre/outer threshold, fitted coefficient, damping, vendor edit, runtime switch, accepted-state mutation, or long campaign.
6. Stop after the contract and focused implementation checks pass. Only then propose the next native port in a separate bounded PR.

## Migration sequence after the contract

1. **CMC stage contract** — establish ownership and compatible substitution points.
2. **Native geometry RHS wiring** — one stage/block at a time; retain the existing adapter path as comparator until exact-head tests pass.
3. **Native radial grid/operator ownership** — compare cell centers, derivative/parity closures, origin behavior, and fixed-areal-radius witnesses before caller migration.
4. **Native constraint helpers** — component-wise comparisons and negative controls; do not conflate a diagnostic helper with the production boundary.
5. **Native characteristic boundary** — only after the exact pinned vendor behavior and production face conventions are traced.
6. **Retire obsolete vendor calls** — remove a borrowed path only when no production or diagnostic caller needs it, and preserve the source pin/provenance needed for reproducibility.

## Acceptance record for each replacement

Each migration PR should identify the exact base/head SHA, the borrowed entry point and its native replacement, relevant existing artifacts, the unchanged archive/source convention, same-input comparison, state immutability and regression checks, exact-head workflow results, and known limitations. Re-run only checks affected by changed sources, interfaces, inputs, or behavior.

A passing implementation contract means the component is ready for the next engineering gate; it does not by itself prove full-trajectory conservation, radiation/cone admissibility, centre/constraint convergence, fixed-areal-radius balance, full-horizon convergence, or scientific admission.

## References

- [Current engineering status](STATUS.md)
- [Research workflow](RESEARCH_WORKFLOW.md)
- [SI-agent communications protocol](SI_AGENT_FIRST_READ_COMMS_PROTOCOL.md)
- [Issue #5 cross-agent coordination](https://github.com/miahpenn/0-Cosmos-engine/issues/5)
- [PR #329 — native sequential PIRK2](https://github.com/miahpenn/0-Cosmos-engine/pull/329)
- [PR #339 — actual native PIRK2 stage snapshots](https://github.com/miahpenn/0-Cosmos-engine/pull/339)
- [PR #340 — frozen-slice RHS/source decomposition](https://github.com/miahpenn/0-Cosmos-engine/pull/340)
- [PR #341 — primary-L2 initialized-slice refinement](https://github.com/miahpenn/0-Cosmos-engine/pull/341)
