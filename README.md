# 0-Cosmos Engine

Canonical repository for the coupled GEAR/COSMOS research engine. The goal is a reproducible, production-quality numerical machine—not a collection of disconnected demonstrations.

## Start here

1. **[Current engineering status](docs/STATUS.md)** — current production/source anchors, software validation, open scientific-admission gates and evidence provenance.
2. **[Research workflow](docs/RESEARCH_WORKFLOW.md)** — canonical claim → overlap check → bounded work → exact-head evidence → handoff procedure.
3. **[SI-agent communications protocol](docs/SI_AGENT_FIRST_READ_COMMS_PROTOCOL.md)** — required onboarding and scalable reading strategy for the long coordination thread.
4. **[Cross-agent board / comms](https://github.com/miahpenn/0-Cosmos-engine/issues/5)** — ownership claims, current coordination, corrections and durable findings.
5. **[Production-readiness task board](https://github.com/miahpenn/0-Cosmos-engine/issues/135)** — canonical task inventory, dependencies, PR state and exit criteria.
6. **[Main-branch governance](https://github.com/miahpenn/0-Cosmos-engine/issues/112)** and [momentum-constraint propagation](https://github.com/miahpenn/0-Cosmos-engine/issues/4) — dedicated blocker threads.

Use the status document for the verified engineering snapshot; use Issue #135 for task tracking; use Issue #5 for team coordination. Avoid creating parallel status lists or duplicate diagnostics when an existing task/artifact covers the same question.

## Physics and evidence guardrails

- Preserve archive-supported GEAR/COSMOS conventions, signs and normalization. Treat unresolved bridges as open questions, not permission to invent a coupling law.
- No fitted feedback coefficients, arbitrary damping, floors/clamps, radiation-cone relaxation, forced bounce, imposed clock conversion or unsupported D-to-matter law.
- The exact spherical branch does not by itself represent the archive's full Bianchi-I shear.
- Keep production physics, diagnostics, test harnesses and campaign orchestration separate.
- A green unit suite, symbolic identity, frozen-state counterfactual or short smoke does not establish trajectory stability or scientific admission.
- Pin every result to the exact source SHA, workflow and artifact. Preserve negative results and update the shared board when evidence is superseded.

## Code map

### Matter, scalar fields and stress-energy
- `engine/scalar_system.py` — S/D/COSMOS scalar evolution and projections.
- `engine/cosmos.py` — archive homogeneous COSMOS lane, radiation and shear semantics.
- `engine/local.py` — GEAR local S/D sector and geometric barrier.
- `engine/valencia.py` — conservative/primitive identities and tensor-source contraction.
- `engine/matter_system.py` — conservative spherical dust/radiation transport.
- `engine/v55_matter.py` — V5.5 matter bundle and operating point.
- `engine/stress_energy.py` — total Einstein-source assembly.

### Geometry, evolution and observables
- `vendor/bb-palatini-unified-r0` — pinned reference-metric spherical BSSN/PIRK kernel.
- `engine/v55_pirk_adapter.py` and `engine/production_kernel.py` — geometry adapter and coupled state evolution.
- `engine/v55_initial.py` — production initial-data mapping.
- `engine/production_contract.py` — production capability contract.
- `engine/worldtube.py` — areal radius, Misner–Sharp mass and current diagnostics.
- `engine/invariant_diagnostics.py` — numerical invariant witnesses.
- `engine/run_controls.py` — low-cost controls.
- `engine/run_production.py` — canonical production campaign runner.

Historical/reference modules remain for provenance; consult the status and test contracts before treating one as an active production entry point.

## Current release posture

**Scientific admission: NOT ADMITTED for production campaigns** until the current gates in [docs/STATUS.md](docs/STATUS.md) are closed. The bounded software checks and one-step smoke are necessary software evidence only; campaign continuation, convergence, centre/constraint propagation and fixed-areal-radius conservation must be evaluated independently.

## Contribution workflow

- Start from current `main`; use a focused branch and pull request.
- Before work, claim the task in Issue #5 with lane, branch, exact source SHA, proposed bounded action and an overlap check.
- Inspect existing PR discussions, checks and artifacts before re-running tests.
- Run the applicable exact-head checks after the final change. Do not report a workflow as passing if it did not run on that SHA.
- Update Issue #135 and the relevant shared thread with result, limits and the next gate.
