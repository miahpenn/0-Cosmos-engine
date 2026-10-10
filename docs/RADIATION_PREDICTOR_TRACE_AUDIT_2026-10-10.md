# Radiation predictor trace audit — 2026-10-10

## Run identity

- Repository: `miahpenn/0-Cosmos-engine`
- Physics branch: `physics/spatial-beta-covariant-source-closure`
- Trace workflow: [Run #7](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38068938182)
- Repository CI: [Run #576](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38068938171)
- Trigger commit: `c8a9967126ba07cd12e4da9c86b7066a703d3419`
- Production-source baseline before trace instrumentation: `9eed1cf95c52f001c2b87afbbe6da37faa56d50e`
- Trace/test instrumentation commit: `527b44931b059f4d118a86c898b91d79ca4d96f2`
- Workflow commit: `ed0bc6f61103f2439556b27ad588116a149f4add`
- Pinned vendor submodule: `d6052d605673ce9d82cdc99b0df79855fa2ea215`
- GitHub artifact: `0star-radiation-stage-trace`, artifact ID `11676486591`
- Artifact ZIP SHA-256: `bef03745ad4ea1ab62d8a42de9b957c57210571d933cf2728c96a6196b970c43`
- `trace.json` SHA-256: `642a82b8b97ddc849b08f70b4379672b08e81f842d3ae8bdeaaa1d5374afd8a3`
- `summary.json` SHA-256: `f36105276c2ef4fc7ee91ed8af5593792cf05ee27324149a3e17e4cca87fee22`

## Outcome

Both GitHub workflows concluded `success`. This means CI passed and the diagnostic successfully captured and uploaded the expected failure; it does **not** mean the radiation trajectory reached its target.

Configuration: N=80, r_max=80, dr=1, CFL=0.03, dt=0.03, target coordinate time 50, radiation enabled, D amplitude 1e-10. The existing production light-constraint boundary and radiation outer closure were left unchanged.

The run reached accepted coordinate time t=46.98000000000102 and central proper time tau=22.016115393070738 after 1,566 completed steps. The next attempted predictor at t=47.01000000000102 failed at stage `explicit_predictor_first_CMC_lapse`, outermost cell i=79, r=79.5.

- Last accepted physical radiation cone ratio: |S|/E = 0.9983867169143378.
- Last accepted cone margin: C = +4.6971156703181716e-7.
- Predictor evaluated with predictor metric: |S|/E = 1.00005699412463422.
- Predictor cone margin: C = -1.6900205171508503e-8.
- The predictor exceeds the cone by only about 2.27e-12 in the reported physical E and |S|, but the realizability condition is genuinely violated.

Do not treat t=47.01 as an accepted trajectory point: the failure occurs while constructing its predictor.

## Ordered cone-margin budget at the terminal cell

The tracer evaluates the flux-only and source-updated conservative states on the accepted metric, then evaluates the same predictor conservative state on the predictor metric.

| Stage | C | Increment |
|---|---:|---:|
| Accepted state, accepted metric | +4.69711567e-7 | — |
| Flux-only update, accepted metric | +1.48122402e-6 | +1.01151245e-6 |
| Then geometric sources, accepted metric | +2.21713013e-6 | +7.35906106e-7 |
| Predictor conservative state, predictor metric | -1.69002052e-8 | -2.23403033e-6 |
| Net accepted-to-predictor change | — | -4.86611772e-7 |

The budget closure error is zero at recorded precision. Predictor conservative variables match the RHS reconstruction exactly at recorded precision; source capture is present for every traced cell, and the instrumentation error count is zero.

**Strongest supported proximate finding:** the flux-only and geometric-source updates, considered sequentially on the accepted metric, increase the cone margin. Evaluating the predictor state on the evolved predictor metric then reduces that margin by 2.2340e-6 and crosses the cone. This associates the immediate crossing with the predictor metric/geometry stage rather than a standalone flux-only overshoot. It is not yet proof of which geometry RHS, outer closure, CMC-lapse interaction, or stage-coupling term is ultimately defective.

At cell 79, accepted-to-predictor values include:
- a: 0.60895534 -> 0.60320350
- b: 1.28146656 -> 1.28756177
- X: 0.94389565 -> 0.94655832
- sqrt(gamma): 7515.5765 -> 7452.3308
- inverse radial metric gamma_rr: 1.46306132 -> 1.48535718

The inversion exception reports alpha=1 and beta=0 at the outer cell. The trace did not alter the boundary normalization or conservative radiation variables.

## Approach to failure

| Accepted t | Central proper time tau | Outer |S|/E | Outer C |
|---:|---:|---:|---:|
| 46.92 | 21.99851536 | 0.99536187 | 1.29832904e-6 |
| 46.95 | 22.00732617 | 0.99684596 | 9.00958790e-7 |
| 46.98 | 22.01611539 | 0.99838672 | 4.69711567e-7 |

The accepted outer-cell margin was already shrinking over the last three accepted steps; the predictor metric stage pushes it just beyond the realizability boundary.

## Integrity and diagnostic coverage

The artifact records the branch, trigger, source baseline, instrumentation commit, workflow commit, vendor submodule SHA, runtime environment, run ID and SHA-256 fingerprints for the production and diagnostic files. It contains 394 stage snapshots, 66 predictor budgets, 1,566 progress samples, and zero instrumentation errors. The downloaded ZIP's SHA-256 matches the digest reported by GitHub; the embedded trace hash matches `SHA256SUMS.txt`.

No clipping, floors, damping, source insertion, alternative metric, projection, or timestep adjustment was applied. No production physics changes were included in the trace instrumentation.

## What remains unproven

This run does not yet establish the ultimate cause of the predictor metric change, its convergence with resolution/domain, whether the covariant-source changes caused or merely coexist with the failure, or any full-cycle/bounce conclusion.

## Next work — algebraic audit before another expensive trajectory

1. Audit the exact predictor metric update at cell 79 and corresponding radial geometry RHS, including outer-boundary closure and CMC solve inputs/outputs.
2. Independently reconstruct the physical radiation cone invariant from conservative variables, metric and lapse at accepted and predictor stages; compare with the production inversion formula and verify metric dependence/sign conventions.
3. Convert those checks into focused, deterministic regression tests using captured terminal values; no new trajectory is needed for these steps.
4. Only after those gates pass, choose one controlled follow-up run to discriminate among surviving causes. Keep it on the patched branch and leave boundary physics unchanged unless a separate derivation supports a change.

## Handoff

- Physics branch remains isolated from `main`.
- Diagnostic capture is complete and source-pinned.
- No second expensive run has been launched.
- Current conclusion: **late outer-cell radiation realizability failure; the terminal predictor's metric update is the strongest measured immediate trigger. Ultimate root cause remains open.**
