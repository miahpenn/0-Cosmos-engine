# Uniform-quadrature CMC centre-row short-trajectory comparison

**Date:** 13 October 2026 UTC  
**Workflow:** [38327936024 — SUCCESS](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38327936024)  
**Diagnostic head:** `c05dfd7c2959dbcf05beb3f56ae49908a278b2cc`  
**Production source:** `6d8b3f2790644591326df57c55908863f994cb0c`  
**Pinned vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`  
**Artifact:** [11781256471](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38327936024/artifacts/11781256471)  
**Artifact ZIP SHA-256:** `92d060f3e44fd87e33ce9b4be9d4cf0c222f349f178429d8e4672521865fa25a`  
**Report JSON SHA-256:** `b715f8e74c8ded8a521ba889ca71fecc87977ae27d356662c468b542c44844b2`

## Test definition

The diagnostic applies the current, second-order even-ghost, and direct-fourth-order even-parity centre rows to identical reconstructed uniform-quadrature initial profiles and runs **four real production steps per arm** at N=80/120/160, with `dt=0.03*dr`. The selected row is reused at initialization and every internal CMC lapse solve in each step. This is 9 arms × 4 ordinary steps = **36 completed steps**, with 17 CMC solves per arm (initializer plus four solves at each of four steps).

The source and vendor identities were verified in CI; the physical source/equations, matter transport, PIRK sequence, outer boundary, and real stage/admissibility checks were unchanged. The current row remains as a counterfactual control. No clipping, damping, projection, fit, retry, source change, or physical parameter edit was introduced.

## Final-step comparison

Values are from the final completed slice of each four-step arm. `Aa/r₀²` and `max |K_rhs−target|` are captured at the final internal CMC solve using the same stage state; the latter is the maximum over cells 0..N−2.

| N | Centre row | Final t | `Aa_rhs/r₀²` at centre | Max `|K_rhs−target|` | Full `|H|` max | Full `|M|` max | Full `|C_lambda|` max | `|S|/E` outer |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 80 | Current | 0.12 | +1.669e−4 | 1.464e−3 | 5.226e−5 | 2.346e−5 | 2.913e−6 | 4.973e−4 |
| 80 | Even ghost | 0.12 | −1.420e−5 | 2.311e−5 | 3.591e−5 | 1.651e−6 | 2.014e−7 | 4.973e−4 |
| 80 | Direct D4 | 0.12 | −1.138e−5 | 2.073e−5 | 3.616e−5 | 1.326e−6 | 1.607e−7 | 4.973e−4 |
| 120 | Current | 0.08 | +3.804e−4 | 1.533e−3 | 2.627e−5 | 2.392e−5 | 1.972e−6 | 3.302e−4 |
| 120 | Even ghost | 0.08 | −4.137e−5 | 1.092e−5 | 1.048e−5 | 1.272e−6 | 1.063e−7 | 3.302e−4 |
| 120 | Direct D4 | 0.08 | −3.839e−5 | 1.053e−5 | 1.038e−5 | 1.422e−6 | 1.188e−7 | 3.302e−4 |
| 160 | Current | 0.06 | +7.073e−4 | 1.558e−3 | 2.075e−5 | 2.393e−5 | 1.477e−6 | 2.471e−4 |
| 160 | Even ghost | 0.06 | −5.152e−5 | 6.272e−6 | 6.820e−6 | 1.134e−6 | 7.085e−8 | 2.471e−4 |
| 160 | Direct D4 | 0.06 | −4.848e−5 | 6.205e−6 | 6.760e−6 | 1.220e−6 | 7.617e−8 | 2.471e−4 |

## What this establishes

1. **The frozen-slice signal persists over four real steps.** At N=160 the current row retains a `K_rhs−target` mismatch of roughly `1.56e−3`; the two alternatives remain near `6.2e−6`. Their centre `Aa_rhs/r₀²` witnesses remain around `−5e−5`, while the current-row witness grows to `+7.1e−4`.
2. **Momentum/connection behaviour improves materially.** At N=160 the final full-domain momentum maximum is `2.39e−5` with the current row, versus `1.13e−6` even-ghost and `1.22e−6` direct-D4. The connection maximum is `1.48e−6` versus `7.1–7.6e−8`. Hamiltonian full-domain maxima also rise in the current control over these four steps but stay near the initial order in both candidates.
3. **No arm encountered a step/stage rejection in this bounded test.** Radiation outer `|S|/E` rose gradually with time (e.g. `2.47e−4` at N=160 after four steps), and all stage admissibility checks passed. Same-resolution Misner–Sharp mass values remained essentially identical between row variants during this short interval; this is not a full flux-balance audit.

## Limits and next gate

These are strong diagnostic results for the centre-row / CMC-series defect, not production admission. Four steps are too short to establish long-horizon stability, convergence, Hamiltonian propagation, complete Misner–Sharp conservation, or behaviour near the known late-time radiation boundary failure. The even-ghost and direct-D4 options both remain candidates; neither is promoted by this run. The direct-D4 centre-row residual closes to roundoff partly by construction, so its apparent row closure is not treated as an independent pass criterion.

Next gate: retain the current-row control and run the same three arms over a longer but still bounded trajectory only after deciding explicit numerical acceptance witnesses in advance. Continue recording H/M/C_lambda/determinant, signed centre-series rates, actual CMC target and row residuals, radiation admissibility at every stage, proper time, and Misner–Sharp current/flux/work closure. Keep production main and physics locks untouched until multi-resolution and conservation evidence supports a change.
