# Cross-Agent Radiation Retry Findings and Handoff

**Updated:** 12 October 2026 UTC
**Classification:** Diagnostic-only; production admission remains open.

## Current disposition

The real adaptive retry path has now been exercised against the archived accepted numerical state immediately preceding a known radiation-cone rejection. The nominal step is rejected by the existing exact admissibility guard; one half-step is accepted without clipping matter, weakening the cone inequality, or changing physics sources. The retry helper also leaves the frozen input arrays unchanged.

This is local software/evolution-step evidence only. It does not establish a stable multi-step trajectory, full-campaign admission, constraint propagation, Misner–Sharp balance, resolution convergence, or a physical bounce/cycle.

## Evidence records

1. **Signed frozen-state timestep budget:** [report](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/radiation-frozen-state-stage-budget-2026-10-12/docs/RADIATION_FROZEN_STATE_DT_STAGE_BUDGET_2026-10-12.md). Workflow [38324794759](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38324794759), 91 repository-wide tests passed. Exact accepted-state NPZ SHA-256: `ccd1fa0fbabe65c51c3dd7ac2c9227d6d25aa680c169eef1afd3fca4326294cf`.
2. **Real retry from frozen failure state:** [report](https://github.com/miahpenn/0-Cosmos-engine/blob/fix/radiation-adaptive-step-on-cone-2026-10-12/docs/RADIATION_FROZEN_STATE_REAL_RETRY_2026-10-13.md). Workflow [38325032498](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38325032498), exact tested candidate head `651dc124484e7c056b784033b69bedba6da0b5a5`; 95 repository-wide tests passed. Artifact [11781015455](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38325032498/artifacts/11781015455).

## Numerical findings

- Frozen accepted state starts at coordinate time approximately `t=47.13`, with `N=80`, `r_max=80`, and requested `dt=0.03`.
- Actual production attempt at `dt=0.03` rejects at `stage='explicit-predictor'`, cell `i=79` / `r=79.5`, with `|S|/E = 1.0004012501659771`.
- The actual retry helper accepts `dt=0.015`; final maximum full-grid `|S|/E = 0.9996001341484262`, minimum densitized cone margin `C = +2.0388462e-9`.
- The numerical field-array digest before and after the retry helper is identical: `6ec904c223640d1b43d84946b54b909032705705727bc395dd1b87c774ccedf9`.
- Independent frozen-state arms at `dt=0.015`, `0.0075`, and `0.00375` each complete one real production step. Signed stage accounting localizes the nominal negative cone margin to evaluating the predicted conservative radiation state on the explicit predictor metric; the inverse metric `gamma^rr=X^2/a` rises through both `a` and `X` updates.
- The adaptive retry underflow logic no longer imposes an unexplained 8-ULP minimum. It stops halving only when the candidate timestep no longer advances stored floating-point coordinate time.

## Guardrails and unresolved work

- Do not alter or relax `E >= |S|`; no floors, clamps, projection rescue, damping, fitted coefficients, or boundary/source edits are justified by these tests.
- Full sourced Hamiltonian, momentum, connection, determinant, and CMC witnesses remain nonzero; Misner–Sharp change over one step is not a conservation-balance proof.
- Keep the adaptive candidate unmerged while resolving overlap with the campaign-admission work tracked by draft PR #3 on `engine/run_production.py`.
- Next bounded engineering gate: verify that the normal production campaign runner writes the actual rejected attempt and accepted retry into its usual retry ledger, starting from the same pinned archived state. Then, only if that passes, consider one preregistered short multi-step retry test. Do not launch a costly full campaign until diagnostics and admission criteria are explicit.
- No production `main`, numerical vendor pin, or canonical GEAR archive was changed by these diagnostic audits.

## Cross-agent communication

Research coordination issue: [#5](https://github.com/miahpenn/0-Cosmos-engine/issues/5). Use this report and both evidence reports before duplicating the retry investigation. The results close only the local real-kernel retry-mechanics gap; they do not close production readiness.
