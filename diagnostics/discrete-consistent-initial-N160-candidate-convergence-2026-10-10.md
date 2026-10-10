# N=160 candidate initial-residual convergence check

Date: 2026-10-10
Branch: `diag/discrete-consistent-initial-review`

## Purpose

Extend the diagnostic-only central Hamiltonian residual check to N=160 after the corrected N=80 baseline/candidate pair passed and confirmed its resolution metadata.

## Run

- Mode: candidate only
- Resolution: N=160
- Domain: r_max=40
- Initial-data parameters, radiation setting, CFL, and final time remain those pinned by the diagnostic workflow.
- Final time: t=3
- Compare the N=160 central residual at t=3 against the already verified N=40 and N=80 candidate results.

## Guardrails

- Diagnostic only; no production defaults or physics equations change.
- No hand tuning, damping, floors, clamps, or fitted coefficients.
- Check report metadata, solver convergence, accepted-step count, and per-step/cumulative attribution closure before interpreting the result.
- Do not admit this candidate based on the residual comparison alone.
