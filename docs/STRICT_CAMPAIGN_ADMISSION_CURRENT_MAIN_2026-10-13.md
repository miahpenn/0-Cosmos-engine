# Strict campaign execution admission — current-main port

**Date:** 13 October 2026 UTC  
**Scope:** campaign execution contract and result serialization only.

## What this gate establishes

- Configuration and production-capability preflight run before evolution.
- The t=0 slice is recorded; evolved state and per-step diagnostics are checked for finite values and positive geometry.
- Endpoint normalization is limited to floating-point roundoff at an already-computed endpoint; it never inserts an extra evolution step.
- Admission fails closed if requested resolutions are missing or misordered, a resolution reports a numerical failure, the requested horizon is not reached, or a computed Misner–Sharp current residual summary is non-finite.
- When fewer than three history samples exist, the residual is explicitly undefined and represented as JSON `null`. Ledger and campaign summary writes prohibit non-standard JSON `NaN`/Infinity tokens.
- A computed non-finite residual remains an admission error even though its serialized representation is normalized to `null`.

## Validation

The pull-request workflow runs the complete repository test suite and a bounded N=40, r_max=40, CFL=0.03, t=1.0 production smoke. It does not launch the t=12/t=160 campaign.

## Explicit non-claim

This is **execution-completion admission, not physical/scientific admission**. Passing software contracts and the short smoke does not establish Hamiltonian/momentum/connection convergence, radiation continuation, integrated Misner–Sharp conservation, turnaround/re-expansion, or production readiness. Existing evidence still shows open numerical gates; no physical equations, source laws, model parameters, or campaign horizon are changed by this port.
