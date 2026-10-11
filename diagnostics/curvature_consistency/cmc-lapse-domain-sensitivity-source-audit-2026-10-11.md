# CMC Lapse Domain Sensitivity — Source-Level Audit

**Date:** 2026-10-11  
**Scope:** Read-only implementation audit after the three-domain residual campaign  
**Branch:** `diag/independent-metric-ricci-audit`

## What the production code actually solves

The current implementation is `engine/cmc_gauge.py`; `V55ProductionKernel._solve_lapse` delegates to `solve_cmc_lapse`. Its interior discrete equation is assembled as

\[
\alpha'' + c(r)\alpha' - m(r)\alpha
= -\frac{\dot K_{\rm CMC}}{X^2/a},
\]

with

\[
c=-\frac12\frac{a'}a+\frac{b'}b-\frac{X'}X+\frac2r,
\qquad
Q=\frac32 A_a^2+\frac13K^2+4\pi(\rho+p_r+2p_t),
\qquad
m=\frac{Q}{X^2/a}.
\]

The solver uses the repository's cell-centered fourth-order/parity derivative rows (including its one-sided outer closure) and a banded linear system. The last cell is a Dirichlet normalization row, **\(\alpha(R)=1\)**. Interior equations are imposed through the penultimate cell.

The CMC target \(\dot K_{\rm CMC}\) is not a fixed input. `target_kdot` recomputes the current K RHS as the sum of the vacuum L2 contribution and the matter-coupled L3 contribution, then takes its proper-volume weighted mean over the **entire current grid**. The code's weight (omitting the common (4\pi) factor) is

\[
w_i=r_i^2\sqrt{a_i}\,b_i/X_i^3\,\Delta r.
\]

Consequently the domain changes the lapse through at least two direct global paths:
1. the outer normalization row \(\alpha(R)=1\), whose radius changes;
2. the full-volume projection that determines \(\dot K_{\rm CMC}\), whose integration domain and boundary-adjacent RHS values change.

There is also an **evolution-history path**: the production kernel applies a separate incoming-light/characteristic boundary condition to the outer geometry before its final stage-aware lapse solve. Independently evolved runs on different domains therefore need not have identical interior geometry/matter states, even where their radial grids overlap.

## What the previous artifacts can and cannot establish

The domain-size profile artifacts retain full radial profiles of (H_{\rm vendor}) and \(\alpha\), plus per-sample summary diagnostics. They do **not** retain the full per-cell arrays for (a,b,X,A_a,K,\Lambda\), scalar fields, conserved matter variables, or the computed \(\dot K_{\rm CMC}\). Therefore the previous saved profiles are sufficient to measure lapse-profile changes but **not sufficient to independently reconstruct the exact lapse operator/source or separate boundary-value, global-target, and evolved-state effects**.

This limits the causal interpretation of the strong domain dependence found in the three-domain audit. It is not evidence that the lapse causes the residual movement, and it does not identify a defective physics layer.

## Required next diagnostic

The next campaign is a **frozen-state crop-resolve audit**, not a physics change and not another three-domain evolution sweep:

1. Evolve one unchanged production trajectory at (N=320, R=160, \Delta r=0.5) to (t=12), saving the full state at the same (0.25)-time intervals.
2. At each saved state, re-solve the current native CMC lapse system on the full domain and on exact crops at ((N,R)=(160,80),(240,120),(320,160)). The metric, scalar, and matter arrays on every shared interior cell are held identical to the largest-domain state.
3. Record each crop's native \(\dot K_{\rm CMC}\), lapse profile, coefficients, and discrete operator residual.
4. Add diagnostic counterfactual solves using the *full-domain* \(\dot K_{\rm CMC}\) with (a) the native crop normalization \(\alpha(R)=1\), and (b) the full-domain solved \(\alpha\) value at that crop radius. These do not alter evolution. Comparing them isolates the global target-projection effect and, with the operator held fixed, the Dirichlet boundary-value effect.
5. Compare these frozen-state crop results with the actual (r_{\max}=80,120) lapse profiles from the previous campaign. The mismatch then includes evolution-history differences, rather than being mistaken for a pure elliptic boundary effect.

The residual of every diagnostic linear solve must be reported against the exact assembled discrete operator. Full-state capture is necessary for this audit; the previous ((H_{\rm vendor},\alpha))-only profiles cannot substitute for it.

## Guardrails

- Do not change the production lapse equation, boundary conditions, CMC target, gauge, projection, matter/geometry equations, or evolution algorithm.
- Do not infer physical D propagation, a bounce, turnaround, or cycle from the residual ridge or lapse response.
- Do not begin another sweep until the single-trajectory frozen-state audit completes and its saved artifact is verified.
