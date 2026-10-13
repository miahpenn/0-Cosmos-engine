# Constraint-preserving (a,b) versus X Hamiltonian refinement

Date: 2026-10-10  
Branch: \`diag/discrete-consistent-initial-review\`

## Motivation

The eight-group spatial-metric report showed that the individual a and b attribution values were larger than X at the later probes. However, production enforces the algebraic relation \(a\,b^2=1\). Hybrid states that switch a and b independently violate that relation, so their split is only localization under an unconstrained accounting construction.

## Method

Repeat the same \(N=160\), \(t=3\) candidate trajectory with one-full, two-half, and four-quarter-step endpoint comparisons at t = 0, 0.75, 1.5, and 2.25. Use seven Shapley groups:
- \((a,b)\) together;
- \(X\) separately;
- lapse/shift/gauge;
- \(A_a\);
- \(K\);
- \(\Lambda\);
- all scalar fields and conservative fluids together.

Every hybrid takes a and b as a pair from the same endpoint state, so \(a\,b^2=1\) remains satisfied to floating-point precision. Full-grid attribution closure is retained. This does not assert any other physical causality.

## Guardrails

No evolution equation, source, gauge choice, production default, damping, floor, clamp, or fitted coefficient is changed. The result is diagnostic-only; no initializer admission follows from this run alone.

## Artifact

\`runs/constraint-preserving-ab-x-h-refinement/report.json\`, uploaded as \`constraint-preserving-ab-x-h-refinement\`.
