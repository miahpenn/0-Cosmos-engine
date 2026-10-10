# Geometry-versus-matter Hamiltonian time-refinement diagnostic

Date: 2026-10-10  
Branch: \`diag/discrete-consistent-initial-review\`

## Why this follows the S/PS projection check

The S/PS density increment remained almost identical under full, half, and quarter stepping, while the central Hamiltonian increment showed larger time-step sensitivity. The earlier 11-group attribution showed cancelling S/PS, metric, Lambda, and K contributions, but did not yet test the entire geometry and matter sectors as two grouped endpoints.

## Method

At \(t=0, 0.75, 1.5, 2.25\), use the same N=160 candidate-initialized trajectory and compare the next interval using:
- one full step (\(\Delta t=0.0075\));
- two half steps;
- four quarter steps.

For each endpoint, evaluate all four geometry/matter hybrid states and use a symmetric two-group Shapley decomposition of the Hamiltonian residual increment. Geometry includes \((a,b,X,\alpha,\beta,A_a,K,\Lambda,B)\); matter includes all scalars and all conservative fluids.

Report component closure over the full grid and compare how the geometry and matter contributions change under time refinement. The split is diagnostic attribution, not proof of physical causality.

## Guardrails

No evolution equation, source, gauge choice, initializer default, damping, floor, clamp, or fitted coefficient is changed. The candidate remains opt-in and not admitted to production.

## Artifact

\`runs/geometry-matter-h-refinement/report.json\`, uploaded as \`geometry-matter-h-refinement\`.
