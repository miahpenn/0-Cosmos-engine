# S/PS scalar stress projection and step-refinement audit

Date: 2026-10-10  
Branch: \`diag/discrete-consistent-initial-review\`

## Purpose

The first-cell RHS diagnostic found large, cancelling S/PS, metric, Lambda, and K contributions to the accepted-step Hamiltonian increment, while D/PD did not contribute at the report's numerical precision at the four sampled times. No equation-level defect was established.

This diagnostic therefore tests numerical time-step sensitivity before changing any equation. At the same four states (t = 0, 0.75, 1.5, 2.25), it compares one full step, two half steps, and four quarter steps over the same time interval. It also computes the S-sector energy density using the existing scalar-projection formula, then symmetrically attributes its change between the S/PS field pair and the metric dependence through a and X.

## Guardrails

- Same N=160 candidate-initialized full-kernel trajectory and parameters as the previous diagnostic.
- No evolution equation, source, gauge choice, or production default is modified.
- Time-step differences are reported, not forced to meet an invented physics threshold.
- The rho_S component split checks accounting and attribution; it is not proof of physical causality.

## Artifact

\`runs/s-ps-projection-refinement/report.json\`, uploaded as \`s-ps-projection-refinement\`.
