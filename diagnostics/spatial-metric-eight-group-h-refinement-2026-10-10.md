# Spatial-metric eight-group Hamiltonian time-refinement diagnostic

Date: 2026-10-10  
Branch: \`diag/discrete-consistent-initial-review\`

## Motivation

The six-layer audit showed that the central-H refinement difference at t = 0.75, 1.5, and 2.25 is dominated by the \((a,b,X)\) spatial-metric group. Gauge, \(A_a\), and \(K\) differences were negligible, and the \(\Lambda\) difference was much smaller than the spatial-metric difference.

Because the algebraic regularity projection couples \(a,b,X\), the next step is to split those three variables symmetrically while keeping the remaining registered state groups unchanged.

## Method

Use the same N=160 candidate-initialized main trajectory. At t = 0, 0.75, 1.5, and 2.25 compare one full step, two half steps, and four quarter steps. For each endpoint evaluate the 256 combinations of eight hybrid groups:
- \(a\), \(b\), \(X\), each separately;
- lapse/shift/gauge \((\alpha,\beta,B)\);
- \(A_a\);
- \(K\);
- \(\Lambda\);
- all scalar fields and all fluid species together.

Report central and first-five-cell component attributions, full-grid closure, and full/half versus half/quarter changes. These variables remain coupled in the real trajectory; Shapley values are numerical attribution, not independent physical forces.

## Guardrails

No evolution equation, source term, gauge choice, production default, damping, floor, clamp, or fitted coefficient is changed. No result from this diagnostic alone admits the candidate initializer to production.

## Artifact

\`runs/spatial-metric-eight-group-h-refinement/report.json\`, uploaded as \`spatial-metric-eight-group-refinement\`.
