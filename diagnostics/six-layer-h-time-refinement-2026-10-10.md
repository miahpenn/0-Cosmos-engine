# Six-layer Hamiltonian time-refinement diagnostic

Date: 2026-10-10  
Branch: \`diag/discrete-consistent-initial-review\`

## Motivation

The geometry-versus-matter refinement result localizes almost all of the central-Hamiltonian step-size sensitivity at later probes to the geometry group. The previous attribution audit identifies four relevant geometry subgroups—spatial metric, lapse/shift/gauge, trace-free curvature \(A_a\), trace curvature \(K\), and connection \(\Lambda\)—alongside matter. This diagnostic resolves those groups under the same full/half/quarter endpoint comparison.

## Method

At t = 0, 0.75, 1.5, and 2.25, start from the matched N=160 candidate trajectory and compare one full step, two half steps, and four quarter steps. Apply the already-reviewed symmetric Shapley decomposition over 64 hybrids to each endpoint, with groups:
- spatial metric \((a,b,X)\);
- lapse/shift/gauge \((\alpha,\beta,B)\);
- \(A_a\);
- \(K\);
- \(\Lambda\);
- all scalar fields and all fluid species.

Each probe records the component increments and the full-grid closure error. This is intended to localize the numerical time-step sensitivity and does not establish a causal or physical defect.

## Guardrails

No evolution equations, source terms, gauge choices, production defaults, damping, floors, clamps, or fitted coefficients are changed. The candidate remains diagnostic-only and not admitted to production.

## Artifact

\`runs/six-layer-h-refinement/report.json\`, uploaded as \`six-layer-h-refinement\`.
