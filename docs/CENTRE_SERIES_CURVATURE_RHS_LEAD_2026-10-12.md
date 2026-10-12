# Centre-series curvature RHS: current lead and next test

**Date:** 2026-10-12  
**Exact production source:** \`6d8b3f2790644591326df57c55908863f994cb0c\`  
**Pinned vendor:** \`d6052d605673ce9d82cdc99b0df79855fa2ea215\`  
**Branch:** \`diag/discrete-consistent-initial-production-test-2026-10-12\`  
**Status:** diagnostic lead only; no production patch authorized.

## Bottom line

A matched initial-data A/B run removes the initial discrete Hamiltonian residual, but the unchanged production RHS starts regenerating H/M/C after one step. The independent manufactured momentum-operator test passes near fourth-order convergence on a fixed annulus, so repeating generic parity/constraint-evaluator checks is not the next step.

The most specific lead is the production trace-free curvature/lapse RHS for \(A_a\) at the first radial cell. In the frozen production RHS audit at \(t=0\), the initial \(A_a/r^2\) is zero, yet its exact semi-discrete rate grows approximately like \(1/\Delta r^2\):

| N | dr | first-cell d(Aa/r²)/dt | first-cell Ricci contribution to Aa_dot | \(R^r{}_r-R/3\) |
|---:|---:|---:|---:|---:|
| 40 | 2.000 | 1.1044e-3 | 1.0776e-3 | 1.0648e-3 |
| 60 | 1.333 | 2.4561e-3 | 1.0567e-3 | 1.0429e-3 |
| 80 | 1.000 | 4.3565e-3 | 1.0511e-3 | 1.0369e-3 |
| 120 | 0.667 | 9.7915e-3 | 1.0475e-3 | 1.0330e-3 |

The first-cell Ricci contribution stays O(1) as the cell radius \(r=dr/2\) shrinks, so dividing by \(r^2\) makes the quotient rate blow up. A regular spherical trace-free tensor requires \(A_a=O(r^2)\); its RHS must preserve that series class. This identifies a concrete RHS/centre-series inconsistency to test, not yet a final root-cause verdict.

### Evidence

- Frozen centre-series RHS-rate audit, no evolution: [workflow 38215595339](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38215595339), artifact ID 11732501863, SHA-256 \`8e249d7cd2f6e1b06e1fab8bc26c872a1bbd0b99eedfacae81929e31ed14821c\`.
- Matched discrete-consistent initial metric A/B through t=2.16: [workflow 38284401548](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38284401548).
- One ordinary production step from frozen t=0 A/B states: [workflow 38285057887](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38285057887).
- Saved checkpoint constraint history: [workflow 38287772729](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38287772729).
- Signed K RHS / recomputed CMC target budget, no evolution: [workflow 38288404354](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38288404354).
- Manufactured momentum operator identity, fixed physical annulus, no evolution: [workflow 38288419513](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38288419513). Discrete-vs-analytic momentum RMS errors fall from 2.192e-6 to 3.009e-8 and connection RMS errors from 2.651e-6 to 3.658e-8 over N=40→120; fitted descriptive order is 3.927/3.923. This argues against a generic derivative/parity or constraint-evaluator bug.

## Code-path observation needing reconciliation

In vendor \`vacuum_bssn.py\` at the pinned commit, \`primary_l2_rhs\` advances Aa using

\[
\dot A_a\vert_{L2} =
-\left(D_rD_r\alpha-\frac13 D^2\alpha\right)
+\alpha\left(R^r{}_r-\frac13R\right).
\]

The vendor's \`geometry_terms\` builds the Ricci components from sampled metric second derivatives, \(D_r\Lambda\), \((1-a/b)/r^2\), and explicit \(1/r\) terms. These terms must combine into the regular \(O(r^2)\) radial trace-free class near the origin.

The vendor also provides \`enforce_algebraic_regularity\`, which imposes \(a b^2=1\) and extrapolates the first-cell quotients \((1-a/b)/r^2\), \(A_a/r^2\), and odd \(\Lambda/r\). A scan of the exact local \`engine/production_kernel.py\` used by the t=12 and frozen-state audits found no invocation of that projector in \`V55ProductionKernel.step()\`. Older generic status text says a centre projector is reapplied during evolution; that statement must be reconciled against the pinned production path rather than assumed. This is a code-path discrepancy and a test target, not yet authorization to add a projector.

## Next diagnostic (do not modify production yet)

Build an independent manufactured regular-centre Ricci/Aa RHS test at N=40/60/80/120. Use smooth analytic even \(a,b,X,\alpha\) profiles with \(a b^2=1\), and a matching analytic odd conformal connection. Use proper radial distance \(dl=\sqrt{a}/X\,dr\) and areal radius \(R=r\sqrt{b}/X\) to independently construct metric-derived curvature:

\[
\frac{d}{dl}=\frac{X}{\sqrt a}\frac{d}{dr},\quad
R^r{}_r=-2\frac{R_{ll}}R,\quad
R=-4\frac{R_{ll}}R+2\frac{1-R_l^2}{R^2}.
\]

Then compare signed production and analytic \(R^r{}_r-R/3\), trace-free lapse Hessian, the full Aa RHS, \(r^2\dot A_a\), and \(\dot A_a/r^2\) in the first four cells and at fixed physical radii. Also record determinant, connection, and full decomposition closure. The analytic side must not reuse the vendor's \`geometry_terms\`.

Predeclared pass criterion: the full production Aa RHS must preserve the \(O(r^2)\) centre series (bounded/convergent first-cell \(\dot A_a/r^2\)) and its curvature/lapse pieces must converge to the independent analytic reference. Fail closed on non-finite quantities or mismatched signed-term closure. Only after a specific term fails this identity should a single bounded frozen-input numerical counterfactual be considered. H/M/C and Misner–Sharp balance remain separate outcomes.

## Guardrails

- Main is untouched; no changes to scalar/matter physics, sources, gauge, EOS, or potentials.
- No damping, floors, clamps, fitted centre law, added physical source, or interpretation as bounce/ejection.
- Do not launch t=160 before the t=12 constraint and balance admission gates are resolved.
- The detailed current continuation handoff is also saved in the local Library at \`/0* Cosmos Engine/CONTINUATION_HANDOFF_2026-10-12_RHS_COMPATIBILITY.md\`.
