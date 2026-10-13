# Diagnostic report: signed PIRK L2 Ricci-bracket budget at the regular centre

**Date:** 2026-10-13 UTC  
**Disposition:** completed offline diagnostic; no evolution or production edit  
**Branch base:** `cb22285df4c578b1cf76931fbba7661b501d86e6`  
**Current production kernel blob:** `2954e7a79656849b5f088053ea0f4ac97d4bb443`  
**Pinned vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`  
**Archived full-field input:** [run 38325853415, artifact 11780592833](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38325853415/artifacts/11780592833)  
**Archived input ZIP SHA-256:** `b99295c8f485e0ca5386b8e1400d89e010e768d16e2984c727e3251c37854c4d`  
**Current-source comparison:** [PR #39 run 38333176932, artifact 11783348300](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38333176932/artifacts/11783348300); ZIP SHA-256 `356887b0c91b25883914902b61a51b558bfcf083850236ab6455726cd212c270`; report SHA-256 `633f903db8b9618a7c72097b58416a316d776c5e1d4401471f4dc77dc3457d44`.

## Question and scope

PR #39's source-pinned signed stage observer found that the primary-L2 trace-free Ricci term dominates the non-decaying initial first-cell \(\dot A_a\) rate. This offline audit expands that term into the individual subterms of the pinned vendor's radial Ricci bracket and scalar-curvature bracket. It uses existing archived geometry profiles only; no kernel/production/vender code was edited and no production step or new numerical evolution was run.

The older full-field artifact records baseline states at N=40/60/80/120/160 and the initial, predictor, and post-step diagnostic stages. Its source commit is older than current main (`6d8b3f2790644591326df57c55908863f994cb0c`), so this report does **not** claim that all old stage arrays are byte-identical to current-main states. For the initial state, however, recomputing the vendor curvature operators from the stored \(a,b,X,\Lambda,\alpha\) arrays agrees with PR #39's newer current-source profiles at roundoff for every resolution (table below); the newer artifact does not expose the primitive geometry arrays needed to compare every bracket subterm directly.

## Exact operator and signed budget

The pinned vendor implements
\[
R^r{}_r=-\frac{X^2}{a}B_r,\qquad R=-\frac{X^2}{a}B_s,
\]
where (B_r) is assembled from \(a'',\ a\Lambda',\ (a'/a)^2,\ (b'/b)^2,\ \Lambda a',\ a'/(rb),\ (1-a/b)/r^2,\ \chi'',\ \chi' a'/a,\ \chi' b'/b,\ \chi'/r\), and (B_s) contains the corresponding scalar-curvature terms including \(b''/b\), \(b'/(rb)\), \((1-a/b)/r^2\), \(\chi''\), \((\chi')^2\), and explicit \(1/r\) terms. Here \(\chi'=-\tfrac12 X'/X\) and \(\chi''=-\tfrac12[X''/X-(X'/X)^2]\), exactly as in the pinned source.

The target is the actual coded BSSN contribution
\[
\mathcal R_{TF}=\alpha(R^r{}_r-R/3)
=-\alpha\frac{X^2}{a}(B_r-B_s/3).
\]
Each table entry below is the contribution to this final signed \(\mathcal R_{TF}\), not a raw bracket value. Shared terms are grouped only where their algebraic meaning is the same (for example, the paired \(-a\Lambda'\) occurrences). Fourth-order first/second derivatives and the vendor's even/odd centre-ghost conventions were reproduced directly from the pinned `r0_operators.py`.

### Current-main cross-check on the same initial data

| N | max absolute ΔR | max absolute ΔRᵣʳ | max absolute ΔTF Ricci |
|---:|---:|---:|---:|
| 40 | 3.24e-16 | 2.92e-16 | 1.84e-16 |
| 60 | 4.51e-15 | 2.17e-15 | 6.70e-16 |
| 80 | 9.64e-16 | 4.34e-16 | 7.56e-16 |
| 120 | 2.81e-15 | 1.51e-15 | 5.71e-16 |
| 160 | 5.78e-15 | 2.61e-15 | 6.81e-16 |

The signs and sums reconstruct pointwise; the maximum sum-of-subterms roundoff for the grouped trace-free budget is below \(1.4\times10^{-18}\) on the displayed centre calculations. This is arithmetic/implementation closure, not a continuum or physical validation.

### Initial first-cell scaling

| N | r₀ | Hessian TF term | Ricci TF term | Combined L2 \(\dot A_a\) | Ricci TF / r₀² |
|---:|---:|---:|---:|---:|---:|
| 40 | 1.0000 | +2.941834e-5 | +1.077623e-3 | +1.107042e-3 | +1.077623e-3 |
| 60 | 0.6667 | +3.607548e-5 | +1.056736e-3 | +1.092812e-3 | +2.377657e-3 |
| 80 | 0.5000 | +3.873749e-5 | +1.051089e-3 | +1.089826e-3 | +4.204356e-3 |
| 120 | 0.3333 | +4.075529e-5 | +1.047505e-3 | +1.088261e-3 | +9.427549e-3 |
| 160 | 0.2500 | +4.148543e-5 | +1.046320e-3 | +1.087805e-3 | +1.674112e-2 |

The first-cell Ricci term remains approximately \(+1.05\times10^{-3}\) while \(r_0\) falls from 1 to 0.25. The quotient grows about 15.5× for the Ricci component alone; the full L2 quotient grows about 15.7×. That is inconsistent with the required regular-centre class \(\dot A_a=O(r^2)\) for this baseline initialized profile.

### Grouped first-cell Ricci contributions

These are signed contributions to \(\alpha(R^r{}_r-R/3)\), reconstructed from the baseline initial profiles:

| Grouped bracket contribution | N=80 | N=160 |
|---|---:|---:|
| Paired \(-a\Lambda'\) terms | +1.363288e-3 | +1.408027e-3 |
| Metric ratio \((1-a/b)/r^2\) terms | +8.409621e-4 | +8.475179e-4 |
| Radial \(a'/(rb)\) | −5.358158e-4 | −5.591145e-4 |
| Radial \(6b'/(rb)\) scalar term | −5.356006e-4 | −5.590669e-4 |
| \(-2(a/b)b'/(rb)\) | +1.785891e-4 | +1.863702e-4 |
| Paired \(a''/a\) terms | −1.748507e-4 | −1.853091e-4 |
| Explicit \(\chi'/r\) terms | +8.926812e-5 | +9.317792e-5 |
| \(\chi''\) terms | −8.741140e-5 | −9.265116e-5 |
| Scalar \(b''/b\) term | −8.740507e-5 | −9.264952e-5 |
| All remaining grouped terms combined | +6.489832e-8 | +1.765928e-8 |
| **Net trace-free Ricci contribution** | **+1.051089e-3** | **+1.046320e-3** |

The first-cell residual is the incomplete cancellation of several similarly sized, centre-sensitive terms rather than an immediately isolated single summand. The two biggest positive groups are the connection derivative and the metric-ratio block, offset largely by separate explicit radial-gradient terms. This ranking is descriptive; it does not establish which continuum term, initial-data relation, or discrete operator is at fault.

## Reconciliation with existing work and limits

- The earlier independent proper-radial/areal-radius Ricci comparison and direct-D2 Hessian control are already recorded at [Issue #4 comment 6135131094](https://github.com/miahpenn/0-Cosmos-engine/issues/4#issuecomment-6135131094), [6135159478](https://github.com/miahpenn/0-Cosmos-engine/issues/4#issuecomment-6135159478), and [6135213322](https://github.com/miahpenn/0-Cosmos-engine/issues/4#issuecomment-6135213322). I did not rerun or duplicate those controls.
- Those records found a large centre-local baseline discrepancy between the BSSN curvature expression and an independent metric reconstruction, while high-order diagnostic origin quadrature suppresses much of that Ricci-side mismatch. The same alternate initialization still leaves the \(\sim4\times10^{-5}\) lapse-Hessian contribution and the \(\dot A_a/r_0^2\) scaling problem.
- The vendor explicitly uses the independently evolved conformal connection \(\Lambda\) in its BSSN Ricci principal part. Replacing that operator with metric-only Ricci is **not** a valid repair off the connection-constraint surface. A pointwise connection-value identity or a small algebraic closure residual is not enough to adjudicate the derivative/regularity issue.
- Smooth manufactured regular-centre Ricci/\(A_a\) controls and the direct-D2 lapse-Hessian check are existing controls; this audit supplies no contrary proof of a generic vendor sign/parity bug.
- The captured newest artifact has no primitive \(a,b,X,\Lambda\) arrays, so exact current-source term-by-term equality is verified for the initial geometry via curvature outputs, not claimed for every later stage.

## Disposition and next gate

**No unique faulty summand has been identified.** The audit makes the centre-local cancellation budget explicit and supports the pre-existing interpretation that the baseline initialization/regular-centre compatibility is implicated, while preserving the unresolved Hessian-series effect. Do not patch a term by fitting the observed numbers; do not replace BSSN Ricci with metric Ricci; do not promote the alternate quadrature initializer from this artifact alone.

Next discriminator should remain source-pinned and bounded: retain the signed grouped and raw (B_r)/(B_s) subterms at the same actual L2 stage inputs on the *current-source* geometry arrays, then compare directly with the existing independent regular-centre Ricci relation and the already-running CMC/lapse-series lane. This should be an offline/profile-reuse task if the required primitive stage arrays can be recovered; otherwise report the artifact limitation before scheduling another evolution. Physics, vendor, production main, and long-campaign settings remain unchanged.
