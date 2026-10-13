# Origin-half-cell-only PIRK stage Ricci comparison — 13 October 2026

**Status:** completed diagnostic-only; no scientific admission.  
**Draft PR:** [#84](https://github.com/miahpenn/0-Cosmos-engine/pull/84)  
**Exact-head run:** [38344327911](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38344327911)  
**Artifact:** [11788516284](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38344327911/artifacts/11788516284)  
**Artifact ZIP SHA-256:** \`130cc945e6fc8232e8a9a2f902afa48b6821be87d6bbc415e6916c8f9642b060\`

## Exact scope and guardrails

- Source baseline: \`c01f2467347e8aebd57d220374fd7bc4e93045bc\`.
- Baseline initializer blob: \`4cfbfe4be6ac3ce6ed492f4e435186997d94116a\`.
- Origin-half-cell-only initializer blob: \`3be4ec52490bcab9f1d71c5b21c347367c0b89f6\`.
- Production kernel blob: \`2954e7a79656849b5f088053ea0f4ac97d4bb443\`.
- Vendor: \`d6052d605673ce9d82cdc99b0df79855fa2ea215\`.
- Only \`engine/v55_initial.py\` differs under \`engine/\`; vendor, submodule and dependency pins match.
- Full software suite: **234 passed**. Compilation, source provenance, candidate-scope and vendor checks passed.
- N=160/240/320, \`r_max=80\`, radiation ON, \`CFL=0.03\`; the baseline and candidate each take exactly one normal production step (dt=0.015/0.01/0.0075 respectively).
- Every arm recorded the four real production L2 inputs; the source state remained unchanged; all signed arrays were finite. Vendor L2 reconstruction closure was exactly zero for every recorded stage.
- No RHS, equation, CMC gauge, physical source, boundary, production kernel, vendor file, or admission condition changed.

The counterfactual only replaces the first origin half-cell integral with the even-series weights exact for \(1,r^2,r^4\). Every inter-centre interval remains the baseline trapezoid.

## The measurement

At each actual PIRK L2 input, the observer records:
1. \(\alpha(R^r{}_r-R/3)\) from the vendor BSSN Ricci brackets;
2. the signed trace-free lapse-Hessian contribution \(-(\partial_r^2\alpha-\Delta\alpha/3)\);
3. their sum and the vendor's returned L2 \(A_a\) rate;
4. an independent proper-radial reconstruction of radial/tangential Ricci from \(R_{\rm areal}=r\sqrt{b}/X\) and \(f=\sqrt a/X\);
5. the signed difference between the two trace-free Ricci representations on the exact same stage metric.

The independent curvature construction does not use evolved \(\Lambda\) or vendor Ricci brackets, but it shares the pinned radial D1/D2 stencils. Its formula/sign is covered by the manufactured constant-curvature regression in this PR.

## First-cell signed result, all four stages

The centre radius is \(r_0=40/N\). Values below are the unscaled first-cell trace-free Ricci representation gap (vendor \(\alpha R_{\rm TF}\) minus proper-radial \(\alpha R_{\rm TF}\)) and the corresponding vendor L2 \(A_a\) RHS numerator.

| N | Actual L2 input | Baseline Ricci gap | Origin-only Ricci gap | Baseline L2 \(A_a\) RHS | Origin-only L2 \(A_a\) RHS |
|---:|---|---:|---:|---:|---:|
| 160 | Initial RHS | +1.614303e-3 | -1.900102e-4 | +1.087805e-3 | -9.097363e-5 |
| 160 | Explicit predictor | +1.613958e-3 | -1.900336e-4 | +1.087599e-3 | -9.097088e-5 |
| 160 | Primary predictor | +1.613958e-3 | -1.900336e-4 | +1.087599e-3 | -9.097063e-5 |
| 160 | Final implicit L2 | +1.603656e-3 | -1.872800e-4 | +1.081038e-3 | -8.905878e-5 |
| 240 | Initial RHS | +1.612622e-3 | -1.956798e-4 | +1.087505e-3 | -9.388192e-5 |
| 240 | Explicit predictor | +1.612395e-3 | -1.956956e-4 | +1.087375e-3 | -9.387459e-5 |
| 240 | Primary predictor | +1.612395e-3 | -1.956956e-4 | +1.087375e-3 | -9.387444e-5 |
| 240 | Final implicit L2 | +1.602113e-3 | -1.928890e-4 | +1.080831e-3 | -9.192439e-5 |
| 320 | Initial RHS | +1.612045e-3 | -1.974897e-4 | +1.087405e-3 | -9.478899e-5 |
| 320 | Explicit predictor | +1.611876e-3 | -1.975016e-4 | +1.087309e-3 | -9.478144e-5 |
| 320 | Primary predictor | +1.611876e-3 | -1.975016e-4 | +1.087309e-3 | -9.478134e-5 |
| 320 | Final implicit L2 | +1.601601e-3 | -1.946775e-4 | +1.080772e-3 | -9.281864e-5 |

The origin-only rule reduces the absolute central Ricci representation gap by about **88%** at these high resolutions and reduces the primary L2 \(A_a\) numerator by about **91%**. However, it does not remove the gap: its signed magnitude is roughly \(1.90\times10^{-4}\) to \(1.97\times10^{-4}\) from N=160 to N=320, while \(r_0\) halves from 0.25 to 0.125. The metric-derived TF Ricci itself also remains nonzero at the first cell (about \(+5.7\times10^{-5}\) to \(+6.0\times10^{-5}\)); the vendor term remains about \(-1.31\times10^{-4}\) to \(-1.35\times10^{-4}\). Both representations therefore fail to exhibit the expected vanishing first-cell TF-centre sequence in this candidate.

The \(A_a\) L2 numerator likewise remains near \(-9.1\times10^{-5}\) to \(-9.5\times10^{-5}\) as the cell radius shrinks. Dividing by \(r_0^2\), the candidate's full initial \(A_a\) rate is approximately \(-1.46\times10^{-3}\), \(-3.38\times10^{-3}\), and \(-6.07\times10^{-3}\) at N=160/240/320. This is still a growing, non-regular centre quotient.

## Where the discrepancy appears

There is no new large jump from initial RHS to explicit or primary predictor: the signed gap changes by around \(10^{-7}\) or less between those first three stage inputs. The final implicit L2 input reduces its absolute gap only slightly. The same shape holds for baseline and origin-only arms.

The lapse-Hessian numerator is approximately \(+4.15\times10^{-5}\) to \(+4.22\times10^{-5}\) in both arms; its baseline-versus-candidate change is below \(3\times10^{-9}\). The central L2 change is dominated by the signed trace-free Ricci term, while the matter/L3 contribution remains small and does not cancel it. Each vendor L2 split closes to the returned operator exactly in the archived arrays, so this is not an instrumentation arithmetic failure.

This sharpens the attribution:
- the origin half-cell correction changes the initial geometry substantially and removes most of the baseline central source anomaly;
- but the origin-only profile still carries a persistent high-resolution first-cell curvature mismatch *before the evolution update*;
- predictor stages preserve it rather than creating a new defect during this first step.

## Interpretation and next gate

The evidence argues against further tuning of the origin half-cell coefficients. It also does not justify editing the vendor trace-free Ricci equation: the independent metric-based Ricci is itself non-vanishing at the first cell, so the candidate's metric profile / near-origin discrete differentiation remains incompatible with the expected regular-centre series. The broader all-interval quadrature counterfactual is not equivalent to this isolated half-cell experiment and already has a distinct mass-current/flux-work trade-off.

Next, reconcile this stage-local result with the existing uniform radial-quadrature and CMC-centre-row audits. Specifically compare the first-cell source integral and the first two inter-centre increments against the regular even Taylor moments on the *same stored production profiles*, and keep the fixed-areal-radius current ledger as a separate admission gate. No production physics change, fitted correction, or promotion is justified here.

Scientific admission remains false.
