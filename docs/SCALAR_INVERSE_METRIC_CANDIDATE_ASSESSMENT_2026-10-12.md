# Scalar lapse-gradient inverse-metric correction — candidate assessment

Date: 2026-10-12

## Scope

Candidate branch: \`fix/scalar-inverse-metric-lapse-gradient-strict-2026-10-12\`  
Strict-admission parent: \`d5a5b97eadcd5c565722e61262a22195baf3c1d5\`  
Production source baseline: \`6d8b3f2790644591326df57c55908863f994cb0c\`  
Pinned numerical kernel: \`d6052d605673ce9d82cdc99b0df79855fa2ea215\`

The only production-expression change is in \`engine/scalar_system.py\`: in the scalar momentum RHS, replace \`geometry.X**-2 * alphap * fp\` with \`invr * alphap * fp\`, where \`invr = X**2/a\` is the already-used inverse radial metric coefficient.

The declared physical spatial metric is

\[
\gamma_{rr} = \frac{a}{X^2}, \qquad
\gamma^{rr} = \frac{X^2}{a}.
\]

The lapse-gradient contraction in the scalar momentum equation is \(\nabla^i\alpha\,\nabla_i f = \gamma^{rr}\,\partial_r\alpha\,\partial_r f\). The scalar Laplacian in this same RHS already uses \(X^2/a\), as do the metric and stress projections. This fixes an inconsistent contraction; it adds no new term or parameter.

A new regression test, \`tests/test_scalar_lapse_gradient_metric_contraction.py\`, checks S, D, and COSMOS \(\phi\) momentum RHS values on a deliberately non-flat smooth metric and lapse. The test verifies the declared inverse metric identity and compares the production RHS against the expected covariant contraction.

## Execution and tests

Initial candidate run: [38282365811](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38282365811)  
Follow-up run after fixing a syntax typo in the report-comparison helper: [38282790322](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38282790322); this rerun completed successfully.

The first run reports:
- full repository regression suite: passed;
- explicit scalar metric-contraction regression: passed;
- candidate production campaign at \(N=40,60,80\), \(r_{\max}=80\), CFL \(=0.03\), with radiation enabled: every requested resolution reached \(t=12\);
- strict campaign-admission result: passed for all three resolutions.

The first run's last step, the JSON comparison helper, initially failed from a closing-bracket syntax typo in \`tools/compare_scalar_metric_candidate_campaign.py\`. This did **not** invalidate the already-completed numerical trajectory; candidate history, checkpoints, and ledgers were uploaded. The reporting syntax was corrected in commit \`f733bbf5e662e8894e673a998ca915e68246bb3a\`, and the follow-up is rechecking the complete workflow. The trajectory comparison was also independently computed directly from the uploaded baseline and candidate ledgers/checkpoints.

## Matched baseline comparison

Baseline campaign: [38208253687](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208253687)  
Candidate campaign/artifact: [38282365811](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38282365811), artifact \`scalar-inverse-metric-candidate-t12\`.

All baseline and candidate runs are time-aligned: ledger histories have 201 / 301 / 401 rows for N=40 / 60 / 80, respectively, and identical stored coordinate-time sequences. Each reached \(t=12\). The correction only produces very small downstream changes at this operating point.

| Resolution | \(\tau_{\rm base}(12)\) | \(\tau_{\rm candidate}(12)\) | Difference |
|---:|---:|---:|---:|
| N=40 | 11.9559993814 | 11.9559982838 | \(-1.098\times10^{-6}\) |
| N=60 | 11.9284259877 | 11.9284246545 | \(-1.333\times10^{-6}\) |
| N=80 | 11.9184114053 | 11.9184099809 | \(-1.424\times10^{-6}\) |

Maximum absolute history differences in the key diagnostics remain small:
- \(H_{\rm eff}\): \(2.6\)–\(3.9\times10^{-11}\);
- exterior Hamiltonian max witness: \(6.1\)–\(8.8\times10^{-8}\);
- admitted/inner Hamiltonian max witness: \(1.9\times10^{-7}\) (N40), \(3.8\times10^{-7}\) (N60), and \(6.8\times10^{-7}\) (N80);
- inner momentum max witness: \(2.0\times10^{-8}\) (N40), \(5.0\times10^{-8}\) (N60), and \(1.5\times10^{-7}\) (N80);
- inner connection max witness: \(9.3\times10^{-8}\) (N40), \(1.3\times10^{-7}\) (N60), and \(1.5\times10^{-7}\) (N80);
- maximum lapse-rate change is around \(1.1\)–\(1.2\times10^{-6}\).

The Misner–Sharp current balance is effectively unchanged. RMS differences at \(t=12\) are approximately \( -4.36\times10^{-8}\) (N40), \(+2.38\times10^{-9}\) (N60), and \(-1.14\times10^{-9}\) (N80) for the reported maximum/RMS pair; the RMS changes alone are far smaller than their existing residual norms. The candidate MS current RMS remains approximately \(1.262\times10^{-3}\), \(2.315\times10^{-4}\), and \(2.144\times10^{-4}\).

Checkpoint-level field differences at \(t=12\) are likewise small: the largest reported geometry-field differences are at most about \(2.6\times10^{-7}\) for \(\Lambda\), while scalar S differs by at most about \(2.5\times10^{-7}\); the rest of the scalar, metric, and conservative matter fields differ less or similarly at sub-micro scale. These are trajectory sensitivities, not acceptance improvements.

## What it does and does not fix

This change corrects a demonstrable inverse-metric inconsistency in the scalar equation. The small trajectory difference is consistent with the term's measured small instantaneous effect on the total matter-current tangent for the tested states.

It does **not** explain or repair the large central constraint layer. At \(t=12\), the candidate full-domain maxima still sit at the first center cells, with approximate magnitudes:
- Hamiltonian: \(4.04\times10^{-2}\) (N40), \(7.97\times10^{-2}\) (N60), \(1.35\times10^{-1}\) (N80);
- momentum: \(3.32\times10^{-3}\), \(3.00\times10^{-3}\), \(8.08\times10^{-3}\);
- connection: \(2.07\times10^{-2}\), \(2.52\times10^{-2}\), \(3.09\times10^{-2}\).

The first two center cells are not the accepted formal Hamiltonian/momentum stencil region; nevertheless, the admitted inner residuals also remain substantial. The CMC center-row and full fourth-order CMC operator counterfactuals were separately tested and rejected as repairs because they made the admitted momentum residual slightly worse while barely changing H and connection errors. The stage-sync and regularized-\(D_rA_a\) counterfactuals likewise failed to make material progress.

## Decision

1. Keep this metric-contraction correction isolated on the candidate branch; it is equation-correctness work, not an approved production promotion.
2. Do not merge the candidate into \`main\` merely because it passes a test suite and reaches \(t=12\). The numerical/scientific admission gate is still open.
3. Do not claim that it fixes the principal constraint problem or justifies a long \(t=160\) campaign.
4. Continue root-cause diagnosis from the established central spatial layer and geometric momentum-constraint budget, using frozen states and stage-level signed residuals. The dominant failure remains orders of magnitude larger than this correction's trajectory effect.

No change to \`main\`, no new source coefficient, damping, floor, projection, imposed bounce, or other physical intervention was made.


The corrected t=12 baseline-comparison workflow completed successfully in [run 38282790322](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38282790322), artifact ID `11762490670`.
