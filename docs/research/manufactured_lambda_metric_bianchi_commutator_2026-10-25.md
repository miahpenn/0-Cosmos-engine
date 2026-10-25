# Manufactured metric control: metric-Lambda Bianchi commutator — 2026-10-25

**Status:** research diagnostic only. No accepted target-state profiles are read, no evolution is run, and no production equation, vendor source, or dependency is changed.

**Purpose:** test whether the same nested-derivative versus primitive chain-rule commutator seen in the accepted-state annular analysis also appears on prescribed, smooth, origin-regular metrics. Decompose its signed contribution into the three declared terms of
\[
\Lambda_{\rm metric}=\frac{D_1a}{2a^2}-\frac{D_1b}{ab}
+\frac{2}{r}\left(\frac1b-\frac1a\right).
\]
The key question is whether the radial reference-connection term is again a leading channel under independent changes to profile amplitude and length scale.

## Method

For each metric, set
\[
q(r)=1+A\left(1-e^{-(r/L)^2}\right),\qquad
a=q^{-4/3},\quad b=q^{2/3},\quad X=q^{-2/3},\quad R_{\rm areal}=rq.
\]
This guarantees \(ab^2=1\), unit proper-radial factor \(\sqrt a/X=1\), a regular origin, and a known areal-radius map. Use \(r_{\max}=80\), cell-centred grids \(N=60,120,240\), and the pinned parity-aware fourth-order D1/D2 stencils (even parity for \(a,b,X\); odd parity for \(\Lambda_{\rm metric}\)). The annular norm uses native cells satisfying \(2\le R_{\rm areal}<10\).

Each of the three assembled terms is differentiated with D1 and compared with its primitive D1/D2 chain derivative. Its resulting \(\Delta\Lambda'_{\rm metric}\) is mapped through the linear Bianchi response
\[
\Delta\mathcal B_h=
\frac12D_1(X^2\Delta\Lambda'_{\rm metric})
+\frac23P X^2\Delta\Lambda'_{\rm metric},\qquad
P=6\chi_r+\frac32\left(\frac2r+\frac{D_1b}{b}\right),\quad
\chi_r=-\frac{D_1X}{2X}.
\]
The three signed channels are independently summed and checked against the unsplit nested-versus-chain reconstruction. The script rejects non-finite values or a decomposition identity gap above \(5\times10^{-13}\); this is an algebra/accounting check, **not** a physical acceptance threshold.

Reproduce with:

\`\`\`bash
python tools/audit_manufactured_lambda_metric_commutator.py \
  --output results/manufactured-lambda-metric-bianchi-commutator-2026-10-25/summary.json
\`\`\`

## Results from the local reproduction

The radial-connection channel is the largest annular RMS channel in all three profile controls. The observed order between \(N=120\) and \(N=240\) is approximately fourth order for the signed sum and for each component.

### Reference smooth metric: \(A=0.1,\ L=7\)

| \(N\) | Annular cells | \(a\)-term RMS | \(b\)-term RMS | radial-connection RMS | signed-sum RMS | \(N^4\) × signed-sum RMS |
|---:|---:|---:|---:|---:|---:|---:|
| 60 | 6 | \(2.8915\times10^{-6}\) | \(2.6589\times10^{-6}\) | \(1.5569\times10^{-5}\) | \(1.5614\times10^{-5}\) | 202.36 |
| 120 | 11 | \(1.9240\times10^{-7}\) | \(1.5916\times10^{-7}\) | \(6.2945\times10^{-7}\) | \(8.2636\times10^{-7}\) | 171.35 |
| 240 | 22 | \(1.2226\times10^{-8}\) | \(1.0146\times10^{-8}\) | \(4.0681\times10^{-8}\) | \(5.2611\times10^{-8}\) | 174.55 |

Observed order of the signed-sum RMS is 4.24 from \(N=60\to120\), then 3.97 from \(N=120\to240\). At \(N=240\), the radial channel RMS is 3.33 times the \(a\)-term RMS and 4.01 times the \(b\)-term RMS.

### Profile controls

| Control \((A,L)\) | Signed-sum RMS at \(N=60,120,240\) | \(N^4E\) at \(N=60,120,240\) | Signed-sum order, coarse then fine | Fine-grid radial/\(a\)-term RMS |
|---|---|---|---|---:|
| Reference \((0.1,7)\) | \(1.5614e{-5},\,8.2636e{-7},\,5.2611e{-8}\) | 202.36, 171.35, 174.55 | 4.24, 3.97 | 3.33 |
| Half amplitude \((0.05,7)\) | \(6.8987e{-6},\,3.8565e{-7},\,2.3969e{-8}\) | 89.41, 79.97, 79.52 | 4.16, 4.01 | 3.33 |
| Shorter scale \((0.1,5)\) | \(1.0032e{-4},\,5.1033e{-6},\,3.3585e{-7}\) | 1300.13, 1058.21, 1114.27 | 4.30, 3.93 | 1.53 |

The shorter-scale profile is less asymptotic on the coarse grid, but the fine-grid order remains close to four. Changing amplitude substantially changes the coefficient without changing the observed order. Shortening the length scale also changes the relative channel strengths, but the radial-connection term remains the largest single channel in the tested annulus.

Identity closure was at most \\(4.4\\times10^{-17}\\) in the verified CI artifact, well below the accounting guard. Fixed-R signed samples at R=0.75, 1.0, 1.25, 1.5, 1.75, 2, 3, 5, and 8 are included in the machine-readable output so the center-to-annulus transition and cancellation signs are visible rather than reporting RMS alone.

## Center and inner-shell readout

The same control now reports the first three native-cell commutator values and region RMS separately for R<1, 1<=R<2, and the annulus 2<=R<10. These are diagnostics of the commutator Delta B_h, not the full target-state propagation residual.

| Profile control | Apparent R<1 RMS order, 60->120 / 120->240 | 1<=R<2 RMS order, 120->240 | Fixed R=0.75 order, 60->120 / 120->240 | Annular signed-sum order, 60->120 / 120->240 |
|---|---:|---:|---:|---:|
| Reference (A=0.1, L=7) | 2.90 / 3.66 | 4.19 | 3.60 / 4.29 | 4.24 / 3.97 |
| Half amplitude (A=0.05, L=7) | 2.95 / 3.67 | 4.15 | 3.63 / 4.28 | 4.16 / 4.01 |
| Shorter scale (A=0.1, L=5) | 2.83 / 3.64 | 4.16 | 3.52 / 4.25 | 4.30 / 3.93 |

At N=60, the 1<=R<2 region contains no cell centres for these profiles, so no order is inferred from that grid. The origin RMS also samples different numbers of cell centres as N changes (1, 1, and 3 across N=60,120,240); its reported orders are therefore **descriptive, not a formal convergence estimate**. At fixed R>=1, the signed samples generally return to approximately fourth-order refinement, while the R=0.75 coarse-to-middle interval is less asymptotic and its fine interval is closer to fourth order. This supports keeping origin native-cell, fixed-R inner-shell, and annular measurements separate rather than merging them into one error norm or order claim.

## Interpretation and limits

1. **Supported:** the nested D1 of the assembled \(\Lambda_{\rm metric}\) and the primitive D1/D2 chain-rule representation produce a nonzero, refinement-convergent commutator on smooth analytic metrics. The radial \(2(1/b-1/a)/r\) term is the largest RMS contribution in these controls, consistent with its prominence in the target-state channel split.
2. **Not established:** this does not prove the target-state residual is fully explained, set an acceptable production-error budget, establish a discrete Bianchi identity that the production scheme is required to satisfy, or show that any candidate Ricci form is physically authoritative.
3. **No operator selection:** component grouping is exact signed bookkeeping but is representation-dependent; no term is an independently established physical cause. The evolved-\(\Lambda\) BSSN Ricci path must remain intact pending a frozen, archive-backed discrete contract.
4. **Keep regions separate:** the control is about the annular stencil commutator. It does not cure or reclassify the first native cells or the inner-shell fixed-\(R\) transition.

**Disposition:** this independent manufactured control strengthens the finite-grid compatibility/truncation interpretation but does not close the matter-coupled propagation identity, total stress-energy conservation, the origin/inner-shell gate, or scientific admission. Those flags remain false.
