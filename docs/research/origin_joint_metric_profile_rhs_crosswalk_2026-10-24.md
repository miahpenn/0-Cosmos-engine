# Initial-origin RHS crosswalk: partial half-cell versus full continuous metric

**Status:** research diagnostic only; not a production change or scientific-admission decision.  
**PR:** [#417](https://github.com/miahpenn/0-Cosmos-engine/pull/417), stacked on the PR #416 research branch.  
**Diagnostic run:** [GitHub Actions run 40155299835](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40155299835), all three resolution jobs succeeded.  
**Head commit:** `7e83c61d7cf56fb2846465be4387fd2b615a70e9`.

## Question

Does the centre defect arise mainly from reinitializing conserved matter on the half-cell candidate, or from the metric profile that feeds the discrete Ricci stencil?

The bounded check evaluates the initial geometry-stage (A_a) RHS for three copied-state arms at (N=60,120,240), then performs one production-kernel step only as a replay discriminator. It does **not** repeat the full (t=0.60) evolutions.

- **Partial half-cell, matter reprojected:** PR #416 candidate. The continuous-source integral corrects only ([0,r_0]); later radial intervals retain the pinned trapezoid recurrence. Matter is reinitialized on the rebuilt metric.
- **Partial half-cell, matter held fixed:** same partial metric construction, but the source-pinned initializer's original conserved matter arrays are restored before the copied CMC solve. This isolates matter-reprojection effects.
- **Full continuous metric:** the continuous-source radial metric ODE is integrated over the full profile using 2048 substeps per coarse cell. As in the original full-metric shadow, (A_a), (K), scalar fields and conserved matter arrays are held fixed while the geometry metric is rebuilt.

All arms use the same copied regular-origin CMC row and pinned production geometry RHS. The full-profile and partial arms differ in the metric construction; the second arm measures whether the matter-state difference could explain the result.

## Main measurements

The quantity (RHS/r_0^2) is the initial combined (L_2[A_a]+L_3[A_a]) RHS at cell 0, divided by (r_0^2). The final column is the accepted value of (A_a/r_0^2) after the single step.

| N | Partial half-cell: initial RHS / (r_0^2) | Partial half-cell, matter fixed | Full continuous metric | Partial accepted (A_a/r_0^2) | Full-profile accepted (A_a/r_0^2) |
|---:|---:|---:|---:|---:|---:|
| 60 | (-2.4278302\times10^{-4}) | (-2.4278287\times10^{-4}) | (-2.0031306\times10^{-5}) | (-9.67399\times10^{-6}) | (-8.05276\times10^{-7}) |
| 120 | (-1.1796777\times10^{-3}) | (-1.1796776\times10^{-3}) | (-2.2807673\times10^{-5}) | (-2.34986\times10^{-5}) | (-4.58790\times10^{-7}) |
| 240 | (-4.9146818\times10^{-3}) | (-4.9146817\times10^{-3}) | (-2.3223439\times10^{-5}) | (-4.89500\times10^{-5}) | (-2.33635\times10^{-7}) |

The key geometric component, the trace-free Ricci contribution to (L_2[A_a]), is:

| N | Partial half-cell | Full continuous metric |
|---:|---:|---:|
| 60 | (-1.0102454\times10^{-4}) | (-2.0257443\times10^{-6}) |
| 120 | (-1.2925520\times10^{-4}) | (-7.1410734\times10^{-7}) |
| 240 | (-1.3605744\times10^{-4}) | (-1.8355045\times10^{-7}) |

The trace-free Ricci term for the full continuous metric falls by approximately a factor of four under each refinement; its (r_0^{-2})-normalized value approaches a finite scale. In contrast, the partial half-cell term remains (O(10^{-4})), and its normalized value grows by roughly a factor of four with each refinement. Its one-step accepted central value likewise worsens.

## Direct metric-profile / Ricci-stencil crosswalk

A follow-up bounded capture records the partial and full-profile metric fields and the exact pinned Ricci operands at cells 0–5, on the same initial geometry stage used for the RHS. The new run is [40155573527](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40155573527); all N=60/120/240 jobs passed. No extra evolution was performed.

| N | Partial-minus-full `B[0]` | Partial-minus-full `B[1]` | `delta B[1]/dr^2` | Fraction of (S_2/9) | Partial-minus-full TF Ricci at cell 0 |
|---:|---:|---:|---:|---:|---:|
| 60 | (+8.57\times10^{-11}) | (-1.5040\times10^{-4}) | (-8.4599\times10^{-5}) | 0.606 | (-9.8999\times10^{-5}) |
| 120 | (+2.20\times10^{-11}) | (-5.5434\times10^{-5}) | (-1.2473\times10^{-4}) | 0.893 | (-1.2854\times10^{-4}) |
| 240 | (+5.54\times10^{-12}) | (-1.5092\times10^{-5}) | (-1.3583\times10^{-4}) | 0.973 | (-1.3587\times10^{-4}) |

Here (B=X^6) and (S_2=-4\pi A^2=-1.2566370614\times10^{-3}) is the leading (r^2) coefficient of the declared Gaussian scalar's regular-origin source for this setup. Thus (S_2/9=-1.3962634016\times10^{-4}). The first-cell metric value (B[0]) already matches the full-profile reference to near roundoff; the first substantive mismatch appears at **cell 1**, after the first centre-to-centre trapezoid interval.

This agrees with the earlier analytic manufactured-source discriminator in [`tools/audit_origin_integral_quadrature.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/d6c6ffd35f5a698558b7f019f3dd5afb2c0b7909/tools/audit_origin_integral_quadrature.py), which derives (delta B_1=S_2\,\Delta r^2/9) for an even quadratic source when the first half-cell is integrated exactly but the next interval is still trapezoidal. That old test is a manufactured analytic control, not the live initializer itself; the present table is the live-profile measurement, and it approaches the same coefficient under refinement.

The cell-0 trace-free Ricci difference tracks (delta B_1/\Delta r^2) closely by N=240. Mechanistically, the (O(\Delta r^2)) metric error at cell 1 is differentiated by the centre Ricci stencil, creating an (O(1)) raw curvature error at cell 0. Dividing the resulting (A_a) regularity witness by (r_0^2\sim\Delta r^2) then explains why the normalized central quantity worsens as the mesh is refined. This is the strongest current explanation of the observed scaling; it is not a complete constraint-propagation proof.

## Matter-reprojection control

The partial-metric RHS with original conserved matter restored agrees with the reprojected-matter arm to within relative differences of (6.1\times10^{-7}), (1.24\times10^{-7}), and (2.93\times10^{-8}) for (N=60,120,240). The trace-free Ricci change is correspondingly negligible. The (L_3) matter-source contributions are also small compared with the problematic partial-metric (L_2) contribution.

This rules out matter reinitialization as the material cause of the gap. The dominant difference is the geometry profile entering the discrete Ricci stencil.

## Integrity checks and limits

The workflow passed the source/vendor pins, checked that all captured arrays are finite, verified the (A_a) (L_2) decomposition into the trace-free lapse Hessian plus trace-free Ricci to the configured closure bound, and confirmed that every copied input state remained unchanged by its one-step call. Protected production/vendor/dependency paths were unchanged.

The result is still a discrete initial-slice diagnostic. It does not establish the complete constraint-propagation identity, total stress-energy conservation, long-time stability, or scientific admission. Production initialization, vendor source, and production equations remain untouched; `scientific_admission=false`.

## Next discriminating step

Keep the full-profile construction as a research comparator and inspect the radial metric-profile difference and the exact Ricci-stencil operands across cells 0–2. Do not patch `v55_initial.py` based on this result alone. Any production-initializer change needs its own explicit approval and a provenance-backed closure argument.
