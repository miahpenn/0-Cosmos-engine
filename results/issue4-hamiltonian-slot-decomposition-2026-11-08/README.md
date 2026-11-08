# Issue #4 Hamiltonian-slot decomposition (2026-11-08)

Read-only diagnostic. `scientific_admission=false`. No production code, tolerance or gauge change.

## Question

Brown Eq. (11b) has two Hamiltonian slots. The earlier Option C run (outputs
`option_c_hamiltonian_slot_N40-60-80.json`, not in the repo) found that putting the BSSN
Hamiltonian in those slots instead of the metric-only one changed H by 3–7% of max |H_BSSN|. A fit
of that difference against X²·D_r(Λ) gave a coefficient of about −0.03 to −0.04, not 1. The
difference was recorded as unexplained.

## Method

The pinned discriminator shows H_BSSN is affine in Λ, with δH_BSSN = X²·D_r(δΛ) (closure about
3e-19). With C = Λ − Λ_metric from the vendor's own constraint routine:

    H_BSSN − H_metric = X²·D_r(C) + remainder

The remainder vanishes in the continuum. `tools/issue4_hamiltonian_slot_decomposition.py` computes
both parts on the tool's manufactured state (r_max = 40) at N = 40, 60, 80, 120, 160, 240, 320,
twice: with C as manufactured, and with Λ set so C = 0. Cells 0/1 and the outer 4 cells are reported
separately from the interior.

Pins: production source `7f47669a1c0edfafef9fc5f27b2db80840d3f035`, vendor
`d6052d605673ce9d82cdc99b0df79855fa2ea215` (checked by the tool), diagnostic commit `b2618daa`.
Python 3.12.15, numpy 2.5.4, scipy 1.18.1. Result: `summary.json`, SHA-256
`df87a695df97042dfdfa0937370cc4728afc34328cc09b1b76c6da380f16252c`.

## Results

| N | interior H_BSSN − H_metric | Λ part X²·D_r(C) | remainder | remainder / max\|H_BSSN\| | remainder, cell 0 |
|---|---|---|---|---|---|
| 40 | 5.15e-06 | 3.82e-06 | 1.33e-06 | 4.5e-03 | 1.27e-05 |
| 80 | 6.25e-06 | 6.09e-06 | 1.60e-07 | 5.4e-04 | 3.09e-06 |
| 160 | 6.77e-06 | 6.74e-06 | 3.25e-08 | 1.1e-04 | 7.67e-07 |
| 320 | 6.92e-06 | 6.91e-06 | 7.76e-09 | 2.6e-05 | 1.92e-07 |

- The remainder is identical to the last digit with C as manufactured and with C = 0. All C
  dependence is in the Λ part.
- Observed order of the interior remainder between successive N: 3.27, 2.77, 2.38, 2.17, 2.08, 2.04.
  Cells 0 and 1 also converge at about second order.
- The earlier whole-grid maxima (1.947e-05, 1.243e-05, 1.002e-05 at N = 40, 60, 80) are the cell-0
  values of H_BSSN − H_metric in this run. They are the Λ part plus the cell-0 remainder.
- The earlier fit used Λ, not C, as the regressor. Λ ≈ Λ_metric is large (max about 6.1e-4)
  and not what H_BSSN − H_metric depends on, which accounts for the fitted coefficient near zero.

## Reading

The difference between the two Hamiltonian slots is X²·D_r(C), the connection constraint's
contribution to H_BSSN, plus a discretization remainder that converges at about second order. On
this manufactured state the Λ part does not vanish under refinement because C is nonzero by
construction. Brown Eq. (11b) carries C explicitly in its other terms, so the metric-only slot
avoids counting that contribution twice; using H_BSSN in the slots requires subtracting X²·D_r(C).

## Where the second order comes from (`order_source.json`)

`tools/issue4_hamiltonian_slot_order_source.py` (diagnostic commit `5f367be2`, SHA-256 of
`order_source.json` `27cf15253e3b9ed6fa1916a50d533da7c6f8083a807788b6a7f860cf471aa3eb`) measures each
piece against the exact continuum values of the analytic state:

- Vendor H_BSSN with exact Λ, and vendor C with exact Λ: fourth order at every cell, including 0 and 1.
- Metric-only H: fourth order at fixed radius (r = 1, 2, 5, 10), second order at cells 0–2. Its error
  behaves like h⁴/r², from the 1/R and 1/R² quotients of the warped-product Ricci near the centre.
- Candidate slot H_BSSN − X²·D_r(C): fourth order at every cell (2.7e-10 at N = 320, against 1.9e-7
  for metric-only H at cell 0), and its response to the Λ-only mutation is at roundoff (≤ 1.5e-18).

The Gate A decision that follows is recorded in `docs/research/issue4_gate_a_conventions_2026-11-04.md`, §11.
