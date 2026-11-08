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

## Brown vacuum residual with each slot (`brown_residual_by_slot.json`)

`tools/issue4_brown_residual_by_hamiltonian_slot.py` (diagnostic commit `ac418a59`, SHA-256 of
`brown_residual_by_slot.json` `41769ab7aa4c6c98f4e54ff187b1eef5cff18813455482333ffc208c023d2287`)
compares the measured two-sided M_r tangent along the pinned vacuum RHS with Brown Eq. (11b),
N = 40…1280, ε = 2e-5 and 1e-5 (largest change between ε: 3.7e-12).

- **Cell 0.** Metric-only slot: residual converges at about first order (−7.5e-8 at N = 640).
  H_BSSN − X²·D_r(C) slot: about third order (3.2e-10 at N = 640).
- **Fixed radius.** With either slot, the residual tends to the same nonzero limit, for example
  +2.38e-9 at r = 2 (local signal −4.24e-6) and −3.44e-10 at r = 5 (signal −1.99e-6). This limit
  does not depend on the slot.
- **The limit is one term.** A least-squares fit of the N = 640 residual (0.5 < r < 30) against
  candidate shapes gives α X² [1.005·D_r(u) + 3.016·u/r], u = (D_r a / 2a)·C, leaving 0.3% of the
  residual norm unexplained. With the coefficients fixed at exactly 1 and 3, the remaining residual
  converges at fourth order at r = 1, 2 and 10 (for example 3.3e-9, 2.1e-10, 1.3e-11 at r = 1 for
  N = 160, 320, 640) and at about third order at cell 0 and in the interior maximum, down to the
  roundoff floor near 1e-12 at N = 1280.

u = (D_r a / 2a)·C is the vendor-specific undifferentiated-C term in Δ_Λ R^r_r recorded in Gate A
§9.1.

## Derivation of the term (`vendor_connection_term_derivation.json`)

From the pinned vendor code: Λ enters α(R^r_r − R/3) in the A_a rate as α X²[⅔·D_r(C) + u]; the K
rate has no Λ; and the vendor's M_r = D_r(A_a) − ⅔D_r(K) + 6A_a·φ′ + (3/2)A_a·q. So the A_a-rate piece
δ = α X² u changes the M_r tangent by exactly

    δṀ_r = D_r(δ) + (6φ′ + (3/2)q)·δ,   φ′ = −D_r(X)/(2X),  q = 2/r + D_r(b)/b,

whose leading part is the fitted α X²[D_r(u) + 3u/r]. Brown Eq. (11b) carries only the ⅔·D_r(C) part.

`tools/issue4_vendor_connection_term_derivation.py` (diagnostic commit `b9b8d76d`, SHA-256 of
`vendor_connection_term_derivation.json`
`83f7db7f3cbb188f23eed091ab0ff95f0d934d3417be850c78c608d661124296`), N = 40…1280:

- Removing δ from the actual A_a rate and re-measuring changes the tangent by exactly δṀ_r, to at most
  2e-13.
- Measured − Brown (Gate A slot) − δṀ_r converges at fourth order at r = 1 and 2 and about third order at
  cell 0 and in the interior maximum, to the roundoff floor at N = 1280 (r = 2: +2.38e-9 before,
  +2.4e-12 after).
- δṀ_r and the fitted form agree to ≤ 1.5e-11.

- The locked evaluator's geometric predictor (`tools/issue4_full_bianchi_evaluator.py`, transcribed with
  ρ = j = 0) closes the same vacuum tangent identically (r = 2: 1.21e-11 at N = 640) and differs from
  Brown + δṀ_r only by discretization (≤ 6.3e-9 at N = 40, decreasing). It already uses the Gate A slot
  and the full W connection term, including u.

The term is derived, not fitted. Gate A §11.1 records it.
