# CMC lapse centre-row repair — source derivation and bounded evidence

**Date:** 2026-10-12  
**Branch:** `diag/cmc-lapse-cellcenter-row-repair-2026-10-12`  
**Parent diagnostic:** `diag/cmc-uniform-quadrature-operator-2026-10-12` at `4715e8aed714909bba5c4adba92edc4ec77c00ee`  
**Status:** isolated numerical-layer candidate; not merged; no evolution or production-physics edits

## Derivation from the actual grid

The pinned grid is cell-centred: (r_i=(i+1/2)\Delta r), and its fourth-order derivatives use even-reflection ghost cells for the lapse. At the first cell, the production fourth-order operators therefore use

- (D_r\alpha_0=(-8\alpha_0+9\alpha_1-\alpha_2)/(12\Delta r))
- (D_r^2\alpha_0=(-14\alpha_0+15\alpha_1-\alpha_2)/(12\Delta r^2))

The old CMC solver instead made its first matrix row `alpha[0] - alpha[1] = 0`. That is not the centre symmetry condition for samples located at (\Delta r/2) and (3\Delta r/2); it forces equal values at two distinct positive radii and omits the lapse PDE at the first cell. The audit confirms this is consequential: on the uniform-quadrature family, the old solver yields `alpha[0] == alpha[1]` and an even-series fit residual that is still (4.32\times10^{-5}) at N=160. Its reported maximum elliptic residual `1.58e-3` is dominated by cell 0, where the solver was applying the equality row rather than the PDE.

## Candidate implementation

The candidate replaces the separate tridiagonal, second-order CMC matrix rows with the same parity-aware fourth-order first/second derivative rows used by the production K RHS. The centre cell is now a PDE row evaluated with even-reflection ghost weights. Cell (N-2) uses the vendor's one-sided fourth-order closure; cell (N-1) remains the existing `alpha=1` outer normalization and is not counted as a PDE row. The CMC target, (Q) definition, physical sources, outer normalization, positivity admission and model parameters are unchanged.

This is a numerical discretization consistency repair candidate, not a new physical law.

## Bounded acceptance

`tools/audit_cmc_lapse_cellcenter_row_repair.py` evaluates the same uniform-quadrature initial profiles at N=40/60/80/120/160 without advancing them. It reports the CMC row residual and production (\dot K-\dot K_{target}) on cells 0..N-2, lapse even-series fit, and signed (\dot A_a) centre/fixed-radius witnesses. The dedicated regression asserts that the same discrete rows close, lapse remains positive and the artificial equality of the first two cell values is removed.

Even if this gate passes, it does **not** admit the evolved machine. The existing regular-centre (\dot A_a=O(r^2)), momentum-constraint propagation, radiation admissibility at N=240/320, Misner–Sharp balance and convergence gates remain independent. `main` and the canonical GEAR archive remain untouched.
