# Trace-free Ricci / proper-radial reconciliation — 15 October 2026

**Disposition:** diagnostic-only; scientific admission remains **NOT ADMITTED**. This refresh reuses the earlier closed, unmerged Ricci observer rather than creating a competing implementation.

## Exact source pins
- Production source baseline: `fcb32b6b4025d5674e0560b8fb24581933f7f23c`; production-kernel blob: `9685e76bb5a6f26d6bf383e240de24f37523f265`.
- Six engine source blobs are checked against both baseline and diagnostic head.
- Vendor submodule expected: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Any engine/vendor/dependency diff from the source pin blocks the diagnostic.

## Measurement
The observer captures exact inputs to four real primary-L2 geometry calls in one bounded production step at (N=40,60,80,120,160), retaining signed profiles for vendor BSSN trace-free Ricci, an independent areal-radius/proper-radial reconstruction, and direct-D2 lapse Hessian. Production operators and outputs are returned unchanged.

With (R_{m areal}=rsqrt b/X) and (q=X/sqrt a), the independent identities are (R^ell{}_ell=-2R_{m areal}^{-1}d^2R_{m areal}/dell^2) and ({}^{(3)}R=-4R_{m areal}^{-1}d^2R_{m areal}/dell^2+2[1-(dR_{m areal}/dell)^2]/R_{m areal}^2). This route does not use evolved `Lambda` in its Ricci formula, but shares the pinned fourth-order derivative family; independent representation, not independent discretization.

## Regional split and negative control
The report separates centre cells (i=0,1) from the fixed **areal-radius** annulus (5le R_{m areal}le8), retaining first-six-cell, (rge5), and fixed-coordinate-radius samples only as supporting comparisons.

A manufactured smooth constant-positive-curvature 3-space ((R_{m areal}=r), (dell^2=dr^2/(1-kr^2)), (k=10^{-3})) is a negative control: equal radial/angular Ricci eigenvalues imply vanishing trace-free Ricci. Regression checks both representations near the centre and confirms that annulus data exist.

The Hessian difference is reported signed, absolute, and relative; no arbitrary absolute cutoff drives acceptance. Only finiteness and same-stage algebraic closure gate diagnostic integrity.

## Limits
One ordinary step per resolution, not a long campaign. A discrepancy localizes representation/operator compatibility; it does not prove which equation is wrong. No production physics, vendor code, gauge, boundary, source normalization, coefficient, campaign default, or admission policy changes. Numeric results/full signed profiles are emitted to `results/tracefree-ricci-proper-radial-2026-10-15/`.
