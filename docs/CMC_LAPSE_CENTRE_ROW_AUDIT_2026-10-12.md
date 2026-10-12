# CMC lapse centre-row and operator-identity audit — 2026-10-12

**Run:** [38296556473](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38296556473) — successful  
**Artifact:** [11768278320 — cmc-lapse-centre-row-audit](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38296556473/artifacts/11768278320)  
**Artifact ZIP SHA-256:** e02c2c0f384d20692fc02bb3ee771409d0743a2479dc39276ee3a14b1156101f  
**Diagnostic commit:** 95606bb4db9260286191aec4a513a0430aed386a  
**Frozen vendor:** d6052d605673ce9d82cdc99b0df79855fa2ea215  
**Production main:** unchanged. No evolution run or production edit.

## Scope

Read-only audit of two archived initial states:
- R80 shared-core state: N=160, rmax=80, dr=0.5, with outer lapse fixed to the archived R160 value at the R80 boundary.
- Native R160 state: N=320, rmax=160, dr=0.5, with alpha(R)=1.

The script solves the pinned CMC lapse equation on each frozen state, evaluates the exact discrete elliptic rows, reconstructs the independent production K RHS, compares the production lapse-geometric K term to the CMC operator's implied term, and fits alpha as an even polynomial in r over the first eight cells. The outer normalization cell is explicitly not treated as an elliptic PDE row.

## Results

### 1. The CMC lapse row closes

- Both lapse profiles remain finite and positive.
- Maximum production K-RHS minus target over interior PDE rows: 6.04e-15.
- Maximum independently reconstructed CMC elliptic residual over interior rows: 7.32e-15.
- Solved alpha differs from the archived alpha by only 8.58e-11 max / 6.16e-11 RMS on boundary-matched R80 and 3.35e-12 max / 2.44e-12 RMS on native R160.

### 2. The centre lapse profile is smooth under an even-series fit

Using alpha(r)=a0+a2*r^2+a4*r^4+a6*r^6 fitted to the first eight cell centres:

- R80 extrapolated alpha(0): 1.0144792546397763.
- R160 extrapolated alpha(0): 1.0148363888527776.
- Maximum fit residual across the eight fitted cells: 4.49e-8 for both.
- First-cell value minus the fitted polynomial at r=0.25: 3.17e-8 for both.
- The first three lapse-Laplacian brackets vary smoothly: about -1.270e-3, -1.244e-3, -1.195e-3.

These are fixed-dr, two-domain checks, not a formal refinement-order proof. They do not certify every evolved state or the separate uniform-quadrature initial arm.

### 3. CMC and production lapse operators agree at the audited cells

The production vendor's lapse-geometric contribution to K RHS agrees with the independently reconstructed CMC operator term to roundoff:

- At the first cell, R80: production 1.269887358833915e-3; CMC reconstruction 1.269887358831954e-3; difference 1.96e-15.
- At the first cell, R160: production 1.270659413056301e-3; CMC reconstruction 1.270659413056116e-3; difference 1.85e-16.
- Across the recorded first cells and fixed radii r≈2,3,4,5, the absolute discrepancy remains at or below about 2e-15.

This is the most important result of this audit: **for these archived initial states, the CMC centre row and the production lapse-geometric K term are mutually consistent to roundoff.** The lapse central row is therefore not presently supported as the source of the earlier centre-series defect on these states.

### 4. The normalization cell remains a distinct boundary/global-average seam

- R80 outer-cell K-RHS residual: -4.588278121699964e-12.
- R160 outer-cell K-RHS residual: -6.934447319747322e-14.

Those cells impose the lapse normalization and are not elliptic PDE rows, yet their K RHS still enters the global target average. This aligns with the separate nine-target rank probe, which found a near-identity target closure with dH/dKdot=-2.4269e-10 and H around -8.58e-14 across the target interval.

## Interpretation

The read-only audit **weakens the hypothesis that the archived-state CMC lapse centre row itself is defective**. The alpha profile fits an even centre series closely, the exact elliptic rows close, and the production lapse term matches the CMC operator. Do not patch the centre row based on the earlier uniform-arm K-RHS jump alone.

The CMC global target closure is a separate concern: on the frozen R80 slice the target average is nearly algebraically redundant with the interior PDE. That requires derivation against the authoritative CMC condition, not a fitted adjustment or silently dropping the outer normalization cell.

The prior Aa/K signed-term audit remains relevant but is not superseded. It used a different uniform-quadrature initial arm and found a persistent trace-free lapse-Hessian contribution and a large first-cell K-RHS jump. This current audit did **not** reconstruct that uniform-quadrature arm; it therefore cannot dismiss the discrepancy across all initializers.

## Next discriminating action

1. Run the same operator comparison on the exact uniform-fourth-order quadrature initial profiles for N=40/60/80/120/160. Compare the CMC lapse row, production lapse RHS term, trace-free Hessian, Ricci term and K source on the *same* frozen state.
2. If those also agree, shift the search to how the uniform initializer constructs geometry/matter and to the centre regularity identities of Aa/K/X/a/b, while keeping the target-closure rank issue separately open.
3. If they disagree, isolate the exact derivative/row difference before any repair.
4. Keep the radiation realizability failures at N=240/320 as a separate blocker. No centre projection, damping, floors/clamps, source/gauge change, or long production campaign until both the central identity and radiation admission gates pass.

## Reproduction

Source: diagnostics/cmc_lapse_centre_row_audit.py  
Workflow: .github/workflows/cmc-lapse-centre-row-audit.yml  
Rank companion: docs/CMC_TARGET_CLOSURE_RANK_AUDIT_2026-10-12.md  
Artifact ZIP is preserved in the local library under ARTIFACTS_2026-10-12/cmc-lapse-centre-row-audit-operator-compare.zip.
