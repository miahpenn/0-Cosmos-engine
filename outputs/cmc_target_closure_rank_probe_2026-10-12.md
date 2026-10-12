# CMC target-closure rank probe — initial slice

**Date:** 2026-10-12  
**Run:** [38294955148](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38294955148)  
**Artifact:** [11768056222 — cmc-target-closure-rank-probe](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38294955148/artifacts/11768056222)  
**Diagnostic source commit:** `1985daaee34e89fe1da68968f0518249f79cd3bc`  
**Workflow commit:** `be7dac470eca22a1267dbbeac7776719afbd88d6`  
**Artifact ZIP SHA-256:** `64cc1666a2010d73b462cac48c78706ac48be1a1d02fea7de343bd09172c4e24`  
**Extracted report SHA-256:** `8ce2d07e90f0a760fd0683efb72e0af6530893b13f03776dc1ed67efb240b1ab`

## Result

Nine fixed Kdot targets were tested on the same frozen R80 (t=0) state, with all non-lapse fields fixed and the outer lapse boundary matched to R160 at `r=79.75`. The target span extends two native target gaps below/above the R80/R160 values, from (6.1294749186\times10^{-5}) to (6.3028914972\times10^{-5}).

The solved target output tracks the input target with slope **0.999999999757**. The closure (H=\mathrm{target\_kdot}(\alpha_{\mathrm{solved}})-\dot K_{\mathrm{trial}}) remains between (-8.6051\times10^{-14}) and (-8.5613\times10^{-14}), range (4.3788\times10^{-16}). Linear-fit slope (dH/d\dot K=-2.4269\times10^{-10}). The lapse profile itself is affine in trial Kdot to RMS error (4.48\times10^{-15}), max error (6.66\times10^{-15}).

PDE-row K-RHS residuals remain below (3.71\times10^{-15}); the outer normalization-cell residual reaches (4.60\times10^{-12}). The scalar closure equals the weighted mean pointwise residual to (1.14\times10^{-20}).

## Interpretation

This materially strengthens the hypothesis that the current global target projection is nearly redundant with the elliptic PDE rows and boundary closure over the tested range: solving the lapse for almost any trial Kdot returns a lapse whose projected target is almost exactly that same Kdot. This supports a near-null direction / family of near-stationary target–lapse candidates, rather than the global projection independently selecting a unique Kdot.

It is **not** proof of a continuum of exact roots or of a physics error. The authoritative CMC derivation must establish whether an independent condition is supposed to select Kdot. Do not patch the target projection or omit the boundary row by intuition.

## Next step

Write the coupled discrete residual explicitly: (1) the (N-1) CMC PDE rows, (2) the outer lapse-normalization row, and (3) the proposed global target condition. Derive its rank/Schur complement analytically or with a small finite-difference Jacobian. Determine whether the target condition adds an independent equation under the authoritative criterion. Keep diagnostic-only; no expensive evolution or production source change.
