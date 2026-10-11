# 0* Cosmos Engine — CMC Lapse Audit Results and Continuation Handoff

**Date:** 2026-10-11  
**Active branch:** `diag/independent-metric-ricci-audit`  
**Working rule:** one expensive trajectory at a time; diagnose before physics; no tuning or silent physics changes.

## What just completed

The frozen-state CMC crop audit completed with success in [workflow run #38117001920](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38117001920). Artifact ID `11694568292`, name `cmc-lapse-frozen-state-crop-audit`, SHA-256 `04949c7e0332c495385f6afef231b54270225892b1d124f2c19973c7d5b4872d` (verified against the GitHub artifact digest). It contains 49 saved full-state slices and the complete lapse-profile crop comparisons.

Detailed report: [cmc-lapse-frozen-state-crop-audit-findings-2026-10-11.md](cmc-lapse-frozen-state-crop-audit-findings-2026-10-11.md).

## Results that matter

1. The domain dependence of alpha on a **fixed evolved state** is almost entirely the direct effect of the existing outer lapse normalization alpha(R)=1. If the cropped solve uses the full-domain lapse value at its crop boundary, the full-domain solution restricted to that crop is reproduced to maximum absolute difference `2.2639e-11` for R=80 and `1.9454e-12` for R=120.
2. The cropped vs full-domain `Kdot_CMC` change is tiny: the maximum relative difference is `8.08e-8` at R=80 and `8.71e-8` at R=120. Solving the same cropped operator with the full-domain target changes alpha by at most `4.18e-9` and `1.13e-8`, respectively.
3. Discrete operator residuals are approximately `1e-14`. The test finds no implementation defect in the assembled CMC solve/crop path and does not justify changing the gauge.
4. The large lapse differences between independently evolved domains remain. At t=12, the original R=80 run differs from R=160 by `0.116139` at r=39.75; at R=120 the independent-run difference is `0.038486` at r=47.25. The current frozen-state crop shift at those points is only `+0.008267` and `-0.010898`, respectively. Most of the discrepancy therefore remains associated with differences in initial/evolved state histories; the prior smaller-domain artifacts did not save the full per-cell state arrays needed to attribute that discrepancy.

## What must not be inferred

- Do not call the residual ridge a physical D-front, causal propagation, bounce, turnaround, or cycle.
- Do not patch the lapse equation or normalization based on this result. alpha(R)=1 is the current gauge normalization; its domain sensitivity is measured, not shown to be a bug.
- Do not attribute the remaining independently evolved state differences to one field yet. The smaller-domain geometry/scalar/matter histories need to be compared directly.

## Current single active run

**Workflow:** [Run #38150925743 — Diagnostic - R80 evolved-state capture](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38150925743).

At launch status-check, compile and configuration validation passed; dependency installation was in progress and numerical evolution had not started. This is the only new expensive run. It reproduces the untouched production settings at N=160, R=80, dr=0.5 through t=12 and saves every cell of geometry, scalar fields, and conservative matter variables at the same 49 sample times. No simultaneous R=120 run should be launched.

## Next analysis once R80 capture completes

Compare its full state at each sample time with the archived R=160 trajectory cropped to the first 160 cells. Compare all arrays independently: geometry (`a,b,X,alpha,beta,Aa,K,Lambda,B`), scalar fields (`S,PS,D,PD,phi,Pi`), and the conserved variables of dark matter, baryons, and radiation (`rest,energy_t,momentum_r`). Record max absolute and RMS difference over the full shared radius and over r<=20, including peak radius/time and first time a difference grows beyond roundoff. Compare H_vendor profiles too.

Because both simulations use the same dr, dt, model, and sampling rule, the timestamps should match. The R=160 full-state snapshots are already present in artifact `11694568292`; they do not need to be regenerated. This one R=80 capture will determine whether the mismatch starts in the initial data or emerges in a particular evolved field. Only after this comparison should we decide whether R=120 needs a full-state replay.

## Fixed production settings for current and next diagnostic

N=160/R=80 for the active run; CFL=0.0075; D amplitude=1e-4; scalar amplitude=0.01; width=7; radiation enabled; final t=12; samples at initial state plus quarter-time targets. Production physics, defaults, gauge, projection, evolution boundary conditions, and algorithm remain unchanged.

## Source and implementation references

- CMC solver: `engine/cmc_gauge.py`.
- Frozen-state crop diagnostic: `diagnostics/curvature_consistency/cmc_lapse_frozen_state_crop_audit.py`.
- Three-domain residual synthesis: `diagnostics/curvature_consistency/strong-D-domain-size-three-domain-synthesis-2026-10-11.md`.
- Full-state R80 capture script: `diagnostics/curvature_consistency/strong_D_evolved_state_r80_capture.py`.
- Diagnostic workflow: `.github/workflows/strong-D-domain-size-residual-audit.yml`.
