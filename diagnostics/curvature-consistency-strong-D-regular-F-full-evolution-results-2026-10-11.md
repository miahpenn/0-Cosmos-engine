# Strong-D Production vs. Regular-F Cubic/Gauss Full-Evolution Diagnostic — Results

**Date:** 2026-10-11  
**Repository:** `miahpenn/0-Cosmos-engine`  
**Branch:** `diag/independent-metric-ricci-audit`  
**Workflow:** [Run #38100254183](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38100254183)  
**Workflow commit:** `3a9c7585d15bdddb49f1313b80e8316ce793a32b`  
**Artifact:** [strong-D-production-vs-regular-F-full-evolution](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38100254183/artifacts/11688328331)  
**Artifact ZIP SHA-256:** `f2416cc062c7cf8348103582a3a4eb75fceae3d18d4130cd6d82a94cfff77966`  
**Extracted `report.json` SHA-256:** `11e3262d64fdafdc8d2210201fb3a21dab019806c23307d27d7d06fac445ca8c`

## Executive result

The paired run completed successfully for both trajectories: 6,000 accepted production-kernel steps per case, from (t=0) to (t=22.5), with no exceptions and no diagnostic-gate failures. The Hamiltonian term-reconstruction and curvature-replacement closure checks passed at all eight samples.

The experiment cleanly separates two effects:

1. **The production initial radial quadrature is strongly implicated in the large persistent central Hamiltonian residual.** Replacing only the initial (B)-integral quadrature in a diagnostic shadow reduces central (H) from (+1.0790\times10^{-3}) to (-2.997\times10^{-7}) at (t=0). The shadow central residual remains far smaller through the full evolution: (-8.017\times10^{-6}) at (t=22.5), versus (+1.1020\times10^{-3}) in production — a factor of about 137.5 in absolute central residual.

2. **That initial improvement does not eliminate the growing off-centre Hamiltonian residual.** The two trajectories' maximum (|H|) outside the first two cells becomes almost identical by (t=12) and remains so through (t=22.5). This points to an additional residual-growth mechanism shared by the unchanged evolution kernel; it does not, by itself, identify that mechanism.

No production code, source equation, projection, gauge, boundary condition, evolution kernel, or default was changed. The cubic/Gauss result is an isolated diagnostic counterfactual, not a production patch or admission.

## Configuration and controls

- (N=160), (r_{\max}=80), (\Delta r=0.5)
- Strong-D amplitude (10^{-4}); scalar amplitude (0.01); width (7)
- Radiation ON; CFL (0.0075); requested final coordinate time (22.5)
- Both cases advanced through the unchanged `V55ProductionKernel.step`
- **Production case:** untouched `V55ProductionKernel.initialize`
- **Shadow case:** only the initial (B)-radial quadrature differs: regular-origin interval plus piecewise cubic interpolation of (F(r^2)) integrated with five-point Gauss-Legendre; then the same existing centre projection, CMC lapse solve and gauge reset

Both trajectories have the same recorded D amplitude evolution to the displayed precision, and their proper-time and e-fold outcomes are extremely close. That supports the conclusion that this counterfactual mainly removes an initial constraint offset rather than manufacturing a different large-scale evolution.

## Matched checkpoint results

(H_c) is central vendor Hamiltonian residual. “Off-centre max” is the maximum absolute vendor (H) over cells 2 and outward (the first two cells are excluded). Times shown are actual sampled times, which can differ slightly from requested labels because the integrator lands on its fixed time step.

| Actual (t) | Production (H_c) | Cubic/Gauss (H_c) | Production off-centre max (|H|) | Cubic/Gauss off-centre max (|H|) |
|---:|---:|---:|---:|---:|
| 0.00000 | +1.079037e-3 | −2.996607e-7 | 1.441033e-5 | 1.906222e-7 |
| 0.01125 | +1.078649e-3 | −2.996693e-7 | 1.440202e-5 | 1.906785e-7 |
| 4.00125 | +1.145248e-3 | −1.928781e-5 | 1.236833e-4 | 4.808250e-6 |
| 7.99875 | +1.081074e-3 | −6.820359e-6 | 1.012021e-4 | 5.944199e-5 |
| 12.00000 | +1.093302e-3 | −7.588420e-6 | 2.829182e-4 | 2.828119e-4 |
| 16.00125 | +1.097745e-3 | −7.812921e-6 | 4.010170e-4 | 4.010076e-4 |
| 19.99875 | +1.100135e-3 | −7.922360e-6 | 4.506958e-4 | 4.506866e-4 |
| 22.50000 | +1.102028e-3 | −8.016607e-6 | 4.638238e-4 | 4.638125e-4 |

At the final checkpoint the off-centre maxima differ by only (1.125\times10^{-8}), around (2.4\times10^{-5}) of their common magnitude. The shadow off-centre maximum grows from (1.906\times10^{-7}) initially to (4.638\times10^{-4}) at the end (about (2.43\times10^3) times larger). This is a separate and material signal: correcting the initial central offset is not sufficient to keep the full-grid constraint residual small.

## Independent-identity and geometry checks

- Every sampled H-term reconstruction closure passed; recorded maximum discrepancy was zero in these evaluations.
- Every curvature replacement identity passed; maximum recorded discrepancy was below (9\times10^{-19}), versus a tolerance of about (4.55\times10^{-13}).
- At (t=22.5), production central curvature gap (R_{vendor}-R_{metric}=-5.1232\times10^{-5}); shadow gap (+2.1454\times10^{-6}).
- At (t=22.5), production centre connection constraint (C_\Lambda=-5.3336\times10^{-5}); shadow (+7.5979\times10^{-7}). Maximum absolute connection constraint is (8.7191\times10^{-5}) in production and (4.0759\times10^{-5}) in the shadow.
- Both trajectories reached (t=22.5) with positive geometry fields and no numerical exception. Both recorded zero cycle events and zero handoffs.
- Both show severe lapse collapse during the run: minimum lapse reaches about (4.6732\times10^{-3}) near (t=20), rising to about (6.3008\times10^{-3}) at (t=22.5). Maximum lapse reaches about (1.2484). This is shared by both cases and deserves continued monitoring; completion is not a claim that the strong-D geometry is physically admitted.
- Final proper time: production (\tau=6.9833262997), shadow (\tau=6.9833400017), difference (1.37\times10^{-5}). Final e-fold counts both equal about (0.1130544).

## Interpretation guardrails

**Supported:** the initial (B)-integral quadrature is a major contributor to the large central Hamiltonian/connection mismatch observed in the production initializer. The independent regular-(F) integral and cubic/Gauss cross-checks agree, and the paired trajectory carries the central improvement forward.

**Also supported:** a substantial off-centre residual develops under the unchanged evolution scheme even when the initial central Hamiltonian residual is approximately (3\times10^{-7}). By (t=12), production and shadow off-centre maxima are essentially equal. That residual-growth question needs its own controlled diagnosis.

**Not supported yet:** that this quadrature change is sufficient to repair the full constraint system; that the off-centre drift is specifically temporal truncation error rather than spatial/evolution/constraint propagation; or that the current strong-D run demonstrates a bounce, turnaround, cycle or production-ready trajectory.

## Recommended next diagnostic

Run a **matched half-CFL refinement** on the same (N=160), (r_{\max}=80) production/shadow pair, with CFL reduced from (0.0075) to (0.00375), and stop at (t=12). Compare the (t\approx4,8,12) off-centre H maxima and central residuals against this baseline. This isolates time-step sensitivity of the shared growth while keeping the same spatial grid, initializer comparison and evolution equations. It is a diagnostic-only test; no physics patch should be made from the current result alone.
