# Flat LCDM Known-Physics Test Results

**Status: PASS for the homogeneous COSMOS background lane only.**

- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38202378125
- Artifact: `flat-lcdm-smoke` (artifact ID `11725582645`)
- Test code: `diagnostics/flat_lcdm_smoke.py`
- Commit containing test/workflow: `c0a8934d8d73a2609596c5196a72487978193ba9`

## Setup

Dimensionless flat dust + cosmological constant model:
- (H_0=1)
- (Omega_m=0.30)
- (Omega_Lambda=0.70)
- radiation = shear = baryonic density = 0
- dark-sector coupling (eta=0)
- scalar potential made constant

Analytic reference:
`H(a)^2 = 0.3/a^3 + 0.7`

At (t=1), the exact solution gives:
- (a=2.444686302915898)
- (H=0.8488421367015576)

## Results

| Steps | dt | Absolute error in a | Absolute error in H | Max Friedmann residual | Max dust-conservation error |
|---:|---:|---:|---:|---:|---:|
| 100 | 0.01 | 1.23866e-10 | 3.67414e-10 | 3.86979e-10 | 1.05151e-8 |
| 200 | 0.005 | 7.60636e-12 | 2.27056e-11 | 2.39182e-11 | 6.49613e-10 |
| 400 | 0.0025 | 4.69402e-13 | 1.40954e-12 | 1.48392e-12 | 4.03673e-11 |

Halving the timestep reduced errors by:
- (a): 16.28 then 16.20
- (H): 16.18 then 16.11

These ratios are consistent with fourth-order RK4 convergence. The finest-step result matches the analytic flat-LCDM solution very closely; Friedmann closure and dust conservation also improve as the timestep is reduced.

## What this proves—and does not prove

It **does prove** that the existing `engine.cosmos.cosmos_rhs` homogeneous background equations, under this controlled dust+Λ limit, reproduce the expected standard flat-LCDM solution with fourth-order convergence under an independent RK4 driver.

It **does not prove** that the full production kernel is validated. This test does not evolve the spatial BSSN geometry, the local S/D fields, the coupled matter fluxes, or the production PIRK2 stepper. Keep those as separate gates.

## Next check

The separate near-vacuum production-geometry RHS test did **not** pass: the geometry blocks were finite, but the (K) primary-l3 RHS reached (7.41827\times10^{-5}) in a nominally flat state. That deserves a targeted audit before claiming the full machine preserves static flat spacetime. This is a real diagnostic discrepancy, not yet a proven production-evolution failure.
