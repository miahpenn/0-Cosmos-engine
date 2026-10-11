# Strong-D production metric-Ricci/H audit — results

Date: 2026-10-11  
Branch: `diag/independent-metric-ricci-audit`  
Workflow: [run 38098552108](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38098552108)  
Commit: `769f0bb45167c4616e2e315fa331395cb3d5f319`  
Artifact: [strong-D-production-metric-ricci-H-audit](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38098552108/artifacts/11686323614)  
Artifact ZIP SHA-256: `499631b43579903740d5e9a4fe5840ad57ae39d01f980b6cd1333447a6eecd1b`

## Configuration

- Unchanged production initializer and existing centre projection ON
- (N=160), (r_{\max}=80), (\Delta r=0.5)
- D amplitude (10^{-4}), local amplitude (0.01), width 7
- Radiation ON, CFL (0.0075), (Delta t=0.00375)
- 6,000 accepted steps to (t=22.5)
- Registered checkpoints: (t=0,0.01125,4.00125,7.99875,12,16.00125,19.99875,22.5)
- 0 cycle events, 0 handoffs; finite-state checks passed

All algebraic diagnostic gates passed: H-term reconstruction, curvature replacement closure, connection identity closure and vendor connection-constraint reconstruction. The report is diagnostic-only; it makes no production change or admission claim.

## Central-cell history

| Sample coordinate time | (H_{vendor}) | (H_{metric}) counterfactual | (R_{vendor}-R_{metric}) | centre (|C_\Lambda|) |
|---:|---:|---:|---:|---:|
| 0 | (+1.07904\times10^{-3}) | (+9.92790\times10^{-4}) | (+8.62473\times10^{-5}) | (3.90041\times10^{-5}) |
| 0.01125 | (+1.07865\times10^{-3}) | (+9.92490\times10^{-4}) | (+8.61590\times10^{-5}) | (3.89730\times10^{-5}) |
| 4.00125 | (+1.14525\times10^{-3}) | (+1.23646\times10^{-3}) | (-9.12090\times10^{-5}) | (5.38515\times10^{-5}) |
| 7.99875 | (+1.08107\times10^{-3}) | (+1.13665\times10^{-3}) | (-5.55787\times10^{-5}) | (4.99537\times10^{-5}) |
| 12 | (+1.09330\times10^{-3}) | (+1.14746\times10^{-3}) | (-5.41549\times10^{-5}) | (5.20569\times10^{-5}) |
| 16.00125 | (+1.09774\times10^{-3}) | (+1.15069\times10^{-3}) | (-5.29449\times10^{-5}) | (5.27124\times10^{-5}) |
| 19.99875 | (+1.10014\times10^{-3}) | (+1.15215\times10^{-3}) | (-5.20145\times10^{-5}) | (5.30615\times10^{-5}) |
| 22.5 | (+1.10203\times10^{-3}) | (+1.15326\times10^{-3}) | (-5.12322\times10^{-5}) | (5.33357\times10^{-5}) |

The central vendor residual begins at (1.079\times10^{-3}) and ends at (1.102\times10^{-3}), a net change of about 2.13%. Thus the large central offset is already present at the production initial slice and stays similar in magnitude over this run. It is not appropriate to describe the whole centre value as evolution-generated drift.

The curvature gap changes sign and ends at roughly 4.65% of the central vendor residual magnitude. The metric-derived counterfactual has its own large residual at both initial and final samples. Curvature-route disagreement is measurable, but it does not explain the central residual by itself.

## Important second finding: off-centre residual grows

Although the global maximum remains at the centre, the maximum absolute vendor H on cells 2 and beyond grows from (1.441\times10^{-5}) initially to (4.638\times10^{-4}) at (t=22.5), about 32 times larger. So separate two phenomena:
- a large, approximately persistent central residual that pre-exists evolution; and
- growing noncentral Hamiltonian residual outside the first two cells.

The connection constraint is also spatially concentrated near the centre: its maximum remains in cell 0 or cell 1, with whole-grid maximum (8.72\times10^{-5}) at the final sample; the central value at the end is (5.33\times10^{-5}). Do not assume the connection residual fully accounts for the H offset; the same-state Ricci comparison did not support that.

## Final central signed terms

At (t=0):
- (R_{vendor}=+3.58617\times10^{-3})
- ((2/3)K^2=+1.92823\times10^{-4})
- (-16\pi\rho= -2.69996\times10^{-3})
- extrinsic-A term is negligible at the centre
- net (H_{vendor}=+1.07904\times10^{-3})

At (t=22.5):
- (R_{vendor}=-4.35326\times10^{-3})
- ((2/3)K^2=+5.98072\times10^{-5})
- (-16\pi\rho=+5.39548\times10^{-3}), because the reported total projected centre density is negative
- extrinsic-A term is negligible at the centre
- net (H_{vendor}=+1.10203\times10^{-3})

This is term accounting, not a causal explanation. The sign change in projected density and curvature is an observed property of this numerical trajectory and must be audited through the existing stress-energy and geometry definitions before physical interpretation.

## Next test

Trace the production initialization pipeline at strong D and the same (N,r_{\max}): record the Hamiltonian terms before centre regularity projection, immediately after it, after CMC-lapse/gauge setup, and from the public `kernel.initialize` entry point. If the (1.08\times10^{-3}) residual is already present before projection, that rules out the projection as the origin of the initial offset, while leaving the discrete radial reconstruction, central stencil and normalization conventions to be tested separately.
