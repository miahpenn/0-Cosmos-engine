# Full Production Homogeneous Cosmology Smoke Test

## Initial full-kernel result: PASS

- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38202631408
- Test: `diagnostics/production_cosmology_smoke.py`
- Initial test commit: `c8b1334ec5202ab413c52b91731df5c58712f5b2`
- Artifact: `production-cosmology-smoke`, ID `11726231253`

This run used the actual `V55ProductionKernel.initialize` and one actual `V55ProductionKernel.step`, at resolution 64, (r_{max}=32), and (dt=10^{-4}). The archive operating point was used with local S/D amplitudes zero and radiation disabled. The native exponential COSMOS scalar potential remained active and dynamic.

### Observed

- Initial (H_{eff}=0.005668707332082492)
- Independently recomputed initializer (H_0=0.005668707332082496)
- Initial absolute difference: (4.34\times10^{-18})
- After one production step, (H_{eff}=0.005668705253181093)
- Post-step Hamiltonian residual max: (7.08\times10^{-13})
- Post-step momentum residual max: (5.01\times10^{-13})
- Post-step connection residual max: (1.89\times10^{-20})
- All evolved fields finite; workflow passed.

This is a positive full-step smoke test, not full validation of the production engine.

## Stronger reference comparison added

The follow-up comparison **passed**:
- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38202689236
- Production (H_{eff}) after one step: (0.005668705253181093)
- Independent homogeneous RK4 reference: (0.0056687051274573585)
- Absolute difference: (1.25724\times10^{-10})
- Predeclared tolerance: (10^{-9})
- Result: PASS

## Scope and caveat

The analytic 30%-matter / 70%-constant-Λ test passed in the homogeneous COSMOS lane, but the production scalar potential is exponential and rolls. The production run here tests the engine's own archive cosmology, not the identical constant-Λ model. Do not conflate these results. The production kernel remains marked `implemented_not_campaign_validated`.
