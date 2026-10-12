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

## Independent reference comparison: PASS

- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38202689236
- Production (H_{eff}) after one step: (0.005668705253181093)
- Independent homogeneous RK4 reference: (0.0056687051274573585)
- Absolute difference: (1.25724\times10^{-10})
- Predeclared tolerance: (10^{-9})
- Result: PASS

## Production timestep convergence: PASS

- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38202746844
- Commit: `8324107f5b2bf3239ccb46db980046d1445b8d0a`
- Test: actual `V55ProductionKernel.step` at three time steps from fresh initial states, compared against independent RK4 integration of the homogeneous COSMOS equations.
- Grid: resolution 64, (r_{max}=32), local S/D amplitudes zero, radiation off; native rolling exponential potential active.

| dt | Production H_eff | Independent RK4 H_eff | Absolute error |
|---:|---:|---:|---:|
| 1.0e-4 | 0.005668705253181093 | 0.005668705127457359 | 1.25724e-10 |
| 5.0e-5 | 0.005668706292631518 | 0.005668706262189917 | 3.04416e-11 |
| 2.5e-5 | 0.005668706812356937 | 0.005668706804870205 | 7.48673e-12 |

When dt is halved, reference error falls by factors **4.130** and **4.066**, consistent with second-order convergence. All three errors are below the predeclared (10^{-9}) tolerance; all fields remain finite; initial Friedmann normalization agrees to (4.34\times10^{-18}).

Constraint maxima also decrease with dt:

| dt | Hamiltonian max | Momentum max | Connection max |
|---:|---:|---:|---:|
| 1.0e-4 | 7.07694e-13 | 5.01072e-13 | 1.88744e-20 |
| 5.0e-5 | 3.54433e-13 | 2.52230e-13 | 9.48790e-21 |
| 2.5e-5 | 1.77301e-13 | 1.26923e-13 | 4.71860e-21 |

The workflow reports `PASS: production kernel matches the homogeneous reference with second-order convergence`.

## Scope and caveat

The analytic 30%-matter / 70%-constant-Λ test passed in the homogeneous COSMOS lane, but the production scalar potential is exponential and rolls. The production run here tests the engine's own archive cosmology, not the identical constant-Λ model. These results validate the homogeneous operating point and a one-step convergence check only; they do **not** establish long-time stability, inhomogeneous geometry correctness, or coupled S/D physics. The production kernel remains marked `implemented_not_campaign_validated`.
