# Strong-D half-CFL refinement — results and domain-test handoff

**Date:** 2026-10-11 UTC  
**Branch:** `diag/independent-metric-ricci-audit`  
**Classification:** diagnostic only. No production physics, evolution equations, projection, gauge, boundary condition, or defaults changed.

## 1. Completed half-CFL run

- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38106791775
- Commit used by run: `ba118fff7e8a3a1353506019e928fb9d58d79186`
- Artifact: `strong-D-regular-F-half-CFL-refinement`, artifact ID `11691305695`
- Artifact digest reported by GitHub: `sha256:909914aa4af0f1b69d8612e06f7e104728f3d9d1b99f02d4b0fc20c6270ca5b7`
- Job conclusion: success; both trajectories completed to (t=12), 6,400 accepted steps each, no evolution failure or diagnostic-gate failure.
- Configuration: (N=160, r_{max}=80, Delta r=0.5, mathrm{CFL}=0.00375, D=10^{-4}), radiation on. The timestep was halved from the original-CFL baseline.
- Modes: untouched production initializer and regular-F cubic/Gauss initial-(B) quadrature shadow. Both use the unchanged production evolution kernel.

The artifact contains `report.json` and the per-mode progress/checkpoint JSON files. The report's constraint reconstruction and curvature-replacement closures pass their stated tolerances. Neither trajectory reports a cycle event or handoff.

## 2. Matched comparison with the original-CFL baseline

Original-CFL baseline: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38100254183  
Baseline artifact ID: `11688328331`, digest `sha256:f2416cc062c7cf8348103582a3a4eb75fceae3d18d4130cd6d82a94cfff77966`.

The tables use the off-centre maximum `max_abs_H_cells_2plus`, excluding cells 0 and 1, so the large production central residual does not mask the migrating off-centre feature. Nominal target labels are matched; actual sampled coordinate times differ by up to roughly one timestep between the runs.

### Production initializer

| Target (t) | Original CFL (max_{i\ge2}|H_i|) | Half CFL (max_{i\ge2}|H_i|) | Relative change |
|---:|---:|---:|---:|
| 4 | (1.2368329\times10^{-4}) | (1.2376563\times10^{-4}) | +0.0666% |
| 8 | (1.0120214\times10^{-4}) | (1.0123311\times10^{-4}) | +0.0306% |
| 12 | (2.8291823\times10^{-4}) | (2.8292911\times10^{-4}) | +0.0038% |

### Regular-F cubic/Gauss shadow

| Target (t) | Original CFL (max_{i\ge2}|H_i|) | Half CFL (max_{i\ge2}|H_i|) | Relative change |
|---:|---:|---:|---:|
| 4 | (4.8082504\times10^{-6}) | (4.7988102\times10^{-6}) | -0.1963% |
| 8 | (5.9441992\times10^{-5}) | (5.9544987\times10^{-5}) | +0.1732% |
| 12 | (2.8281186\times10^{-4}) | (2.8282273\times10^{-4}) | +0.0038% |

For the shadow trajectory the full-domain maximum is off-centre by (t=8) (reported radius 6.25) and (t=12) (radius 9.75) in both CFL runs. The quarter-time full-profile domain test remains necessary to determine whether that migration is continuous or an argmax switch between local peaks.

At (t=12), the production central (H) is (1.0934304\times10^{-3}), while the shadow central (H) is (-7.4770\times10^{-6}). Nevertheless, their off-centre maxima are nearly equal: (2.8292911\times10^{-4}) and (2.8282273\times10^{-4}), respectively. Thus removing the large initial central residual does not remove the later off-centre growth.

### Lapse comparison

The first-cell lapse is also closely matched between original and half CFL:
- Production, target (t=12): (0.0165558762) versus (0.0165580700).
- Shadow, target (t=12): (0.0165561401) versus (0.0165583339).

At target (t=8), the relative lapse difference is about 0.145%; at (t=4), about 0.0105%. This does not support a material timestep effect over the sampled interval. The late lapse decline remains severe and is a reason to retain lapse in the domain-size test.

## 3. Conclusion allowed by the evidence

**The off-centre residual growth is not materially reduced by halving the timestep at the sampled targets.** Differences from the original-CFL baseline remain below 0.2% for both trajectories at (t\approx4,8,12); at (t\approx12), the mismatch is about 0.004%. The initial production-vs-shadow residual difference is therefore not a sufficient explanation of the later off-centre residual.

This is not yet proof of a continuum, physical characteristic or a specific evolution-equation defect. The present half-CFL report samples only (t=0,0.01,4,8,12) and does not store full radial (H)/lapse profiles. Spatial and domain effects remain unresolved.

## 4. Domain-size experiment launched

- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38108719482
- Trigger commit: `9087dbea4c8507894088cc2603c6608b6cfc8623`
- The domain audit script and fixed-spacing/quarter-time configuration validation compile and pass.
- At time of this note: dependency installation is in progress; numerical trajectories have not yet started and there is no result artifact yet.
- Four cases run sequentially: production at ((N,r_{max})=(160,80)), production at ((320,160)), then the corresponding regular-F cubic/Gauss shadows.
- Both domains hold (Delta r=0.5), CFL (=0.0075), (D=10^{-4}), scalar profile, radiation, and production evolution fixed; target interval is (0ldots12).
- Full radial Hamiltonian-residual and lapse profiles are written at (t=0) and every 0.25 time units, with first-cell lapse, off-centre peak radius/magnitude, minimum lapse, geometry health and accepted-step checkpoints.

The launch marker was changed in the audit branch only after the half-CFL assessment. No other in-progress Actions runs were present immediately before launch, so no expensive campaigns overlapped.

## 5. Next analysis, once the artifact exists

1. Compare the full (H(r,t)) profiles at both domain sizes, not only (arg\max |H|). Establish whether the off-centre feature moves continuously or the maximum switches between peaks.
2. Track (\alpha(0,t)), (\alpha_{min}(t)), peak radius and numerical health together. Because the lapse is from a domain-wide elliptic solve, domain changes may alter the central lapse itself.
3. Treat domain-size dependence as evidence about boundary/lapse coupling, not proof by itself.
4. Only after the domain comparison, stage a separate (Delta r) refinement to test whether the feature's radius-time relation survives a smaller grid spacing.
5. Do not patch physics, assert a bounce/turnaround, or admit production based on this diagnostic alone.
