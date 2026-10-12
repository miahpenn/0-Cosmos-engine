# V5.5 production constraint causal audit — 2026-10-12

**Disposition:** production horizon completion is verified; numerical-constraint and balance admission remains open. This report consolidates the signed source budget, saved-state variable attribution, and independent Ricci-route comparison. All calculations below are offline post-processing of the archived campaign states unless a different run is explicitly identified. No production trajectory was added by these diagnostics.

## 1. Pinned campaign provenance

- Production campaign: [run 38208253687](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208253687)
- Production kernel/source commit: `6d8b3f2790644591326df57c55908863f994cb0c`
- Campaign artifact: `v55-production-t12-n40-n60-n80`, SHA-256 `dbbaa3be94242cdeecd437d0a41c00d4cf3f2d625e1f3aaaad9e1b45833f2047`
- Configuration: N=(40,60,80), r_max=80, CFL=0.03, radiation enabled, requested t_final=12, saved checkpoints around t=(2,4,6,8,10,12), plus explicit t=0 initialization.
- Every resolution completed to t=12. This is **execution/horizon admission only**. The scalar H/M constraints and balance criteria have not passed scientific admission.

The initial slice is already imperfect in the first central cell. On the production operator, the full-domain H residual at t=0 is approximately +3.769e-3, +3.932e-3, +3.989e-3 for N=40,60,80. This central stencil is outside the vendor operator's admitted center region. The defect does not disappear at the first admitted cell: by t=12 the H residual at i=2 is -4.7781e-3 at r=5.0 for N=40, -1.07212e-2 at r=3.3333 for N=60, and -1.96250e-2 at r=2.5 for N=80.

### End-state source budget

Using the production definitions H = H_geometry - 16*pi*rho_total and M = M_geometry - 8*pi*j_total:

| Resolution | Full-domain H max at t=12 | H residual at i=2 (first admitted cell) | H geometry at i=2 | 16*pi*rho_total at i=2 | Full-domain M max at t=12 |
|---:|---:|---:|---:|---:|---:|
| N=40 | 4.04166e-2 at r=1.0 | -4.77807e-3 at r=5.0 | -3.90011e-3 | +8.77958e-4 | 3.32492e-3 |
| N=60 | 7.96984e-2 at r=0.6667 | -1.07212e-2 at r=3.3333 | -9.41176e-3 | +1.30942e-3 | 3.00203e-3 |
| N=80 | 1.34942e-1 at r=0.5 | -1.96250e-2 at r=2.5 | -1.81070e-2 | +1.51798e-3 | 8.08143e-3 |

The central-cell decomposition at N=80, t=12 is H_geometry=-1.33108e-1, 16*pi*rho_total=+1.83396e-3, and H=-1.34942e-1. The central M residual is likewise almost entirely geometric (M_geometry=+8.08156e-3; matter source=+1.27e-7). So this is not explained by a matter-source magnitude error alone.

Location matters. For the same fixed exterior region r>=5, the t=12 H max is 4.778e-3 at N=40 but only 5.077e-5 at N=60 and 6.963e-5 at N=80. The central and exterior behaviors differ; do not use a moving “exclude first two cells” norm as a common physical-radius convergence claim.

## 2. Saved-interval Shapley attribution by evolved variable

- Offline diagnostic workflow: [run 38209482354](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209482354)
- Artifact `v55-production-hamiltonian-source-budget`, SHA-256 `4e104c5709d5a377cabc05b39ac830fde24224b176d95294bd42d3db1d106c06`
- The diagnostic compares each saved pair across six intervals for N=40/60/80. It evaluates 256 pre/post hybrid states per interval and attributes the change in H to a, b, X, lapse/shift variables, A_a, K, Lambda, and matter. No production step is taken.
- Max absolute Shapley closure across all 18 intervals and cells: 5.47e-17.

The split revises the coarser “metric group dominates” statement:

1. In the first central cell (which is not an admitted point of the formal constraint stencil), X is the largest consistently negative contributor to the interval change in H at N=60 and N=80; Lambda often supplies a positive offset. At N=80 from t≈8 to t≈10, the cell-0 change is -8.1314e-2, attributed approximately X=-6.839e-2, a=-2.483e-2, b=-2.47e-3, Lambda=+1.426e-2, matter=+1.09e-4.
2. At the first admitted cell i=2, the story differs: the Lambda contribution is often the largest negative increment. At N=80 from scheduled t=10 to t=12, H at r=2.5 changes by -6.820e-3, with contributions approximately Lambda=-6.69e-3, a=+2.30e-3, b=-2.33e-3, X=-3.29e-4, matter=+2.34e-4.
3. The pointwise H formula contains no direct alpha/beta/B term, so the instantaneous lapse/shift hybrid contribution to H is zero. This does **not** mean gauge is uninvolved in time evolution.
4. At the center, matter and K/A_a changes are much smaller than the dominant geometry-variable terms. At i=2, the signed interaction of a, b, X and Lambda varies with time; don't summarize the admitted region as an X-only problem.

Shapley values are endpoint-state counterfactual allocations for a nonlinear function, not proof that an individual evolution equation is defective. They exclude intrastep variation. They are useful because they identify exactly which variables to audit next without changing the evolution.

## 3. Independent metric-derived Ricci comparison on the exact saved states

- Clean, closure-gated workflow: [run 38209752406](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209752406)
- Artifact `v55-production-hamiltonian-source-budget`, SHA-256 `358b40fd3c9337ee9010752ed2386a1496b84fb3b6ccfc0ad70f68f86dae9565`
- 21 states audited: the production-initialized t=0 slice plus six saved checkpoint states for each resolution. The same state, metric, extrinsic-curvature fields, and matter source are used; only the diagnostic 3-Ricci scalar is replaced by a metric-derived expression on that same slice.
- All algebraic/reconstruction gates pass across all 21 snapshots. Maximum absolute closure errors: H reconstruction 0; curvature replacement identity 1.39e-17; Ricci route split 8.67e-19; evolved-Lambda-to-Ricci identity 2.78e-17; connection constraint reconstruction 0.

The independent metric-derived curvature route changes the Hamiltonian residual materially but does not eliminate it:

| Resolution | t=12 central H using production Ricci | Same-state H using metric-derived Ricci | H at first admitted cell, production | H at first admitted cell, metric-derived |
|---:|---:|---:|---:|---:|
| N=40 | -4.04166e-2 | -1.16464e-2 | -4.77807e-3 | -3.24848e-3 |
| N=60 | -7.96984e-2 | -2.39553e-2 | -1.07212e-2 | -7.08700e-3 |
| N=80 | -1.34942e-1 | -3.82023e-2 | -1.96250e-2 | -1.27665e-2 |

At N=80, t=12, the metric-derived route reduces the magnitude of the central H residual by about 71.7%, but the central residual remains -3.82e-2. At the first admitted cell r=2.5, the magnitude falls by about 34.9%, but a residual of -1.277e-2 remains.

The Ricci difference is exactly split into:
- explicit evolved-connection contribution: R_vendor(Lambda_evolved) - R_vendor(Lambda_metric);
- remaining metric-route/discretization difference: R_vendor(Lambda_metric) - R_metric.

At N=80, t=12, cell 0: total curvature gap=-9.674e-2, explicit connection contribution=+1.800e-2 (an offset), remaining route difference=-1.147e-1. At the first admitted cell r=2.5: total gap=-6.859e-3, connection contribution=+9.574e-3, route difference=-1.643e-2. Thus, in these samples, the remaining metric-route difference dominates the net curvature gap; the evolved-Lambda contribution partly offsets it rather than producing it alone.

Critical caveat: the two Ricci routes already differ on the t=0 production initialization. At N=80, the central t=0 H is +3.989e-3 with the production operator and +1.585e-3 with the metric-derived diagnostic. The initializer was not solved against this alternate discrete curvature route. Therefore this counterfactual is a geometric consistency discriminator, **not authorization to replace production Ricci**. The remaining defect and the differing discretizations must be explained before any numerical intervention.

## 4. Parallel research tracks

### CMC target-history sensitivity

- Full R80 evolving target/boundary counterfactual: [run 38207554072](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38207554072).
- Coarsened archived R160 Kdot history closes over 99.98% of the baseline final proper-time gap in the reported controlled setup. The every-other-sample target history gives final matched-minus-R160 tau about +1.985e-6; the full-sample target history had about +6.30e-7. This is causal evidence of target/boundary-history sensitivity, not proof that the archived R160 target is the unique physical or self-consistent R80 target.
- Initial-slice target/lapse fixed-point diagnostic: [run 38209662328](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209662328), with an earlier successful report at [run 38209333687](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209333687). Initial R80 vs R160 targets differ by only -3.47e-7 on the archived initial states, and the fixed-target elliptic solve repeats identically to 1e-12. The undamped alpha-dependent map did not meet its configured convergence criterion within 20 iterations despite tiny late alpha changes; record that as a non-converged map probe, not a failed production solve.

### Strong-D ridge and cancelled N=640 line

- Fixed-domain N=160/N=320 strong-D ridge refinement: [run 38165315602](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38165315602). Late shell-window peak amplitudes are close (N=160 2.829182e-4 at r=9.75; N=320 2.827630e-4 at r=9.875; ratio 0.99945), but the middle-time peak identity switches. Treat this as late-time agreement, not a blanket convergence proof.
- N=640 signed-pole run [38183314088](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38183314088) was cancelled at the time limit before t=12; it has no admitted full-horizon result.

### Stage-by-stage geometry trace

The parallel one-step production constraint stage-trace line is running on `diag/independent-metric-ricci-audit` (latest observed run [38209707396](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38209707396)). Do not duplicate it while in progress. Reconcile its within-step attribution with the saved-endpoint findings above when its artifact is available.

## 5. Decision and next steps

**Current decision:** no change to the production equations, no promotion to `main`, and no t=160 campaign yet. The t=12 engine executes, but constraint/balance admission remains open.

1. Reconcile the one-step stage trace with this saved-state decomposition, paying particular attention to X, Lambda, and the vendor-Ricci versus metric-route gap. The important distinction is between excluded central cells and admitted (r\geq2\,dr) cells.
2. If the stage trace identifies a specific stage/operator, design one bounded numerical-only counterfactual using frozen identical input states, a matching initialization/discretization, and predeclared pass/fail criteria.
3. Retain both full-domain and fixed-physical-radius (r\geq5) witnesses. Do not change thresholds after seeing results.
4. Keep matter source/flux balance and the CMC boundary/target path as parallel audits. They are not evidence that the central constraint issue has been fixed.
5. Do not add artificial damping, floors, clamps, fitted feedback, a bounce law, or a silent physics change. Re-run the full bounded t=12 N=40/60/80 campaign only after the isolated diagnostic has justified a numerical-layer change.

**Guardrail:** all commits for this audit are isolated to `fix/strict-campaign-admission-2026-10-12` / draft PR #3. The default `main` branch has not been modified by this audit.
