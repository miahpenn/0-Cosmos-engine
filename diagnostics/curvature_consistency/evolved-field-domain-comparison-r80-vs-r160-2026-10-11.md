# Evolved-field domain comparison: R=80 vs R=160 — 2026-10-11

## Admission and inputs

- Diagnostic-only comparison; no engine evolution or production physics was changed by this comparison.
- R=80 capture workflow: [run 38150925743](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38150925743), commit `f7be6306a0f7845f06d02dbb098ae869eec9759c`, successful.
- R=80 artifact: `strong-D-evolved-state-r80-capture`, artifact ID `11706797173`, SHA-256 `8b118d1f81e01f2a5d291d160b8bac42abc6a11051d8453efd9c19d0826d6d8b`.
- R=80 settings: N=160, R=80, dr=0.5, CFL=0.0075, D amplitude=1e-4, scalar amplitude=0.01, width=7, radiation enabled, final t=12, 3,200 accepted steps, 49 saved slices. No cycle events, handoffs, failures, gauge changes, boundary-condition changes, projection changes, or production-physics changes.
- R=160 reference: the 49-slice full-state snapshot in the CMC frozen-state crop artifact from [the crop-audit run 38117001920](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38117001920). The reference state is restricted to the identical R=80 cell centers. Both time arrays match exactly; both radial grids match exactly over the shared domain.
- This is an independent-evolution comparison, not the frozen-state crop operator test. The comparison uses the archived evolved geometry, scalar fields, and conservative matter fields at all saved times.

## Principal result

The evolved fields do diverge substantially despite an almost identical initial physical-state snapshot. The initial alpha field differs by at most 5.66858e-4 (domain-dependent lapse normalization); the initial D, PD, a, b, X, scalar, and matter arrays match at the stored precision, apart from roundoff-scale initial differences in Lambda/Aa/K. Differences then grow during the independent evolutions.

At t=12:
- Maximum lapse difference is 0.116139 at r=39.75. At that location alpha_R80=0.9538121810 and alpha_R160=0.8376731134.
- Maximum difference in a is 0.00700752 at r=17.75; b differs by 0.00371757 at r=17.25.
- D differs by 0.00421259 at r=13.25; PD differs by 0.00399110 at r=13.25.
- X differs by 0.00172424 at r=35.75; Lambda by 0.00226774 at r=13.75; Aa by 0.00190403 at r=18.25.
- Conservative matter differences are smaller in the core but notable near the shared outer edge: dark-matter radial momentum differs by 0.00214747 at r=79.75, while its core maximum difference (r<=20) is 1.83969e-4.
- The lapse's first growth-threshold crossing is at t=5.74875, r=21.75. The growth threshold is only a reporting aid: max(|difference|) > max(10 times the initial max difference, 1e-6 times the field scale at that slice). It is not a causal or physics-admission criterion.
- D and PD pass that reporting threshold at the first post-initial saved slice (t=0.25125, r=0.25). This says the initially tiny lapse normalization difference is followed by an early scalar-field history difference; it does not by itself identify the causal pathway.

## Central proper-time accumulator

The saved `tau` accumulator differs between the two independent trajectories:
- t=4.00125: delta tau (R80 - R160) = +0.00101611.
- t=7.99875: delta tau = +0.01563555.
- t=12: tau_R80=6.91852193294; tau_R160=6.90677928325; delta=+0.01174264969.
- Maximum absolute saved-slice delta tau is 0.01570218.

These are differences in the engine's saved `tau` accumulator. Preserve that wording until its precise convention is confirmed; do not silently relabel it as a formal central proper-time measurement.

## Field-level results

Max absolute difference over the shared grid at t=12. “Core max” means r<=20. Growth times use the threshold defined above. A blank growth time means the threshold was not crossed.

| Field | Initial max difference | t=12 max difference | Radius at t=12 max | t=12 core max | First growth crossing (t, r) |
|---|---:|---:|---:|---:|---:|
| geometry.alpha | 5.66858e-4 | 1.16139e-1 | 39.75 | 5.32211e-2 | 5.74875, 21.75 |
| geometry.a | 0 | 7.00752e-3 | 17.75 | 7.00752e-3 | 4.75125, 6.75 |
| scalars.D | 0 | 4.21259e-3 | 13.25 | 4.21259e-3 | 0.25125, 0.25 |
| scalars.PD | 0 | 3.99110e-3 | 13.25 | 3.99110e-3 | 0.25125, 0.25 |
| geometry.b | 0 | 3.71757e-3 | 17.25 | 3.71757e-3 | 4.99875, 6.75 |
| geometry.Lambda | 2.50354e-13 | 2.26774e-3 | 13.75 | 2.26774e-3 | 0.25125, 3.25 |
| matter.dark_matter.momentum_r | 0 | 2.14747e-3 | 79.75 | 1.83969e-4 | 0.25125, 79.75 |
| geometry.Aa | 1.61329e-16 | 1.90403e-3 | 18.25 | 1.90403e-3 | 0.25125, 6.25 |
| geometry.X | 0 | 1.72424e-3 | 35.75 | 1.21210e-3 | 0.49875, 0.25 |
| scalars.phi | 0 | 1.38020e-3 | 36.75 | 9.68403e-4 | 0.25125, 0.25 |
| geometry.K | 2.42861e-16 | 4.90541e-4 | 79.75 | 4.90527e-4 | 0.25125, 79.75 |
| matter.baryons.momentum_r | 0 | 3.34655e-4 | 79.75 | 2.86706e-5 | 0.25125, 79.75 |
| matter.dark_matter.rest | 0 | 1.56237e-4 | 79.25 | 3.58193e-5 | 3.49875, 78.75 |
| matter.dark_matter.energy_t | 0 | 1.52799e-4 | 79.25 | 6.49114e-5 | 3.49875, 78.75 |
| scalars.S | 0 | 1.13453e-4 | 5.75 | 1.13453e-4 | 0.25125, 0.25 |
| scalars.PS | 0 | 6.73951e-5 | 0.25 | 6.73951e-5 | 0.25125, 0.25 |
| scalars.Pi | 0 | 3.19319e-5 | 12.75 | 3.19319e-5 | 1.50000, 0.25 |
| matter.baryons.rest | 0 | 2.43490e-5 | 79.25 | 5.58240e-6 | 3.49875, 78.75 |
| matter.baryons.energy_t | 0 | 2.40948e-5 | 79.25 | 1.00956e-5 | 3.49875, 78.75 |
| matter.radiation.momentum_r | 0 | 2.27264e-6 | 79.75 | 7.85687e-8 | 0.25125, 79.75 |
| matter.radiation.energy_t | 0 | 2.13791e-6 | 79.25 | 2.65996e-8 | 0.25125, 79.25 |
| geometry.beta | 0 | 0 | 0.25 | 0 | — |
| geometry.B | 0 | 0 | 0.25 | 0 | — |
| matter.radiation.rest | 0 | 0 | 0.25 | 0 | — |

## Interpretation limits and next step

1. The crop audit establishes that the frozen-state CMC crop operator reproduces the full-domain solve. It does not explain this independent-evolution divergence.
2. The new pair establishes a domain-associated evolved-state mismatch, including substantial lapse divergence and growing scalar/geometry differences. It does **not** yet establish whether the lapse normalization, elliptic boundary influence, or subsequent coupled evolution is the dominant cause.
3. Do not infer a physical front, propagation speed, bounce, turnaround, or causal ordering from these threshold locations. Any speed statement requires checking characteristic speeds for the actual evolved metric and matter system.
4. Next, audit the exact initializer and the lapse/CMC solve at t=0, then trace the first saved-slice changes by field using the same state arrays. Keep this diagnostic-only: no gauge/physics patch and no new expensive trajectory until the existing evidence is reconciled.
