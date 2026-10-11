# Boundary-matched CMC R=80 experiment — results

- Status: **completed**
- Commit under test: `35072830f9e4d7029bb10243a0520b9a273a17dd`
- Steps: 3200; saved states: 49; final t: 12
- CMC-alpha reference SHA-256: `d5e80a52ff03848b72545b6ca695c4ceb1b4e4e69a50a43508555c849001a710`
- Matter boundary inputs: archived R=80 last-cell a,b,X,alpha,beta, interpolated in time and frozen for matter RHS evaluation.
- CMC Dirichlet: R=160 alpha(r=79.75), interpolated in time.
- Interpolation caveat: source histories are sampled about every 0.25 time units, much coarser than dt=0.00375; stage-level interpolation error is not rigorously bounded.

## t=12 mismatch versus archived R=160 state

| Field | R80 baseline max error | Matched R80 max error | Max-error reduction | Baseline core max | Matched core max |
|---|---:|---:|---:|---:|---:|
| geometry.alpha | 0.116139 | 0.106702 | 8.13% | 0.0532211 | 0.0521883 |
| geometry.a | 0.00700752 | 0.00604251 | 13.77% | 0.00700752 | 0.00604251 |
| geometry.b | 0.00371757 | 0.00320811 | 13.70% | 0.00371757 | 0.00320811 |
| geometry.X | 0.00172424 | 0.0014107 | 18.18% | 0.0012121 | 0.00102151 |
| geometry.Aa | 0.00190403 | 0.00165661 | 12.99% | 0.00190403 | 0.00165661 |
| geometry.K | 0.000490541 | 0.000478226 | 2.51% | 0.000490527 | 0.000478207 |
| geometry.Lambda | 0.00226774 | 0.00196329 | 13.43% | 0.00226774 | 0.00196329 |
| scalars.D | 0.00421259 | 0.00365105 | 13.33% | 0.00421259 | 0.00365105 |
| scalars.PD | 0.0039911 | 0.00345949 | 13.32% | 0.0039911 | 0.00345949 |
| scalars.phi | 0.0013802 | 0.00115549 | 16.28% | 0.000968403 | 0.000833993 |
| scalars.PS | 6.73951e-05 | 6.10017e-05 | 9.49% | 6.73951e-05 | 6.10017e-05 |
| matter.dark_matter.momentum_r | 0.00214747 | 0.00213567 | 0.55% | 0.000183969 | 0.000155067 |
| matter.dark_matter.energy_t | 0.000152799 | 0.000173894 | -13.81% | 6.49114e-05 | 5.48216e-05 |
| matter.baryons.momentum_r | 0.000334655 | 0.000332815 | 0.55% | 2.86706e-05 | 2.41666e-05 |
| matter.radiation.momentum_r | 2.27264e-06 | 1.4408e-06 | 36.60% | 7.85687e-08 | 6.62526e-08 |

## Proper-time accumulator (retain the engine's exact label)

- Baseline R80 minus R160 at final saved slice: 0.0117426496868
- Matched R80 minus R160 at final saved slice: 0.0106210937742

## Interpretation limits

This is a controlled numerical counterfactual, not proof that a gauge choice is physically irrelevant. Positive error reduction means the boundary-matched trajectory is numerically closer to the archived R=160 trajectory for that field and metric; it does not by itself prove causation. Negative reduction means the mismatch increased. The outer-boundary histories are interpolated from quarter-time snapshots, so small changes must be interpreted against that limitation. No physical front, propagation speed, bounce, or turnaround is inferred here.

