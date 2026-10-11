# Reference-boundary interpolation audit and confounds — 2026-10-11

## Scope

Offline analysis of the admitted R=160 CMC frozen-state artifact. No engine trajectory was run and no production code was changed for this analysis.

Input:
- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38117001920
- Artifact: `cmc-lapse-frozen-state-crop-audit`, ID `11694568292`
- Source file: `frozen_state_crop_profiles.jsonl`, 49 slices from t=0 to t=12.

## Correct shared-boundary index

The R=80 outer cell center is r=79.75. In the R=160 grid with dr=0.5, this is index 159 (zero-based). The proposed boundary time series must use the R=160 full-domain solved lapse at that index, not the R=160 outermost cell at r=159.75, whose native normalization is alpha=1.

## Held-out temporal interpolation check

For an empirical check, every other interior saved slice was held out and linearly predicted from its adjacent retained slices. This tests interpolation over a 0.5 coordinate-time span using 24 held-out values. The intended trajectory's saved reference slices are spaced about 0.25 apart, and its internal stage spacing is dt=0.00375.

| Metric | Held-out result |
|---|---:|
| Boundary alpha range at r=79.75 | 0.9860461 to 1.0002207 |
| Held-out points | 24 |
| Maximum absolute interpolation error | 1.48936e-4 |
| RMS absolute interpolation error | 5.62025e-5 |
| Median absolute interpolation error | 3.54937e-5 |
| 95th percentile absolute error | 1.03596e-4 |
| Largest error | t=5.74875, prediction low by 1.48936e-4 |

This is a held-out interpolation diagnostic, **not a rigorous error bound**. It shows that linear interpolation is not exact and has errors of order 1e-4 on this coarser-than-target holdout. It does not establish a reliable stage-level error bound at dt=0.00375; the 49 saved samples cannot resolve arbitrary between-slice temporal structure. The proposed counterfactual must report this limitation rather than claim exact boundary matching.

## Code-path confound

In `engine/cmc_gauge.py`, `solve_cmc_lapse` constructs the CMC banded operator and sets the last row to the Dirichlet normalization `alpha(R)=1`; `_cmc_operator_banded` installs that final row. A diagnostic override can change the final RHS entry, but it must be passed at every solve performed during all stages, not only at initialization.

In `engine/production_kernel.py`, the solved `geom.alpha[-1]` is also read by the outer radiation characteristic-speed calculation (`s_minus_rad`, `s_plus_rad`) and the radiation outer-face flux metric. The same boundary metric array is also used in the outer dust flux/source diagnostics. Therefore changing the CMC boundary lapse can change more than the elliptic normalization: it changes boundary metric inputs used in conservative matter closure/flux evaluation. This is a real confound for interpreting a boundary-matched trajectory as a pure gauge-normalization experiment.

## Decision before any trajectory

Do not launch the simple boundary-matched trajectory as if it isolated the CMC normalization. The next implementation review must choose and document one of these diagnostic questions:

1. **Joint boundary intervention:** change CMC boundary normalization at every solve and explicitly accept/report that boundary matter flux/closure inputs change too. This can test whether the combined intervention suppresses divergence, but cannot isolate normalization from the boundary closure.
2. **Separated intervention:** add an opt-in diagnostic-only path that changes the CMC outer Dirichlet value while explicitly retaining the baseline R=80 boundary metric used by the relevant matter boundary flux/characteristic closure. This is a deliberately split counterfactual, not production physics; it must be narrowly scoped, have assertions proving the override touches only the declared boundary calculation, and be compared with the joint intervention if attribution is needed.
3. **Higher-resolution reference boundary history:** if stage-level attribution is necessary, produce a reference history with sufficient temporal output or validated dense-output/replay from the R=160 trajectory. Do not assume interpolation of 49 saved slices provides a bound at every RK/PIRK stage.

Before running, verify that the R=80 and R=160 initialization algorithms are identical except for domain size, and preserve the current production trajectory and source equations unchanged. The one-step diagnostic script is now present in the audit branch, but the new one-step CI workflow still needs a visible run and result before that diagnostic is independently admitted.

## Status

- R=80 full-state trajectory: completed successfully.
- Frozen-state crop audit: completed successfully.
- One-step diagnostic source: committed in audit branch; its independent CI reproduction has not yet been observed.
- Boundary-matched trajectory: not launched.
- No physics, production gauge, or matter equations changed.
