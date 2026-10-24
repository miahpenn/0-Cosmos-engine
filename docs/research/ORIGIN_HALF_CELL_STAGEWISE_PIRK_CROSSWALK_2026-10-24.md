# Saved-stage origin Ricci/source join — 2026-10-24

**Disposition:** offline research diagnostic only. No evolution was rerun, and no production equation, initializer, coefficient, gauge, boundary, vendor pin, threshold, or scientific-admission flag was changed. The baseline and archived origin-half-cell candidate were joined across every saved accepted state and captured PIRK stage.

## Provenance

- Production source anchor: `4aec88d0e3ff028e0a65a8ac36b30948625df792`.
- Vendor pin: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Baseline artifacts: workflow [#40127716287](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40127716287), N=120 artifact [#12573419382](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40127716287/artifacts/12573419382), N=240 artifact [#12572764318](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40127716287/artifacts/12572764318).
- Candidate artifacts: workflow [#40128647367](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40128647367), N=120 artifact [#12574010545](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40128647367/artifacts/12574010545), N=240 artifact [#12574710341](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40128647367/artifacts/12574710341).
- The tool checks the three internal artifact checksums, source/vendor anchors, matching radial grids, matching stage labels, and exact matched start times. It recomputes proper-radial trace-free Ricci from each saved `(a,b,X)) profile; it does not change or advance either state.

## First captured turn in the source and gap histories

Ratios below are candidate/baseline. The fixed-support source RMS uses every cell centre with r≤10, while the first-eight-cell statistic is retained for continuity with the earlier crosswalk. Times label each saved stage's accepted-step start time.

| Resolution | First explicit stage with fixed-support (r≤10) source RMS ratio ≥ 1 | First explicit stage with first-8-cell source RMS ratio ≥ 1 | First explicit stage with cell-0 source absolute ratio ≥ 1 | First explicit stage with Ricci-gap max / RMS ratios ≥ 0.5 |
|---:|---|---|---|---|
| N=120 | step011_explicit_predictor, t=0.22, ratio 1.033 | step011_explicit_predictor, t=0.22, ratio 1.089 | step012_explicit_predictor, t=0.24, ratio 1.012 | step012_explicit_predictor, t=0.24, max 0.539 / RMS 0.509 |
| N=240 | step010_explicit_predictor, t=0.10, ratio 1.005 | step010_explicit_predictor, t=0.10, ratio 1.044 | step012_explicit_predictor, t=0.12, ratio 1.042 | step012_explicit_predictor, t=0.12, max 0.556 / RMS 0.526 |

The fixed physical support makes the resolution comparison cleaner than eight cells alone: the first source-RMS crossing occurs at t=0.22 for N=120 (ratio 1.033) and t=0.10 for N=240 (ratio 1.005), at the same respective stages as the first-eight-cell statistic. At those stages the Ricci gap is still reduced versus baseline: first-eight-cell max/RMS ratios are 0.279/0.240 and 0.103/0.098, respectively. One step later for N=120 and two steps later for N=240, the gap has rebounded to roughly half of baseline. Thus the multi-cell source channel signals the turn earlier than a centre-only source test.

## Non-monotone behaviour and resolution sensitivity

| N | Saved stage/time | First-8-cell Ricci-gap max ratio | Ricci-gap RMS ratio | Source RMS ratio, first 8 cells | Source RMS ratio, fixed r≤10 | Cell-0 source absolute ratio |
|---:|---|---:|---:|---:|---:|---:|
| 120 | initial accepted state, t=0.00 | 0.112 | 0.124 | 0.540 | 0.677 | 0.079 |
| 120 | explicit predictor, t=0.16 | 0.021 | 0.029 | 0.763 | 0.873 | 0.025 |
| 120 | explicit predictor, t=0.20 | 0.103 | 0.098 | 0.958 | 0.982 | 0.158 |
| 120 | explicit predictor, t=0.22 | 0.279 | 0.240 | 1.089 | 1.033 | 0.339 |
| 120 | explicit predictor, t=0.24 | 0.539 | 0.509 | 1.201 | 1.067 | 1.012 |
| 120 | accepted state, t=0.60 | 0.195 | 0.295 | 0.312 | 0.337 | 0.082 |
| 240 | initial accepted state, t=0.00 | 0.121 | 0.133 | 0.463 | 0.776 | 0.086 |
| 240 | explicit predictor, t=0.08 | 0.027 | 0.034 | 0.721 | 0.941 | 0.028 |
| 240 | explicit predictor, t=0.10 | 0.103 | 0.098 | 1.044 | 1.005 | 0.174 |
| 240 | explicit predictor, t=0.11 | 0.288 | 0.247 | 1.376 | 1.030 | 0.369 |
| 240 | explicit predictor, t=0.12 | 0.556 | 0.526 | 1.895 | 1.046 | 1.042 |
| 240 | accepted state, t=0.60 | 0.168 | 0.169 | 0.506 | 0.552 | 0.433 |

Both resolutions show a deep initial reduction followed by a rebound and later decline. But the rebound is not aligned in physical time: the first fixed-support source-RMS crossing occurs around t=0.10 at N=240 and t=0.22 at N=120. At t=0.60 the Ricci gap is smaller than baseline in the first eight cells, but source and cell-0 ratios remain resolution-dependent; the proper-radial representation itself is not thereby proved regular.

The per-stage JSON report also preserves fixed-support (r≤10) source and Ricci-gap max/RMS ratios, signed first-four-cell Ricci gaps and source values, the location of each first-eight-cell maximum, the saved Ricci and lapse-Hessian source terms at cell 0, and the matching accepted/PIRK-stage identity. The fixed-radius ratios confirm the first source-RMS crossing at t=0.22 (N=120) and t=0.10 (N=240); this is not an artifact of comparing different physical radii via the first-eight-cell window. This guards against overinterpreting a cell-0 norm or a single accepted time.

## Interpretation

**Supported:** the archived half-cell initialization materially changes early near-origin geometry and the (A_a) source; the effect survives multiple PIRK stages. The source begins to lose its multi-cell advantage before the cell-0 ratio alone reveals it.

**Not supported:** a production fix; continuum convergence; a diagnosis that uniquely blames the vendor Ricci formula; or promotion of scientific admission. The Ricci crosswalk uses the same pinned D1/D2 spatial operators as the vendor representation, so it is a different geometric representation—not an independent discretization.

## Next gate

The captured history is now exhausted at the existing recorded stage granularity. The next useful comparison is not another same-output rerun: it is a source-pinned, explicitly matched physical-time and stage-state convergence analysis that accounts for the CFL-dependent timestep, and checks whether the near-origin geometric/source profiles approach a common limit. Keep any candidate initializer counterfactual offline until that comparison distinguishes quadrature truncation from evolution/stencil response. Scientific admission remains false.

## Reproducibility

- Offline stage join: [`tools/audit_origin_half_cell_stage_crosswalk.py`](../../tools/audit_origin_half_cell_stage_crosswalk.py)
- Focused regression tests: [`tests/test_origin_half_cell_stage_crosswalk.py`](../../tests/test_origin_half_cell_stage_crosswalk.py)
- CI workflow downloads the existing artifacts and recomputes the join; it does not launch the production integrator.
