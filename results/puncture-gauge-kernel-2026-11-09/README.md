# Research puncture-gauge kernel: validation (2026-11-09)

`research/puncture_gauge_kernel_v1.py` (`PunctureGaugeKernel`) evolves the lapse (advective 1+log) and shift
(β_t = B, damped Gamma-driver B_t = ¾Λ_t − ηB) instead of solving CMC, stage for stage as the vendor's
`damped_moving_puncture_pirk2_step`, with the production kernel's own scalar/fluid RHS and stress-energy
sources. Research only; no production file changed. Vendor `d6052d60`; Python 3.12.15, numpy 2.5.4, scipy 1.18.1.

## 1. Vacuum: bit-identical to the vendor step

Schwarzschild puncture (M = 1, ψ = 1 + M/2r, α = 1), r_max = 30, η = 2, CFL 0.1, fluids at 1e-30, V0 = 0:
after 100 steps at N = 150 (KO ε = 0.1) every geometry field (a, b, X, α, β, Aa, K, Λ, B) equals the vendor
step's to 0.0; minimum lapse 0.72, max |β| 0.30. Locked by `tests/test_research_puncture_gauge_kernel.py`
(N = 60, 20 steps, ε = 0 and 0.1). Tool: `tools/puncture_gauge_vacuum_check.py`.

## 2. With matter: the Misner–Sharp balance at r = 10 converges

`tools/puncture_gauge_ledger_run.py`: self-gravitating massive S pulse (A = 0.05, Gaussian width 3, time-symmetric,
areal reconstruction), V0 = 0, r_max = 60, η = 2, CFL 0.1, to t = 14; diagnostics every 0.1. Balance residual
max |dM/dt − (flux_T + work_pR + flux_shift)| over 1 ≤ t ≤ 13:

| N | residual | max flux_T | max |flux_shift| | max |β(10)| |
|---|---|---|---|---|
| 100 | 3.54e-3 | 2.61e-3 | 1.1e-6 | 1.41e-3 |
| 200 | 2.48e-4 | 2.00e-3 | 9.4e-7 | 1.40e-3 |
| 400 | 2.42e-5 | 2.08e-3 | 9.8e-7 | 1.42e-3 |

Factors 14 and 10 per doubling (order ≈ 3.4–3.8); 1% of the flux at N = 400, part of which is the 0.1 sampling of
dM/dt. KO ε = 0.1 at N = 200 changes nothing material (2.41e-4, `pl200ko.json`). The N = 100 mismatch was
under-resolution (dr = 0.6 for a width-3 pulse), not the ledger.

## Not shown here

The shift term itself end to end: β at the worldtube is ~1.4e-3, so flux_shift (~1e-6) is below the residual.
Its algebra is tested exactly in `tests/test_misner_sharp_shift_flux.py` (#454). A run with a larger shift at the
worldtube would test it in evolution. The radial Maxwell sector is not supported by this kernel.
