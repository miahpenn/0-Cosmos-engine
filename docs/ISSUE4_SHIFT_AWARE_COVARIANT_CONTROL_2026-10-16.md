# Issue #4 — shift-aware covariant-divergence manufactured control

**Date:** 16 October 2026 (UTC)  
**Status:** diagnostic-only formula control; not production-source evidence and not admission.

## Motivation

The frozen-state total-stress-energy artifact from run [38939174153](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38939174153) already contains baseline, conditional rest-current-only, and joint conditional-source arms, but its configuration sets `nonzero_shift_at_frozen_slice=false`. Before extending that production-source-pinned observer, the coordinate mixed-tensor identity must pass an analytic nonzero-shift control.

## Identity and implementation

For the spherical ADM metric
`ds²=(-α²+γ_rr β²)dt² + 2γ_rr β dt dr + γ_rr dr² + γ_θθ dΩ²`,
the independent coordinate components used by the test are

- `T^{tt}=E/α²`
- `T^{tr}=−βE/α² + J_r/(αγ_rr)`
- `T^{rr}=Eβ²/α² − 2βJ_r/(αγ_rr) + p_r/γ_rr`
- `T^{θθ}=p_t/γ_θθ`

with mixed components `T^t_t=−E+βJ_r/α`, `T^t_r=J_r/α`, `T^r_t=β(E+p_r)−αJ_r/γ_rr−β²J_r/α`, and `T^r_r=p_r−βJ_r/α`. Both metric-source contractions retain the `2 T^{tr} ∂g_tr` cross term.

## Manufactured control

The test represents flat Minkowski space in the regular time-dependent radial chart `R=r+t v(r)`, `T=t`, where `v(r)=0.035 r exp[−(r/4.3)²]`. On the measured slice, `α=γ_rr=1`, `γ_θθ=r²`, and the shift is nonzero but odd/regular at the origin. A homogeneous isotropic stress with Eulerian `E=2.3`, `p=0.7`, and `J_r=0` has exactly zero covariant divergence in this flat chart.

The local analytic harness produced annulus RMS residuals of order (10^{-17}) for both `∇_μT^μ{}_t` and `∇_μT^μ{}_r`, including maximum absolute residuals at or below (1.2×10^{-16}) on the tested grids. Deliberately reversing the time-component flux sign produces annulus RMS near (1.44×10^{-1}); omitting the spatial `g_tr` metric-source term produces RMS near (9.69×10^{-4}); omitting the time `g_tr` term produces RMS near (4.93×10^{-5}). The zero-shift limit reproduces the prior Issue #4 component contract.

## Limits and next gate

These are analytic derivative inputs, so the control tests tensor-index algebra and shift cross-term signs—not a production radial stencil, a semidiscrete RHS, the physical exchange law, or a full evolved-state conservation law. No engine/vendor equation, coupling, state, or production trajectory changed.

Next: integrate this helper into a source-pinned copy of the existing observer, preserve the zero-shift artifact as the comparator, and run the baseline/rest-current/joint-source arms on a nonzero-shift state. Report signed first-two-cell and fixed-annulus budgets. Keep scientific admission false until the conditional premise, full RHS/source accounting, stencil/PIRK closure, and convergence gates independently pass.
