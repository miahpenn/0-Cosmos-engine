# Gear 69/70B nonlinear cavity in the production engine (2026-11-09)

The production engine (`V55ProductionKernel`, CMC slicing, β = 0) reproduces the gear's centre-proper-time
growth coefficient of the uncoupled mode D to about 1e-4, from a three-resolution study. Research only:
`research/gear69_cavity_engine_harness_v1.py` changes nothing in production; it monkeypatches the D right-hand
side and stress projection for the cavity.

## Setup

- Pure D cavity, V_D = −D²/2, μ = 1, G = 1 normalization (the engine's D sector is G = 1; φ is the one at 8πG = 1).
  V0 = 0 (the COSMOS dark-energy background), φ = Π = S = PS = 0, fluid densities 1e-30 (the fluid solvers need a
  positive rest current).
- Dirichlet wall at **fixed areal radius** R = 10, re-located on every RHS call; D's derivatives see odd-reflected
  ghost values about the exact wall position, and D's stress-energy and RHS are zero beyond it (vacuum exterior).
- Data: the gear's growing mode to O(A³), D = A u + A³ w, ∂ₜD = γ(A u + 3A³ w), PD = ∂ₜD / (1 + A² α₂(r)) (normal
  derivative on the gear's polar slice, K_θ = 0, wall-normalized lapse); areal reconstruction as in
  `engine/v55_initial`. Initial momentum constraint at round-off; Hamiltonian constraint 5e-6 (N = 100), 5e-8 (N = 200).
- Observable: Γ_τ = PD(0)/D(0) along the centre worldline (β = 0), ratio (Γ_τ − γ_num)/γ_num/D_c², γ_num from an
  A = 1e-5 run at the same N corrected for its own D_c² shift. Degree-2 fit over D_c² ∈ [5e-4, 5e-3], t ≥ 1.5.
- r_max = 20, dt = 0.06 dr, T = 5.5; N = 100, 150, 200 (h = 0.2, 0.133, 0.1).

## Result (`intercepts.txt`)

| N | intercept A = 0.001 | intercept A = 0.002 | O(D_c⁴) slope |
|---|---|---|---|
| 100 | 4.58131 | 4.58137 | 11.26–11.30 |
| 150 | 4.58204 | 4.58209 | 11.04–11.07 |
| 200 | 4.58229 | 4.58234 | 10.96–10.99 |
| Richardson (observed p ≈ 2.6) | 4.58251 | 4.58256 | |

Target (GEAR-86, AMENDMENT-3-GEAR-87, polar-areal): intercept 4.58259 (prediction 4.58260), slope ≈ 10.9.

- The gap to the target converges: −1.25e-3, −5.3e-4, −2.8e-4 (order ≈ 2.6); the slope converges to 10.9. Against
  the polar-areal table (`research/gear69_polar_areal_table_N100.csv`) the N = 200 point gaps are −1.2e-4 to −3.4e-4.
- Pre-registered criterion (|I∞ − 4.58259| ≤ 5e-5 pass; ≤ 2e-4 consistent): A = 0.002 passes (−3e-5); A = 0.001 is
  consistent (−8e-5). With p = 2 instead of 2.6 the A = 0.002 extrapolation is 4.58266, so the bracket is
  4.5825–4.5827, i.e. agreement to about ±8e-5 (2e-5 relative).
- Amplitudes 0.001 and 0.002 collapse to about 5e-5 at every N. A = 0.004 (N = 100, 200 runs) carries the O(A⁴)
  initial-data transient at small D_c², the same pattern as the polar-areal A = 0.004 column.
- C = 1.443222 would shift the coefficient by 2.1e-3, about 25× this uncertainty: the engine, independently of the
  polar-areal code, supports C = 1.4432004.

## Limits

R = 10, μ = 1, G = 1 only; fixed CFL 0.06 (time-step dependence not separately bounded); the wall and centre
extrapolation are second-order harness pieces, which set the observed order. Not a test of the moving-puncture gauge.
