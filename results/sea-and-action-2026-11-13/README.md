# Sea, action and settled D (HANDOFF to projects Nov 13 2026: GEAR-100, AMENDMENT-2-GEAR-100)

RESEARCH ONLY. No production code changed. Scripts: `research/sd_sea_v1.py` (flat pocket, optional sea field χ),
`research/sea_action_analysis_v1.py` (S's action, D split, profiles), `research/dm_one_budget_v1.py` (gravity runs; now
with a time-step option and γ_rr recorded). Full outputs: `flat_dr004.txt`, `flat_dr002.txt`, `gravity_*.txt`.

## Measures

- **S's action** J = Σ_pocket √(⟨S_τ²⟩⟨S²⟩) dV over 10-unit windows (8 in the gravity runs), S_τ in local proper time.
  For a harmonic oscillation this is E/ω with E the full oscillation energy (interaction share included), so no split
  convention enters. Constant J: retuning; falling: drain; rising: real gain.
- **D's split.** δ = D − ⟨D⟩ (oscillation about the window mean), separated at spatial wavenumber k_c = 0.5 (group
  velocity < 0.33) into (a) settled (k < k_c) and (b) travelling (k ≥ k_c). Energy and isotropic pressure to quadratic
  order with the local curvature m² = −1 + (3λ/2)(⟨D⟩² + ⟨S²⟩) (= 2 at D*). Rulers: k_c = 0.3 and 1.0.

## Unpaid inputs

| | flat (Tasks 1, 2) | gravity (Task 0) |
|---|---|---|
| λ | 5000 | 5000 |
| floor / V0 | none ("no φ / no vacuum") | one-budget: V0 and densities × 304.92, D's ejected vacuum is the floor |
| φ_today, β_D | — | 0.179055; β_D = 0, −0.04 (labelled) |
| sea mass m | √2 (D's settled frequency); 1.0 and 2.0 off-tune | — |
| sea amplitude | uniform displacement 0.1 D*, at rest (tapered to zero 5 units before the sponge) | — |
| seed, start | 1e-10, cold start | 1e-10, cold start |
| pocket | x_in = 2.0, 3.266, R = 8, edge width 1 | x_c = 1.5, 2.0 |
| rulers | dr = 0.04 / 0.02; absorbing / reflecting edge (and row / periodic row) | dr = 0.1 / 0.05; r_max = 24 / 36 |

G_N = 1 (gravity runs only); G_reset not involved.

## Task 2, first: what the sea model does to the vacuum (derived)

With Ψ = D + χ inside the quartic and +½m²χ², minimising over χ at fixed Ψ gives χ = −Ψ/(m² − 1), and Ψ feels
V(Ψ) = −½μ²Ψ² + (λ/8)Ψ⁴ with μ² = m²/(m² − 1). So the sea is not neutral:

| m | μ² | vacuum depth (× 1/(2λ)) | D, χ at the vacuum | normal-mode frequencies | S's mass² outside matter | matter keeps Ψ out for |
|---|---|---|---|---|---|---|
| no sea | 1 | 1 | D*, 0 | √2 | 4 | x > 4/3 |
| √2 | 2 | 4 | 2√2 D*, −√2 D* | 0.562, 3.562 | 7 | x > 8/3 |
| 2.0 | 4/3 | 1.78 | 1.54 D*, −0.38 D* | 0.885, 3.196 | 5 | x > 16/9 |
| 1.0 | — | unbounded below | — | — | — | — |

At m = √2 the sea does not ring at √2: once coupled through the quartic, D and χ mix into modes at 0.56 and 3.56.
At m = 1, V is linear in χ at fixed Ψ ≠ 0, so it has no floor; that control ran away in the cold start in every case
(both edges, both x_in, both dr). At m = √2 the x_in = 2.0 pocket no longer excludes Ψ (2.0 < 8/3).
If GEAR-100 means the sea to leave D's vacuum unchanged, the potential needs another term; that is a derivation, not
a setting, so it is left to chat.

RESULTS_PLACEHOLDER
