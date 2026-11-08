# Issue #4 centre defect: attribution to the CMC lapse solver's centre row (2026-11-08)

Read-only diagnostic. `scientific_admission=false`, `acceptance_oracle_satisfied=false`. No production
code, vendor code, archived state or tolerance changed. The ghost-row lapse below exists only as a
diagnostic copy inside the driver.

## Question

The matter-violation control (`results/issue4-gate-c-matter-violation-control-2026-11-08/`) found that a
fresh CMC lapse adds a centre defect to the M_r tangent that grows with resolution and that no predictor
term represents. Two hypotheses: (H1) the vendor operators (parity ghosts, lapse Hessian/Laplacian)
mishandle a lapse at the origin; (H2) the CMC solver's centre row makes the lapse non-smooth at the origin
and correct operators then expose it.

## What was run

`tools/issue4_gate_c_matter_violation_control.py --study lapse`, workflow run
https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42507992361 (head
`099c4a13f363fed1c4f634bdab5c1fd35bbaa064`; production engine/vendor equal to `4aec88d0`; same pins and
environment as the matter control). Matter unchanged; on each archived t = 0.60 state (N = 60, 120, 240)
the lapse is replaced by:

| Variant | Lapse used |
|---|---|
| `stored` (baseline) | archived lapse |
| `lapse_smooth_bump` | stored + 1e-4·exp(−(r/3)²) |
| `lapse_row_imposed` | stored, with only α₀ := α₁ |
| `lapse_cmc` | fresh solve, production `engine/cmc_gauge.py:solve_cmc_lapse` |
| `lapse_cmc_centre_even` | fresh solve with α₀ replaced by its even quartic extrapolation from cells 1–3 |
| `lapse_cmc_ghost_row` | fresh solve from a diagnostic copy of `solve_cmc_lapse`, built from its own source with only the centre row changed (below) |

Each state goes through the unmodified exact-profile tangent pipeline; the change in the measured M_r
tangent relative to `stored` is compared with an **independent prediction**: the lapse-derivative terms of
the A_a and K rates (−(D_rD_r − ∇²/3)δα and −∇²δα, transcribed from the vendor formulas) carried through
the vendor M_r chain rule, D₁(δL^{A_a}) − ⅔D₁(δL^K) + (6φ′ + 1.5q)δL^{A_a}, using fourth-order
parity-ghost stencils implemented separately in the driver (unit-checked: errors fall 16× per doubling).

The production centre row (`engine/cmc_gauge.py`, lines 215–217):

    # Regular center: alpha'(0)=0.
    band[1, 0] = 1.0
    band[0, 1] = -1.0          (right-hand side 0: alpha0 - alpha1 = 0)

The diagnostic ghost row applies the same three-point equation at cell 0 with the even reflection
α₋₁ = α₀ (the convention every other parity-ghost operator in the engine uses):

    band[1, 0] = -1/h^2 - c[0]/(2h) - m[0];  band[0, 1] = 1/h^2 + c[0]/(2h);  bvec[0] = rhs[0]

Provenance: artifact `issue4-gate-c-matter-violation-control-42507992361` (id 13627827839), published
digest `sha256:66181672354fd8f1cc12fd8f392d5648e59364760d15a5e8685633276c04f777`, JSON inside SHA-256
`4f427a66763a844d6d0c51eed6849226ee6022b91aeadd1953f4eea98988c077`;
`centre_lapse_study_from_run_log.json` here is the same report from the raw job log, SHA-256
`3180807ea45e007752dd42f9e48c8e83ed233ac5fd1f9e2c4ddfdb582b05d703`. Every build passes the pipeline's
two-epsilon tangent closures.

## Results (ε = 1e-5; change in measured M_r tangent relative to `stored`)

| Variant | max, r ≤ 10, N = 60 | N = 120 | N = 240 | trend | measured − predicted, max (N = 240) |
|---|---|---|---|---|---|
| `lapse_cmc` (production row) | 2.21e-5 | 4.50e-5 | 9.05e-5 | doubles per doubling | 1.3e-8 |
| `lapse_row_imposed` | 2.68e-5 | 5.47e-5 | 1.10e-4 | doubles per doubling | 1.2e-8 |
| `lapse_cmc_centre_even` | 9.2e-7 | 2.87e-6 | 5.98e-6 | grows | 3.5e-9 |
| `lapse_smooth_bump` | 1.45e-6 | 3.46e-7 | 5.1e-8 | falls | 2.1e-8 |
| `lapse_cmc_ghost_row` | 1.37e-9 | 7.2e-10 | 3.6e-10 | falls | 4.9e-13 |

All variants peak in cells 0–2. Lapse centre structure:

| | (α₁ − α₀), N = 60 / 120 / 240 | α₀ − even extrapolation from cells 1–3 |
|---|---|---|
| stored | 1.10e-4 / 2.79e-5 / 7.01e-6 (∝ h², smooth) | −4.6e-6 / −8.8e-8 / −1.5e-9 |
| production CMC | 0 / 0 / 1e-15 (imposed) | 8.9e-5 / 2.4e-5 / 6.0e-6 (∝ h²) |
| ghost-row CMC | 1.10e-4 / 2.79e-5 / 7.01e-6; (α₁ − α₀)/h² = 6.17e-5, 6.28e-5, 6.31e-5 | −4.6e-6 / −8.8e-8 / −1.5e-9 |

The ghost-row solution differs from the production CMC solution by 1.63e-4, 4.17e-5, 1.05e-5, the same as
the stored-minus-production difference, and it reproduces the stored lapse's centre structure to the digits
shown.

## 0* audit (Gate A §6)

| Point | Finding |
|---|---|
| (i) Source and convention pin | Production solver `engine/cmc_gauge.py` at `4aec88d0`, centre row lines 215–217; vendor `d6052d60`; pipeline helpers `3f00a4e0`, `c10f9e8b`; run 42507992361. |
| (ii) Independent derivation / witness | Separately implemented stencils predict the measured centre change for every variant at every N to ≤ 4.3e-8 absolute (≤ 2e-3 relative; ≤ 5e-13 for the ghost row). The vendor operators reproduce standard fourth-order stencil algebra; the smooth bump's effect falls under refinement. **H1 is not supported.** |
| (iii) Explicit discrete cost | The row α₀ = α₁ moves α₀ off the smooth even profile by c·h² (8.9e-5, 2.4e-5, 6.0e-6). The fourth-order lapse Hessian turns that into an O(1) error in the A_a rate at cells 0–1, and the M_r chain rule's derivative into an O(1/h) error: 2.2e-5, 4.5e-5, 9.1e-5. |
| (iv) No inferred term | No fit and no new coefficient; the prediction uses only the vendor formulas and the stencil definitions. |
| (v) Falsifiable same-state comparison | Imposing α₀ = α₁ alone on the smooth stored lapse reproduces the defect; replacing only the centre row with the parity-ghost row removes it (250,000× smaller at N = 240 and falling). **H2 holds.** |

## Reading

The growing centre defect in the M_r tangent is caused by the CMC lapse solver's centre row, which imposes
α₀ = α₁ on a cell-centred grid where a smooth even lapse has α₁ − α₀ = α″(0)·h². It is not a vendor-operator
fault: the vendor's parity-ghost operators behave exactly as the stencil algebra predicts, and they are what
exposes the row. A row built from the same equation with the even ghost α₋₁ = α₀ gives a smooth lapse and no
growing defect on the same states.

Replacing α₀ alone after the production solve (`lapse_cmc_centre_even`) is not enough: the row's influence
reaches cells 1–3 through the banded solve, and the defect still grows, at about 15× smaller size.

This is the recorded, signed, localized mismatch that Gate D asks for. It does not by itself authorize the
change: replacing the production centre row is a production edit for owner review, and it would need its
own pinned evolution comparison (the solver also refreshes the lapse at each PIRK stage). The solver's
interior rows are three-point (second order); that is a separate property and was not changed here.
