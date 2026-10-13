# Archived Accepted-State H/M/Connection Constraint Audit

**Date:** 2026-10-13 UTC  
**Scope:** Offline recomputation from the existing radiation-replay artifact; no numerical evolution and no engine/vendor/physics edits.  
**Coordination claim:** [Issue #5 comment 6137339970](https://github.com/miahpenn/0-Cosmos-engine/issues/5#issuecomment-6137339970)  
**Production source baseline:** `29d828b4d72b017c12a20b6dc90549f2511f025c`  
**Pinned geometry/PIRK vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`  
**Current accepted replay observer head:** `a854c3d6eafde1a1bae8066b7cc212bac7511ffd`

## Why this audit

The outer radiation cone is the immediate reason the saved strict probe stops, but constraint propagation is a separate production-admission requirement. This audit recomputes production Hamiltonian, momentum, conformal-connection, and determinant witnesses on the **27 full-field accepted snapshots already saved near the terminal radiation failure**, to make the central-vs-exterior constraint picture explicit without duplicating the active CMC target/operator audit.

## Provenance and implementation cross-checks

- Reused [PR #61 frozen replay run 38337821862](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38337821862), [artifact 11785239489](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38337821862/artifacts/11785239489).
- Artifact ZIP SHA-256: `5c3dee4cbe62fb410ed7bcdf97707385019e00bae4f925ab51d4242b42933d09`. All five entries listed in the embedded `SHA256SUMS.txt` were recomputed and matched.
- Exact source/vendor pins and the capture schema were checked. The saved history contains 27 accepted states from `t=47.130000000001026` through `t=47.15300258518951` at N=80, `r_max=80`.
- Stress was reconstructed from the saved scalar profiles and recovered dark matter, baryon, and radiation states with the same production primitive/projection conventions. Against the captured total-stress inputs for the same physical slice, the maximum absolute differences were `2.1684e-19` (`rho`), `3.4694e-18` (`pr`), and `4.6640e-20` (`pt`). The saved metric profiles match the accepted terminal slice bitwise in `a,b,X,Aa,K,Lambda`; only the lapse changes during the stage CMC solve.
- The reproduced pinned Ricci expression `alpha*(R^r_r - R/3)` agrees with the archived L2 trace-free Ricci witness to maximum absolute difference `1.4633e-15` (maximum relative difference to the output scale `4.2871e-13`). This is an implementation consistency check, not an independent continuum-curvature validation.

The verifier is [`tools/audit_radiation_accepted_constraint_history.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/radiation-accepted-constraint-history-2026-10-13/tools/audit_radiation_accepted_constraint_history.py). It validates the exact ZIP digest, internal manifest, provenance, state shapes, stress-source comparison, and finite outputs, then prints per-snapshot signed witnesses in JSON. It imports no engine code and advances no state.

## Constraint definitions used

The expressions match the production diagnostic assembly:

- `H = R - (Aa² + 2 Ab²) + (2/3)K² - 16π rho_total`, with `Ab=-Aa/2`.
- `M = D_r Aa - (2/3)D_r K + 6 Aa chi_r + (Aa-Ab)(2/r + D_r b/b) - 8π j_total`, with `chi_r=-D_r X/(2X)`.
- `C_Lambda = Lambda - [D_r a/(2a²) - D_r b/(ab) + 2/r*(1/b - 1/a)]`.
- Conformal determinant witness: `a*b² - 1`.

All derivatives use the pinned vendor's cell-centred fourth-order first/second derivative and its one-sided closure for the last two cells. Residuals are kept signed; the `r>=5` witnesses align with the production fixed-physical inner-exclusion convention. The terminal outer band is shown separately.

## Findings: first and last saved states

Values below are maxima by absolute value, with the signed value and cell/radius identified.

| Witness | First saved state, t=47.130000000001 | Terminal saved state, t=47.153002585190 |
|---|---:|---:|
| Full-domain H | `-0.756017` at r=0.5 | `-0.759857` at r=0.5 |
| H on r>=5 | `-0.0109945` at r=63.5 | `-0.0110070` at r=63.5 |
| H on outer 5 cells | `-0.00194176` at r=79.5 | `-0.00198971` at r=79.5 |
| Full-domain M | `+0.0618083` at r=0.5 | `+0.0622176` at r=0.5 |
| M on r>=5 | `+0.00118709` at r=79.5 | `+0.00123655` at r=79.5 |
| M on outer 5 cells | `+0.00118709` at r=79.5 | `+0.00123655` at r=79.5 |
| Full-domain C_Lambda | `+0.164105` at r=0.5 | `+0.164781` at r=0.5 |
| C_Lambda on r>=5 | `-0.00249855` at r=64.5 | `-0.00260206` at r=64.5 |
| C_Lambda on outer 5 cells | `+0.000797839` at r=79.5 | `+0.000860939` at r=79.5 |
| `a*b²-1` | `+1.63893e-5` at r=70.5 | `+1.64512e-5` at r=70.5 |

The saved central-cell term budget at the terminal state is:

- Hamiltonian: Ricci `R=-0.8039106`; trace-free-curvature term `-2.51045e-5`; `+(2/3)K²=+0.1015356`; matter `-16πrho=-0.0574566`; total `H=-0.7598567`.
- Momentum: `D_r Aa=-0.00583734`; `-(2/3)D_r K=+0.0515255`; `6Aa chi_r=+0.00076549`; metric-ratio term `+0.0243859`; matter-current term `-8πj=-0.00862193`; total `M=+0.0622176`.
- Connection: `Lambda=+0.2537241`; subtracted metric-connection expression `-0.08894335`; total `C_Lambda=+0.1647807`.

## Evolution across the saved terminal window

Across these 27 accepted snapshots, the full-domain peaks move only modestly: `|H|` at the first cell grows about 0.51%, `|M|` about 0.66%, and `|C_Lambda|` about 0.41%. On `r>=5`, the maximum exterior `|H|` grows about 0.11%, exterior `|M|` about 4.17%, and exterior `|C_Lambda|` about 4.14%. The outer five-cell connection maximum increases about 7.91%. The determinant maximum is around `1.64e-5` throughout.

Thus the saved terminal window contains a large central constraint residual and smaller but nonzero exterior/outer residuals. Importantly, the central residuals are already large in its **first saved state**, at `t=47.13`; this artifact has no earlier accepted snapshots, so it cannot tell when during the full trajectory they first developed. The radiation cone approaches its realizability boundary at the outermost cell in this same window, but co-observation does not establish that either defect causes the other.

## Disposition and next gate

This result strengthens the case that the production machine has **more than one unadmitted numerical witness**: a central Hamiltonian/momentum/connection constraint problem is present alongside the late outer radiation-predictor failure. Do not attribute the central residual to the CMC operator mismatch without a matched term-level test, and do not attribute the radiation crossing to the central residual based on correlation alone.

The fixed-slice CMC K-operator / centre-row review remains with its current owner. A useful next independent discriminator is a source-pinned comparison of these same constraint definitions on an earlier complete saved state or existing t<=12 campaign artifact, with identical operators, fixed radii, and signed term budgets. No new long run is justified until the team identifies a specific failing discrete identity and a bounded counterfactual that improves the relevant admission witnesses.
