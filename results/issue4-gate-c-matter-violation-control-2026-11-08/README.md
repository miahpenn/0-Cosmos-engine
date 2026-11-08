# Issue #4 Gate C, matter side: matter-violation control and lapse re-solve (2026-11-08)

Read-only. `scientific_admission=false`, `acceptance_oracle_satisfied=false`. No production code,
evaluator logic, tolerance or archived state changed; zero production steps.

## Why

On the accepted states every predictor operand is truncation-sized
(`results/issue4-gate-c-accepted-convergence-2026-11-08/`), so terms proportional to matter-sourced
constraint violations were tested only weakly. This control gives H_total and M_total a matter-sourced,
resolution-independent violation and asks whether the locked predictor captures the resulting change in
the measured M_r tangent.

## What was run

`tools/issue4_gate_c_matter_violation_control.py`, workflow run
https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42506718478 (head
`69874450ef4ac0c71f80356546df48ca5f70f10d`; production engine/vendor checked equal to `4aec88d0`;
helper blobs `3f00a4e0`, `c10f9e8b`; Python 3.12.15, numpy 2.5.4, scipy 1.18.1).

1. **Build** (one process): restore the exact #420 t = 0.60 profiles (N = 60, 120, 240; artifacts of run
   40156226086, checksum-verified by the stage-source helper) and build four states per resolution:
   - `baseline`: unchanged;
   - `baseline_lapse_resolved`: unchanged matter, lapse replaced by a fresh CMC solve (kernel's own solver);
   - `matter_stored_lapse`: dark matter rest, energy_t, momentum_r scaled by (1 + 0.5 g) and baryon
     momentum_r += 0.2·energy_t·(r/3)·g, g = exp(−(r/3)²), geometry and stored lapse unchanged;
   - `matter_lapse_resolved`: the same matter change with a fresh CMC lapse.
   Each runs through the unmodified exact-profile tangent pipeline (`_run` of
   `tools/audit_exact_profile_momentum_tangent.py`).
2. **Evaluate** (fresh process, vendor from the pinned diagnostic-source checkout): the unchanged locked
   evaluator via the accepted-fixture loader, with every loader check active except the archive identity
   pins, which are replaced by each fixture's own computed digests. The N = 60 reference regression is not
   executed for these control states.

Checks: the `baseline` build reproduces the archived accepted-state and tangent-helper digests at all
three resolutions (N = 60 arrays byte-identical to the archive). Every build passes the tangent pipeline's
two-epsilon closures and total-T projection cross-checks. ε sensitivity ≤ 7e-13.

- Artifact `issue4-gate-c-matter-violation-control-42506718478` (id 13626809765), published digest
  `sha256:73f10a8102cf3d32c108725474a60bb4ef7afd8341fb6848db20cc908dbcaa24`; JSON inside SHA-256
  `73102b4c7619b39d1ed6c580e7732e05a99596b989857d73856f14c29c86fdfd` (run log).
- `gate_c_matter_violation_control_from_run_log.json`: the same report from the raw job log, SHA-256
  `55960b5566a1203eef40128857c30fc0ed530c586aecfc21ab5ed3e93824469d`.

## Results (ε = 1e-5; differences are each state minus `baseline` at the same N)

**Injected violation** (`matter_stored_lapse`): max ΔH_total = 2.44e-5, 2.53e-5, 2.55e-5 and
max ΔM_total = 3.4e-7 at N = 60, 120, 240, against baseline |H_total| 8.4e-8 and |M_total| 3.7e-8 at
N = 240. The measured tangent changes by up to 1.2e-6 at every N.

**Matter side, stored lapse: the predictor captures the change.** Change in the residual (measured −
predictor) against change in the measured tangent:

| N | Δ tangent, r = 2 | Δ residual, r = 2 (pw) | Δ residual, r = 1 (pw) | Δ residual, r = 3 (pw) | Δ residual, r = 2 (FV) |
|---|---|---|---|---|---|
| 60 | +1.19e-6 | −1.22e-8 | +6.0e-8 | −2.9e-9 | +2.4e-9 |
| 120 | +1.19e-6 | −1.31e-9 | +4.3e-9 | −8.8e-11 | −1.5e-10 |
| 240 | +1.21e-6 | −1.19e-10 | +3.7e-10 | −4.4e-12 | +2.9e-11 |

The uncaptured part falls at about third to fifth order while the change it sits on stays fixed, to 1e-4
of it at r = 2 by N = 240. At cells 0–2 the uncaptured part is ≤ 6e-9 at N = 240 (tangent change up to
7.3e-7).

**Lapse re-solve: the predictor captures none of it.** For `baseline_lapse_resolved` (no matter change),
Δ residual equals Δ tangent to the digits shown at every radius and cell. The change is confined to the
first cells and grows with resolution: cells 0–2 at N = 60, 120, 240 ≈ (+2.0, −2.2, +1.4)e-5,
(+4.2, −4.5, +2.8)e-5, (+8.4, −9.1, +5.7)e-5, doubling per doubling, while at r = 3 it falls from 4.9e-6
to 6.2e-8 to 8.0e-10. `matter_lapse_resolved` shows the sum of the two effects.

## Reading

- The matter-side coupling of the locked projected-Bianchi predictor (H_total, M_total with the matter
  substitution, total-stress source) is confirmed on the production RHS: a resolution-independent,
  matter-sourced constraint violation changes the measured M_r tangent by ~1e-6, and the predictor matches
  that change to within a residual that converges away under refinement.
- A fresh CMC lapse produces a centre defect in the M_r tangent that grows about linearly in N and that no
  predictor term represents. Its pattern (doubling per doubling, first cells only) matches the centre-row
  localization recorded in Gate A §5, where the CMC solver's first row imposes α₀ = α₁; a lapse with that
  row is not smooth at the origin, so the discrete lapse-Hessian contributions to the tangent no longer
  cancel there. The archived t = 0.60 endpoints carry the unrefreshed PIRK lapse
  (stored − fresh CMC lapse = 1.6e-4, 4.2e-5, 1.0e-5 at N = 60, 120, 240), which is why the accepted-state
  test did not show it.
- This is a diagnostic localization, consistent with Gate A §5. It does not authorize any change to the
  CMC solver, its centre row, the lapse, or the derivative stencils.

**Not settled here:** the analytic derivation of the centre defect from the α₀ = α₁ row (its coefficient
and the exact N scaling), and how much of it enters the actual PIRK evolution, where the start stage
refreshes the lapse.
