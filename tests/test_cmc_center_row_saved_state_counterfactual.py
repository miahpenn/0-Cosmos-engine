from __future__ import annotations

import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.production_kernel import V55ProductionKernel
from tools import audit_evolved_origin_aa_regularity as source_audit
from tools.audit_cmc_center_row_saved_state_counterfactual import (
    _operator_compatible_alpha,
    _replace_only_lapse,
)


def test_copied_cmc_center_pde_row_closes_without_mutating_frozen_state():
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=60, r_max=80.0, include_radiation=True)
    before = source_audit._state_digest(state)

    alpha, solve_report = _operator_compatible_alpha(state)
    assert alpha.shape == state.grid.centers.shape
    assert np.all(np.isfinite(alpha))
    assert np.all(alpha > 0.0)
    assert alpha[-1] == 1.0
    assert solve_report["raw_operand_scaled_rows_pass"] is True
    assert state.geometry.alpha[0] == state.geometry.alpha[1]

    counterfactual = _replace_only_lapse(state, alpha)
    assert source_audit._state_digest(state) == before

    _, vacuum, _ = adapter.vendor_modules()
    baseline_constraints = vacuum.constraints(state.grid, state.geometry)
    counterfactual_constraints = vacuum.constraints(
        counterfactual.grid, counterfactual.geometry
    )
    for name in baseline_constraints:
        assert np.array_equal(baseline_constraints[name], counterfactual_constraints[name])

    l2 = vacuum.primary_l2_rhs(state.grid, counterfactual.geometry)
    l3 = adapter.primary_l3_with_matter(
        state.grid, counterfactual.geometry,
        counterfactual.scalars, counterfactual.matter,
    )
    kdot = l2["K"] + l3["K"]
    target = solve_report["target_kdot_from_existing_CMC_projection"]
    assert np.all(np.isfinite(kdot))
    assert np.max(np.abs(kdot[:-1] - target)) < 2.0e-11
    assert source_audit._state_digest(state) == before
