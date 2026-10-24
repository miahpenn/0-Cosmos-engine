from __future__ import annotations

import numpy as np

from engine.production_kernel import V55ProductionKernel
from tools import audit_current_source_momentum_constraint_tangent as tangent
from tools.audit_saved_state_momentum_tangent import _compact_total_t_divergence


def test_saved_state_observer_accepts_exact_kernel_geometry_terms_without_mutation():
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=40, r_max=80.0, amplitude=0.01,
        width=7.0, D_amplitude=1.0e-10, include_radiation=True,
    )
    digest_before = tangent.state_digest(state)
    terms = kernel._geometry_stage_terms(
        state.grid, state.geometry, state.scalars, state.matter,
        state.electromagnetic,
    )
    rates = tangent.rhs_direction(
        kernel, state, geometry_terms_override=terms
    )
    assert np.array_equal(rates["geometry"]["a"], terms["explicit"]["a"])
    assert np.array_equal(rates["geometry"]["b"], terms["explicit"]["b"])
    assert np.array_equal(rates["geometry"]["X"], terms["explicit"]["X"])
    assert np.array_equal(
        rates["geometry"]["Aa"],
        terms["primary_l2"]["Aa"] + terms["primary_l3"]["Aa"],
    )
    assert np.array_equal(
        rates["geometry"]["Lambda"],
        terms["lambda_l2"] + terms["lambda_l3"],
    )

    initial = tangent.observe_constraint(state)
    finite_difference = [
        tangent.tangent_for_epsilon(kernel, state, rates, eps)
        for eps in (2.0e-5, 1.0e-5)
    ]
    assert initial["M_geom_closure_ok"] is True
    assert initial["J_closure_ok"] is True
    total_t = _compact_total_t_divergence(kernel, state, rates, 1.0e-5)
    assert total_t["independent_total_projection_crosscheck_passed"] is True
    assert set(total_t["total_stress_energy_divergence"]) == {
        "div_T_t", "div_T_r", "div_T_t_fv", "div_T_r_fv"
    }
    assert set(total_t["sector_resolved_divergence"]) == {
        "dark_matter", "baryons", "radiation", "S", "D", "phi"
    }
    assert "time_density_component_accounting_passed" in total_t["sector_reconstruction"]
    for item in finite_difference:
        assert item["closure_ok"] is True
        assert item["Mdot_geom_tangent_closure_ok"] is True
        assert item["plus_minus_current_projection_closure_ok"] is True
        assert np.all(np.isfinite(item["Mdot"]))
        assert np.all(np.isfinite(item["Mdot_geom"]))
    assert tangent.state_digest(state) == digest_before
