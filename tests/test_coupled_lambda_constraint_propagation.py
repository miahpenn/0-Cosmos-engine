import numpy as np
import pytest

from engine.production_kernel import V55ProductionKernel
from engine.v55_pirk_adapter import vendor_modules
from tools.audit_coupled_lambda_constraint_propagation import _run_pair
from tools.audit_current_main_centre_proper_radial_ricci import _state_digest


def _state():
    kernel = V55ProductionKernel()
    state = kernel.initialize(
        resolution=40, r_max=80.0, amplitude=0.01,
        width=7.0, D_amplitude=1.0e-10, include_radiation=True,
    )
    return kernel, state


def test_one_step_lambda_perturbation_keeps_constraints_signed_and_finite():
    _, vacuum, _ = vendor_modules()
    kernel, state = _state()
    before = _state_digest(state)
    dt = 0.01 * float(state.grid.dr)
    case, profiles = _run_pair(kernel, state, vacuum, wavelength=8.0, dt=dt)

    assert _state_digest(state) == before
    assert case["original_baseline_state_unchanged"] is True
    assert case["baseline_input_state_unchanged"] is True
    assert case["off_mode_input_state_unchanged"] is True
    assert case["finite"] is True
    assert case["delta_Aa_final_max_abs"] > 0.0
    assert case["actual_Aa_rate_delta_vs_frozen_L2_prediction_relative_rms_error"] < 0.25

    # Initial H/M/determinant should not change when only Lambda is perturbed;
    # the connection residual should change by precisely that odd mode.
    for name in ("hamiltonian", "momentum", "determinant"):
        off = profiles[f"initial_constraint_off_mode_{name}"]
        base = profiles[f"initial_constraint_baseline_{name}"]
        assert np.max(np.abs(off - base)) < 2.0e-13
    connection_delta = (
        profiles["initial_constraint_off_mode_connection"]
        - profiles["initial_constraint_baseline_connection"]
    )
    assert np.allclose(
        connection_delta, profiles["delta_lambda_input"],
        rtol=0.0, atol=5.0e-13,
    )
    for key in (
        "hamiltonian", "momentum", "connection", "determinant"
    ):
        assert np.all(np.isfinite(profiles[f"constraint_delta_{key}"]))


def test_smaller_dt_does_not_degrade_frozen_l2_response_match_materially():
    _, vacuum, _ = vendor_modules()
    kernel, state = _state()
    dr = float(state.grid.dr)
    coarse, _ = _run_pair(kernel, state, vacuum, wavelength=8.0, dt=0.01*dr)
    fine, _ = _run_pair(kernel, state, vacuum, wavelength=8.0, dt=0.005*dr)
    e0 = coarse["actual_Aa_rate_delta_vs_frozen_L2_prediction_relative_rms_error"]
    e1 = fine["actual_Aa_rate_delta_vs_frozen_L2_prediction_relative_rms_error"]
    assert np.isfinite(e0) and np.isfinite(e1)
    assert e1 <= 1.25 * e0 + 1.0e-6
