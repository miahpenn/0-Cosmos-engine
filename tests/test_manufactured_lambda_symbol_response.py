import numpy as np
import pytest

from tools.audit_manufactured_lambda_symbol_response import _manufactured_test_case


def test_manufactured_odd_lambda_sweep_closes_response_and_chain_rule_split():
    report, profiles = _manufactured_test_case(160)

    assert report["all_profiles_finite"] is True
    assert report["all_response_identities_close"] is True
    modes = report["mode_results"]
    assert [x["wavelength"] for x in modes] == [32.0, 16.0, 8.0]
    assert all(x["max_abs_response_closure_error"] < 5.0e-12 for x in modes)
    assert all(x["max_abs_chain_rule_algebraic_only_closure_error"] < 5.0e-12 for x in modes)
    assert all(x["max_abs_Aa_L2_RHS_response_identity_closure"] < 5.0e-12 for x in modes)
    assert all(x["max_abs_chainrule_Aa_L2_algebraic_only_closure_error"] < 5.0e-12 for x in modes)
    assert all(np.isfinite(x["max_abs_Lambda_L3_RHS_delta"]) for x in modes)

    # The centred D4 operator's measured response matches its known modified
    # wavenumber symbol on interior samples; error grows as the mode is shortened.
    for row in modes:
        assert row["measured_D4_to_analytic_derivative_norm_ratio"] == pytest.approx(
            row["fourth_order_modified_wavenumber_ratio_theory"], rel=2.0e-12, abs=2.0e-12
        )
    errors = [x["D4_analytic_derivative_relative_rms_error_fixed_interior"] for x in modes]
    assert errors[0] < errors[1] < errors[2]

    # A smooth odd connection perturbation is O(r) at the origin. The
    # chain-rule arm removes its D4 derivative response but retains algebraic
    # Lambda*D1(a) response; therefore its change must equal the algebraic term.
    assert np.all(np.isfinite(profiles["lambda32_delta_lambda"]))
    assert np.all(np.isfinite(profiles["lambda32_response_closure_error"]))


def test_manufactured_symbol_response_is_linear_in_off_constraint_amplitude():
    from tools.audit_manufactured_lambda_symbol_response import (
        _evaluate_modes,
        _modified_wavenumber_ratio,
    )
    from tools.audit_current_main_centre_proper_radial_ricci import proper_radial_ricci
    from engine.v55_pirk_adapter import vendor_modules

    grid_ops, vacuum, _ = vendor_modules()
    grid = grid_ops.SphericalCellGrid(160, 80.0)
    r = np.asarray(grid.centers, dtype=float)
    q = 1.0 + 1.0e-5 * r**2 + 2.0e-8 * r**4
    qp = 2.0e-5 * r + 8.0e-8 * r**3
    a, b, X = q**(-4.0/3.0), q**(2.0/3.0), q**(1.0/3.0)
    lam = (
        (-4.0/3.0) * a * qp / q / (2.0 * a**2)
        - ((2.0/3.0) * b * qp / q) / (a*b)
        + 2.0/r * (1.0/b - 1.0/a)
    )
    small, _ = _evaluate_modes(
        grid, vacuum, a, b, X, lam, wavelengths=(8.0,), amplitude=1.0e-4
    )
    large, _ = _evaluate_modes(
        grid, vacuum, a, b, X, lam, wavelengths=(8.0,), amplitude=2.0e-4
    )
    sm = small["mode_results"][0]
    lg = large["mode_results"][0]
    assert lg["max_abs_total_observed_Ricci_gap_delta"] == pytest.approx(
        2.0 * sm["max_abs_total_observed_Ricci_gap_delta"], rel=2.0e-12, abs=2.0e-15
    )
    assert lg["max_abs_derivative_Ricci_response"] == pytest.approx(
        2.0 * sm["max_abs_derivative_Ricci_response"], rel=2.0e-12, abs=2.0e-15
    )
    assert lg["max_abs_algebraic_Ricci_response"] == pytest.approx(
        2.0 * sm["max_abs_algebraic_Ricci_response"], rel=2.0e-12, abs=2.0e-15
    )
    assert _modified_wavenumber_ratio(2.0*np.pi*grid.dr/8.0) > 0.99
