import numpy as np

from engine.v55_pirk_adapter import vendor_modules
from tools.audit_current_main_centre_proper_radial_ricci import proper_radial_ricci


def test_proper_radial_ricci_reconstructs_constant_curvature_geometry():
    grid_ops, _, _ = vendor_modules()
    grid = grid_ops.SphericalCellGrid(160, 8.0)
    r = np.asarray(grid.centers, dtype=float)
    curvature = 1.0e-3

    # A constant-positive-curvature 3-space in areal-radius coordinates:
    # dl^2 = dr^2/(1-k r^2), so R(l)=r has R_ll = -k R.
    areal_radius = r.copy()
    proper_factor = 1.0 / np.sqrt(1.0 - curvature * r * r)

    result = proper_radial_ricci(grid, areal_radius, proper_factor)
    interior = r <= 2.0

    np.testing.assert_allclose(
        result["R_l"][interior],
        np.sqrt(1.0 - curvature * r[interior] ** 2),
        rtol=2.0e-6,
        atol=2.0e-8,
    )
    np.testing.assert_allclose(
        result["R_ll"][interior],
        -curvature * r[interior],
        rtol=2.0e-5,
        atol=2.0e-8,
    )
    np.testing.assert_allclose(
        result["Ric_l_l"][interior],
        2.0 * curvature,
        rtol=2.0e-5,
        atol=2.0e-8,
    )
    np.testing.assert_allclose(
        result["Ric_theta_theta"][interior],
        2.0 * curvature,
        rtol=2.0e-5,
        atol=2.0e-8,
    )
    assert np.max(np.abs(result["tf_ricci"][interior])) < 5.0e-8


def test_vendor_tf_ricci_matches_proper_radial_on_constant_curvature_slice():
    grid_ops, vacuum, _ = vendor_modules()
    grid = grid_ops.SphericalCellGrid(160, 8.0)
    r = np.asarray(grid.centers, dtype=float)
    curvature = 1.0e-3
    proper_factor = 1.0 / np.sqrt(1.0 - curvature * r * r)

    # Conformal variables for gamma = f(r)^2 dr^2 + r^2 dOmega^2
    # with det(conformal gamma)=1 and areal radius R=r.
    state = vacuum.flat_state(grid)
    state.a[:] = proper_factor ** (4.0 / 3.0)
    state.b[:] = proper_factor ** (-2.0 / 3.0)
    state.X[:] = proper_factor ** (-1.0 / 3.0)
    state.alpha[:] = 1.0
    state.beta[:] = 0.0
    d1 = lambda values, parity: grid.cell_derivative_fourth(values, parity=parity)
    state.Lambda[:] = (
        d1(state.a, 1) / (2.0 * state.a * state.a)
        - d1(state.b, 1) / (state.a * state.b)
        + 2.0 / r * (1.0 / state.b - 1.0 / state.a)
    )

    direct = proper_radial_ricci(grid, r, proper_factor)
    terms = vacuum.geometry_terms(grid, state)
    vendor_tf = np.asarray(terms["Rr_r"] - terms["R"] / 3.0, dtype=float)
    interior = r <= 2.0

    # A 3-space of constant curvature has equal radial and tangential Ricci
    # eigenvalues, so its trace-free Ricci vanishes. Both representations must
    # reproduce that regular limit on the same pinned radial grid.
    assert np.max(np.abs(direct["tf_ricci"][interior])) < 5.0e-6
    assert np.max(np.abs(vendor_tf[interior])) < 5.0e-6
    assert np.max(np.abs(vendor_tf[interior] - direct["tf_ricci"][interior])) < 5.0e-6


def test_regular_smooth_profile_and_parity_negative_control():
    """A regular odd areal radius has bounded TF-Ricci/r^2; even contamination does not."""
    grid_ops, _, _ = vendor_modules()
    grid = grid_ops.SphericalCellGrid(160, 8.0)
    r = np.asarray(grid.centers, dtype=float)
    areal_radius = r * (1.0 + 1.0e-2 * r**2 + 1.0e-3 * r**4)
    proper_factor = 1.0 + 2.0e-2 * r**2 + 1.0e-3 * r**4

    regular = proper_radial_ricci(grid, areal_radius, proper_factor)
    quotient = np.abs(regular["tf_ricci"][:6]) / (r[:6] ** 2)
    assert np.all(np.isfinite(quotient))
    assert float(np.max(quotient)) < 0.2

    # The radius here includes an even-in-r term and violates regular-centre
    # odd parity; the same observer should expose a clear counterexample.
    nonregular_radius = areal_radius + 1.0e-3 * r**2
    nonregular = proper_radial_ricci(grid, nonregular_radius, proper_factor)
    regular_peak = float(np.max(np.abs(regular["tf_ricci"][:6])))
    nonregular_peak = float(np.max(np.abs(nonregular["tf_ricci"][:6])))
    assert nonregular_peak > 10.0 * regular_peak


def test_even_half_cell_integral_uses_regular_source_moments():
    """The origin half-cell integral must integrate an even source, not sample it at r0."""
    from tools.audit_current_main_centre_proper_radial_ricci import _even_half_cell_integral

    coefficients = np.array([1.25, -0.75, 0.2, -0.015], dtype=float)
    r0 = 0.125
    expected = sum(
        coefficients[k] * r0 ** (2 * k + 1) / (2 * k + 1)
        for k in range(len(coefficients))
    )
    actual = _even_half_cell_integral(coefficients, r0)
    sampled_trapezoid = r0 * sum(
        coefficients[k] * r0 ** (2 * k) for k in range(len(coefficients))
    )

    assert abs(actual - expected) < 1.0e-16
    expected_sampling_bias = sum(
        (2.0 * k / (2.0 * k + 1.0)) * coefficients[k] * r0 ** (2 * k + 1)
        for k in range(1, len(coefficients))
    )
    assert abs((sampled_trapezoid - actual) - expected_sampling_bias) < 1.0e-16


def test_lambda_chain_rule_counterfactual_reports_only_offline_shift():
    """Check the algebraic TF-Ricci effect of an offline Lambda-prime substitution."""
    from tools.audit_current_main_centre_proper_radial_ricci import (
        _lambda_chain_rule_counterfactual,
    )

    grid_ops, _, _ = vendor_modules()
    grid = grid_ops.SphericalCellGrid(48, 8.0)
    r = np.asarray(grid.centers, dtype=float)
    d1 = lambda values, parity=1: grid.cell_derivative_fourth(values, parity=parity)
    a = 1.0 + 1.0e-3 * r**2 + 1.0e-5 * r**4
    b = 1.0 - 5.0e-4 * r**2 + 2.0e-6 * r**4
    X = 1.0 + 2.0e-4 * r**2
    Lambda = (
        d1(a, 1) / (2.0 * a * a)
        - d1(b, 1) / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )
    vendor_r = 2.0e-3 + 1.0e-4 * r**2
    vendor_scalar = -1.0e-3 + 2.0e-4 * r**2
    direct_tf = 3.0e-4 - 1.0e-5 * r**2
    vendor_terms = {"Rr_r": vendor_r, "R": vendor_scalar}

    arrays, summary = _lambda_chain_rule_counterfactual(
        grid, r, a, b, X, Lambda, vendor_terms, direct_tf
    )
    delta = arrays["lambda_prime_chain_minus_D4"]
    expected_shift = (2.0 / 3.0) * X**2 * delta
    base_tf = vendor_r - vendor_scalar / 3.0

    np.testing.assert_allclose(
        arrays["lambda_chain_rule_induced_tf_ricci_shift"],
        expected_shift,
        rtol=1.0e-14,
        atol=1.0e-14,
    )
    np.testing.assert_allclose(
        arrays["lambda_chain_rule_vendor_tf_ricci"],
        base_tf + expected_shift,
        rtol=1.0e-14,
        atol=1.0e-14,
    )
    np.testing.assert_allclose(
        arrays["lambda_chain_rule_vendor_minus_direct_tf_ricci"],
        base_tf + expected_shift - direct_tf,
        rtol=1.0e-14,
        atol=1.0e-14,
    )
    assert summary["production_or_vendor_code_changed"] is False
    assert summary["evolution_steps"] == 0
    assert summary["counterfactual_profiles_finite"] is True


def test_combined_origin_and_lambda_counterfactual_adds_only_chain_rule_shift():
    """Check combined frozen-metric counterfactual arithmetic without a state edit."""
    from tools.audit_current_main_centre_proper_radial_ricci import (
        _combined_origin_and_lambda_counterfactual,
    )

    grid_ops, _, _ = vendor_modules()
    grid = grid_ops.SphericalCellGrid(48, 8.0)
    r = np.asarray(grid.centers, dtype=float)
    d1 = lambda values, parity=1: grid.cell_derivative_fourth(values, parity=parity)
    a = 1.0 + 1.0e-3 * r**2 + 1.0e-5 * r**4
    b = 1.0 - 5.0e-4 * r**2 + 2.0e-6 * r**4
    X = 1.0 + 2.0e-4 * r**2
    Lambda = (
        d1(a, 1) / (2.0 * a * a)
        - d1(b, 1) / (a * b)
        + 2.0 / r * (1.0 / b - 1.0 / a)
    )
    direct_tf = 3.0e-4 - 1.0e-5 * r**2
    origin_gap = 1.0e-3 * np.exp(-r)
    vendor_tf = direct_tf + origin_gap
    origin_arrays = {
        "origin_cf_a": a.copy(),
        "origin_cf_b": b.copy(),
        "origin_cf_X": X.copy(),
        "origin_cf_Lambda": Lambda.copy(),
        "origin_cf_vendor_tf_ricci": vendor_tf.copy(),
        "origin_cf_direct_tf_ricci": direct_tf.copy(),
        "origin_cf_tf_ricci_difference": origin_gap.copy(),
    }
    before = {key: value.copy() for key, value in origin_arrays.items()}

    arrays, summary = _combined_origin_and_lambda_counterfactual(grid, r, origin_arrays)
    expected_shift = (2.0 / 3.0) * X**2 * (
        arrays["combined_origin_lambda_cf_lambda_prime_chain_rule"]
        - arrays["combined_origin_lambda_cf_lambda_prime_discrete_D4"]
    )
    np.testing.assert_allclose(
        arrays["combined_origin_lambda_cf_induced_tf_ricci_shift"],
        expected_shift,
        rtol=1.0e-14,
        atol=1.0e-14,
    )
    np.testing.assert_allclose(
        arrays["combined_origin_lambda_cf_vendor_minus_direct_tf_ricci"],
        origin_gap + expected_shift,
        rtol=1.0e-14,
        atol=1.0e-14,
    )
    for key in before:
        np.testing.assert_array_equal(origin_arrays[key], before[key])
    assert summary["production_or_vendor_code_changed"] is False
    assert summary["evolution_steps"] == 0
    assert summary["counterfactual_profiles_finite"] is True
