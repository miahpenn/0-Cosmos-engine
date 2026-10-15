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


def test_lambda_chain_rule_witness_reports_signed_operand_and_exact_tf_correction():
    """The read-only Lambda derivative witness must close algebraically with signed outputs."""
    from tools.audit_current_main_centre_proper_radial_ricci import (
        _lambda_chain_rule_witness,
        proper_radial_ricci,
    )

    grid_ops, vacuum, _ = vendor_modules()
    grid = grid_ops.SphericalCellGrid(160, 8.0)
    r = np.asarray(grid.centers, dtype=float)
    c3, c5 = 1.0e-2, -2.0e-4
    areal = r * (1.0 + c3 * r**2 + c5 * r**4)
    q = areal / r
    state = vacuum.flat_state(grid)
    state.a[:] = q ** (-4.0 / 3.0)
    state.b[:] = q ** (2.0 / 3.0)
    state.X[:] = q ** (-2.0 / 3.0)
    state.alpha[:] = 1.0
    state.Lambda[:] = (
        grid.cell_derivative_fourth(state.a, parity=1) / (2.0 * state.a**2)
        - grid.cell_derivative_fourth(state.b, parity=1) / (state.a * state.b)
        + 2.0 / r * (1.0 / state.b - 1.0 / state.a)
    )

    direct = proper_radial_ricci(grid, areal, np.ones_like(r))
    terms = vacuum.geometry_terms(grid, state)
    vendor_tf = np.asarray(terms["Rr_r"] - terms["R"] / 3.0, dtype=float)
    arrays, summary = _lambda_chain_rule_witness(
        grid, state.a, state.b, state.X, state.Lambda, vendor_tf, direct["tf_ricci"]
    )

    expected = (2.0 / 3.0) * state.X**2 * arrays["lambda_prime_chain_minus_D4"]
    np.testing.assert_allclose(arrays["lambda_chain_tf_ricci_correction"], expected, rtol=0.0, atol=1.0e-15)
    np.testing.assert_allclose(
        arrays["tf_ricci_gap_chain_rule_counterfactual"],
        arrays["vendor_tf_ricci_chain_rule_counterfactual"] - direct["tf_ricci"],
        rtol=0.0,
        atol=1.0e-15,
    )
    assert summary["evolution_steps"] == 0
    assert summary["all_counterfactual_profiles_finite"] is True
    assert len(summary["first_two_cells"]) == 2
    assert summary["first_two_cells"][0]["gap_original"] != summary["first_two_cells"][0]["gap_chain_rule_counterfactual"]
    assert summary["fixed_annulus_r5_to10"]["cells"] > 0
