"""Independent analytic and differential controls for native spherical geometry."""

from types import SimpleNamespace

import numpy as np
import pytest

from diagnostics.native_spherical_operators import NativeSphericalCellGrid
from diagnostics.native_bssn_geometry import (
    metric_connection_lambda,
    native_bssn_constraints,
    native_bssn_geometry_terms,
)
from engine.v55_pirk_adapter import vendor_modules


@pytest.fixture(scope="module")
def pinned_vacuum():
    _, vacuum, _ = vendor_modules()
    return vacuum


def _state(grid, *, flat_chart=False, alpha_profile=None, lambda_offset=0.0):
    r = grid.centers
    if flat_chart:
        # A regular flat 3-geometry written in a nonlinear areal-radius chart:
        # R=r(1+c r^2), gamma_rr=(dR/dr)^2, gamma_OmegaOmega=R^2.
        # X is chosen so the conformal determinant a*b^2 equals one.
        c = 0.01
        areal = r * (1.0 + c * r**2)
        d_areal = 1.0 + 3.0 * c * r**2
        X = (d_areal**2 * (areal / r)**4) ** (-1.0 / 6.0)
        a = X**2 * d_areal**2
        b = X**2 * (areal / r)**2
    else:
        # A smooth non-flat BSSN profile that preserves the conformal
        # determinant a*b^2=1, so the vendor comparison tests the intended
        # reference-metric constraint class.
        f = 0.001 * r**2 + 2.0e-5 * r**4
        a = np.exp(2.0 * f)
        b = np.exp(-f)
        X = 1.0 + 0.0005 * r**2 + 1.0e-5 * r**4
        areal = r * np.sqrt(b) / X

    alpha = 1.0 + 0.001 * r**2 + 2.0e-5 * r**4 if alpha_profile is None else alpha_profile(areal)
    state = SimpleNamespace(
        a=np.asarray(a, dtype=float),
        b=np.asarray(b, dtype=float),
        X=np.asarray(X, dtype=float),
        alpha=np.asarray(alpha, dtype=float),
        beta=np.zeros_like(r),
        Aa=np.zeros_like(r),
        K=np.zeros_like(r),
        Lambda=np.zeros_like(r),
        B=np.zeros_like(r),
    )
    state.Lambda = metric_connection_lambda(grid, state) + lambda_offset
    return state


def _annulus(grid):
    r = grid.centers
    return (r >= 2.0) & (r <= 6.0)


def test_native_geometry_flat_spherical_chart_has_zero_curvature_and_constraints():
    grid = NativeSphericalCellGrid(80, 12.0)
    state = SimpleNamespace(
        a=np.ones(grid.n),
        b=np.ones(grid.n),
        X=np.ones(grid.n),
        alpha=np.ones(grid.n),
        beta=np.zeros(grid.n),
        Aa=np.zeros(grid.n),
        K=np.zeros(grid.n),
        Lambda=np.zeros(grid.n),
        B=np.zeros(grid.n),
    )
    geometry = native_bssn_geometry_terms(grid, state)
    constraints = native_bssn_constraints(grid, state)
    for name in ("R", "R_metric", "Rr_r", "Rtheta_theta_metric", "DrDr_alpha", "lap_alpha"):
        assert np.max(np.abs(geometry[name])) < 2e-10, name
    for name in ("H", "M", "connection", "determinant"):
        assert np.max(np.abs(constraints[name])) < 2e-10, name


@pytest.mark.parametrize("n", [48, 96, 192])
def test_nonlinear_areal_chart_flat_metric_and_lapse_match_analytic_geometry(n):
    grid = NativeSphericalCellGrid(n, 12.0)
    q = 3.0e-3
    s = 2.0e-5
    state = _state(
        grid,
        flat_chart=True,
        alpha_profile=lambda R: 1.0 + q * R**2 + s * R**4,
    )
    geometry = native_bssn_geometry_terms(grid, state)
    constraints = native_bssn_constraints(grid, state)
    mask = _annulus(grid)
    R = geometry["areal_radius"]

    assert np.max(np.abs(geometry["R"][mask])) < 2e-7
    assert np.max(np.abs(geometry["Rr_r"][mask])) < 2e-7
    assert np.max(np.abs(geometry["Rtheta_theta_metric"][mask])) < 2e-7
    np.testing.assert_allclose(
        geometry["DrDr_alpha"][mask],
        2.0 * q + 12.0 * s * R[mask]**2,
        rtol=4e-4,
        atol=2e-7,
    )
    np.testing.assert_allclose(
        geometry["lap_alpha"][mask],
        6.0 * q + 20.0 * s * R[mask]**2,
        rtol=4e-4,
        atol=2e-7,
    )
    assert np.max(np.abs(constraints["connection"][mask])) < 2e-10
    assert np.max(np.abs(constraints["determinant"][mask])) < 3e-14


def test_flat_chart_lapse_hessian_errors_refine_on_fixed_annulus():
    errors = []
    for n in (48, 96, 192):
        grid = NativeSphericalCellGrid(n, 12.0)
        q = 3.0e-3
        s = 2.0e-5
        state = _state(
            grid,
            flat_chart=True,
            alpha_profile=lambda R: 1.0 + q * R**2 + s * R**4,
        )
        geometry = native_bssn_geometry_terms(grid, state)
        mask = _annulus(grid)
        R = geometry["areal_radius"]
        expected_radial = 2.0 * q + 12.0 * s * R**2
        expected_laplace = 6.0 * q + 20.0 * s * R**2
        err = max(
            float(np.max(np.abs(geometry["DrDr_alpha"][mask] - expected_radial[mask]))),
            float(np.max(np.abs(geometry["lap_alpha"][mask] - expected_laplace[mask]))),
        )
        errors.append(err)
    assert errors[0] / errors[1] > 5.0
    assert errors[1] / errors[2] > 5.0


def test_connection_offset_correction_sign_matches_derived_identity():
    grid = NativeSphericalCellGrid(80, 12.0)
    base = _state(grid, flat_chart=True, alpha_profile=lambda R: np.ones_like(R))
    reference = native_bssn_geometry_terms(grid, base)
    perturbed = SimpleNamespace(**vars(base))
    epsilon = 1.25e-3
    delta = epsilon * grid.centers
    perturbed.Lambda = base.Lambda + delta
    actual = native_bssn_geometry_terms(grid, perturbed)
    mask = _annulus(grid)

    da = grid.cell_derivative_fourth(base.a, parity=1)
    db = grid.cell_derivative_fourth(base.b, parity=1)
    expected_rr = base.X**2 * (
        grid.cell_derivative_fourth(delta, parity=-1)
        + 0.5 * (da / base.a) * delta
    )
    expected_scalar = base.X**2 * grid.cell_derivative_fourth(delta, parity=-1)
    np.testing.assert_allclose(
        (actual["R"] - reference["R"])[mask],
        expected_scalar[mask],
        rtol=2e-11,
        atol=2e-12,
    )
    np.testing.assert_allclose(
        (actual["Rr_r"] - reference["Rr_r"])[mask],
        expected_rr[mask],
        rtol=2e-11,
        atol=2e-12,
    )
    assert np.max(np.abs(expected_rr[mask])) > 1.0e-4
    # Wrong-sign controls: neither correction may agree with its negation.
    assert np.max(np.abs((actual["Rr_r"] - reference["Rr_r"])[mask] + expected_rr[mask])) > 1.0e-4
    assert np.max(np.abs(expected_scalar[mask])) > 1.0e-4
    assert np.max(np.abs((actual["R"] - reference["R"])[mask] + expected_scalar[mask])) > 1.0e-4


def test_matter_sources_are_explicit_and_do_not_mutate_frozen_state():
    grid = NativeSphericalCellGrid(48, 12.0)
    state = _state(grid, flat_chart=True, alpha_profile=lambda R: np.ones_like(R))
    arrays_before = {
        name: np.asarray(getattr(state, name)).copy()
        for name in ("a", "b", "X", "alpha", "Aa", "K", "Lambda")
    }
    zero = native_bssn_constraints(grid, state)
    sourced = native_bssn_constraints(
        grid,
        state,
        rho=np.full(grid.n, 1.0e-8),
        radial_momentum_source=np.full(grid.n, -2.0e-9),
    )
    np.testing.assert_allclose(
        sourced["H"] - zero["H"], -16.0 * np.pi * 1.0e-8,
        rtol=1e-12, atol=1e-15,
    )
    np.testing.assert_allclose(
        sourced["M"] - zero["M"], 8.0 * np.pi * 2.0e-9,
        rtol=1e-12, atol=1e-15,
    )
    for name, before in arrays_before.items():
        np.testing.assert_array_equal(getattr(state, name), before)




def test_constant_positive_curvature_warped_product_has_correct_ricci_sign_and_refines():
    """Independent sign witness: R(l)=L sin(l/L) has scalar curvature +6/L^2.

    Use proper radial distance l=r (A=gamma_rr=1). Choosing X and a,b so that
    a/X^2=1, r*sqrt(b)/X=L*sin(r/L), and a*b^2=1 exercises the native metric
    representation without using vendor outputs as the oracle.
    """
    curvature_radius = 5.0
    errors = []
    for n in (48, 96, 192):
        grid = NativeSphericalCellGrid(n, 3.0)
        r = grid.centers
        areal = curvature_radius * np.sin(r / curvature_radius)
        ratio = areal / r

        # These choices give A=a/X^2=1 and R_areal=r*sqrt(b)/X=areal,
        # while preserving the candidate conformal determinant a*b^2=1.
        X = (r / areal) ** (2.0 / 3.0)
        a = X**2
        b = ratio ** (2.0 / 3.0)
        state = SimpleNamespace(
            a=a,
            b=b,
            X=X,
            alpha=np.ones_like(r),
            beta=np.zeros_like(r),
            Aa=np.zeros_like(r),
            K=np.zeros_like(r),
            Lambda=np.zeros_like(r),
            B=np.zeros_like(r),
        )
        state.Lambda = metric_connection_lambda(grid, state)

        geometry = native_bssn_geometry_terms(grid, state)
        mask = (r >= 0.5) & (r <= 2.5)
        expected_radial = 2.0 / curvature_radius**2
        expected_scalar = 6.0 / curvature_radius**2

        assert np.max(np.abs(a * b**2 - 1.0)) < 3e-14
        assert np.all(geometry["Rr_r_metric"][mask] > 0.0)
        assert np.all(geometry["Rtheta_theta_metric"][mask] > 0.0)
        assert np.all(geometry["R_metric"][mask] > 0.0)

        errors.append(max(
            float(np.max(np.abs(geometry["Rr_r_metric"][mask] - expected_radial))),
            float(np.max(np.abs(geometry["Rtheta_theta_metric"][mask] - expected_radial))),
            float(np.max(np.abs(geometry["R_metric"][mask] - expected_scalar))),
        ))

    observed_order = np.log2(np.asarray(errors[:-1]) / np.asarray(errors[1:]))
    assert np.all(observed_order > 3.5), (
        f"positive-curvature Ricci observed orders: {observed_order}; errors={errors}"
    )

def test_native_geometry_matches_pinned_vendor_away_from_boundaries_under_refinement(pinned_vacuum):
    errors = []
    for n in (48, 96, 192):
        grid = NativeSphericalCellGrid(n, 12.0)
        state = _state(
            grid,
            flat_chart=False,
            lambda_offset=1.0e-3 * grid.centers * (1.0 + 2.0e-3 * grid.centers**2),
        )
        native = native_bssn_geometry_terms(grid, state)
        reference = pinned_vacuum.geometry_terms(grid, state)
        mask = _annulus(grid)
        names = (
            ("R", "R"),
            ("Rr_r", "Rr_r"),
            ("DrDr_alpha", "DrDr_alpha"),
            ("lap_alpha", "lap_alpha"),
        )
        err = max(
            float(np.max(np.abs(native[key][mask] - reference[ref_key][mask])))
            for key, ref_key in names
        )
        errors.append(err)
    # Distinct discrete realizations of the same continuum identities should
    # converge together; this is not a request for bitwise equality.
    assert errors[0] / errors[1] > 2.0
    assert errors[1] / errors[2] > 2.0


@pytest.mark.parametrize("bad", ["a", "b", "X", "alpha", "Lambda"])
def test_native_geometry_rejects_invalid_states(bad):
    grid = NativeSphericalCellGrid(16, 12.0)
    state = _state(grid)
    setattr(state, bad, np.full(grid.n - 1, 1.0))
    with pytest.raises(ValueError, match="shape"):
        native_bssn_geometry_terms(grid, state)

    state = _state(grid)
    setattr(state, bad, np.full(grid.n, np.nan))
    with pytest.raises(ValueError, match="finite"):
        native_bssn_geometry_terms(grid, state)
