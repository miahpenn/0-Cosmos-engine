"""Exact geometric controls for the independent covariant stress divergence."""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("COSMOS_PINNED_SOURCE_ROOT", str(ROOT))
SPEC = importlib.util.spec_from_file_location(
    "_issue4_covariant_total_t_audit",
    ROOT / "tools" / "audit_total_stress_energy_divergence.py",
)
if SPEC is None or SPEC.loader is None:
    raise ImportError("cannot load independent total-T divergence auditor")
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def _d1(values: np.ndarray, parity: int) -> np.ndarray:
    """Polynomial-exact five-point derivative for the analytic test profiles."""
    del parity  # The test functions are explicit polynomials on the positive half-line.
    values = np.asarray(values, dtype=np.float64)
    n = values.size
    dr = 0.25
    if values.ndim != 1 or n < 5:
        raise ValueError("test derivative expects a one-dimensional profile with n>=5")
    result = np.empty_like(values)
    target = np.zeros(5, dtype=np.float64)
    target[1] = 1.0
    for i in range(n):
        start = min(max(i - 2, 0), n - 5)
        indices = np.arange(start, start + 5)
        offsets = indices.astype(np.float64) - float(i)
        vandermonde = np.vstack([offsets**power for power in range(5)])
        weights = np.linalg.solve(vandermonde, target)
        result[i] = float(weights @ values[indices]) / dr
    return result


def _zero(n: int) -> np.ndarray:
    return np.zeros(n, dtype=np.float64)


def test_flat_static_isotropic_stress_has_zero_covariant_divergence():
    n = 128
    dr = 0.25
    r = (np.arange(n, dtype=float) + 0.5) * dr
    alpha = np.ones(n)
    gamma_rr = np.ones(n)
    gamma_thth = r**2
    rho = np.full(n, 2.3)
    pressure = np.full(n, 0.7)
    j = _zero(n)
    zeros = _zero(n)
    div_t, div_r = AUDIT.covariant_divergence_from_profiles(
        r, _d1, alpha, gamma_rr, gamma_thth,
        rho, j, pressure, pressure,
        zeros, zeros, zeros, zeros, zeros,
    )
    assert float(np.max(np.abs(div_t))) < 2.0e-12
    assert float(np.max(np.abs(div_r))) < 2.0e-11


def test_flat_static_anisotropic_stress_matches_independent_radial_identity():
    n = 128
    dr = 0.25
    r = (np.arange(n, dtype=float) + 0.5) * dr
    alpha = np.ones(n)
    gamma_rr = np.ones(n)
    gamma_thth = r**2
    rho = np.full(n, 0.2)
    pr = 0.4 + 0.03 * r**2
    pt = 0.4 + 0.02 * r**2
    j = _zero(n)
    zeros = _zero(n)
    div_t, div_r = AUDIT.covariant_divergence_from_profiles(
        r, _d1, alpha, gamma_rr, gamma_thth,
        rho, j, pr, pt,
        zeros, zeros, zeros, zeros, zeros,
    )
    # In flat static spherical coordinates: D_r = p_r' + 2(p_r-p_t)/r.
    expected = 0.06 * r + 2.0 * (pr - pt) / r
    assert np.allclose(div_t, 0.0, rtol=0.0, atol=2.0e-12)
    assert np.allclose(div_r, expected, rtol=0.0, atol=2.0e-11)


def test_finite_volume_spherical_face_flux_reconstruction_is_exact_for_quadratic_flux():
    n = 64
    dr = 0.25
    r = (np.arange(n, dtype=np.float64) + 0.5) * dr
    flux = 0.7 * r**2
    derivative = AUDIT.finite_volume_radial_flux_derivative(flux, dr)
    expected = 1.4 * r
    assert np.allclose(derivative, expected, rtol=0.0, atol=2.0e-12)


def test_flat_static_isotropic_stress_has_zero_finite_volume_covariant_divergence():
    n = 128
    dr = 0.25
    r = (np.arange(n, dtype=float) + 0.5) * dr
    alpha = np.ones(n)
    gamma_rr = np.ones(n)
    gamma_thth = r**2
    rho = np.full(n, 2.3)
    pressure = np.full(n, 0.7)
    j = _zero(n)
    zeros = _zero(n)
    div_t, div_r = AUDIT.covariant_divergence_fv_from_profiles(
        r, _d1, alpha, gamma_rr, gamma_thth,
        rho, j, pressure, pressure,
        zeros, zeros, zeros, zeros, zeros,
    )
    assert float(np.max(np.abs(div_t))) < 2.0e-12
    assert float(np.max(np.abs(div_r))) < 2.0e-11

def test_finite_volume_face_reconstruction_respects_regular_origin_parity():
    n = 32
    dr = 0.25
    r = (np.arange(n, dtype=np.float64) + 0.5) * dr

    # F_r = r^2 g(r^2), with quadratic g, is reconstructed exactly at r=dr.
    even_flux = 0.7 * r**2 * (1.0 + 0.03 * r**2 + 0.002 * r**4)
    even_faces = AUDIT.spherical_flux_face_values(even_flux, dr, origin_power=2)
    expected_even = 0.7 * dr**2 * (1.0 + 0.03 * dr**2 + 0.002 * dr**4)
    assert np.isclose(even_faces[1], expected_even, rtol=0.0, atol=2.0e-14)

    # F_t = r^3 g(r^2) uses the odd regular-origin power.
    odd_flux = 0.4 * r**3 * (1.0 + 0.02 * r**2 + 0.001 * r**4)
    odd_faces = AUDIT.spherical_flux_face_values(odd_flux, dr, origin_power=3)
    expected_odd = 0.4 * dr**3 * (1.0 + 0.02 * dr**2 + 0.001 * dr**4)
    assert np.isclose(odd_faces[1], expected_odd, rtol=0.0, atol=2.0e-14)

    # The radial metric-source density scales as r at a regular origin.
    source_density = 0.6 * r * (1.0 + 0.015 * r**2 + 0.0005 * r**4)
    source_faces = AUDIT.spherical_flux_face_values(source_density, dr, origin_power=1)
    expected_source = 0.6 * dr * (1.0 + 0.015 * dr**2 + 0.0005 * dr**4)
    assert np.isclose(source_faces[1], expected_source, rtol=0.0, atol=2.0e-14)

def test_finite_volume_flat_anisotropic_stress_matches_volume_averaged_identity():
    n = 128
    dr = 0.25
    r = (np.arange(n, dtype=np.float64) + 0.5) * dr
    alpha = np.ones(n)
    gamma_rr = np.ones(n)
    gamma_thth = r**2
    rho = np.full(n, 0.2)
    pr = 0.4 + 0.03 * r**2
    pt = 0.4 + 0.02 * r**2
    j = _zero(n)
    zeros = _zero(n)
    div_t, div_r = AUDIT.covariant_divergence_fv_from_profiles(
        r, _d1, alpha, gamma_rr, gamma_thth,
        rho, j, pr, pt,
        zeros, zeros, zeros, zeros, zeros,
    )

    # Pointwise flat spherical identity is p_r' + 2(p_r-p_t)/r = 0.08 r.
    # The FV operator returns its exact spherical-volume average per cell.
    left = np.arange(n, dtype=np.float64) * dr
    right = left + dr
    expected = 0.06 * (right**4 - left**4) / (right**3 - left**3)
    assert np.allclose(div_t, 0.0, rtol=0.0, atol=2.0e-12)
    assert np.allclose(div_r, expected, rtol=0.0, atol=2.0e-11)

def test_finite_volume_curved_regular_flat_chart_has_expected_refinement():
    # Flat Minkowski space in a nonlinear areal-radius coordinate R(r)=r(1+c r^2).
    # The spatial metric varies, but a static isotropic stress has exact div(T)=0.
    # Cellwise FV residuals have O(dr^3) origin behavior and O(dr^4) annulus
    # behavior because the per-cell source-density integral uses Simpson quadrature.
    def d1_seven(values, _parity):
        values = np.asarray(values, dtype=np.float64)
        n = values.size
        out = np.empty_like(values)
        target = np.zeros(7, dtype=np.float64)
        target[1] = 1.0
        for i in range(n):
            start = min(max(i - 3, 0), n - 7)
            indices = np.arange(start, start + 7, dtype=np.float64)
            offsets = indices - float(i)
            vandermonde = np.vstack([offsets**power for power in range(7)])
            weights = np.linalg.solve(vandermonde, target)
            out[i] = (weights @ values[start:start + 7]) / dr
        return out

    center_errors = []
    annulus_errors = []
    for n in (40, 80, 160):
        dr = 8.0 / n
        r = (np.arange(n, dtype=np.float64) + 0.5) * dr
        c = 1.0e-3
        areal_radius = r * (1.0 + c * r**2)
        d_areal_dr = 1.0 + 3.0 * c * r**2
        alpha = np.ones(n)
        gamma_rr = d_areal_dr**2
        gamma_thth = areal_radius**2
        rho = np.full(n, 2.3)
        pressure = np.full(n, 0.7)
        j = _zero(n)
        zeros = _zero(n)
        div_t, div_r = AUDIT.covariant_divergence_fv_from_profiles(
            r, d1_seven, alpha, gamma_rr, gamma_thth,
            rho, j, pressure, pressure,
            zeros, zeros, zeros, zeros, zeros,
        )
        assert float(np.max(np.abs(div_t))) < 2.0e-12
        center_errors.append(float(abs(div_r[0])))
        annulus = (r >= 2.0) & (r <= 6.0)
        annulus_errors.append(float(np.sqrt(np.mean(div_r[annulus]**2))))

    # Halving dr should reduce the first-cell Simpson-origin error by ~8 and
    # the fixed-physical-annulus error by ~16; do not fit a production order.
    for coarse, fine in zip(center_errors[:-1], center_errors[1:]):
        assert 7.5 < coarse / fine < 8.5
    for coarse, fine in zip(annulus_errors[:-1], annulus_errors[1:]):
        assert 14.0 < coarse / fine < 18.0



def test_scalar_rhs_term_budget_reconstructs_pinned_scalar_rhs_and_labels_conditional_phi_source():
    """Lock RHS term bookkeeping and expose the implemented phi source separately."""
    kernel = AUDIT.V55ProductionKernel()
    state = kernel.initialize(resolution=40, r_max=24.0, include_radiation=True)
    before = AUDIT.state_digest(state)

    uncoupled = AUDIT.scalar_rhs_term_budget(state, 0.0)
    coupled = AUDIT.scalar_rhs_term_budget(state, -0.04)

    for budget in (uncoupled, coupled):
        assert budget["reported"] is True
        assert budget["all_profiles_finite"] is True
        assert budget["reconstruction_passed"] is True
        assert budget["conditional_exchange_model_law_accepted"] is False
        assert set(budget["sectors"]) == {"S", "D", "phi"}
        for sector in budget["sectors"].values():
            assert sector["reconstruction"]["field_rate_passed"] is True
            assert sector["reconstruction"]["momentum_rate_passed"] is True
            np.testing.assert_allclose(
                sector["reconstructed_rates"]["field_rate"],
                sector["implementation_rates"]["field_rate"],
                rtol=5.0e-12, atol=5.0e-18,
            )
            np.testing.assert_allclose(
                sector["reconstructed_rates"]["momentum_rate"],
                sector["implementation_rates"]["momentum_rate"],
                rtol=5.0e-12, atol=5.0e-18,
            )

    rho_dm = AUDIT.dm_density(
        AUDIT.metric_slice_from_q(state.grid, state.geometry), state.matter
    )
    phi_coupling = np.asarray(
        coupled["sectors"]["phi"]["momentum_rate_terms"]["raw_signed_arrays"][
            "conditional_dm_exchange"
        ],
        dtype=np.float64,
    )
    np.testing.assert_allclose(
        phi_coupling, -0.04 * np.asarray(state.geometry.alpha, dtype=np.float64) * np.asarray(rho_dm, dtype=np.float64),
        rtol=5.0e-12, atol=5.0e-18,
    )

    increment = AUDIT.scalar_rhs_term_coupling_increment(
        uncoupled, coupled, state.grid.centers
    )
    assert increment["reported"] is True
    assert increment["same_base_state"] is True
    assert increment["conditional_exchange_model_law_accepted"] is False
    assert increment["phi_exchange_increment_matches_declared_rhs_term"] is True
    assert increment["only_declared_phi_exchange_term_changed"] is True
    assert AUDIT.state_digest(state) == before


def test_scalar_sector_time_density_budget_splits_rhs_and_geometry_terms_without_admission():
    """Report termwise scalar d_t(sqrt(-g) T^t_t) on the same frozen state pair."""
    kernel = AUDIT.V55ProductionKernel()
    state = kernel.initialize(resolution=40, r_max=24.0, include_radiation=True)
    before = AUDIT.state_digest(state)
    old_beta = AUDIT.production_kernel_module.BETA_DM
    try:
        AUDIT.production_kernel_module.BETA_DM = -0.04
        rates = AUDIT.rhs_direction(kernel, state)
    finally:
        AUDIT.production_kernel_module.BETA_DM = old_beta

    epsilon = 2.0e-5
    plus_state = AUDIT.shifted_state(kernel, state, rates, +epsilon)
    minus_state = AUDIT.shifted_state(kernel, state, rates, -epsilon)
    rhs_budget = AUDIT.scalar_rhs_term_budget(state, -0.04)
    report = AUDIT.scalar_sector_energy_time_density_budget(
        state, plus_state, minus_state, epsilon, rhs_budget
    )

    assert report["reported"] is True
    assert report["all_profiles_finite"] is True
    assert report["conditional_exchange_model_law_accepted"] is False
    assert report["status"] == "scalar_sector_time_density_rhs_ledger_only"
    assert set(report["sectors"]) == {"S", "D", "phi"}
    assert "without a closure tolerance" in report["interpretation"]

    for sector_name, sector in report["sectors"].items():
        assert sector["state_pairing"]["base_state_sha256"] == before
        assert sector["state_pairing"]["epsilon"] == epsilon
        assert set(sector["field_rate_terms"]) == set(
            rhs_budget["sectors"][sector_name]["field_rate_terms"]["raw_signed_arrays"]
        )
        assert set(sector["momentum_rate_terms"]) == set(
            rhs_budget["sectors"][sector_name]["momentum_rate_terms"]["raw_signed_arrays"]
        )
        for key in (
            "measured_density_t_dot",
            "predicted_density_t_dot",
            "measured_minus_predicted_residual",
        ):
            values = np.asarray(sector[key]["raw_signed_array"], dtype=np.float64)
            assert values.shape == state.grid.centers.shape
            assert np.all(np.isfinite(values))
        assert "density_t_dot_volume_factor_contribution" in sector["geometry_terms"]
        assert "density_t_dot_metric_prefactor_contribution" in sector["geometry_terms"]

    # Algebraic implementation accounting for the conditional phi Pi RHS term:
    # it is visible in the diagnostic ledger but remains unaccepted as a physical law.
    phi = report["sectors"]["phi"]
    observed = np.asarray(
        phi["conditional_phi_exchange"]["density_t_dot_contribution"]["raw_signed_array"],
        dtype=np.float64,
    )
    rho_dm = AUDIT.dm_density(
        AUDIT.metric_slice_from_q(state.grid, state.geometry), state.matter
    )
    r = np.asarray(state.grid.centers, dtype=np.float64)
    gamma_rr = np.asarray(state.geometry.a) / np.asarray(state.geometry.X)**2
    gamma_thth = np.asarray(state.geometry.b) * r**2 / np.asarray(state.geometry.X)**2
    sqrt_minus_g = np.asarray(state.geometry.alpha) * np.sqrt(gamma_rr) * gamma_thth
    expected = (
        -sqrt_minus_g * np.asarray(state.scalars.Pi) * (-0.04 * np.asarray(rho_dm))
        / float(AUDIT.KAPPA)
    )
    np.testing.assert_allclose(observed, expected, rtol=5.0e-12, atol=5.0e-18)
    assert AUDIT.state_digest(state) == before


def test_scalar_sector_divergence_operator_components_reconstruct_existing_observers():
    """Split existing scalar divergence values without changing the operator result."""
    kernel = AUDIT.V55ProductionKernel()
    state = kernel.initialize(resolution=40, r_max=24.0, include_radiation=True)
    before = AUDIT.state_digest(state)
    old_beta = AUDIT.production_kernel_module.BETA_DM
    try:
        AUDIT.production_kernel_module.BETA_DM = -0.04
        rates = AUDIT.rhs_direction(kernel, state)
    finally:
        AUDIT.production_kernel_module.BETA_DM = old_beta

    epsilon = 2.0e-5
    plus_state = AUDIT.shifted_state(kernel, state, rates, +epsilon)
    minus_state = AUDIT.shifted_state(kernel, state, rates, -epsilon)
    paired = AUDIT.sector_resolved_divergence_tangent_from_states(
        state, plus_state, minus_state, epsilon
    )
    sectors = paired["sector_decomposition"]["sectors"]
    assert set(sectors) >= {"S", "D", "phi"}

    for name in ("S", "D", "phi"):
        sector = sectors[name]
        ledger = sector["operator_component_decomposition"]
        assert ledger["reported"] is True
        assert ledger["all_profiles_finite"] is True
        assert ledger["component_reconstruction_passed"] is True
        assert "no full conservation" in ledger["interpretation"]
        for observer_name, total_suffix in (
            ("pointwise", ""),
            ("finite_volume", "_fv"),
        ):
            observer = ledger[observer_name]
            assert observer["reported"] is True
            assert observer["all_profiles_finite"] is True
            terms = observer["raw_signed_arrays"]
            for component_name in (
                "time_density_div_T_t",
                "radial_flux_div_T_t",
                "metric_source_div_T_t_gtt",
                "metric_source_div_T_t_gamma_rr",
                "metric_source_div_T_t_gamma_theta_theta",
                "component_sum_div_T_t",
                "component_sum_gap_div_T_t",
                "time_density_div_T_r",
                "radial_flux_div_T_r",
                "metric_source_div_T_r_gtt",
                "metric_source_div_T_r_gamma_rr",
                "metric_source_div_T_r_gamma_theta_theta",
                "component_sum_div_T_r",
                "component_sum_gap_div_T_r",
            ):
                values = np.asarray(terms[component_name], dtype=np.float64)
                assert values.shape == state.grid.centers.shape
                assert np.all(np.isfinite(values))
            for channel in ("div_T_t", "div_T_r"):
                expected = np.asarray(
                    sector["raw_signed_arrays"][channel + total_suffix],
                    dtype=np.float64,
                )
                reconstructed = np.asarray(
                    terms[f"component_sum_{channel}"], dtype=np.float64
                )
                gap = np.asarray(
                    terms[f"component_sum_gap_{channel}"], dtype=np.float64
                )
                np.testing.assert_allclose(
                    reconstructed, expected, rtol=5.0e-12, atol=5.0e-18
                )
                np.testing.assert_allclose(
                    gap, 0.0, rtol=5.0e-12, atol=5.0e-18
                )
    assert AUDIT.state_digest(state) == before
