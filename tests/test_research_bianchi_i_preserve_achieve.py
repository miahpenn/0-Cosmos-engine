"""Regression tests for the research-only Bianchi-I A6 numerical oracle."""
from __future__ import annotations

import math

import pytest

from tools.audit_bianchi_i_preserve_achieve import (
    run_benchmark_suite,
    run_unsourced_case,
    gear45_shear_norm_source,
    slope_domega_dcontraction_efold,
    slope_domega_dln_a,
    slope_with_sources_dln_a,
)


@pytest.mark.parametrize("w", [-1.0 / 3.0, 1.0, 2.0])
@pytest.mark.parametrize("omega0", [0.1, 0.5, 0.9])
def test_unsourced_integrator_matches_exact_scaling_and_branch_trend(w: float, omega0: float) -> None:
    result = run_unsourced_case(w=w, omega_sigma_initial=omega0)

    assert result["solver_success"] is True
    assert result["max_component_relative_error"] < 2.0e-9
    assert result["max_omega_abs_error"] < 2.0e-9
    assert result["max_interior_finite_difference_slope_abs_error"] < 1.0e-4

    if w < 1.0:
        assert result["expected_contraction_trend"] == "increases_during_contraction"
        assert result["omega_sigma_delta"] > 0.0
    elif w == 1.0:
        assert result["expected_contraction_trend"] == "preserved"
        assert abs(result["omega_sigma_delta"]) < 2.0e-9
    else:
        assert result["expected_contraction_trend"] == "decreases_during_contraction"
        assert result["omega_sigma_delta"] < 0.0


@pytest.mark.parametrize("w", [-1.0 / 3.0, 1.0, 2.0])
@pytest.mark.parametrize("omega", [0.0, 0.2, 0.5, 0.9, 1.0])
def test_ln_a_and_contraction_efold_slopes_have_opposite_signs(w: float, omega: float) -> None:
    assert math.isclose(
        slope_domega_dcontraction_efold(omega, w),
        -slope_domega_dln_a(omega, w),
        rel_tol=0.0,
        abs_tol=1.0e-14,
    )


def test_free_identity_is_zero_at_marginal_stiff_case_and_zero_shear() -> None:
    assert slope_domega_dln_a(0.2, 1.0) == 0.0
    assert slope_domega_dln_a(0.0, 2.0) == 0.0
    assert slope_domega_dln_a(1.0, 2.0) == 0.0


def test_source_aware_formula_reduces_to_free_identity() -> None:
    sigma2, rho = 0.2, 0.8
    H = -math.sqrt((sigma2 + rho) / 3.0)
    got = slope_with_sources_dln_a(
        sigma2=sigma2,
        rho=rho,
        w=2.0,
        H=H,
    )
    assert math.isclose(got, slope_domega_dln_a(0.2, 2.0), rel_tol=0.0, abs_tol=1.0e-14)


def test_gear45_source_and_opposite_exchange_derive_a6_slope() -> None:
    # GEAR 45 supplies this gamma_i row (gamma_star=17) and isotropic gamma=5.
    # It does not pin a unique numerical mu/sigma_i initial state, so this is an
    # explicit constructed equation-level fixture, not an archived trajectory.
    gamma_i = (12.0, 18.0, 21.0)
    gamma_star = sum(gamma_i) / 3.0
    anisotropic_pressures = tuple(value - gamma_star for value in gamma_i)
    sigma2 = 0.2
    sigma_scale = math.sqrt(2.0 * sigma2 / sum(value**2 for value in anisotropic_pressures))
    sigma_i = tuple(sigma_scale * value for value in anisotropic_pressures)
    mu = 0.1
    rho_isotropic = 0.7
    rho_total = rho_isotropic + mu
    gamma_isotropic = 5.0
    w_eff = (
        (gamma_isotropic - 1.0) * rho_isotropic
        + (gamma_star - 1.0) * mu
    ) / rho_total
    H = -math.sqrt((sigma2 + rho_total) / 3.0)

    sigma2_source = gear45_shear_norm_source(
        mu=mu, sigma_i=sigma_i, gamma_i=gamma_i
    )
    expected_source = mu * sum(
        sigma * pressure for sigma, pressure in zip(sigma_i, anisotropic_pressures)
    )
    assert math.isclose(sum(sigma_i), 0.0, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(0.5 * sum(sigma * sigma for sigma in sigma_i), sigma2, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(sigma2_source, expected_source, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(sigma2_source, mu * math.sqrt(2.0 * sigma2 * sum(p * p for p in anisotropic_pressures)), rel_tol=0.0, abs_tol=1.0e-14)
    assert sigma2_source > 0.0

    got = slope_with_sources_dln_a(
        sigma2=sigma2,
        rho=rho_total,
        w=w_eff,
        H=H,
        sigma2_source=sigma2_source,
        rho_exchange=-sigma2_source,
    )

    # Independent quotient-rule check using the GEAR-45 component ODEs:
    # dot(sigma2) = -6 H sigma2 + S_sigma; dot(rho+mu) = -3 H(1+w_eff)(rho+mu) - S_sigma.
    omega_denominator = sigma2 + rho_total
    sigma2_dot = -6.0 * H * sigma2 + expected_source
    rho_total_dot = -3.0 * H * (1.0 + w_eff) * rho_total - expected_source
    omega_dot = (
        sigma2_dot * rho_total - sigma2 * rho_total_dot
    ) / (omega_denominator * omega_denominator)
    expected_slope = omega_dot / H
    assert math.isclose(got, expected_slope, rel_tol=0.0, abs_tol=1.0e-14)

    wrong_sign = slope_with_sources_dln_a(
        sigma2=sigma2,
        rho=rho_total,
        w=w_eff,
        H=H,
        sigma2_source=-sigma2_source,
        rho_exchange=sigma2_source,
    )
    assert not math.isclose(got, wrong_sign, rel_tol=0.0, abs_tol=1.0e-8)


def test_energy_exchange_changes_marginal_preserve_case() -> None:
    sigma2, rho = 0.2, 0.8
    H = -math.sqrt((sigma2 + rho) / 3.0)
    unsourced = slope_with_sources_dln_a(
        sigma2=sigma2, rho=rho, w=1.0, H=H
    )
    exchanged = slope_with_sources_dln_a(
        sigma2=sigma2,
        rho=rho,
        w=1.0,
        H=H,
        rho_exchange=0.04,
    )
    assert unsourced == 0.0
    assert math.isclose(exchanged, 0.013856406460551019, rel_tol=0.0, abs_tol=1.0e-14)


@pytest.mark.parametrize("H", [0.0, -1.0e-8])
def test_slope_gate_rejects_turnaround_or_excluded_hubble_band(H: float) -> None:
    kwargs = {
        "sigma2": 0.2,
        "rho": 0.8,
        "w": 1.0,
        "H": H,
        "min_abs_H": 1.0e-6,
    }
    with pytest.raises(ValueError, match="H is zero or inside"):
        slope_with_sources_dln_a(**kwargs)


def test_slope_gate_rejects_constraint_inconsistent_sample() -> None:
    with pytest.raises(ValueError, match="constraint"):
        slope_with_sources_dln_a(
            sigma2=0.2,
            rho=0.8,
            w=1.0,
            H=-1.0,
        )


def test_slope_gate_is_scale_invariant_for_small_densities() -> None:
    sigma2, rho = 2.0e-20, 8.0e-20
    consistent_H = -math.sqrt((sigma2 + rho) / 3.0)
    got = slope_with_sources_dln_a(
        sigma2=sigma2,
        rho=rho,
        w=0.5,
        H=consistent_H,
    )
    omega = sigma2 / (sigma2 + rho)
    expected = 3.0 * (0.5 - 1.0) * omega * (1.0 - omega)
    assert math.isclose(got, expected, rel_tol=0.0, abs_tol=1.0e-14)

    # The same small density inputs with an inconsistent Hubble value must fail.
    with pytest.raises(ValueError, match="constraint"):
        slope_with_sources_dln_a(
            sigma2=sigma2,
            rho=rho,
            w=0.5,
            H=-1.0e-6,
        )


def test_slope_gate_rejects_nonfinite_inputs_and_bad_component_domains() -> None:
    with pytest.raises(ValueError, match="finite"):
        slope_domega_dln_a(float("nan"), 1.0)
    with pytest.raises(ValueError, match="strictly positive"):
        slope_with_sources_dln_a(
            sigma2=0.0,
            rho=1.0,
            w=1.0,
            H=-math.sqrt(1.0 / 3.0),
        )
    with pytest.raises(ValueError, match="strictly between"):
        run_unsourced_case(w=1.0, omega_sigma_initial=1.0)


def test_deterministic_benchmark_report_is_bounded_and_admission_gated() -> None:
    report = run_benchmark_suite()

    assert report["case_count"] == 9
    assert report["all_checks_passed"] is True
    assert report["scientific_admission"] == "NOT ADMITTED"
    assert "no production-engine calls" in report["model_scope"]

def test_gear45_source_changes_sign_when_shear_orientation_is_reversed() -> None:
    gamma_i = (12.0, 18.0, 21.0)
    sigma_i = (-0.5, 0.1, 0.4)  # trace-free, deliberately not tied to archive initial data
    source = gear45_shear_norm_source(mu=0.1, sigma_i=sigma_i, gamma_i=gamma_i)
    reversed_source = gear45_shear_norm_source(
        mu=0.1,
        sigma_i=tuple(-value for value in sigma_i),
        gamma_i=gamma_i,
    )

    assert source > 0.0
    assert math.isclose(reversed_source, -source, rel_tol=0.0, abs_tol=1.0e-14)


def test_gear45_source_vanishes_for_isotropic_directional_pressure() -> None:
    sigma_i = (-0.5, 0.1, 0.4)
    source = gear45_shear_norm_source(
        mu=0.1,
        sigma_i=sigma_i,
        gamma_i=(5.0, 5.0, 5.0),
    )

    assert source == 0.0

