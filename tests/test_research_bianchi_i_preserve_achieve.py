"""Regression tests for the research-only Bianchi-I A6 numerical oracle."""
from __future__ import annotations

import math

import pytest

from tools.audit_bianchi_i_preserve_achieve import (
    run_benchmark_suite,
    run_unsourced_case,
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


def test_anisotropic_shear_source_changes_marginal_preserve_case() -> None:
    sigma2, rho = 0.2, 0.8
    H = -math.sqrt((sigma2 + rho) / 3.0)
    unsourced = slope_with_sources_dln_a(
        sigma2=sigma2, rho=rho, w=1.0, H=H
    )
    sourced = slope_with_sources_dln_a(
        sigma2=sigma2,
        rho=rho,
        w=1.0,
        H=H,
        sigma2_source=0.03,
    )
    assert unsourced == 0.0
    assert math.isfinite(sourced)
    assert sourced != 0.0


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
    assert math.isfinite(exchanged)
    assert exchanged != 0.0


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
