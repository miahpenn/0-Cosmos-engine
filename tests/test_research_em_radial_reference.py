"""Contract tests for the research-only radial Coulomb reference oracle."""
import math

import pytest

from engine.valencia import spherical_metric_from_bssn
from tools.research_em_radial_reference import (
    CurvedRadialCoulombReference,
    RadialCoulombReference,
)


@pytest.fixture
def reference():
    # Non-unit scales make accidental normalization omissions visible.
    return RadialCoulombReference(rho0=0.75, source_radius=2.0)


def test_total_charge_matches_closed_form(reference):
    expected = (32.0 * math.pi / 105.0) * 0.75 * 2.0**3
    assert reference.total_charge == pytest.approx(expected, rel=0.0, abs=1e-15)


def test_charge_density_is_regular_and_compactly_supported(reference):
    assert reference.charge_density(0.0) == pytest.approx(0.75)
    assert reference.charge_density(reference.source_radius) == 0.0
    assert reference.charge_density(1.5 * reference.source_radius) == 0.0
    assert reference.charge_density(0.5 * reference.source_radius) == pytest.approx(
        0.75 * (1.0 - 0.25) ** 2
    )


def test_central_field_has_the_required_linear_limit(reference):
    r = 1.0e-5 * reference.source_radius
    expected_slope = 4.0 * math.pi * reference.rho0 / 3.0
    assert reference.electric_field(0.0) == 0.0
    assert reference.electric_field(r) / r == pytest.approx(
        expected_slope, rel=3.0e-10
    )


def test_exterior_field_is_coulomb_and_enclosed_charge_is_constant(reference):
    for r in (1.1 * reference.source_radius, 2.0 * reference.source_radius, 8.0):
        assert reference.enclosed_charge(r) == pytest.approx(
            reference.total_charge, rel=0.0, abs=1e-15
        )
        assert r * r * reference.electric_field(r) == pytest.approx(
            reference.total_charge, rel=2.0e-15, abs=1e-15
        )


def test_differential_gauss_identity_from_independent_profile_formula(reference):
    # Compare a numerical derivative of Q against the separately evaluated
    # profile 4*pi*rho_e*r**2, away from the compact-support boundary.
    r = 0.37 * reference.source_radius
    h = 1.0e-5 * reference.source_radius
    dQ_dr = (
        reference.enclosed_charge(r + h) - reference.enclosed_charge(r - h)
    ) / (2.0 * h)
    source = 4.0 * math.pi * reference.charge_density(r) * r * r
    assert dQ_dr == pytest.approx(source, rel=2.0e-9, abs=1.0e-12)


def test_pure_radial_electric_stress_projection_signs(reference):
    projection = reference.stress_energy(0.61 * reference.source_radius)
    assert projection.energy_density > 0.0
    assert projection.radial_pressure == -projection.energy_density
    assert projection.tangential_pressure == projection.energy_density
    assert projection.radial_momentum_density == 0.0
    assert (
        projection.radial_pressure
        + 2.0 * projection.tangential_pressure
        - projection.energy_density
    ) == pytest.approx(0.0, abs=1.0e-15)


def test_finite_volume_gauss_residual_vanishes_for_exact_field(reference):
    faces = (
        0.0,
        0.25 * reference.source_radius,
        0.5 * reference.source_radius,
        reference.source_radius,
        1.5 * reference.source_radius,
        3.0 * reference.source_radius,
    )
    tolerance = 5.0e-14 * max(1.0, abs(reference.total_charge))
    residuals = []
    for inner, outer in zip(faces, faces[1:]):
        residual = reference.gauss_shell_residual(
            inner,
            outer,
            reference.electric_field(inner),
            reference.electric_field(outer),
        )
        residuals.append(residual)
    assert max(abs(value) for value in residuals) <= tolerance


def test_finite_volume_gauss_residual_detects_perturbed_field(reference):
    inner = 0.25 * reference.source_radius
    outer = 0.5 * reference.source_radius
    perturbed_outer_field = reference.electric_field(outer) + 1.0e-5
    residual = reference.gauss_shell_residual(
        inner,
        outer,
        reference.electric_field(inner),
        perturbed_outer_field,
    )
    assert abs(residual) > 1.0e-6


@pytest.mark.parametrize(
    "rho0, source_radius",
    [
        (math.nan, 1.0),
        (1.0, math.inf),
        (1.0, 0.0),
        (True, 1.0),
    ],
)
def test_invalid_fixture_parameters_fail_closed(rho0, source_radius):
    with pytest.raises(ValueError):
        RadialCoulombReference(rho0=rho0, source_radius=source_radius)


@pytest.mark.parametrize("radius", [-1.0, math.nan, math.inf, True])
def test_invalid_radius_fails_closed(reference, radius):
    with pytest.raises(ValueError):
        reference.electric_field(radius)


def test_invalid_shell_order_and_nonfinite_field_fail_closed(reference):
    with pytest.raises(ValueError, match="r_outer"):
        reference.gauss_shell_residual(1.0, 0.5, 0.0, 0.0)
    with pytest.raises(ValueError, match="E_outer"):
        reference.gauss_shell_residual(0.0, 1.0, 0.0, math.nan)



@pytest.fixture
def curved_reference():
    # psi=1+eta*(r/R0)^2 gives a regular but nontrivial frozen metric.
    return CurvedRadialCoulombReference(
        rho0=0.75, source_radius=2.0, metric_strength=0.4
    )


def test_curved_metric_factors_match_production_bssn_mapping(curved_reference):
    for r in (0.2, 0.74, 1.4, 3.0):
        psi = curved_reference.conformal_factor(r)
        metric = spherical_metric_from_bssn(
            r=r, a=1.0, b=1.0, X=1.0 / psi, alpha=1.0, beta=0.0
        )
        assert metric.gamma_rr == pytest.approx(psi**2, rel=2.0e-15)
        assert metric.gamma_thth == pytest.approx(
            curved_reference.areal_radius(r) ** 2, rel=2.0e-15
        )
        assert metric.sqrt_gamma == pytest.approx(
            curved_reference.sqrt_gamma(r), rel=2.0e-15
        )


def test_curved_gauss_flux_uses_areal_radius_and_proper_volume(curved_reference):
    for r in (0.2, 0.6, 1.4, 2.8):
        q = curved_reference.enclosed_charge(r)
        field_hat = curved_reference.electric_field(r)
        field_coordinate = curved_reference.coordinate_electric_field(r)
        assert curved_reference.areal_radius(r) ** 2 * field_hat == pytest.approx(
            q, rel=2.0e-15, abs=1.0e-14
        )
        assert curved_reference.sqrt_gamma(r) * field_coordinate == pytest.approx(
            q, rel=2.0e-15, abs=1.0e-14
        )


def test_curved_differential_gauss_identity_includes_metric_volume(curved_reference):
    r = 0.37 * curved_reference.source_radius
    h = 1.0e-5 * curved_reference.source_radius
    dQ_dr = (
        curved_reference.enclosed_charge(r + h)
        - curved_reference.enclosed_charge(r - h)
    ) / (2.0 * h)
    proper_volume_source = (
        4.0
        * math.pi
        * curved_reference.charge_density(r)
        * curved_reference.sqrt_gamma(r)
    )
    assert dQ_dr == pytest.approx(
        proper_volume_source, rel=3.0e-9, abs=1.0e-12
    )


def test_curved_field_is_regular_at_center_and_coulomb_in_areal_radius(
    curved_reference,
):
    r = 1.0e-5 * curved_reference.source_radius
    expected_slope = 4.0 * math.pi * curved_reference.rho0 / 3.0
    assert curved_reference.electric_field(0.0) == 0.0
    assert curved_reference.electric_field(r) / r == pytest.approx(
        expected_slope, rel=3.0e-10
    )
    for r_outer in (1.2 * curved_reference.source_radius, 3.0 * curved_reference.source_radius):
        assert (
            curved_reference.areal_radius(r_outer) ** 2
            * curved_reference.electric_field(r_outer)
        ) == pytest.approx(
            curved_reference.total_charge, rel=2.0e-15, abs=1.0e-14
        )


def test_curved_finite_volume_gauss_residual_and_negative_control(
    curved_reference,
):
    faces = (
        0.0,
        0.25 * curved_reference.source_radius,
        0.5 * curved_reference.source_radius,
        curved_reference.source_radius,
        1.5 * curved_reference.source_radius,
        3.0 * curved_reference.source_radius,
    )
    tolerance = 8.0e-14 * max(1.0, abs(curved_reference.total_charge))
    for inner, outer in zip(faces, faces[1:]):
        residual = curved_reference.gauss_shell_residual(
            inner,
            outer,
            curved_reference.coordinate_electric_field(inner),
            curved_reference.coordinate_electric_field(outer),
        )
        assert abs(residual) <= tolerance

    inner = 0.25 * curved_reference.source_radius
    outer = 0.5 * curved_reference.source_radius
    residual = curved_reference.gauss_shell_residual(
        inner,
        outer,
        curved_reference.coordinate_electric_field(inner),
        curved_reference.coordinate_electric_field(outer) + 1.0e-5,
    )
    assert abs(residual) > 1.0e-6


@pytest.mark.parametrize("metric_strength", [0.0, -0.1, math.nan, math.inf])
def test_curved_reference_rejects_nonpositive_or_nonfinite_metric_strength(
    metric_strength,
):
    with pytest.raises(ValueError):
        CurvedRadialCoulombReference(
            rho0=1.0, source_radius=1.0, metric_strength=metric_strength
        )
