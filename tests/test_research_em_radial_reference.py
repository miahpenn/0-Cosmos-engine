"""Contract tests for the research-only radial Coulomb reference oracle."""
import math

import pytest

from tools.research_em_radial_reference import RadialCoulombReference


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
        0.75 * (1.0 - 0.25**2) ** 2
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
