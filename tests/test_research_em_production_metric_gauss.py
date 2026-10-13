"""Research-only tests on a frozen production-initialized spherical metric."""
from dataclasses import FrozenInstanceError, replace
import math

import numpy as np
import pytest

from engine.production_kernel import V55ProductionKernel
from tools.research_em_radial_reference import RadialCoulombReference
from tools.research_em_production_metric_gauss import (
    FrozenProductionMetric,
    solve_discrete_production_metric_gauss,
)


@pytest.fixture(scope="module")
def production_metric():
    # This is the real current V5.5 initializer and CMC gauge, but no step is
    # advanced and no electromagnetic source is inserted into the state.
    state = V55ProductionKernel().initialize(
        resolution=128, r_max=40.0, include_radiation=True
    )
    frozen = FrozenProductionMetric.from_production_state(state)
    # Prevent accidentally substituting the manufactured conformal metric
    # used by the separate analytic-oracle tests.
    assert len(frozen.centers) == 128
    assert max(
        max(abs(value - 1.0) for value in component)
        for component in (frozen.a, frozen.b, frozen.X)
    ) > 1.0e-12
    return frozen


@pytest.fixture
def source():
    # A manufactured diagnostic source only; not a charge assignment to matter.
    return RadialCoulombReference(rho0=0.75, source_radius=10.0)


@pytest.mark.parametrize("resolution", [32, 64, 128])
def test_production_metric_gauss_shell_closure_and_flux_identity(
    production_metric, source, resolution
):
    result = solve_discrete_production_metric_gauss(
        production_metric, source, resolution=resolution
    )
    scale = max(1.0, abs(result.total_charge))
    assert len(result.faces) == resolution + 1
    assert len(result.cell_charges) == resolution
    assert result.face_charges[0] == 0.0
    assert result.face_electric_field_hat[0] == 0.0
    assert result.face_electric_field_r[0] == 0.0
    assert result.max_abs_cell_residual <= 5.0e-14 * scale

    for i, radius in enumerate(result.faces):
        if radius == 0.0:
            continue
        flux_coordinate = (
            production_metric.sqrt_gamma(radius)
            * result.face_electric_field_r[i]
        )
        flux_areal = (
            production_metric.areal_radius(radius) ** 2
            * result.face_electric_field_hat[i]
        )
        assert flux_coordinate == pytest.approx(
            result.face_charges[i], rel=2.0e-13, abs=2.0e-13 * scale
        )
        assert flux_areal == pytest.approx(
            result.face_charges[i], rel=2.0e-13, abs=2.0e-13 * scale
        )


def test_production_metric_gauss_source_budget_matches_integral_on_one_frozen_geometry(
    production_metric, source
):
    # All three target grids use the exact same frozen production metric.
    # Cell charge integrates rho_e*sqrt(gamma) rather than midpoint-sampling
    # rho_e, and is checked against a separate higher-order integral.
    solutions = [
        solve_discrete_production_metric_gauss(
            production_metric, source, resolution=n
        )
        for n in (32, 64, 128)
    ]
    for result in solutions:
        reference_charge = tuple(
            production_metric.reference_enclosed_charge(source, radius)
            for radius in result.faces
        )
        scale = max(1.0, abs(reference_charge[-1]))
        tolerance = 5.0e-11 * scale
        assert all(math.isfinite(value) for value in result.cell_charges)
        assert max(
            abs(numerical - exact)
            for numerical, exact in zip(result.face_charges, reference_charge)
        ) <= tolerance
        assert result.total_charge == pytest.approx(
            reference_charge[-1], rel=5.0e-11, abs=tolerance
        )
        for i, radius in enumerate(result.faces):
            if radius == 0.0:
                continue
            reference_field = (
                reference_charge[i]
                / production_metric.areal_radius(radius) ** 2
            )
            assert result.face_electric_field_hat[i] == pytest.approx(
                reference_field,
                rel=5.0e-11,
                abs=5.0e-11 * max(1.0, abs(reference_field)),
            )


def test_production_metric_gauss_has_a_regular_centre_and_detects_field_error(
    production_metric, source
):
    result = solve_discrete_production_metric_gauss(
        production_metric, source, resolution=64
    )
    assert all(math.isfinite(value) for value in result.face_electric_field_hat)
    assert all(math.isfinite(value) for value in result.face_electric_field_r)
    assert result.face_electric_field_hat[1] > 0.0
    assert math.isfinite(
        result.face_electric_field_hat[1] / result.faces[1]
    )

    # Perturb a face field deliberately; the signed finite-volume Gauss
    # residual must detect it rather than merely checking that the solve ran.
    cell = 2
    left = result.faces[cell]
    right = result.faces[cell + 1]
    perturbed = result.face_electric_field_r[cell + 1] + 1.0e-5
    residual = (
        production_metric.sqrt_gamma(right) * perturbed
        - production_metric.sqrt_gamma(left)
        * result.face_electric_field_r[cell]
        - result.cell_charges[cell]
    )
    assert abs(residual) > 1.0e-6


def test_production_metric_gauss_respects_polarity_and_neutral_limit(
    production_metric
):
    positive = RadialCoulombReference(rho0=0.75, source_radius=10.0)
    negative = RadialCoulombReference(rho0=-0.75, source_radius=10.0)
    neutral = RadialCoulombReference(rho0=0.0, source_radius=10.0)
    plus = solve_discrete_production_metric_gauss(
        production_metric, positive, resolution=64
    )
    minus = solve_discrete_production_metric_gauss(
        production_metric, negative, resolution=64
    )
    zero = solve_discrete_production_metric_gauss(
        production_metric, neutral, resolution=64
    )

    for q_plus, q_minus in zip(plus.face_charges, minus.face_charges):
        assert q_minus == pytest.approx(-q_plus, rel=2.0e-14, abs=2.0e-14)
    for e_plus, e_minus in zip(
        plus.face_electric_field_r, minus.face_electric_field_r
    ):
        assert e_minus == pytest.approx(-e_plus, rel=2.0e-14, abs=2.0e-14)
    assert all(value == 0.0 for value in zero.cell_charges)
    assert all(value == 0.0 for value in zero.face_charges)
    assert all(value == 0.0 for value in zero.face_electric_field_hat)
    assert all(value == 0.0 for value in zero.face_electric_field_r)


def test_production_metric_gauss_fails_closed_on_bad_geometry_and_source_edge(
    production_metric, source
):
    bad_x = list(production_metric.X)
    bad_x[4] = 0.0
    with pytest.raises(ValueError, match="positive"):
        replace(production_metric, X=tuple(bad_x))

    bad_a = list(production_metric.a)
    bad_a[7] = math.nan
    with pytest.raises(ValueError, match="finite"):
        replace(production_metric, a=tuple(bad_a))

    with pytest.raises(ValueError, match="align"):
        solve_discrete_production_metric_gauss(
            production_metric,
            RadialCoulombReference(rho0=0.75, source_radius=9.5),
            resolution=32,
        )
    with pytest.raises(ValueError, match="integer"):
        solve_discrete_production_metric_gauss(
            production_metric, source, resolution=True
        )


def test_production_metric_snapshot_is_immutable(production_metric):
    assert isinstance(production_metric.a, tuple)
    with pytest.raises(FrozenInstanceError):
        production_metric.a = production_metric.b



def test_tiny_production_metric_rejects_scale_mismatched_source_edge():
    """The face-alignment tolerance must scale with tiny coordinate radii."""
    r_max = 1.0e-18
    metric_resolution = 64
    dr = r_max / metric_resolution
    centers = tuple((i + 0.5) * dr for i in range(metric_resolution))
    unit = (1.0,) * metric_resolution
    metric = FrozenProductionMetric(
        centers=centers,
        a=unit,
        b=unit,
        X=unit,
        dr=dr,
        r_max=r_max,
    )
    source = RadialCoulombReference(rho0=0.75, source_radius=1.0e-20)
    # At target N=67, source radius is 0.67 cells from the origin; the old
    # ulp(max(1, radius)) tolerance accepted a face mismatch ~5e-21.
    with pytest.raises(ValueError, match="align"):
        solve_discrete_production_metric_gauss(metric, source, resolution=67)


def test_tiny_production_metric_snapshot_rejects_misplaced_center():
    r_max = 1.0e-18
    resolution = 64
    dr = r_max / resolution
    centers = [(i + 0.5) * dr for i in range(resolution)]
    centers[5] += 1.0e-20
    unit = (1.0,) * resolution
    with pytest.raises(ValueError, match="uniform"):
        FrozenProductionMetric(
            centers=tuple(centers),
            a=unit,
            b=unit,
            X=unit,
            dr=dr,
            r_max=r_max,
        )


def test_production_metric_gauss_rejects_overflowed_areal_radius_square():
    """Finite positive geometry must not turn an overflowed R**2 into E_hat=0."""
    resolution = 4
    dr = 10.0
    metric = FrozenProductionMetric(
        centers=tuple((i + 0.5) * dr for i in range(resolution)),
        a=(1.0e-320,) * resolution,
        b=(1.0e308,) * resolution,
        X=(1.0,) * resolution,
        dr=dr,
        r_max=40.0,
    )
    # Source edge is an exact face; the test isolates only the areal-square guard.
    source = RadialCoulombReference(rho0=0.75, source_radius=20.0)
    with pytest.raises(ValueError, match="face geometry must remain finite"):
        solve_discrete_production_metric_gauss(
            metric, source, resolution=resolution
        )

