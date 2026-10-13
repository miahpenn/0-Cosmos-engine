"""Boundary contract between radiation admissibility and primitive inversion."""
from __future__ import annotations

import numpy as np
import pytest

from engine.matter_system import (
    BSSNMetricSlice,
    ConservedSpecies,
    RadiationAdmissibilityError,
    Species,
    _metric_arrays,
    primitives,
    validate_species_admissibility,
)
from engine.valencia import (
    FluidConserved,
    RadiationAdmissibilityError as RecoveryAdmissibilityError,
    SphericalMetric,
    recover_radiation,
)


def _flat_bssn_slice() -> BSSNMetricSlice:
    one = np.array([1.0])
    zero = np.array([0.0])
    return BSSNMetricSlice(
        r=one.copy(),
        a=one.copy(),
        b=one.copy(),
        X=one.copy(),
        alpha=one.copy(),
        beta=zero.copy(),
        Aa=zero.copy(),
        K=zero.copy(),
        Lambda=zero.copy(),
        B=zero.copy(),
    )


def _flat_spherical_metric() -> SphericalMetric:
    return SphericalMetric(
        alpha=1.0,
        beta=0.0,
        gamma_rr=1.0,
        gamma_rr_inv=1.0,
        gamma_thth=1.0,
        gamma_thth_inv=1.0,
        sqrt_gamma=1.0,
    )


def _conserved(E: float, S: float) -> FluidConserved:
    return FluidConserved(rest=0.0, energy_t=E, momentum_r=S)


def _state(E: float, S: float) -> ConservedSpecies:
    return ConservedSpecies(
        rest=np.array([0.0]),
        energy_t=np.array([E]),
        momentum_r=np.array([S]),
    )


def test_validator_rejects_exact_E_equals_absS_with_nonzero_flux():
    metric = _flat_bssn_slice()
    with pytest.raises(
        RadiationAdmissibilityError,
        match=r"stage='exact-boundary'.*ratio=1.00000000000000000e\+00",
    ):
        validate_species_admissibility(
            metric, _state(E=1.0, S=1.0), Species.RADIATION,
            stage="exact-boundary",
        )


def test_recovery_and_vector_primitive_path_preserve_typed_boundary_failure():
    metric = _flat_spherical_metric()
    state = _conserved(E=1.0, S=1.0)
    assert RecoveryAdmissibilityError is RadiationAdmissibilityError

    with pytest.raises(RadiationAdmissibilityError, match="strict interior"):
        recover_radiation(metric, state)

    # The vectorized wrapper adds cell/neighbour diagnostics but must preserve
    # the retryable error type rather than downgrade it to generic ValueError.
    with pytest.raises(RadiationAdmissibilityError, match="strict interior"):
        primitives(
            _metric_arrays(_flat_bssn_slice()),
            _state(E=1.0, S=1.0),
            Species.RADIATION,
        )


def test_strict_interior_E_greater_than_absS_is_admitted_and_invertible():
    validate_species_admissibility(
        _flat_bssn_slice(), _state(E=1.0, S=0.5), Species.RADIATION,
        stage="strict-interior",
    )
    primitive = recover_radiation(_flat_spherical_metric(), _conserved(E=1.0, S=0.5))
    assert np.isfinite(primitive.rho)
    assert np.isfinite(primitive.pressure)
    assert np.isfinite(primitive.v_r)
    assert primitive.rho > 0.0
    assert primitive.pressure > 0.0
    assert abs(primitive.v_r) < 1.0


@pytest.mark.parametrize("energy", [0.0, 1.0])
def test_zero_momentum_branch_remains_valid_including_vacuum(energy):
    validate_species_admissibility(
        _flat_bssn_slice(), _state(E=energy, S=0.0), Species.RADIATION,
        stage="zero-momentum",
    )
    primitive = recover_radiation(
        _flat_spherical_metric(), _conserved(E=energy, S=0.0)
    )
    assert primitive.rho == energy
    assert primitive.pressure == pytest.approx(energy / 3.0)
    assert primitive.v_r == 0.0

@pytest.mark.parametrize("ratio", [0.5, 0.9999, 0.999999])
@pytest.mark.parametrize("energy", [3.4091116247062796e-8, 1.0e-10])
def test_radiation_pressure_recovery_is_invariant_under_density_scaling(ratio, energy):
    """Scaling E and |S| together must not change the dimensionless pressure."""
    metric = _flat_spherical_metric()
    unit_scale = recover_radiation(metric, _conserved(E=1.0, S=ratio))
    scaled = recover_radiation(
        metric, _conserved(E=energy, S=energy * ratio)
    )

    assert scaled.pressure / energy == pytest.approx(
        unit_scale.pressure, rel=2.0e-10, abs=1.0e-14
    )
    assert scaled.rho / energy == pytest.approx(
        unit_scale.rho, rel=2.0e-10, abs=1.0e-14
    )


@pytest.mark.parametrize("ratio", [0.5, 0.9999, 0.999999])
def test_low_energy_recovery_meets_relative_pressure_residual(ratio):
    """A tiny absolute energy scale must not turn the residual test into a loose one."""
    energy = 3.4091116247062796e-8
    momentum = energy * ratio
    pressure = recover_radiation(
        _flat_spherical_metric(), _conserved(E=energy, S=momentum)
    ).pressure
    v2 = (momentum / (energy + pressure)) ** 2
    residual = (3.0 * pressure - energy) + 4.0 * pressure * v2 / (1.0 - v2)
    assert abs(residual) / energy <= 2.0e-12
