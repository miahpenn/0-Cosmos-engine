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
