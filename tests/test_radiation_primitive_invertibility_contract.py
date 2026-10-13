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
    _radiation_conserved,
    _radiation_flux,
    _valencia_source,
    KAPPA,
    primitives,
    species_projection,
    validate_species_admissibility,
)
from engine.valencia import (
    FluidConserved,
    RadiationAdmissibilityError as RecoveryAdmissibilityError,
    SphericalMetric,
    recover_radiation,
    spherical_metric_from_bssn,
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
    # Recomputing v^2 through the ratio form is mildly ill-conditioned near |S|/E=1.
    # This bound admits the observed float64 rounding floor but not the old scale-loose solve.
    assert abs(residual) / energy <= 5.0e-11


def test_near_cone_recovery_preserves_velocity_norm_and_outer_flux_identity():
    """A captured null-cone-near state must survive primitive/flux round-trip.

    State and metric are from the SHA-verified current-main radiation replay
    (#61, workflow 38337821862, artifact 11785239489), outer cell r=79.5.
    """
    metric = spherical_metric_from_bssn(
        r=79.5,
        a=0.5883984147590013,
        b=1.3036606161142992,
        X=0.9545747827302261,
        alpha=1.0,
        beta=0.0,
    )
    conserved = _conserved(
        E=2.4771114076139816e-4,
        S=1.990541924362441e-4,
    )

    primitive = recover_radiation(metric, conserved)
    rebuilt = _radiation_conserved(metric, primitive)
    np.testing.assert_allclose(
        rebuilt[1:],
        np.asarray([conserved.energy_t, conserved.momentum_r]),
        rtol=2.0e-14,
        atol=0.0,
    )

    # With U_E=sqrt(gamma)E and U_r=sqrt(gamma)S_r, the physical
    # Valencia energy flux is alpha*gamma^rr*U_r - beta*U_E.
    face_flux = _radiation_flux(metric, primitive)
    expected_energy_flux = (
        metric.alpha * metric.gamma_rr_inv * conserved.momentum_r
        - metric.beta * conserved.energy_t
    )
    assert face_flux[1] == pytest.approx(expected_energy_flux, rel=2.0e-14, abs=0.0)


def _captured_near_cone_metric_and_primitive():
    metric = spherical_metric_from_bssn(
        r=79.5,
        a=0.5883984147590013,
        b=1.3036606161142992,
        X=0.9545747827302261,
        alpha=1.0,
        beta=0.0,
    )
    primitive = recover_radiation(
        metric,
        _conserved(
            E=2.4771114076139816e-4,
            S=1.990541924362441e-4,
        ),
    )
    assert primitive.velocity_norm_squared is not None
    # This is the recorded one-ulp discrepancy that invalidated recomputing v^2.
    assert primitive.v2 != metric.gamma_rr * primitive.v_r * primitive.v_r
    return metric, primitive


def test_near_cone_einstein_projection_uses_preserved_velocity_norm():
    metric, primitive = _captured_near_cone_metric_and_primitive()
    projected = species_projection([metric], [primitive])
    h = primitive.rho + primitive.pressure
    W = primitive.lorentz()
    expected_pr = (h * W * W * primitive.v2 + primitive.pressure) / KAPPA
    assert projected["pr"][0] == pytest.approx(expected_pr, rel=2.0e-15, abs=0.0)


def test_near_cone_active_rhs_source_uses_preserved_velocity_norm():
    """The source consumed by _species_rhs_components must share the solved v^2."""
    metric, primitive = _captured_near_cone_metric_and_primitive()
    _, source_s = _valencia_source(
        metric, primitive,
        K=0.0, Aa=0.0,
        dr_alpha=0.0, dr_beta=0.0,
        dr_grr=1.0, dr_gthth=0.0,
    )
    h = primitive.rho + primitive.pressure
    W = primitive.lorentz()
    expected_Srr = (
        h * W * W * primitive.v2 * metric.gamma_rr_inv
        + primitive.pressure * metric.gamma_rr_inv
    )
    expected_source_s = metric.sqrt_gamma * 0.5 * metric.alpha * expected_Srr
    assert source_s == pytest.approx(expected_source_s, rel=2.0e-15, abs=0.0)
