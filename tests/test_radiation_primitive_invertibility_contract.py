"""Boundary contract between radiation admissibility and primitive inversion."""
from __future__ import annotations

import numpy as np
import pytest

from engine.matter_system import (
    BSSNMetricSlice,
    ConservedSpecies,
    MetricDerivativeSet,
    RadiationAdmissibilityError,
    Species,
    _metric_arrays,
    _radiation_conserved,
    _radiation_flux,
    _valencia_flux,
    _species_rhs_components,
    _valencia_source,
    KAPPA,
    primitives,
    species_projection,
    validate_species_admissibility,
)
from engine.valencia import (
    FluidConserved,
    FluidPrimitive,
    RadiationAdmissibilityError as RecoveryAdmissibilityError,
    SphericalMetric,
    recover_radiation,
    stress_tensor_contravariant,
    mixed_flux,
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


def test_near_cone_rhs_assembly_uses_preserved_source_norm_and_outer_flux():
    """Exercise the shared conservative RHS entry with the captured outer state."""
    g = BSSNMetricSlice(
        r=np.asarray([79.5]),
        a=np.asarray([0.5883984147590013]),
        b=np.asarray([1.3036606161142992]),
        X=np.asarray([0.9545747827302261]),
        alpha=np.asarray([1.0]),
        beta=np.asarray([0.0]),
        Aa=np.asarray([0.0]),
        K=np.asarray([0.0]),
        Lambda=np.asarray([0.0]),
        B=np.asarray([0.0]),
    )
    metric = _metric_arrays(g)[0]
    state = ConservedSpecies(
        rest=np.asarray([0.0]),
        energy_t=np.asarray([2.4771114076139816e-4]),
        momentum_r=np.asarray([1.990541924362441e-4]),
    )
    derivatives = MetricDerivativeSet(
        time={},
        radial={
            "alpha": np.asarray([0.0]),
            "beta": np.asarray([0.0]),
            "rr": np.asarray([1.0]),
            "thth": np.asarray([0.0]),
        },
    )

    face_flux, _source_e, source_s, _q_e, _q_s = _species_rhs_components(
        g, derivatives, state, Species.RADIATION
    )
    primitive = recover_radiation(
        metric,
        _conserved(
            E=2.4771114076139816e-4,
            S=1.990541924362441e-4,
        ),
    )
    h = primitive.rho + primitive.pressure
    W = primitive.lorentz()
    expected_Srr = (
        h * W * W * primitive.v2 * metric.gamma_rr_inv
        + primitive.pressure * metric.gamma_rr_inv
    )
    expected_source_s = metric.sqrt_gamma * 0.5 * metric.alpha * expected_Srr
    assert source_s[0] == pytest.approx(expected_source_s, rel=2.0e-14, abs=0.0)

    expected_outer_flux = (
        metric.alpha * metric.gamma_rr_inv * state.momentum_r[0]
        - metric.beta * state.energy_t[0]
    )
    assert face_flux[-1, 1] == pytest.approx(expected_outer_flux, rel=2.0e-14, abs=0.0)

@pytest.mark.parametrize("shift", [0.0, 0.17])
def test_near_cone_radial_momentum_flux_uses_preserved_norm_with_shift(shift):
    """The direct cell flux must not reconstruct v^2 from rounded coordinate v."""
    metric = spherical_metric_from_bssn(
        r=79.5,
        a=0.5883984147590013,
        b=1.3036606161142992,
        X=0.9545747827302261,
        alpha=1.0,
        beta=shift,
    )
    primitive = recover_radiation(
        metric,
        _conserved(
            E=2.4771114076139816e-4,
            S=1.990541924362441e-4,
        ),
    )
    assert primitive.velocity_norm_squared is not None
    assert primitive.v2 != metric.gamma_rr * primitive.v_r * primitive.v_r

    W = primitive.lorentz()
    h = primitive.rho + primitive.pressure
    S_r = h * W * W * metric.gamma_rr * primitive.v_r
    expected = metric.sqrt_gamma * (
        metric.alpha * h * W * W * primitive.v2
        + metric.alpha * primitive.pressure
        - metric.beta * S_r
    )

    # Radiation's active direct-cell outer flux and the shared Valencia
    # fluid flux use the same radial stress identity. The reconstructed HLL
    # face primitives have no saved norm override, so their own face norm is
    # still derived from the face metric and reconstructed coordinate speed.
    assert _radiation_flux(metric, primitive)[2] == pytest.approx(
        expected, rel=2.0e-15, abs=0.0
    )
    assert _valencia_flux(metric, primitive)[2] == pytest.approx(
        expected, rel=2.0e-15, abs=0.0
    )

@pytest.mark.parametrize("shift", [0.0, 0.17])
def test_near_cone_mixed_tensor_helpers_preserve_solved_norm_with_shift(shift):
    """Shifted mixed-tensor helpers retain the primitive's solved spatial v^2."""
    metric = spherical_metric_from_bssn(
        r=79.5,
        a=0.5883984147590013,
        b=1.3036606161142992,
        X=0.9545747827302261,
        alpha=1.0,
        beta=shift,
    )
    primitive = recover_radiation(
        metric,
        _conserved(
            E=2.4771114076139816e-4,
            S=1.990541924362441e-4,
        ),
    )
    assert primitive.velocity_norm_squared is not None

    W = primitive.lorentz()
    h = primitive.rho + primitive.pressure
    shift_over_alpha = metric.beta / metric.alpha
    gamma_rr_contravariant = (
        metric.gamma_rr_inv - shift_over_alpha * shift_over_alpha
    )
    expected_Trr = h * W * W * (
        primitive.v2 * metric.gamma_rr_inv
        - 2.0 * shift_over_alpha * primitive.v_r
        + shift_over_alpha * shift_over_alpha
    ) + primitive.pressure * gamma_rr_contravariant
    T = stress_tensor_contravariant(metric, primitive)
    assert T["rr"] == pytest.approx(expected_Trr, rel=2.0e-15, abs=0.0)

    expected_Tr_t = h * W * W * (
        -metric.alpha * primitive.v_r
        + metric.beta
        + metric.beta * primitive.v2
        - shift_over_alpha * metric.shift_cov * primitive.v_r
    )
    expected_Tr_r = h * W * W * (
        primitive.v2
        - shift_over_alpha * metric.gamma_rr * primitive.v_r
    ) + primitive.pressure
    flux = mixed_flux(metric, primitive)
    assert flux.energy_t == pytest.approx(
        metric.sqrt_minus_g * expected_Tr_t, rel=2.0e-15, abs=0.0
    )
    assert flux.momentum_r == pytest.approx(
        metric.sqrt_minus_g * expected_Tr_r, rel=2.0e-15, abs=0.0
    )

def test_near_cone_flux_regression_detects_rounded_v_squared_product():
    """Use a homogeneous density rescaling where the legacy product differs by one ULP."""
    metric = spherical_metric_from_bssn(
        r=79.5,
        a=0.5883984147590013,
        b=1.3036606161142992,
        X=0.9545747827302261,
        alpha=1.0,
        beta=0.0,
    )
    recovered = recover_radiation(
        metric,
        _conserved(
            E=2.4771114076139816e-4,
            S=1.990541924362441e-4,
        ),
    )
    # Density-rescale the same near-null velocity while retaining rho=3p.
    # This keeps the algebraic state admissible and makes the old/new product
    # difference survive the final binary64 flux rounding.
    h = 1.1853981428841001e-25
    pressure = h / 4.0
    primitive = FluidPrimitive(
        rho=3.0 * pressure,
        pressure=pressure,
        v_r=recovered.v_r,
        gamma_rr=recovered.gamma_rr,
        velocity_norm_squared=recovered.v2,
    )
    W = primitive.lorentz()
    S_r = h * W * W * metric.gamma_rr * primitive.v_r
    legacy = metric.sqrt_gamma * (
        metric.alpha * S_r * primitive.v_r
        + metric.alpha * primitive.pressure
        - metric.beta * S_r
    )
    expected = metric.sqrt_gamma * (
        metric.alpha * h * W * W * primitive.v2
        + metric.alpha * primitive.pressure
        - metric.beta * S_r
    )
    assert legacy != expected, "fixture must detect the rounded v_r^2 product"
    assert _radiation_flux(metric, primitive)[2] == pytest.approx(
        expected, rel=0.0, abs=0.0
    )
    assert _valencia_flux(metric, primitive)[2] == pytest.approx(
        expected, rel=0.0, abs=0.0
    )

