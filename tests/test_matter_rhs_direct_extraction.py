"""Regression tests for direct semi-discrete matter RHS extraction."""
import numpy as np
import pytest

import engine.matter_rhs as matter_rhs_module
import engine.matter_system as matter_system_module
from engine.matter_rhs import species_rhs
from engine.matter_system import (
    BSSNMetricSlice,
    ConservedSpecies,
    MetricDerivativeSet,
    Species,
    evolve_species,
)


def flat_metric(n=2):
    r = np.arange(n, dtype=float) + 0.5
    one = np.ones(n)
    zero = np.zeros(n)
    return BSSNMetricSlice(
        r=r,
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


def derivatives(n=2):
    zero = np.zeros(n)
    return MetricDerivativeSet(
        time={},
        radial={
            "alpha": zero.copy(),
            "beta": zero.copy(),
            "rr": zero.copy(),
            "thth": zero.copy(),
        },
    )


def components(n=2):
    # Nontrivial face fluxes exercise each conservative channel.
    face_flux = np.zeros((n + 1, 3), dtype=float)
    face_flux[:, 0] = [0.0, 0.13, -0.04][:n + 1]
    face_flux[:, 1] = [0.0, -0.2, 0.05][:n + 1]
    face_flux[:, 2] = [0.0, 0.07, -0.03][:n + 1]
    source_e = np.linspace(0.11, -0.09, n)
    source_s = np.linspace(-0.03, 0.06, n)
    q_e = np.linspace(0.02, -0.01, n)
    q_s = np.linspace(-0.015, 0.025, n)
    return face_flux, source_e, source_s, q_e, q_s


def test_species_rhs_uses_components_not_a_fictitious_finite_step(monkeypatch):
    metric = flat_metric()
    terms = components()
    monkeypatch.setattr(
        matter_rhs_module,
        "_species_rhs_components",
        lambda *args, **kwargs: terms,
    )

    def forbidden(*args, **kwargs):
        raise AssertionError("species_rhs must not call evolve_species")

    # Protect against accidentally reintroducing the old unit-step wrapper.
    monkeypatch.setattr(matter_rhs_module, "evolve_species", forbidden, raising=False)
    got = species_rhs(metric, derivatives(), object(), Species.RADIATION)

    face, se, ss, qe, qs = terms
    inv_dr = 1.0 / (metric.r[1] - metric.r[0])
    np.testing.assert_array_equal(
        got.rest, -inv_dr * (face[1:, 0] - face[:-1, 0])
    )
    np.testing.assert_array_equal(
        got.energy_t,
        -inv_dr * (face[1:, 1] - face[:-1, 1]) + se + qe,
    )
    np.testing.assert_array_equal(
        got.momentum_r,
        -inv_dr * (face[1:, 2] - face[:-1, 2]) + ss + qs,
    )


def test_finite_update_retains_componentwise_flux_then_source_order(monkeypatch):
    metric = flat_metric()
    state = ConservedSpecies(
        rest=np.array([0.2, 0.3]),
        # Large conserved values make floating-point operation order observable.
        energy_t=np.array([1.0e12, 1.0e12]),
        momentum_r=np.array([0.1, -0.08]),
    )
    terms = list(components())
    face, se, ss, qe, qs = terms
    # Cell 0: separate subtract/add rounds to 1000000000000.0299, while
    # grouping the complete RHS before adding rounds to 1000000000000.03.
    face[:, 1] = np.array([0.0, -1.9, -0.7])
    se[0] = -2.2
    qe[0] = 0.6
    monkeypatch.setattr(
        matter_system_module,
        "_species_rhs_components",
        lambda *args, **kwargs: tuple(terms),
    )
    dt = 0.1
    expected = state.copy()
    inv_dr = 1.0 / (metric.r[1] - metric.r[0])
    for i in range(len(metric.r)):
        expected.rest[i] -= dt * inv_dr * (face[i + 1, 0] - face[i, 0])
        expected.energy_t[i] -= dt * inv_dr * (face[i + 1, 1] - face[i, 1])
        expected.momentum_r[i] -= dt * inv_dr * (face[i + 1, 2] - face[i, 2])
        expected.energy_t[i] += dt * (se[i] + qe[i])
        expected.momentum_r[i] += dt * (ss[i] + qs[i])

    got = evolve_species(
        metric, derivatives(), state, Species.DARK_MATTER, dt=dt,
        dphi_t=np.zeros(2), dphi_r=np.zeros(2),
    )
    np.testing.assert_array_equal(got.rest, expected.rest)
    np.testing.assert_array_equal(got.energy_t, expected.energy_t)
    np.testing.assert_array_equal(got.momentum_r, expected.momentum_r)

    regrouped_energy = state.energy_t.copy()
    for i in range(len(metric.r)):
        flux_rate = -inv_dr * (face[i + 1, 1] - face[i, 1])
        regrouped_energy[i] += dt * (flux_rate + (se[i] + qe[i]))
    assert got.energy_t[0] != regrouped_energy[0], (
        "the fixture must detect regrouping of finite-step flux/source arithmetic"
    )


def test_actual_finite_time_radiation_update_still_checks_cone(monkeypatch):
    metric = flat_metric()
    state = ConservedSpecies(
        rest=np.zeros(2),
        energy_t=np.ones(2),
        momentum_r=np.full(2, 0.9),
    )
    face = np.zeros((3, 3), dtype=float)
    source_e = np.full(2, -1.0)
    source_s = np.full(2, 0.2)
    q_e = np.zeros(2)
    q_s = np.zeros(2)
    monkeypatch.setattr(
        matter_system_module,
        "_species_rhs_components",
        lambda *args, **kwargs: (
            face, source_e, source_s, q_e, q_s
        ),
    )

    small = evolve_species(
        metric, derivatives(), state, Species.RADIATION, dt=0.01
    )
    assert np.all(small.energy_t > 0.0)
    assert np.all(np.abs(small.momentum_r) < small.energy_t)

    with pytest.raises(ValueError, match="radiation state inadmissible"):
        evolve_species(
            metric, derivatives(), state, Species.RADIATION, dt=1.0
        )


def test_components_helper_preserves_radiation_zero_rest_current_flux():
    metric = flat_metric()
    metric_arrays = matter_system_module._metric_arrays(metric)
    # A unit Eulerian radiation state with exactly zero radial momentum.
    state = ConservedSpecies(
        rest=np.zeros(2),
        energy_t=np.asarray([m.sqrt_gamma for m in metric_arrays]),
        momentum_r=np.zeros(2),
    )
    face_flux, _, _, _, _ = matter_system_module._species_rhs_components(
        metric, derivatives(), state, Species.RADIATION
    )
    np.testing.assert_array_equal(face_flux[:, 0], np.zeros(3))


def test_flat_homogeneous_radiation_pressure_flux_uses_true_face_area():
    """The origin pressure flux and spherical source must cancel in flat FLRW."""
    metric = flat_metric(n=3)
    metric_arrays = matter_system_module._metric_arrays(metric)
    state = ConservedSpecies(
        rest=np.zeros(3),
        energy_t=np.asarray([m.sqrt_gamma for m in metric_arrays]),
        momentum_r=np.zeros(3),
    )
    zero = np.zeros(3)
    terms = MetricDerivativeSet(
        time={},
        radial={
            "alpha": zero.copy(),
            "beta": zero.copy(),
            "rr": zero.copy(),
            "thth": 2.0 * metric.r,
        },
    )

    face_flux, _source_e, source_s, _q_e, _q_s = (
        matter_system_module._species_rhs_components(
            metric, terms, state, Species.RADIATION
        )
    )

    # dr=1: the first interior face is at r=1, so sqrt(gamma)=r^2=1.
    # The old arithmetic average of centre volume factors gave 1.25 here.
    expected_face_pressure_flux = 1.0 / 3.0
    assert face_flux[1, 2] == pytest.approx(
        expected_face_pressure_flux, rel=2.0e-14, abs=0.0
    )
    rhs = -(face_flux[1:, 2] - face_flux[:-1, 2]) + source_s
    np.testing.assert_allclose(rhs, np.zeros(3), rtol=0.0, atol=5.0e-15)
    # The outer boundary at r=3 also carries the true spherical area.
    assert face_flux[-1, 2] == pytest.approx(9.0 / 3.0, rel=2.0e-14, abs=0.0)


def test_nontrivial_metric_radiation_flux_uses_interpolated_face_fields():
    """Flux geometry is evaluated at the face radius, not averaged after derivation."""
    metric = BSSNMetricSlice(
        r=np.asarray([0.5, 1.5, 2.5]),
        a=np.asarray([0.90, 1.00, 1.10]),
        b=np.asarray([1.10, 1.20, 1.30]),
        X=np.asarray([0.90, 0.95, 1.00]),
        alpha=np.asarray([0.80, 0.90, 1.00]),
        beta=np.zeros(3),
        Aa=np.zeros(3),
        K=np.zeros(3),
        Lambda=np.zeros(3),
        B=np.zeros(3),
    )
    metric_arrays = matter_system_module._metric_arrays(metric)
    state = ConservedSpecies(
        rest=np.zeros(3),
        energy_t=np.asarray([m.sqrt_gamma for m in metric_arrays]),
        momentum_r=np.zeros(3),
    )

    face_flux, _source_e, _source_s, _q_e, _q_s = (
        matter_system_module._species_rhs_components(
            metric, derivatives(n=3), state, Species.RADIATION
        )
    )

    r_face = 0.5 * (metric.r[0] + metric.r[1])
    a_face = 0.5 * (metric.a[0] + metric.a[1])
    b_face = 0.5 * (metric.b[0] + metric.b[1])
    X_face = 0.5 * (metric.X[0] + metric.X[1])
    alpha_face = 0.5 * (metric.alpha[0] + metric.alpha[1])
    expected_sqrt_gamma = np.sqrt(a_face) * b_face * r_face**2 / X_face**3
    expected_flux = expected_sqrt_gamma * alpha_face / 3.0
    assert face_flux[1, 2] == pytest.approx(
        expected_flux, rel=3.0e-13, abs=0.0
    )


def test_nonzero_shift_radiation_flux_uses_face_metric():
    """The shared face metric must preserve shift terms with curved geometry."""
    metric = BSSNMetricSlice(
        r=np.asarray([0.5, 1.5, 2.5]),
        a=np.asarray([0.90, 1.00, 1.10]),
        b=np.asarray([1.10, 1.20, 1.30]),
        X=np.asarray([0.90, 0.95, 1.00]),
        alpha=np.asarray([0.80, 0.90, 1.00]),
        beta=np.asarray([0.02, 0.06, 0.10]),
        Aa=np.zeros(3),
        K=np.zeros(3),
        Lambda=np.zeros(3),
        B=np.zeros(3),
    )
    metric_arrays = matter_system_module._metric_arrays(metric)
    state = ConservedSpecies(
        rest=np.zeros(3),
        energy_t=np.asarray([m.sqrt_gamma for m in metric_arrays]),
        momentum_r=np.zeros(3),
    )
    face_flux, _source_e, _source_s, _q_e, _q_s = (
        matter_system_module._species_rhs_components(
            metric, derivatives(n=3), state, Species.RADIATION
        )
    )

    r_face = 0.5 * (metric.r[0] + metric.r[1])
    a_face = 0.5 * (metric.a[0] + metric.a[1])
    b_face = 0.5 * (metric.b[0] + metric.b[1])
    X_face = 0.5 * (metric.X[0] + metric.X[1])
    alpha_face = 0.5 * (metric.alpha[0] + metric.alpha[1])
    beta_face = 0.5 * (metric.beta[0] + metric.beta[1])
    sqrt_gamma_face = np.sqrt(a_face) * b_face * r_face**2 / X_face**3
    # The Eulerian radiation state has E=1 and S_r=0; the energy flux is
    # -sqrt(gamma)*beta*E, while the radial pressure flux is sqrt(gamma)*alpha*p.
    assert face_flux[1, 1] == pytest.approx(
        -sqrt_gamma_face * beta_face, rel=3.0e-13, abs=0.0
    )
    assert face_flux[1, 2] == pytest.approx(
        sqrt_gamma_face * alpha_face / 3.0, rel=3.0e-13, abs=0.0
    )
