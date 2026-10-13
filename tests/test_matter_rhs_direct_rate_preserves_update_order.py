"""Direct matter RHS extraction must not alter finite-dt update ordering."""
import numpy as np

import engine.matter_system as matter_system
from engine.matter_rhs import species_rhs
from engine.matter_system import (
    BSSNMetricSlice,
    ConservedSpecies,
    MetricDerivativeSet,
    Species,
    evolve_species,
)


def _two_cell_metric():
    one = np.ones(2, dtype=float)
    zero = np.zeros(2, dtype=float)
    return BSSNMetricSlice(
        r=np.array([0.5, 1.5], dtype=float),
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


def _operator_components():
    face_flux = np.array([
        [-9.148061683015676e-08, 1.5804727032989443e-06, -8.298641512462363e-07],
        [-2.8179407189044416e-07, 1.5244525849834649e-06, 4.826129610674838e-07],
        [4.948221662761881e-07, -3.8358587868706044e-07, -1.2299037470675234e-06],
    ], dtype=float)
    source_e = np.array([2.1810060918682255, -24.555196128011747], dtype=float)
    source_s = np.array([1.1241769638378333e-06, 4.189773720805749e-07], dtype=float)
    coupling_e = np.array([-0.09784490776549772, -1.283098472616614], dtype=float)
    coupling_s = np.array([6208.109864789437, 9317.861603178651], dtype=float)
    inv_dr = 17.3
    return face_flux, source_e, source_s, coupling_e, coupling_s, inv_dr


def test_species_rhs_returns_direct_rates_without_calling_finite_update(monkeypatch):
    components = _operator_components()

    def fake_components(*args, **kwargs):
        return components

    def forbidden_update(*args, **kwargs):
        raise AssertionError("species_rhs must not advance a fictitious dt=1 state")

    monkeypatch.setattr(matter_system, "_species_operator_components", fake_components)
    monkeypatch.setattr(matter_system, "evolve_species", forbidden_update)

    result = species_rhs(object(), object(), object(), Species.BARYON)

    face, source_e, source_s, coupling_e, coupling_s, inv_dr = components
    expected_rest = -inv_dr * (face[1:, 0] - face[:-1, 0])
    expected_energy = -inv_dr * (face[1:, 1] - face[:-1, 1])
    expected_energy += source_e + coupling_e
    expected_momentum = -inv_dr * (face[1:, 2] - face[:-1, 2])
    expected_momentum += source_s + coupling_s
    np.testing.assert_array_equal(result.rest, expected_rest)
    np.testing.assert_array_equal(result.energy_t, expected_energy)
    np.testing.assert_array_equal(result.momentum_r, expected_momentum)


def test_evolve_species_retains_legacy_componentwise_finite_dt_order(monkeypatch):
    """Pin flux subtraction first, then source addition; do not regroup as U+dt*RHS."""
    components = _operator_components()

    def fake_components(*args, **kwargs):
        return components

    monkeypatch.setattr(matter_system, "_species_operator_components", fake_components)
    state_values = np.array([-1601268.464806936, 810297.4043016228], dtype=float)
    state = ConservedSpecies(
        rest=state_values.copy(),
        energy_t=state_values.copy(),
        momentum_r=state_values.copy(),
    )
    dt = 0.03
    updated = evolve_species(
        _two_cell_metric(),
        MetricDerivativeSet(time={}, radial={}),
        state,
        Species.BARYON,
        dt,
    )

    face, source_e, source_s, coupling_e, coupling_s, inv_dr = components
    expected = state.copy()
    for i in range(len(state.rest)):
        expected.rest[i] -= dt * inv_dr * (face[i + 1, 0] - face[i, 0])
        expected.energy_t[i] -= dt * inv_dr * (face[i + 1, 1] - face[i, 1])
        expected.momentum_r[i] -= dt * inv_dr * (face[i + 1, 2] - face[i, 2])
        expected.energy_t[i] += dt * (source_e[i] + coupling_e[i])
        expected.momentum_r[i] += dt * (source_s[i] + coupling_s[i])

    np.testing.assert_array_equal(updated.rest, expected.rest)
    np.testing.assert_array_equal(updated.energy_t, expected.energy_t)
    np.testing.assert_array_equal(updated.momentum_r, expected.momentum_r)

    # This fixture would distinguish the preceding sequence from PR #8's
    # grouped direct-rate update: U + dt * (-div(F) + source).
    grouped_energy_rhs = -inv_dr * (face[1:, 1] - face[:-1, 1])
    grouped_energy_rhs += source_e + coupling_e
    grouped_energy = state.energy_t + dt * grouped_energy_rhs
    assert not np.array_equal(grouped_energy, expected.energy_t)
