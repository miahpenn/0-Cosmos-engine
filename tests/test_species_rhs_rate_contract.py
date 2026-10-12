"""Regression tests for differential matter-rate extraction and cone checks."""
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


def test_species_rhs_returns_rates_without_a_unit_time_update(monkeypatch):
    rates = ConservedSpecies(
        rest=np.asarray([0.1, 0.2]),
        energy_t=np.asarray([0.3, 0.4]),
        momentum_r=np.asarray([0.5, 0.6]),
    )

    def rate_components(*args, **kwargs):
        zeros = np.zeros(2)
        return rates, zeros, zeros, zeros, zeros

    def forbidden_update(*args, **kwargs):
        raise AssertionError(
            "species_rhs must not emulate a unit-time finite-volume update"
        )

    # This fails against the old wrapper, which calls evolve_species(..., 1.0).
    monkeypatch.setattr(
        matter_rhs_module, "_species_rhs_components", rate_components,
        raising=False,
    )
    monkeypatch.setattr(
        matter_rhs_module, "evolve_species", forbidden_update, raising=False,
    )

    result = species_rhs(
        object(), object(), object(), Species.RADIATION
    )
    np.testing.assert_array_equal(result.rest, rates.rest)
    np.testing.assert_array_equal(result.energy_t, rates.energy_t)
    np.testing.assert_array_equal(result.momentum_r, rates.momentum_r)


def one_cell_flat_metric():
    one = np.asarray([1.0])
    zero = np.asarray([0.0])
    return BSSNMetricSlice(
        r=one.copy(), a=one.copy(), b=one.copy(), X=one.copy(),
        alpha=one.copy(), beta=zero.copy(), Aa=zero.copy(), K=zero.copy(),
        Lambda=zero.copy(), B=zero.copy(),
    )


def test_radiation_cone_check_uses_the_requested_finite_dt(monkeypatch):
    metric = one_cell_flat_metric()
    derivatives = MetricDerivativeSet(time={}, radial={})
    state = ConservedSpecies(
        rest=np.asarray([0.0]),
        energy_t=np.asarray([1.0]),
        momentum_r=np.asarray([0.9]),
    )
    rates = ConservedSpecies(
        rest=np.asarray([0.0]),
        energy_t=np.asarray([-1.0]),
        momentum_r=np.asarray([0.2]),
    )

    def fixed_rates(*args, **kwargs):
        return (
            rates,
            np.asarray([-1.0]),
            np.asarray([0.2]),
            np.asarray([0.0]),
            np.asarray([0.0]),
        )

    monkeypatch.setattr(
        matter_system_module, "_species_rhs_components", fixed_rates
    )

    # The same RHS preserves admissibility over a real small step.
    updated = evolve_species(
        metric, derivatives, state, Species.RADIATION, dt=0.01
    )
    assert updated.energy_t[0] > 0.0
    assert abs(updated.momentum_r[0]) < updated.energy_t[0]

    # A unit-size update of that RHS is inadmissible and must still be rejected.
    with pytest.raises(ValueError, match="finite-time update became inadmissible"):
        evolve_species(
            metric, derivatives, state, Species.RADIATION, dt=1.0
        )
