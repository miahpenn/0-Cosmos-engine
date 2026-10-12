"""Small deterministic checks for the one-step radiation/MS diagnostic."""
from types import SimpleNamespace

import numpy as np
import pytest

from tools.audit_cmc_lapse_radiation_ms_one_step_ab import (
    misner_sharp_interval_budget,
    radiation_cone_snapshot,
)


def test_cone_snapshot_reconstructs_known_predictor_crossing():
    grid = SimpleNamespace(
        centers=np.asarray([79.5], dtype=float),
        r_max=80.0,
    )
    geometry = SimpleNamespace(
        a=np.asarray([0.6032035002114748], dtype=float),
        X=np.asarray([0.946558319470313], dtype=float),
    )
    radiation = SimpleNamespace(
        energy_t=np.asarray([0.0002965253923984548], dtype=float),
        momentum_r=np.asarray([0.0002433162925154048], dtype=float),
    )
    matter = SimpleNamespace(radiation=radiation)

    result = radiation_cone_snapshot(grid, geometry, matter, "known-crossing")
    full = result["regions"]["full_domain"]
    assert full["min_cone_margin"] == pytest.approx(-1.6900205171508503e-8, abs=2e-15)
    assert full["max_absS_over_E"] == pytest.approx(1.0000569941246342, abs=2e-12)
    assert full["admissible"] is False
    assert result["top_five_ratio_cells"][0]["i"] == 0



def test_species_rhs_does_not_emulate_a_unit_time_update(monkeypatch):
    """A rate query must not run finite-time radiation cone checks."""
    import engine.matter_rhs as rhs_module
    from engine.matter_system import ConservedSpecies, Species

    rates = ConservedSpecies(
        rest=np.asarray([0.1, 0.2]),
        energy_t=np.asarray([0.3, 0.4]),
        momentum_r=np.asarray([0.5, 0.6]),
    )

    def rate_components(*args, **kwargs):
        return rates, np.zeros(2), np.zeros(2), np.zeros(2), np.zeros(2)

    def forbidden_finite_update(*args, **kwargs):
        raise AssertionError(
            "species_rhs must not emulate an arbitrary unit-time update"
        )

    # The first patch returns an exact rate sentinel. The second makes the old
    # evolve-by-one-then-subtract implementation fail this regression test.
    monkeypatch.setattr(
        rhs_module, "_species_rhs_components", rate_components, raising=False
    )
    monkeypatch.setattr(
        rhs_module, "evolve_species", forbidden_finite_update, raising=False
    )
    result = rhs_module.species_rhs(
        object(), object(), object(), Species.RADIATION
    )
    np.testing.assert_array_equal(result.rest, rates.rest)
    np.testing.assert_array_equal(result.energy_t, rates.energy_t)
    np.testing.assert_array_equal(result.momentum_r, rates.momentum_r)

def test_ms_one_step_budget_uses_existing_flux_plus_work_identity():
    initial = {"worldtube_r10": {"M_MS": 1.0, "R_sigma": 10.0, "flux_T": 0.07, "work_pR": 0.03}}
    final = {"worldtube_r10": {"M_MS": 1.02, "R_sigma": 10.01, "flux_T": 0.20, "work_pR": 0.10}}

    result = misner_sharp_interval_budget(initial, final, dt=0.1)
    assert result["trapezoidal_flux_work_integral"] == pytest.approx(0.02)
    assert result["integrated_residual"] == pytest.approx(0.0, abs=1e-15)
    assert result["scaled_integrated_residual"] == pytest.approx(0.0, abs=1e-14)
    assert "no tolerance or pass/fail threshold" in result["admission_note"]
