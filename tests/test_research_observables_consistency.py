from types import SimpleNamespace

import numpy as np
import pytest

from engine import bidirectional_production, coupled_engine
from engine.bidirectional_production import save_ledger_archive
from engine.cosmology_observables import effective_hubble


def test_effective_hubble_uses_volume_weighted_mean_curvature():
    geometry = SimpleNamespace(K=np.asarray([-1.0, -3.0]))
    volumes = np.asarray([1.0, 9.0])

    observed = effective_hubble(geometry, volumes)
    expected = -float(np.sum(volumes * geometry.K) / np.sum(volumes)) / 3.0

    assert np.isclose(observed, expected, rtol=0.0, atol=1.0e-15)
    assert not np.isclose(observed, -float(np.mean(geometry.K)) / 3.0)


def test_research_ledger_archive_uses_requested_output_directory(tmp_path):
    output_dir = tmp_path / "nested" / "artifacts"
    history = {"t": np.asarray([0.0, 0.1]), "H_eff": np.asarray([0.2, 0.1])}

    output_path = save_ledger_archive(history, 40, output_dir)

    assert output_path.parent == output_dir
    assert output_path.name == (
        "GEAR_V5_5_BIDIRECTIONAL_PRODUCTION_LEDGER_N40_T0p1.npz"
    )
    with np.load(output_path, allow_pickle=False) as archive:
        np.testing.assert_array_equal(archive["t"], history["t"])
        np.testing.assert_array_equal(archive["H_eff"], history["H_eff"])


def test_ledger_observables_uses_shared_volume_weighted_hubble(monkeypatch):
    grid = SimpleNamespace(
        centers=np.asarray([1.0, 2.0, 3.0]),
        volumes=np.asarray([1.0, 9.0, 90.0]),
        r_max=3.0,
    )
    geometry = SimpleNamespace(
        a=np.ones(3),
        b=np.ones(3),
        X=np.ones(3),
        alpha=np.ones(3),
        Aa=np.zeros(3),
        K=np.asarray([-1.0, -3.0, -5.0]),
    )
    fields = [np.zeros(3) for _ in range(6)]
    zeros = np.zeros(3)
    monkeypatch.setattr(
        bidirectional_production.q,
        "matter_projection",
        lambda *_: (zeros, zeros, zeros, zeros),
    )
    monkeypatch.setattr(
        bidirectional_production.q,
        "D",
        lambda _grid, _values, _parity: np.ones(3),
    )

    observed = bidirectional_production.ledger_observables(
        grid, geometry, fields, Rcoord=2.0
    )

    assert observed["H_eff"] == effective_hubble(geometry, grid.volumes)
    assert not np.isclose(observed["H_eff"], -float(np.mean(geometry.K)) / 3.0)


def test_coupled_driver_uses_shared_hubble_for_cycle_observation(monkeypatch):
    old_state = SimpleNamespace(K=np.asarray([-1.0, -3.0]), alpha=np.ones(2))
    new_state = SimpleNamespace(
        K=np.asarray([-2.0, -4.0]), alpha=np.asarray([0.9, 0.9])
    )
    grid = SimpleNamespace(volumes=np.asarray([1.0, 9.0]))
    calls = []
    previous_rate = 0.123
    current_rate = 0.456

    def fake_effective_hubble(state, volumes):
        calls.append((state, volumes))
        return previous_rate

    captured = []
    cycle = SimpleNamespace(
        observe=lambda *args: (captured.append(args) or "cycle-event")
    )
    fields = tuple(np.asarray([float(i)]) for i in range(6))
    run = SimpleNamespace(
        grid=grid,
        state=old_state,
        fields=fields,
        t=0.0,
        tau=0.0,
        cycle=cycle,
        handoffs=[],
        history=[],
    )

    monkeypatch.setattr(coupled_engine, "effective_hubble", fake_effective_hubble)
    monkeypatch.setattr(
        coupled_engine,
        "step_true_cmc",
        lambda _g, _s, current_fields, _dt: (new_state, current_fields),
    )
    monkeypatch.setattr(
        coupled_engine,
        "ledger_observables",
        lambda *_args: {"H_eff": current_rate, "tau_rate": 1.0},
    )
    monkeypatch.setattr(coupled_engine, "handoff_from_ledger", lambda *_args: "handoff")

    event = coupled_engine.advance(run, 0.1)

    assert event == "cycle-event"
    assert len(calls) == 1 and calls[0][0] is old_state
    assert calls[0][1] is grid.volumes
    assert captured and captured[0][-2:] == (previous_rate, current_rate)
    assert run.history == [{"H_eff": current_rate, "tau_rate": 1.0}]
    assert run.handoffs == ["handoff"]


def test_ledger_archive_refuses_time_unstamped_history(tmp_path):
    with pytest.raises(ValueError, match="sampled 't' history"):
        save_ledger_archive({"H_eff": np.asarray([0.2])}, 40, tmp_path)
