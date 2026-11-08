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


def test_drivers_run_on_the_production_kernel_not_the_reference_stack():
    import ast
    import inspect

    for module in (bidirectional_production, coupled_engine):
        tree = ast.parse(inspect.getsource(module))
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom):
                imported.add(node.module or "")
                imported |= {alias.name for alias in node.names}
            elif isinstance(node, ast.Import):
                imported |= {alias.name for alias in node.names}
        assert not any("reference_pirk_unified" in name or "true_cmc_reference" in name
                       for name in imported), module.__name__


def test_ledger_observables_are_the_production_diagnostics():
    from engine.production_kernel import V55ProductionKernel

    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=24, r_max=40.0)
    row = bidirectional_production.ledger_observables(kernel, state)
    diagnostics = kernel.diagnostics(state, profiles=True)
    for key in bidirectional_production.LEDGER_KEYS:
        assert row[key] == float(diagnostics[key])
    assert row["H_eff"] == effective_hubble(state.geometry, state.grid.volumes)
    # The COSMOS scalar starts at the archive operating point phi = 0.179055,
    # Pi = 0.00341913 (the old reference initializer had phi = 33.8983, Pi = 0.179055).
    assert row["phi_outer"] == pytest.approx(0.179055, abs=1e-12)
    assert row["Pi_outer"] == pytest.approx(0.00341913, rel=1e-6)
    with pytest.raises(ValueError, match="fixed at coordinate r = 10"):
        bidirectional_production.ledger_observables(kernel, state, Rcoord=12.0)


def test_coupled_driver_steps_are_production_kernel_steps():
    from engine.production_kernel import V55ProductionKernel

    run = coupled_engine.initialize(n=24)
    kernel = V55ProductionKernel()
    reference = kernel.initialize(resolution=24, r_max=40.0)
    dt = 0.06 * float(run.grid.dr)
    for _ in range(3):
        coupled_engine.advance(run, dt)
        reference = kernel.step(reference, dt)
    assert run.t == reference.t and run.tau == reference.tau
    for name in ("a", "b", "X", "alpha", "K", "Aa", "Lambda"):
        np.testing.assert_array_equal(getattr(run.state.geometry, name),
                                      getattr(reference.geometry, name))
    assert len(run.history) == 3 and len(run.handoffs) == 3
    traj = coupled_engine.trajectory(run)
    np.testing.assert_array_equal(traj["t"], [row["t"] for row in reference.history])


def test_ledger_archive_refuses_time_unstamped_history(tmp_path):
    with pytest.raises(ValueError, match="sampled 't' history"):
        save_ledger_archive({"H_eff": np.asarray([0.2])}, 40, tmp_path)
