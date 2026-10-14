from types import SimpleNamespace

import numpy as np

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
        "GEAR_V5_5_BIDIRECTIONAL_PRODUCTION_LEDGER_N40_T24.npz"
    )
    with np.load(output_path, allow_pickle=False) as archive:
        np.testing.assert_array_equal(archive["t"], history["t"])
        np.testing.assert_array_equal(archive["H_eff"], history["H_eff"])
