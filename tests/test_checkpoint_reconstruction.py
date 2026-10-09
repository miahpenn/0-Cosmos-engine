"""Tests for the archived N=320 checkpoint reconstruction path.

These tests do not evolve the numerical machine. The integration replay is
performed separately from the frozen archived checkpoints.
"""
import math

import numpy as np
import pytest

from engine import checkpoint_reconstruction as cr
from engine.production_kernel import V55ProductionKernel, write_checkpoint_npz


@pytest.fixture(scope="module")
def small_state():
    return V55ProductionKernel().initialize(
        resolution=32, r_max=32.0, D_amplitude=1.0e-4,
        include_radiation=True,
    )


def _write(state, path):
    write_checkpoint_npz(state, path)
    return path


def test_loader_round_trips_checkpoint_fields(small_state, tmp_path):
    path = _write(small_state, tmp_path / "state.npz")
    loaded = cr.load_checkpoint_npz(path, r_max=32.0)
    assert np.array_equal(loaded.grid.centers, small_state.grid.centers)
    for key in ("a", "b", "X", "alpha", "beta", "Aa", "K", "Lambda", "B"):
        assert np.array_equal(getattr(loaded.geometry, key), getattr(small_state.geometry, key))
    for key in ("S", "PS", "D", "PD", "phi", "Pi"):
        assert np.array_equal(getattr(loaded.scalars, key), getattr(small_state.scalars, key))
    for species in ("dark_matter", "baryons", "radiation"):
        for key in ("rest", "energy_t", "momentum_r"):
            assert np.array_equal(
                getattr(getattr(loaded.matter, species), key),
                getattr(getattr(small_state.matter, species), key),
            )
    assert loaded.t == small_state.t
    assert loaded.tau == small_state.tau
    assert loaded.e_folds == small_state.e_folds


def test_reconstructed_diagnostics_equal_state_diagnostics(small_state, tmp_path):
    path = _write(small_state, tmp_path / "state.npz")
    loaded = cr.load_checkpoint_npz(path, r_max=32.0)
    before = cr.reconstruct_l2(small_state)
    after = cr.reconstruct_l2(loaded)
    assert before == after
    assert all(math.isfinite(value) for value in after.values())


def test_loader_rejects_wrong_domain(small_state, tmp_path):
    path = _write(small_state, tmp_path / "state.npz")
    with pytest.raises(ValueError, match="radial centers"):
        cr.load_checkpoint_npz(path, r_max=35.0)


def test_loader_rejects_non_finite_checkpoint_field(small_state, tmp_path):
    path = _write(small_state, tmp_path / "state.npz")
    with np.load(path, allow_pickle=False) as source:
        arrays = {key: np.array(source[key], copy=True) for key in source.files}
    arrays["D"][2] = np.nan
    bad = tmp_path / "bad.npz"
    np.savez_compressed(bad, **arrays)
    with pytest.raises(FloatingPointError, match="non-finite"):
        cr.load_checkpoint_npz(bad, r_max=32.0)


def test_exact_checkpoint_time_pairing_selects_its_own_step():
    checkpoint_t = 22.499999999994646
    rows = [{"t": checkpoint_t, "step": "checkpoint"}, {"t": 22.5, "step": "final"}]
    result = cr.match_ledger_row_by_checkpoint_time(checkpoint_t, rows)
    assert result["status"] == cr.MATCHED
    assert result["row_index"] == 0
    assert result["row"]["step"] == "checkpoint"
    assert result["offset"] == 0.0
    assert result["candidates"] == [checkpoint_t]


def test_near_but_unequal_timestamp_is_not_a_match():
    result = cr.match_ledger_row_by_checkpoint_time(
        22.499999999994646, [{"t": 22.5}]
    )
    assert result["status"] == cr.UNMATCHED


def test_duplicate_exact_timestamps_are_ambiguous():
    result = cr.match_ledger_row_by_checkpoint_time(4.0, [{"t": 4.0}, {"t": 4.0}])
    assert result["status"] == cr.AMBIGUOUS


def test_exact_comparison_passes_and_nonzero_difference_fails():
    assert cr.compare_field(1.25e-4, 1.25e-4)["outcome"] == "PASS"
    result = cr.compare_field(1.25e-4 + 1.0e-16, 1.25e-4)
    assert result["outcome"] == "FAIL"
    assert result["abs_err"] > 0.0
    assert result["rel_err"] > 0.0


def test_zero_or_missing_reference_is_undefined():
    assert cr.compare_field(1.0e-6, 0.0)["outcome"] == "UNDEFINED"
    assert cr.compare_field(1.0e-6, None)["outcome"] == "UNDEFINED"


def test_non_finite_value_fails():
    assert cr.compare_field(float("nan"), 1.0)["outcome"] == "FAIL"
