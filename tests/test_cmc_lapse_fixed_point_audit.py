"""Unit tests for the read-only CMC fixed-point audit helpers."""
import numpy as np
import pytest

from tools.audit_cmc_lapse_fixed_point import profile_stats, weighted_mean


def test_weighted_mean_uses_supplied_volume_weights():
    values = np.asarray([1.0, 3.0, 7.0])
    weights = np.asarray([1.0, 2.0, 1.0])
    assert weighted_mean(values, weights) == pytest.approx(3.5)


def test_profile_stats_reports_slice_and_weighted_rms():
    values = np.asarray([100.0, 100.0, 1.0, 3.0, 99.0])
    weights = np.asarray([1.0, 1.0, 1.0, 3.0, 1.0])
    stats = profile_stats(values, weights, start=2, stop=4)
    assert stats["count"] == 2
    assert stats["cell_start_inclusive"] == 2
    assert stats["cell_stop_exclusive"] == 4
    assert stats["max_abs"] == 3.0
    assert stats["rms"] == pytest.approx(np.sqrt(5.0))
    assert stats["proper_volume_weighted_rms"] == pytest.approx(np.sqrt(7.0))


@pytest.mark.parametrize(
    "values,weights",
    [
        (np.asarray([1.0, np.nan]), np.asarray([1.0, 1.0])),
        (np.asarray([1.0, 2.0]), np.asarray([0.0, 0.0])),
        (np.asarray([1.0]), np.asarray([1.0, 2.0])),
    ],
)
def test_weighted_mean_fails_closed_on_invalid_inputs(values, weights):
    with pytest.raises(ValueError):
        weighted_mean(values, weights)
