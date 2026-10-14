"""Contracts for the reconstructed, nonphysical GEAR cycle operator."""
import mpmath as mp
import pytest

from research.gear_cycle_operator_v1 import (
    DEFAULT_A0,
    DEFAULT_A_T,
    DEFAULT_H0,
    DEFAULT_HDOT0,
    DEFAULT_PRECISION_DIGITS,
    pole_swap,
    run_closure_regression,
)


def test_archive_seed_closes_exactly_after_27_cycles_and_54_swaps():
    result = run_closure_regression()

    assert result["status"] == "reconstructed_numerical_operator_regression"
    assert result["source_status"] == "reconstruction_not_recovered_original_source"
    assert result["physical_bounce_claim"] is False
    assert result["precision_digits"] == 80
    assert result["full_cycles"] == 27
    assert result["pole_swaps"] == 54
    assert result["exact_closure"] is True
    assert result["exact_zero_deltas"] == {"a": True, "H": True, "Hdot": True}
    assert result["deltas"] == {"a": "0.0", "H": "0.0", "Hdot": "0.0"}


def test_one_full_cycle_is_two_pole_swaps_and_returns_archive_seed():
    with mp.workdps(80):
        initial = (
            mp.mpf(DEFAULT_A0),
            mp.mpf(DEFAULT_H0),
            mp.mpf(DEFAULT_HDOT0),
        )
        turnaround = mp.mpf(DEFAULT_A_T)
        once = pole_swap(*initial, turnaround)
        twice = pole_swap(*once, turnaround)
        assert twice == initial


def test_zero_cycles_returns_seed_without_a_swap():
    result = run_closure_regression(full_cycles=0)
    assert result["pole_swaps"] == 0
    assert result["exact_closure"] is True
    assert result["initial"] == result["final"]


def test_custom_cycle_count_preserves_two_swaps_per_cycle():
    result = run_closure_regression(full_cycles=3, precision_digits=60)
    assert result["full_cycles"] == 3
    assert result["pole_swaps"] == 6
    assert result["exact_closure"] is True


@pytest.mark.parametrize(
    "kwargs",
    [
        {"full_cycles": -1},
        {"full_cycles": True},
        {"full_cycles": 2.5},
        {"precision_digits": 14},
        {"precision_digits": True},
        {"a_T": "0"},
        {"a0": "-1"},
        {"H0": "nan"},
    ],
)
def test_closure_runner_rejects_invalid_parameters(kwargs):
    with pytest.raises(ValueError):
        run_closure_regression(**kwargs)


@pytest.mark.parametrize(
    "args",
    [
        ("0", "0", "0", "1"),
        ("-1", "0", "0", "1"),
        ("1", "nan", "0", "1"),
        ("1", "0", "0", "0"),
    ],
)
def test_single_swap_rejects_invalid_domain_or_nonfinite_values(args):
    with pytest.raises(ValueError):
        pole_swap(*args)


def test_runner_does_not_change_callers_global_mpmath_precision():
    before = mp.mp.dps
    result = run_closure_regression(precision_digits=70)
    assert result["precision_digits"] == 70
    assert mp.mp.dps == before


def test_decimal_seed_is_retained_as_high_precision_text():
    result = run_closure_regression(full_cycles=0, precision_digits=80)
    assert result["initial"]["a"].startswith("3.271828182845904523536028747135")
    assert DEFAULT_PRECISION_DIGITS == 80
