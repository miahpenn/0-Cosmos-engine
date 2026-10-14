"""High-precision reconstruction of the GEAR cyclic pole-swap operator.

This is a numerical involution/regression reconstructed from the preserved
archive. It is not recovered source for the original notebook and is not a
physical bounce law or a cosmological time integrator.
"""
from __future__ import annotations

from typing import Any

import mpmath as mp


DEFAULT_PRECISION_DIGITS = 80
DEFAULT_FULL_CYCLES = 27
DEFAULT_A_T = "9.80859460229"
DEFAULT_A0 = "3.2718281828459045235360287471352662497757"
DEFAULT_H0 = "0.123456789012345678901234567890123456789"
DEFAULT_HDOT0 = "-0.1131778"


def _to_finite_mpf(value: Any, name: str) -> mp.mpf:
    try:
        result = mp.mpf(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a valid real number") from exc
    if not mp.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def pole_swap(
    a: Any,
    H: Any,
    Hdot: Any,
    a_T: Any,
) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    """Apply one archive-reconstructed pole swap to (a, H, Hdot).

    The caller controls mpmath's active precision. Use decimal strings to seed
    high-precision values; avoid first converting long decimals to Python float.
    """
    a_mpf = _to_finite_mpf(a, "a")
    H_mpf = _to_finite_mpf(H, "H")
    Hdot_mpf = _to_finite_mpf(Hdot, "Hdot")
    a_T_mpf = _to_finite_mpf(a_T, "a_T")
    if a_mpf <= 0:
        raise ValueError("a must be strictly positive")
    if a_T_mpf <= 0:
        raise ValueError("a_T must be strictly positive")
    return (a_T_mpf * a_T_mpf / a_mpf, -H_mpf, -Hdot_mpf)


def run_closure_regression(
    *,
    full_cycles: int = DEFAULT_FULL_CYCLES,
    precision_digits: int = DEFAULT_PRECISION_DIGITS,
    a_T: str = DEFAULT_A_T,
    a0: str = DEFAULT_A0,
    H0: str = DEFAULT_H0,
    Hdot0: str = DEFAULT_HDOT0,
) -> dict[str, object]:
    """Run two pole swaps per cycle and return a JSON-friendly closure report.

    The calculation is performed inside an isolated mpmath work-precision
    context; it does not permanently modify the caller's global precision.
    """
    if isinstance(full_cycles, bool) or not isinstance(full_cycles, int) or full_cycles < 0:
        raise ValueError("full_cycles must be an integer >= 0")
    if (
        isinstance(precision_digits, bool)
        or not isinstance(precision_digits, int)
        or precision_digits < 15
    ):
        raise ValueError("precision_digits must be an integer >= 15")

    with mp.workdps(precision_digits):
        turn = _to_finite_mpf(a_T, "a_T")
        initial = (
            _to_finite_mpf(a0, "a0"),
            _to_finite_mpf(H0, "H0"),
            _to_finite_mpf(Hdot0, "Hdot0"),
        )
        if turn <= 0:
            raise ValueError("a_T must be strictly positive")
        if initial[0] <= 0:
            raise ValueError("a0 must be strictly positive")

        state = initial
        swaps = 0
        for _ in range(full_cycles):
            state = pole_swap(*state, turn)
            swaps += 1
            state = pole_swap(*state, turn)
            swaps += 1

        deltas = tuple(end - start for start, end in zip(initial, state, strict=True))
        exact = tuple(delta == 0 for delta in deltas)

        def fmt(value: mp.mpf) -> str:
            return mp.nstr(value, n=precision_digits)

        names = ("a", "H", "Hdot")
        initial_out = {name: fmt(value) for name, value in zip(names, initial, strict=True)}
        final_out = {name: fmt(value) for name, value in zip(names, state, strict=True)}
        delta_out = {name: fmt(value) for name, value in zip(names, deltas, strict=True)}
        exact_out = {name: value for name, value in zip(names, exact, strict=True)}

        return {
            "schema_version": 1,
            "status": "reconstructed_numerical_operator_regression",
            "physical_bounce_claim": False,
            "source_status": "reconstruction_not_recovered_original_source",
            "precision_digits": precision_digits,
            "full_cycles": full_cycles,
            "pole_swaps": swaps,
            "turnaround_scale_factor": fmt(turn),
            "initial": initial_out,
            "final": final_out,
            "deltas": delta_out,
            "exact_zero_deltas": exact_out,
            "exact_closure": all(exact),
        }


__all__ = [
    "DEFAULT_PRECISION_DIGITS",
    "DEFAULT_FULL_CYCLES",
    "DEFAULT_A_T",
    "DEFAULT_A0",
    "DEFAULT_H0",
    "DEFAULT_HDOT0",
    "pole_swap",
    "run_closure_regression",
]
