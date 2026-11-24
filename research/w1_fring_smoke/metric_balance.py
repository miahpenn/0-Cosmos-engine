"""Audit the stored and widened metric updates at the first shared W1 endpoint."""
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

from engine.native_geometry.evolution import native_cmc_explicit_metric_rhs
from engine.production_kernel import NativeSphericalCellGrid
from research.w1_fring_smoke import areal_trace, origin_trace, quarter_stage_trace

METRIC = ("a", "b", "X")
STAGES = {
    "predictor": "explicit_predictor",
    "corrector": "corrected_before_boundary",
}
WIDE = np.longdouble


def metric_arrays(metric, dtype):
    return tuple(np.asarray(getattr(metric, name), dtype=dtype) for name in METRIC)


def curvature(grid, arrays, width):
    a, b, X = arrays
    x2_offset = (X - 1) * (X + 1)
    q = (x2_offset - (a - 1)) / a
    u = ((b - 1) - x2_offset) / (X * (np.sqrt(b) + X))
    return areal_trace.areal_scalar(grid, q, u, width)[:6]


def split_fields(grid, left, right, width):
    # Average each field's marginal contribution over all six swap orders.
    values = {
        choice: curvature(grid, tuple(pair[side] for pair, side in zip(zip(left, right), choice)), width)
        for choice in itertools.product((0, 1), repeat=3)
    }
    shares = [np.zeros(6, dtype=WIDE) for _ in METRIC]
    for order in itertools.permutations(range(3)):
        choice = [1, 1, 1]
        for i in order:
            before = values[tuple(choice)]
            choice[i] = 0
            shares[i] += values[tuple(choice)] - before
    result = {name: share/6 for name, share in zip(METRIC, shares)}
    np.testing.assert_allclose(sum(result.values()), values[0, 0, 0] - values[1, 1, 1],
                               atol=2e-19, rtol=0)
    return result


def collect(divisor):
    captured = {}
    grids = []

    def record(name, grid, geometry):
        captured[name] = geometry
        if not grids:
            grids.append(grid)

    trace = origin_trace.run(0.025, "A", divisor, end_time=quarter_stage_trace.END,
                             stage_observer=record, stage_step=divisor-1)
    grid = grids[0]
    native = NativeSphericalCellGrid(grid.n, grid.r_max)
    initial = captured["after_initial_lapse_refresh"]
    primary = captured["primary_predictor_after_lapse"]
    initial_rates = native_cmc_explicit_metric_rhs(native, initial)
    primary_rates = native_cmc_explicit_metric_rhs(native, primary)
    dt = trace["dt"]
    initial_arrays = metric_arrays(initial, np.float64)
    increments = {
        "predictor": tuple(dt * initial_rates[name] for name in METRIC),
        "corrector": tuple(0.5 * dt * (initial_rates[name] + primary_rates[name]) for name in METRIC),
    }
    states = {}
    for phase, name in STAGES.items():
        stored = metric_arrays(captured[name], np.float64)
        for actual, y0, increment in zip(stored, initial_arrays, increments[phase]):
            np.testing.assert_array_equal(actual, y0 + increment)
        prewrite = tuple(y0.astype(WIDE) + delta.astype(WIDE)
                         for y0, delta in zip(initial_arrays, increments[phase]))
        states[phase] = {
            "binary64": stored,
            "wide_stored": metric_arrays(captured[name], WIDE),
            "wide_prewrite": prewrite,
        }

    rate_readback = {}
    for name, metric in (("start", initial), ("primary", primary)):
        shear = metric.alpha * metric.Aa
        expansion = metric.alpha * metric.K / 3
        rate_readback[name] = {
            "alpha_Aa": shear[:6].tolist(), "alpha_K_over_3": expansion[:6].tolist(),
            "dlog_inverse_radial": (2 * (shear + expansion))[:6].tolist(),
            "dlog_areal_over_r": (0.5 * shear - expansion)[:6].tolist(),
        }
    return grid, states, rate_readback, trace["attempts"]


def numbers(values):
    return np.asarray(values, dtype=float).tolist()


def main():
    if np.finfo(WIDE).eps >= np.finfo(float).eps:
        raise RuntimeError("extended precision is needed for the pre-write replay")
    reference_path = quarter_stage_trace.RESULTS / "stages.json"
    reference = json.loads(reference_path.read_text())
    replays = {divisor: collect(divisor) for divisor in (1, 2, 4)}
    ledger = {}
    for phase, stage in STAGES.items():
        ledger[phase] = {}
        for width in areal_trace.WIDTHS:
            pair_results = {}
            for left, right in ((1, 2), (2, 4)):
                grid, lstate, _, _ = replays[left]
                _, rstate, _, _ = replays[right]
                lp, rp = lstate[phase], rstate[phase]
                signed = {key: curvature(grid, lp[key], width) - curvature(grid, rp[key], width)
                          for key in lp}
                pair = f"{left}-{right}"
                np.testing.assert_array_equal(
                    signed["binary64"], reference["stages"][stage][str(width)][pair]["total"]
                )
                rounded = signed["wide_stored"] - signed["wide_prewrite"]
                shares = split_fields(grid, lp["wide_stored"], rp["wide_stored"], width)
                write_shares = split_fields(grid, lp["wide_stored"], lp["wide_prewrite"], width)
                right_shares = split_fields(grid, rp["wide_stored"], rp["wide_prewrite"], width)
                pair_results[pair] = {
                    **{key: numbers(value) for key, value in signed.items()},
                    "write_rounding": numbers(rounded),
                    "stored_field_shares": {key: numbers(value) for key, value in shares.items()},
                    "write_rounding_field_shares": {
                        name: numbers(write_shares[name] - right_shares[name]) for name in METRIC
                    },
                }
            ledger[phase][width] = pair_results
    provenance = [Path(__file__), Path(origin_trace.__file__), Path(areal_trace.__file__),
                  reference_path]
    print(json.dumps({
        "time": quarter_stage_trace.END, "dr": 0.025,
        "widened_significand_bits": np.finfo(WIDE).nmant + 1,
        "source": origin_trace.SOURCE, "vendor": origin_trace.VENDOR,
        "input_sha256": origin_trace.DIGEST,
        "sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in provenance},
        "attempts": {key: v[3] for key, v in replays.items()},
        "rates": {key: v[2] for key, v in replays.items()}, "curvature": ledger,
    }, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
