"""Measure the determinant and tracefree parts of the W1 conformal metric updates."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from engine.native_geometry.evolution import native_cmc_explicit_metric_rhs
from engine.production_kernel import NativeSphericalCellGrid
from engine.v55_pirk_adapter import vendor_modules
from research.w1_fring_smoke import areal_trace, metric_balance, origin_trace, quarter_stage_trace, quarter_trace, x_endpoint

WIDE = metric_balance.WIDE
FIELDS = ("a", "b")
TERMS = ("tangent", "integration", "increment_arithmetic", "write")
LABELS = ("determinant", "tracefree")
REFERENCE = Path("results/research/w1_fring_x_endpoint/ledger.json")


def readback(values):
    return metric_balance.numbers(values[:6])


def coords(log_fields):
    a, b = log_fields
    return np.stack((a + 2*b, (b - a)/3))


def source_hashes():
    paths = [Path(__file__), Path(origin_trace.__file__), Path(metric_balance.__file__),
             Path(areal_trace.__file__), Path(x_endpoint.__file__), REFERENCE]
    paths.extend(quarter_trace.PREVIOUS / f"A-0.025{suffix}.json" for suffix in ("", "-half"))
    paths.append(quarter_stage_trace.RESULTS / "quarter.json")
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def run(divisor, directory):
    hashes = source_hashes()
    reference = json.loads(REFERENCE.read_text())["runs"][str(divisor)]["readback"]
    archived = (quarter_stage_trace.RESULTS / "quarter.json" if divisor == 4 else
                quarter_trace.PREVIOUS / ("A-0.025-half.json" if divisor == 2 else "A-0.025.json"))
    prior = json.loads(archived.read_text())["quarter_rows" if divisor == 4 else "rows"]
    n = round(60/x_endpoint.DR)
    native = NativeSphericalCellGrid(n, 60.)
    dt = .03*x_endpoint.DR/divisor
    steps = 20*divisor
    totals = {name: np.zeros((2, n), dtype=WIDE) for name in TERMS}
    history, seen = [], []
    stages = {}
    previous = None
    initial = None
    endpoint = None

    def record(name, grid, geometry):
        nonlocal previous, initial, endpoint
        if name == "after_initial_lapse_refresh":
            stages.clear()
        stages[name] = geometry
        if name != "corrected_after_boundary":
            return
        start = stages["after_initial_lapse_refresh"]
        primary = stages["primary_predictor_after_lapse"]
        predictor = stages["explicit_predictor"]
        stored = stages["corrected_before_boundary"]
        start_fields = np.stack([getattr(start, key) for key in FIELDS])
        if initial is None:
            initial = start_fields.astype(WIDE)
        else:
            np.testing.assert_array_equal(start_fields, previous)
        end_fields = np.stack([getattr(stored, key) for key in FIELDS])
        endpoint = metric_balance.metric_arrays(geometry, WIDE)
        np.testing.assert_array_equal(end_fields, endpoint[:2])
        rates0 = native_cmc_explicit_metric_rhs(native, start)
        rates1 = native_cmc_explicit_metric_rhs(native, primary)
        predicted = [getattr(start, key) + dt*rates0[key] for key in FIELDS]
        np.testing.assert_array_equal(predicted, [getattr(predictor, key) for key in FIELDS])
        increment = np.stack([(dt/2)*(rates0[key]+rates1[key]) for key in FIELDS])
        np.testing.assert_array_equal(start_fields + increment, end_fields)

        wide_start = start_fields.astype(WIDE)
        primary_fields = np.stack([getattr(primary, key) for key in FIELDS]).astype(WIDE)
        f0 = start.alpha.astype(WIDE) * start.Aa.astype(WIDE)
        f1 = primary.alpha.astype(WIDE) * primary.Aa.astype(WIDE)
        rates = np.array([[-2], [1]], dtype=WIDE)
        tangent = (WIDE(dt)/2) * rates * (f0+f1)
        ideal_increment = (WIDE(dt)/2)*rates*(wide_start*f0 + primary_fields*f1)
        ideal_log = np.log1p(ideal_increment/wide_start)
        delta_increment = increment.astype(WIDE) - ideal_increment
        write = end_fields.astype(WIDE) - wide_start - increment.astype(WIDE)
        log_parts = {
            "tangent": tangent, "integration": ideal_log-tangent,
            "increment_arithmetic": np.log1p(delta_increment/(wide_start+ideal_increment)),
            "write": np.log1p(write/(wide_start+increment.astype(WIDE))),
        }
        for key, value in log_parts.items():
            totals[key] += coords(value)
        history.append({key: [readback(v) for v in coords(value)] for key, value in log_parts.items()})
        previous = end_fields

    def accepted(state, vacuum, row):
        index = len(seen)
        if index < len(prior) or index == steps:
            actual = areal_trace.observe(state, vacuum, row)
            expected = prior[index] if index < len(prior) else reference
            if actual != expected:
                raise ValueError(f"archived readback disagrees at step {index}")
        seen.append(state.t)

    result = origin_trace.run(x_endpoint.DR, "A", divisor, end_time=x_endpoint.END,
                              accepted_observer=accepted, stage_observer=record, stage_step=None)
    if len(history) != steps or source_hashes() != hashes:
        raise ValueError("unexpected number of steps or changed inputs")
    delta = coords(np.log1p((endpoint[:2]-initial)/initial))
    np.testing.assert_allclose(delta, sum(totals.values()), atol=4e-20, rtol=0)
    metadata = {"attempts": result["attempts"], "time": seen[-1], "sha256": hashes, "steps": history}
    np.savez(directory/f"{divisor}.npz", endpoint=endpoint, initial=initial,
             metadata=json.dumps(metadata, allow_nan=False), **totals)


def combine(directory):
    hashes = source_hashes()
    prior = json.loads(REFERENCE.read_text())
    grid = vendor_modules()[0].SphericalCellGrid(2400, 60.)
    endpoints, totals, runs = {}, {}, {}
    for divisor in (1, 2, 4):
        with np.load(directory/f"{divisor}.npz", allow_pickle=False) as data:
            runs[divisor] = json.loads(data["metadata"].item())
            if runs[divisor].pop("sha256") != hashes or runs[divisor]["attempts"] != 20*divisor:
                raise ValueError("inconsistent replay provenance")
            endpoints[divisor] = data["endpoint"]
            totals[divisor] = {name: data[name] for name in TERMS}
            if divisor == 1:
                initial = data["initial"]
            else:
                np.testing.assert_array_equal(initial, data["initial"])
            drift = sum(totals[divisor].values())
            runs[divisor]["drift"] = [readback(v) for v in drift]
            runs[divisor]["centre_max_by_term"] = {
                name: metric_balance.numbers(np.max(abs(term[:, :6]), axis=1))
                for name, term in totals[divisor].items()
            }
            runs[divisor]["centre_max_total"] = metric_balance.numbers(np.max(abs(drift[:, :6]), axis=1))
            runs[divisor]["grid_max_total"] = metric_balance.numbers(np.max(abs(drift), axis=1))

    pairs = {}
    for left, right in ((1, 2), (2, 4)):
        lv, rv = endpoints[left], endpoints[right]
        x = (lv[2]+rv[2])/2
        delta = {name: totals[left][name]-totals[right][name] for name in TERMS}
        np.testing.assert_allclose(coords(np.log1p((lv[:2]-rv[:2])/rv[:2])),
                                   sum(delta.values()), atol=6e-20, rtol=0)
        pair = f"{left}-{right}"
        widths = {}
        for width in areal_trace.WIDTHS:
            base = metric_balance.curvature(grid, (*rv[:2], x), width)
            total = metric_balance.curvature(grid, (*lv[:2], x), width)-base
            def response(determinant, tracefree):
                a = rv[0]*np.exp(determinant/3 - 2*tracefree)
                b = rv[1]*np.exp(determinant/3 + tracefree)
                return metric_balance.curvature(grid, (a,b,x), width)-base

            contributions = {
                name: {
                    key: response(value[0], 0) if i == 0 else response(0, value[1])
                    for key, value in delta.items()
                }
                for i, name in enumerate(LABELS)
            }
            closure = total - sum(v for part in contributions.values() for v in part.values())
            np.testing.assert_allclose(closure, 0, atol=8e-15, rtol=0)
            x_ledger = prior["pairs"][pair]["curvature"][str(width)]
            remainder = np.asarray(x_ledger["remaining_metric"], dtype=WIDE)
            np.testing.assert_allclose(total, remainder, atol=3e-17, rtol=0)
            write = sum(contributions[label]["write"] for label in LABELS)
            combined_write = write + np.asarray(x_ledger["write"], dtype=WIDE)
            widths[width] = {
                "ab_total": readback(total), "remainder_from_x_ledger": readback(remainder),
                "contributions": {label: {key: readback(v) for key,v in part.items()}
                                  for label,part in contributions.items()},
                "closure": readback(closure), "ab_write": readback(write),
                "combined_metric_write": readback(combined_write),
                "absolute_component_sum": readback(
                    abs(np.asarray(x_ledger["write"], dtype=WIDE))
                    + sum(abs(contributions[label]["write"]) for label in LABELS)
                ),
            }
        pairs[pair] = {"log_delta": {key: [readback(v) for v in value] for key,value in delta.items()},
                       "curvature": widths}
    print(json.dumps({"t": x_endpoint.END, "dr": x_endpoint.DR, "coords": LABELS,
                      "sha256": hashes, "source": origin_trace.SOURCE, "vendor": origin_trace.VENDOR,
                      "input_sha256": origin_trace.DIGEST, "runs": runs, "pairs": pairs},
                     indent=2, allow_nan=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--divisor", type=int, choices=(1,2,4))
    args = parser.parse_args()
    if np.finfo(WIDE).nmant < 63:
        raise RuntimeError("at least 64 significand bits required")
    if args.divisor is None:
        combine(args.directory)
    else:
        run(args.divisor, args.directory)
