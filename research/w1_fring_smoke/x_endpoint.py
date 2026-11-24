"""Run and combine the three fine-grid X ledgers at t=0.015."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from engine.v55_pirk_adapter import vendor_modules
from research.w1_fring_smoke import areal_trace, metric_balance, origin_trace, x_update_ledger as ledger

END = 0.015
DR = 0.025


def provenance():
    modules = (origin_trace, areal_trace, metric_balance, ledger)
    paths = [Path(__file__), *(Path(module.__file__) for module in modules)]
    references = Path("results/research/w1_fring_areal_trace")
    paths.extend(references / f"A-0.025{suffix}.json" for suffix in ("", "-half"))
    paths.append(Path("results/research/w1_fring_quarter_trace/quarter.json"))
    return {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}


def run(divisor, directory):
    hashes = provenance()
    _, endpoint, totals, steps, attempts, baseline, readback = ledger.collect(divisor, end_time=END)
    if hashes != provenance():
        raise ValueError("replay sources changed during execution")
    metadata = {"attempts": attempts, "steps": steps, "readback": readback, "sha256": hashes}
    np.savez(directory / f"{divisor}.npz", endpoint=endpoint, baseline=baseline,
             metadata=json.dumps(metadata, allow_nan=False), **totals)


def combine(directory):
    hashes = provenance()
    grid_ops, _, _ = vendor_modules()
    grid = grid_ops.SphericalCellGrid(round(60 / DR), 60.)
    endpoints, runs, totals = {}, {}, {}
    for divisor in ledger.LABELS:
        with np.load(directory / f"{divisor}.npz", allow_pickle=False) as data:
            runs[divisor] = json.loads(data["metadata"].item())
            if runs[divisor].pop("sha256") != hashes or runs[divisor]["attempts"] != 20 * divisor:
                raise ValueError("replay metadata mismatch")
            endpoints[divisor] = data["endpoint"]
            totals[divisor] = {name: data[name] for name in ledger.COMPONENTS}
            if divisor == 1:
                baseline = data["baseline"]
            else:
                np.testing.assert_array_equal(baseline, data["baseline"])
    pairs = {}
    for left, right in ((1, 2), (2, 4)):
        lv, rv = endpoints[left], endpoints[right]
        delta = {key: totals[left][key] - totals[right][key] for key in ledger.COMPONENTS}
        np.testing.assert_allclose(lv[2] - rv[2], sum(delta.values()), atol=4e-20, rtol=0)
        a, b = [(l + r) / 2 for l, r in zip(lv[:2], rv[:2])]
        widths = {}
        for width in areal_trace.WIDTHS:
            base = metric_balance.curvature(grid, (a, b, rv[2]), width)
            x_only = metric_balance.curvature(grid, (a, b, lv[2]), width) - base
            full = metric_balance.curvature(grid, lv, width) - metric_balance.curvature(grid, rv, width)
            components = {
                name: metric_balance.curvature(grid, (a, b, rv[2] + dx), width) - base
                for name, dx in delta.items()
            }
            rate_dx = sum(delta[name] for name in (*ledger.FACTORS, "increment_arithmetic"))
            write_dx = delta["earlier_corrector_writes"] + delta["last_corrector_write"]
            rate = metric_balance.curvature(grid, (a, b, rv[2] + rate_dx), width) - base
            write = metric_balance.curvature(grid, (a, b, rv[2] + write_dx), width) - base
            closure = x_only - rate - write
            np.testing.assert_allclose(closure, 0, atol=5e-15, rtol=0)
            widths[width] = {
                "full_metric": metric_balance.numbers(full), "X_only": metric_balance.numbers(x_only),
                "remaining_metric": metric_balance.numbers(full - x_only),
                "rate": metric_balance.numbers(rate), "write": metric_balance.numbers(write),
                "components": {key: metric_balance.numbers(value) for key, value in components.items()},
                "closure": metric_balance.numbers(closure),
            }
        pairs[f"{left}-{right}"] = {"X_differences": ledger.samples(delta), "curvature": widths}
    print(json.dumps({"t": END, "dr": DR, "source": origin_trace.SOURCE, "vendor": origin_trace.VENDOR,
                      "input_sha256": origin_trace.DIGEST, "sha256": hashes,
                      "runs": runs, "pairs": pairs}, indent=2, allow_nan=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--divisor", type=int, choices=tuple(ledger.LABELS))
    args = parser.parse_args()
    if np.finfo(ledger.WIDE).nmant < 63:
        raise RuntimeError("expected at least 64 significand bits")
    if args.divisor is None:
        combine(args.directory)
    else:
        run(args.divisor, args.directory)


if __name__ == "__main__":
    main()
