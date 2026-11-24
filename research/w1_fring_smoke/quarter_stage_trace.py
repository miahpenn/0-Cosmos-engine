"""Trace metric curvature at stages approaching the first matched W1 endpoint."""
import hashlib
import json
from pathlib import Path

import numpy as np

from research.w1_fring_smoke import areal_trace, origin_trace, quarter_trace

RESULTS = Path("results/research/w1_fring_quarter_trace")
END = 0.00075


def metric_offsets(metric):
    x2_offset = (metric.X - 1) * (metric.X + 1)
    q = (x2_offset - (metric.a - 1)) / metric.a
    u = ((metric.b - 1) - x2_offset) / (metric.X * (np.sqrt(metric.b) + metric.X))
    return q, u


def run_stages(divisor, reference):
    stages = {}
    grids = []
    accepted = []

    def record_stage(name, grid, metric):
        if not grids:
            grids.append(grid)
        stages[name] = metric_offsets(metric)

    def record_accepted(state, vacuum, row):
        accepted.append(areal_trace.observe(state, vacuum, row))

    meta = origin_trace.run(
        0.025, "A", divisor, end_time=END, accepted_observer=record_accepted,
        stage_observer=record_stage, stage_step=divisor - 1,
    )
    if accepted[-1] != reference:
        raise ValueError("matched endpoint does not reproduce the recorded run")
    return grids[0], stages, meta["attempts"]


def main():
    sources = {
        1: quarter_trace.PREVIOUS / "A-0.025.json",
        2: quarter_trace.PREVIOUS / "A-0.025-half.json",
        4: RESULTS / "quarter.json",
    }
    histories = {}
    attempts = {}
    for divisor, path in sources.items():
        report = json.loads(path.read_text())
        rows = report["quarter_rows" if divisor == 4 else "rows"]
        grid, histories[divisor], attempts[divisor] = run_stages(divisor, rows[divisor])

    comparisons = {}
    # The explicit predictor and subsequent stages target the common endpoint.
    # The lapse refresh is at the start of each run's last substep.
    for stage in histories[1]:
        if stage == "after_initial_lapse_refresh":
            continue
        comparisons[stage] = {}
        for width in areal_trace.WIDTHS:
            pairs = {}
            for left, right in ((1, 2), (2, 4)):
                ql, ul = histories[left][stage]
                qr, ur = histories[right][stage]
                offsets = ((ql, ul), (ql, ur), (qr, ul), (qr, ur))
                ll, lr, rl, rr = [
                    areal_trace.areal_scalar(grid, q, u, width)[:6] for q, u in offsets
                ]
                radial = ((ll - rl) + (lr - rr)) / 2
                angular = ((ll - lr) + (rl - rr)) / 2
                np.testing.assert_allclose(radial + angular, ll - rr, atol=1e-19, rtol=0)
                pairs[f"{left}-{right}"] = {
                    "total": (ll - rr).tolist(), "radial_metric": radial.tolist(),
                    "areal_radius": angular.tolist(),
                }
            comparisons[stage][width] = pairs

    paths = [Path(__file__), Path(origin_trace.__file__), Path(areal_trace.__file__), *sources.values()]
    hashes = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    print(json.dumps({"t": END, "attempts": attempts, "stages": comparisons, "sha256": hashes},
                     indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
