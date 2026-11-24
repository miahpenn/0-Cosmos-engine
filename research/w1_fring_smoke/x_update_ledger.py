"""Account for the W1 endpoint X differences using every stored metric write."""
import hashlib
import json
from pathlib import Path

import numpy as np

from engine.native_geometry.evolution import native_cmc_explicit_metric_rhs
from engine.production_kernel import NativeSphericalCellGrid
from research.w1_fring_smoke import areal_trace, metric_balance, origin_trace, quarter_stage_trace, quarter_trace

WIDE = metric_balance.WIDE
FACTORS = ("X", "alpha", "K")
COMPONENTS = (*FACTORS, "increment_arithmetic", "earlier_corrector_writes", "last_corrector_write")
REFERENCE = Path("results/research/w1_fring_metric_balance/balance.json")
LABELS = {1: "nominal", 2: "half", 4: "quarter"}


def factors(metric):
    return tuple(np.asarray(getattr(metric, key), dtype=WIDE) for key in FACTORS)


def product_components(profile, baseline):
    changes = tuple(current - initial for current, initial in zip(profile, baseline))
    components = {}
    for i, name in enumerate(FACTORS):
        j, k = [idx for idx in range(3) if idx != i]
        # Split each two-factor interaction equally and each three-factor interaction in thirds.
        components[name] = changes[i] / 3 * (
            baseline[j] * baseline[k]
            + (changes[j] * baseline[k] + baseline[j] * changes[k]) / 2
            + changes[j] * changes[k] / 3
        )
    return components


def samples(values):
    return {name: metric_balance.numbers(value[:6]) for name, value in values.items()}


def collect(divisor, *, end_time=quarter_stage_trace.END):
    step_stages = []
    accepted_rows = []
    grids = []

    def record(stage, grid, metric):
        if not grids:
            grids.append(grid)
        if stage == "after_initial_lapse_refresh":
            step_stages.append({})
        step_stages[-1][stage] = metric

    def accepted(state, vacuum, row):
        accepted_rows.append(areal_trace.observe(state, vacuum, row))

    steps = round(end_time / (0.03 * 0.025 / divisor))
    replay = origin_trace.run(
        0.025, "A", divisor, end_time=end_time,
        stage_observer=record, accepted_observer=accepted, stage_step=None,
    )
    suffix = "-half" if divisor == 2 else ""
    path = (quarter_stage_trace.RESULTS / "quarter.json" if divisor == 4
            else quarter_trace.PREVIOUS / f"A-0.025{suffix}.json")
    reference = json.loads(path.read_text())["quarter_rows" if divisor == 4 else "rows"]
    prefix_length = min(steps + 1, len(reference))
    if accepted_rows[:prefix_length] != reference[:prefix_length] or len(step_stages) != steps:
        raise ValueError("replay disagrees with recorded accepted states")

    grid = grids[0]
    native = NativeSphericalCellGrid(grid.n, grid.r_max)
    dt = replay["dt"]
    dt_half = WIDE(dt) / 2
    baseline = factors(step_stages[0]["after_initial_lapse_refresh"])
    baseline_rate = np.prod(baseline, axis=0) / 3
    totals = {name: np.zeros(grid.n, dtype=WIDE) for name in COMPONENTS}
    per_step = []
    previous = baseline[0]
    for i, stages in enumerate(step_stages):
        start = stages["after_initial_lapse_refresh"]
        primary = stages["primary_predictor_after_lapse"]
        predicted = stages["explicit_predictor"].X
        corrected = stages["corrected_before_boundary"].X
        np.testing.assert_array_equal(start.X, previous)
        np.testing.assert_array_equal(corrected, stages["corrected_after_boundary"].X)
        native_start = native_cmc_explicit_metric_rhs(native, start)["X"]
        native_primary = native_cmc_explicit_metric_rhs(native, primary)["X"]
        euler = dt * native_start
        trapezoid = (dt / 2) * (native_start + native_primary)
        np.testing.assert_array_equal(start.X + euler, predicted)
        np.testing.assert_array_equal(start.X + trapezoid, corrected)

        start_factors, primary_factors = factors(start), factors(primary)
        start_parts = product_components(start_factors, baseline)
        primary_parts = product_components(primary_factors, baseline)
        parts = {name: dt_half * (start_parts[name] + primary_parts[name]) for name in FACTORS}
        wide_rate_increment = dt_half * (
            np.prod(start_factors, axis=0) + np.prod(primary_factors, axis=0)
        ) / 3
        parts["increment_arithmetic"] = trapezoid.astype(WIDE) - wide_rate_increment
        parts["corrector_write"] = corrected.astype(WIDE) - start.X.astype(WIDE) - trapezoid.astype(WIDE)
        predictor_write = predicted.astype(WIDE) - start.X.astype(WIDE) - euler.astype(WIDE)
        write_key = "last_corrector_write" if i == steps - 1 else "earlier_corrector_writes"
        totals[write_key] += parts["corrector_write"]
        for name in (*FACTORS, "increment_arithmetic"):
            totals[name] += parts[name]
        per_step.append({"end_t": replay["rows"][i+1]["t"], **samples(parts),
                         "predictor_write": metric_balance.numbers(predictor_write[:6])})
        previous = corrected
    endpoint = metric_balance.metric_arrays(step_stages[-1]["corrected_after_boundary"], WIDE)
    reference_integral = WIDE(dt) * steps * baseline_rate
    np.testing.assert_allclose(
        endpoint[2] - baseline[0], reference_integral + sum(totals.values()), atol=3e-20, rtol=0,
    )
    return grid, endpoint, totals, per_step, replay["attempts"], baseline, accepted_rows[-1]


def main():
    if np.finfo(WIDE).nmant < 63:
        raise RuntimeError("expected at least 64 significand bits for the write ledger")
    replay = {divisor: collect(divisor) for divisor in LABELS}
    old = json.loads(REFERENCE.read_text())
    pairs = {}
    for left, right in ((1, 2), (2, 4)):
        grid, lv, lt, _, _, lb, _ = replay[left]
        _, rv, rt, _, _, rb, _ = replay[right]
        np.testing.assert_array_equal(lb, rb)
        delta = {key: lt[key] - rt[key] for key in COMPONENTS}
        np.testing.assert_allclose(lv[2] - rv[2], sum(delta.values()), atol=4e-20, rtol=0)
        widths = {}
        for width in areal_trace.WIDTHS:
            # Use the arithmetic mean of the two frozen endpoint a,b backgrounds.
            # It supplies one common metric on which to read the effect of each X component.
            a, b = [(l + r) / 2 for l, r in zip(lv[:2], rv[:2])]
            base = metric_balance.curvature(grid, (a, b, rv[2]), width)
            x_response = metric_balance.curvature(grid, (a, b, lv[2]), width) - base
            full = metric_balance.curvature(grid, lv, width) - metric_balance.curvature(grid, rv, width)
            values = {
                key: metric_balance.curvature(grid, (a, b, rv[2] + value), width) - base
                for key, value in delta.items()
            }
            closure = x_response - sum(values.values())
            np.testing.assert_array_equal(
                metric_balance.numbers(full),
                old["curvature"]["corrector"][str(width)][f"{left}-{right}"]["wide_stored"],
            )
            widths[width] = {
                "full_metric": metric_balance.numbers(full),
                "X_only": metric_balance.numbers(x_response),
                "contributions": {key: metric_balance.numbers(value) for key, value in values.items()},
                "closure": metric_balance.numbers(closure),
            }
        pairs[f"{left}-{right}"] = {"X_differences": samples(delta), "curvature": widths}
    paths = [Path(__file__), Path(origin_trace.__file__), Path(metric_balance.__file__),
             Path(areal_trace.__file__), REFERENCE]
    metadata = {
        "dr": 0.025, "t": quarter_stage_trace.END,
        "source": origin_trace.SOURCE, "vendor": origin_trace.VENDOR,
        "input_sha256": origin_trace.DIGEST, "significand_bits": np.finfo(WIDE).nmant + 1,
        "sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
    }
    print(json.dumps({
        "metadata": metadata,
        "runs": {LABELS[i]: {"attempts": r[4], "steps": r[3]} for i, r in replay.items()},
        "pairs": pairs,
    }, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
