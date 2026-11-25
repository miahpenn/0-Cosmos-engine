"""Carry stable primary-L2 Ricci arithmetic through the pinned W1 window."""
import argparse
from copy import deepcopy
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np

from research.w1_fring_smoke import live_lapse, stage_momentum as replay
from research.w1_fring_smoke.first_step_rate_precision import FirstStepComplete
from research.w1_fring_smoke.live_ricci_precision import ARMS, StableRicci, stats

ROOT = Path("results/research")
STAGES = ("explicit_predictor", "primary_before_boundary", "primary_after_lapse",
          "corrector_before_boundary", "accepted_end")


def initial_state(case, go, vacuum, meta):
    captured = []

    def take_checkpoint(kernel, state, dt, watch):
        captured.append((state, dt))
        raise FirstStepComplete

    with patch.object(replay, "advance", take_checkpoint):
        try:
            replay.fine_replay(float(case), go, vacuum, meta)
        except FirstStepComplete:
            pass
    (checkpoint,) = captured
    return checkpoint


def reduced_step(r, step):
    stages = step["stages"]
    results = {}
    for name in STAGES:
        row = stages[name]
        change = row["delta_from_refreshed_start"]
        profiles = {
            "M": np.asarray(row["rows"]["M"]),
            "delta_M": np.asarray(change["M"]),
            "delta_shear": sum(np.asarray(change[key]) for key in replay.PARTS[:2]),
            "delta_K": np.asarray(change["K_gradient"]),
            "delta_metric_coupling": np.asarray(change["metric_coupling"]),
            "delta_matter": sum(np.asarray(change[key]) for key in replay.PARTS[4:]),
        }
        results[name] = {key: {**stats(r, values), "signed": values.tolist()}
                         for key, values in profiles.items()}
        closure = profiles["delta_M"] - sum(v for k, v in profiles.items() if k.startswith("delta_") and k != "delta_M")
        results[name]["momentum_closure_max"] = float(np.max(abs(closure)))
    return {"t_start": step["t_start"], "t_end": step["t_end"], "stages": results}


def run(state, dt, observer, steps, first_reference=None, show_progress=False):
    r = state.grid.centers[:replay.LOCAL_CELLS]
    rows = []
    for index in range(steps):
        state, ledger = replay.advance(observer, state, dt, True)
        if index == 0 and first_reference is not None:
            for name, stage in ledger["stages"].items():
                expected = {key: value for key, value in first_reference["stages"][name].items()
                            if key != "M_statistics"}
                if stage != expected:
                    raise ValueError("first step differs from the prior branch")
            if live_lapse.state_digest(state) != first_reference["end_sha256"]:
                raise ValueError("first-step endpoint digest differs")
        rows.append(reduced_step(r, ledger))
        if show_progress and (index+1) % 10 == 0:
            print(f"dr={state.grid.dr:g} accepted {index+1}/{steps}, t={state.t:.6f}",flush=True)
    return {"dt": dt, "rows": rows, "constraints": live_lapse.end_constraints(observer.vacuum,state),
            "end_sha256": live_lapse.state_digest(state)}


def trajectory(case, state, dt, vacuum, first):
    r = state.grid.centers[:replay.LOCAL_CELLS]
    archived = ROOT/"w1_fring_momentum_history"/f"{case}.json"
    baseline = json.loads(archived.read_text())
    if baseline["dt"] != dt:
        raise ValueError("baseline timestep differs")
    np.testing.assert_array_equal(baseline["r"], r[:len(baseline["r"])])
    observer = StableRicci(vacuum, ARMS["stable_metric"])
    result = run(state, dt, observer, len(baseline["rows"]), first["stable_metric"], True)
    comparisons = []
    for old, new in zip(baseline["rows"], result["rows"], strict=True):
        if old["t"] != new["t_end"]:
            raise ValueError("expected matched accepted times")
        reference = np.asarray(old["M"]["signed"])
        updated = np.asarray(new["stages"]["accepted_end"]["M"]["signed"][:len(reference)])
        comparisons.append({
            "t": old["t"], "production": {**stats(r[:len(reference)],reference), "signed": reference.tolist()},
            "stable_minus_production": {**stats(r[:len(reference)],updated-reference),
                                        "signed": (updated-reference).tolist()},
        })
    return {"stable": result, "comparisons": comparisons, "history_sha256": replay.sha(archived),
            "L2_evaluations": observer.calls}


def timing(state, dt, vacuum, first):
    result = {}
    for arm in ("production", "stable_metric"):
        for divisor in (1,2):
            observer = (replay.Observe(vacuum) if arm == "production"
                        else StableRicci(vacuum, ARMS[arm]))
            result[f"{arm}_dt_over_{divisor}"] = run(
                deepcopy(state), dt/divisor, observer, divisor, first[arm] if divisor==1 else None)
    times = {row["rows"][-1]["t_end"] for row in result.values()}
    if times != {dt}:
        raise ValueError("timing arms differ in endpoint")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("0.025","0.0125"))
    parser.add_argument("mode", choices=("trajectory","timing"))
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    meta = replay.provenance()
    go,vacuum,_ = replay.pins.ad.vendor_modules()
    if vacuum.KO_EPSILON != 0:
        raise ValueError("expected KO off")
    first_path = ROOT/"w1_fring_live_ricci_precision"/f"{args.case}.json"
    first = json.loads(first_path.read_text())
    state,dt = initial_state(args.case,go,vacuum,meta)
    digest = live_lapse.state_digest(state)
    if digest != first["start_sha256"]:
        raise ValueError("initializer differs from previous first-step arms")
    r = state.grid.centers[:replay.LOCAL_CELLS].tolist()
    if args.mode == "trajectory":
        result = trajectory(args.case,state,dt,vacuum,first)
    else:
        result = timing(state,dt,vacuum,first)
    report = {"case": args.case, "mode": args.mode, "r": r, "start_sha256": digest,
              "engine": replay.pins.SOURCE, "vendor": replay.pins.VENDOR,
              "input_sha256": first["input_sha256"], "script_sha256": replay.sha(Path(__file__)),
              "first_sha256": replay.sha(first_path), "result": result}
    args.output.write_text(json.dumps(report,allow_nan=False,separators=(",", ":"))+"\n")


if __name__ == "__main__":
    main()
