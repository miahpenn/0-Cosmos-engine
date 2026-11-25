"""Read shell signs and last-corrector Aa/K rounding from existing W1 operands."""
import argparse
from decimal import localcontext
import json
from pathlib import Path

import numpy as np

from research.w1_fring_smoke import stage_momentum as replay
from research.w1_fring_smoke.momentum_precision import derivative, wide


ARMS = ("production", "three_point_control", "fourth_order_match")
EDGES = np.arange(2, 11) * .025


def signs(r, values):
    values = np.asarray(values)
    shells = [values[(r >= a) & (r < b)] for a, b in zip(EDGES[:-1], EDGES[1:])]
    occupied = [len(v) for v in shells]
    if 0 in occupied:
        raise ValueError("expected occupied common shells")
    mean = np.asarray([
        np.average(values[(r >= a) & (r < b)], weights=r[(r >= a) & (r < b)]**2)
        for a, b in zip(EDGES[:-1], EDGES[1:])
    ])
    def flips(v):
        nonzero = v[v != 0]
        return int(np.count_nonzero(np.signbit(nonzero[1:]) != np.signbit(nonzero[:-1])))
    return {"shell_means": mean.tolist(), "across_shells": flips(mean),
            "inside_each_shell": [flips(v) for v in shells],
            "native_band": flips(np.concatenate(shells))}


def last_write(snapshot):
    dr, dt = map(float, (snapshot["dr"], snapshot["dt"]))
    n = replay.LOCAL_CELLS
    r = (np.arange(n) + .5) * dr
    go, vacuum, _ = replay.pins.ad.vendor_modules()
    grid = go.SphericalCellGrid(len(snapshot["start_K"]), 60.)
    final = {key: snapshot[f"corrector_copy_{key}"] for key in ("a", "b", "X")}
    dx = vacuum._d(grid, final["X"], 1)[:n]
    db = vacuum._d(grid, final["b"], 1)[:n]
    shear_factor = 3/r - 3*dx/final["X"][:n] + 1.5*db/final["b"][:n]
    result = {}
    profiles = {}
    with localcontext() as ctx:
        ctx.prec = 60
        for field in ("Aa", "K"):
            names = ("start_primary_l2", "corrector_copy_primary_l2",
                     "start_primary_l3", "primary_primary_l3")
            rhs = [snapshot[f"{name}_{field}"][:n] for name in names]
            initial = snapshot[f"start_{field}"][:n]
            final = snapshot[f"corrector_before_boundary_{field}"]
            increment = .5 * dt * (rhs[0] + rhs[1] + rhs[2] + rhs[3])
            np.testing.assert_array_equal(final, initial + increment)
            stored = [wide(b) - wide(a) for a, b in zip(initial, final)]
            computed = [wide(dt)*sum(map(wide, parts))/2 for parts in zip(*rhs)]
            rounded = list(map(wide, increment))
            steps = {"stored": stored, "computed": computed,
                     "increment_arithmetic": [a-b for a,b in zip(rounded, computed)],
                     "write": [a-b for a,b in zip(stored, rounded)]}
            profiles[field] = {}
            for label, values in steps.items():
                profile = [derivative(values, i, wide(dr)) for i in range(n-2)]
                if field == "K":
                    profile = [-2*v/3 for v in profile]
                else:
                    profile = [v + wide(c)*w for v,c,w in zip(profile, shear_factor, values)]
                profiles[field][label] = np.array(profile, dtype=float)
            result[field] = {label: replay.local_rms(r[:n-2], v)
                             for label,v in profiles[field].items()}
            result[field]["write_over_stored"] = result[field]["write"] / result[field]["stored"]
    result["joint_last_write"] = replay.local_rms(r[:n-2], profiles["Aa"]["write"] + profiles["K"]["write"])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("live_dir", type=Path)
    parser.add_argument("snapshot_dir", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    replay.provenance()
    result = {}
    for case in ("0.025", "0.0125"):
        live_path = args.live_dir/f"{case}.json"
        snapshot_path = args.snapshot_dir/f"{case}.npz"
        live = json.loads(live_path.read_text())
        r = np.asarray(live["r"])
        signs_by_arm = {}
        for arm in ARMS:
            stages = live[arm]["stages"]
            signs_by_arm[arm] = {
                "accepted_start": signs(r, stages["accepted_start"]["rows"]["M"]),
                "accepted_end": signs(r, stages["accepted_end"]["rows"]["M"]),
                "step_delta": signs(r, stages["accepted_end"]["delta_from_refreshed_start"]["M"]),
            }
            start = np.asarray(stages["accepted_start"]["rows"]["M"])
            end = np.asarray(stages["accepted_end"]["rows"]["M"])
            native_reversal = (start != 0) & (end != 0) & (np.signbit(start) != np.signbit(end))
            start_shell = np.asarray(signs_by_arm[arm]["accepted_start"]["shell_means"])
            end_shell = np.asarray(signs_by_arm[arm]["accepted_end"]["shell_means"])
            signs_by_arm[arm]["start_to_end_reversals"] = {
                "shell_means": ((start_shell != 0) & (end_shell != 0) &
                                (np.signbit(start_shell) != np.signbit(end_shell))).tolist(),
                "native_cells_per_shell": [
                    int(np.count_nonzero(native_reversal & (r >= a) & (r < b)))
                    for a, b in zip(EDGES[:-1], EDGES[1:])
                ],
            }
            if arm != "production":
                delta = live[arm]["changes_vs_production"]["accepted_end"]["signed"]["M"]
                signs_by_arm[arm]["minus_production"] = signs(r, delta)
        with np.load(snapshot_path, allow_pickle=False) as snapshot:
            np.testing.assert_array_equal(snapshot["corrector_before_boundary_M"],
                                          live["production"]["stages"]["corrector_before_boundary"]["rows"]["M"])
            rounding = last_write(snapshot)
        rounding["joint_write_over_step_M"] = (rounding["joint_last_write"] /
            live["production"]["stages"]["accepted_end"]["band_delta_rms"]["M"])
        result[case] = {"sign_counts": signs_by_arm, "production_rounding": rounding,
                        "live_sha256": replay.sha(live_path), "snapshot_sha256": replay.sha(snapshot_path)}
    args.output.write_text(json.dumps(result, allow_nan=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
