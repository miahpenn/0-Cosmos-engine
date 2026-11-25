"""Stage 0: explicit off/on eight-node control, research only; no acceptance gates.

Run in separate processes, from the same clean exact head and environment as
the control suite.  Stage 1 results and gates are not rescored.  The flat
arms reproduce #606's construction, including its tiny positive dust and
its V0/lambda overrides; native arms retain the production initial state
and original couplings.  The switch is initialization-only.

  python -m research.control_tests.stage0 preflight OUT
  python -m research.control_tests.stage0 arm NAME flat off OUT
  python -m research.control_tests.stage0 arm NAME flat on OUT
  python -m research.control_tests.stage0 arm NAME native off OUT
  python -m research.control_tests.stage0 arm NAME native on OUT

NAME is any of test1.ARMS with kind 1a. Native evolution is a separate
comparison, not flat vacuum and not a substitute for the flat pair.
No centre clamp, damping, altered integrator, or tolerance adjustment.
"""
import argparse
import contextlib
import hashlib
import io
import json
import os
import platform
import subprocess
from pathlib import Path

import numpy as np
import scipy

from engine.production_kernel import V55ProductionKernel
from research.control_tests import ctl, observer, test1


def provenance():
    def git(*args):
        return subprocess.check_output(["git", *args], text=True).strip()
    runtime = io.StringIO()
    with contextlib.redirect_stdout(runtime):
        np.show_runtime()
    return {
        "head": git("rev-parse", "HEAD"),
        "submodules": git("submodule", "status", "--recursive"),
        "tracked_status": git("status", "--porcelain", "--untracked-files=no"),
        "python": platform.python_version(),
        "numpy": np.__version__, "scipy": scipy.__version__,
        "numpy_runtime": runtime.getvalue(),
        "threads": {k: os.environ.get(k) for k in
                    ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")},
    }


def build(dr, r_max, mode, on):
    if mode == "flat":
        return ctl.make(dr, r_max, eight_node_areal=on)
    ctl.ss.V0, ctl.ss.SD_QUARTIC_LAMBDA = ctl._ORIG
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=round(r_max / dr), r_max=r_max,
                              outer_edge="open_edge", initial_structure="none",
                              eight_node_areal=on)
    return kernel, state


def snapshot(state):
    assert state.electromagnetic is None
    arrays = {"r": np.asarray(state.grid.centers).copy()}
    for prefix, obj, fields in (
        ("g", state.geometry, ctl.VAC.ALL_FIELDS),
        ("s", state.scalars, ctl.FIELDS_S),
    ):
        for field in fields:
            arrays[f"{prefix}_{field}"] = np.asarray(getattr(obj, field)).copy()
    for species in ("dark_matter", "baryons", "radiation"):
        fluid = getattr(state.matter, species)
        for field in ("rest", "energy_t", "momentum_r"):
            arrays[f"m_{species}_{field}"] = np.asarray(getattr(fluid, field)).copy()
    arrays["clocks"] = np.asarray([state.t, state.tau, state.e_folds])
    return arrays


def hashes(arrays):
    return {k: hashlib.sha256(str(v.dtype).encode() + str(v.shape).encode()
                               + np.ascontiguousarray(v).tobytes()).hexdigest()
            for k, v in arrays.items()}


def reading(state):
    H, M = ctl.constraints(state)
    return {"H_cell0": float(H[0]), "M_cell0": float(M[0]),
            "H_max": float(np.max(np.abs(H))), "M_max": float(np.max(np.abs(M))),
            "a_minus_b": float(state.geometry.a[0] - state.geometry.b[0]),
            "Aa": float(state.geometry.Aa[0])}


def preflight(out):
    out.mkdir(parents=True, exist_ok=True)
    report = {"provenance": provenance(), "pairs": {}}
    for dr, rm in ((0.05, 8.0), (0.025, 8.0), (0.05, 12.0), (0.0125, 8.0)):
        for mode in ("flat", "native"):
            name = f"{mode}_dr{dr}_r{rm:g}"
            pair, reads = [], []
            for on in (False, True):
                _, state = build(dr, rm, mode, on)
                arrays = snapshot(state)
                pair.append(arrays)
                reads.append({"eight_node_areal": on, "hashes": hashes(arrays),
                              "reading": reading(state),
                              "V0": float(ctl.ss.V0),
                              "SD_QUARTIC_LAMBDA": float(ctl.ss.SD_QUARTIC_LAMBDA)})
                np.savez_compressed(out / f"{name}_{'on' if on else 'off'}.npz", **arrays)
            left, right = pair
            if left.keys() != right.keys():
                raise RuntimeError("state fields differ")
            report["pairs"][name] = {"arms": reads, "fields": {
                k: {"bitwise_equal": reads[0]["hashes"][k] == reads[1]["hashes"][k],
                    "max_abs_diff": float(np.max(np.abs(left[k] - right[k])))}
                for k in left}}
    ctl.ss.V0, ctl.ss.SD_QUARTIC_LAMBDA = ctl._ORIG
    (out / "stage0_preflight.json").write_text(json.dumps(report, indent=2, allow_nan=False))


def arm(name, mode, on, out):
    kind, dr, dt, rm = test1.ARMS[name]
    if kind != "1a":
        raise ValueError("Stage 0 only compares 1a configurations")
    out.mkdir(parents=True, exist_ok=True)
    label = f"stage0_{mode}_{'on' if on else 'off'}_{name}"
    kernel, state = build(dr, rm, mode, on)
    meta = {"provenance": provenance(), "eight_node_areal": on, "mode": mode,
            "V0": float(ctl.ss.V0), "SD_QUARTIC_LAMBDA": float(ctl.ss.SD_QUARTIC_LAMBDA),
            "initial_hashes": hashes(snapshot(state)), "initial_reading": reading(state)}
    # Same production-step/recorder equivalence contract as #606's flat arms.
    observer.check_bitwise(kernel, observer.Recorder(), state, dt)

    def on_save(st, g0, kdot, R):
        return {"state_hashes": hashes(snapshot(st)), **reading(st)}

    def on_step(n, st, per_save):
        # Every step, no extra RHS evaluations. Cells are at r=(i+1/2)dr.
        g = st.geometry
        return {"step": n, "t": n * dt,
                "a_minus_b": (g.a[:3] - g.b[:3]).tolist(),
                "Aa": g.Aa[:3].tolist(), "X": g.X[:3].tolist(),
                "Lambda": g.Lambda[:3].tolist()}

    _, rec = ctl.run(label, kernel, state, dt, out, with_cdot=False,
                     t_end=test1.T_END, on_save=on_save, on_step=on_step)
    rec["stage0"] = meta
    rec["centre_trace"] = rec.pop("near_save")
    (out / f"{label}.json").write_text(json.dumps(rec, allow_nan=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("preflight")
    p.add_argument("out", type=Path)
    p = sub.add_parser("arm")
    p.add_argument("name", choices=[k for k, v in test1.ARMS.items() if v[0] == "1a"])
    p.add_argument("mode", choices=("flat", "native"))
    p.add_argument("switch", choices=("off", "on"))
    p.add_argument("out", type=Path)
    args = parser.parse_args()
    if args.cmd == "preflight":
        preflight(args.out)
    else:
        arm(args.name, args.mode, args.switch == "on", args.out)


if __name__ == "__main__":
    main()
