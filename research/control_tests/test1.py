"""Stage 1 runner (design 1BK-Xj3r...; arm list per Amendment 13, wording per Amendment 15). Research only.

Code only until Jeremiah's go for stage 1. Nothing here reads a gate: arms write JSON/npz readings, and the 60-digit
readback (mp_readback, under python3 with mpmath) scores them against the predeclared gates afterwards.

  python -m research.control_tests.test1 inputs OUT     t = 0 input-check slices (zero steps; A13 section 3)
  python -m research.control_tests.test1 arm NAME OUT   one evolved arm (names in ARMS)

Each evolved arm: saves at t = 0.3k to t_end 2.1; at every save it also takes one recorded step (the R-stage
observer, A4 section 3) from the saved state and keeps its stage arrays. The recorded step is checked bitwise
against the production step on the arm's first save before any reading is kept; the run then steps with the
production kernel as before, so the recording never feeds the evolution.
"""
import json
import sys
from pathlib import Path

import numpy as np

from research.control_tests import ctl, observer, t0_checks

DT0, DT1 = 3.75e-4, 1.875e-4
T_END = 2.1
K0_1B, C0_1B = -0.1, 0.5

# name: (kind, dr, dt, r_max)
ARMS = {
    "1a_dr05_dt0_r8": ("1a", 0.05, DT0, 8.0),
    "1a_dr025_dt0_r8": ("1a", 0.025, DT0, 8.0),
    "1a_dr05_dt1_r8": ("1a", 0.05, DT1, 8.0),
    "1a_dr05_dt0_r12": ("1a", 0.05, DT0, 12.0),
    "1b_K-0.1_c0.5_dr05_dt0_r8": ("1b", 0.05, DT0, 8.0),
}
INPUTS_1B = [(K0, c0, dr) for K0 in (-0.1, -0.017) for c0 in (0.5, 0.9) for dr in (0.05, 0.025, 0.0125)]


def build(kind, dr, r_max):
    if kind == "1a":
        k, s = ctl.make(dr, r_max)
        return k, s, {}
    return t0_checks.build_1b(dr, K0_1B, C0_1B, r_max)


def inputs(out):
    """Zero-step slices for the t = 0 input checks (B0-B6, A13 section 1), 1b inputs on all grids plus 1a arms."""
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    meta = {}
    cases = [(f"1b_K{K0}_c{c0}_dr{dr}", lambda K0=K0, c0=c0, dr=dr: t0_checks.build_1b(dr, K0, c0, 8.0))
             for K0, c0, dr in INPUTS_1B]
    cases += [(f"1a_dr{dr}_r{rm:g}", lambda dr=dr, rm=rm: (*ctl.make(dr, rm), {}))
              for dr, rm in ((0.05, 8.0), (0.025, 8.0), (0.05, 12.0))]
    for name, mk in cases:
        k, s, p = mk()
        g0, kdot, R, _ = ctl.semi(k, s)
        H, M = ctl.constraints(s)
        tot = observer.matter_totals(s.grid, s.geometry, s.scalars, s.matter)
        extra = {"H": H, "M": M, **{f"m_{x}": v for x, v in tot.items()}}
        np.savez_compressed(out / f"{name}.npz", **ctl.save_arrays(s, g0, kdot, R, extra))
        meta[name] = {"dr": s.grid.dr, "n": s.grid.n, "r_max": s.grid.r_max, **p}
    (out / "inputs.json").write_text(json.dumps(meta, indent=1))


def arm(name, out):
    kind, dr, dt, r_max = ARMS[name]
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    kernel, state, p = build(kind, dr, r_max)
    rec_kernel = observer.Recorder()
    rstage, checked = [], []

    def on_save(st, g0, kdot, R):
        if not checked:
            observer.check_bitwise(kernel, rec_kernel, st, dt)
            checked.append(True)
        new, log = rec_kernel.recorded_step(st, dt)
        reading, arrays = observer.r_stage(st.grid, st, new, log, dr)
        k = len(rstage)
        np.savez_compressed(out / f"{name}_rstage{k}.npz", **arrays)
        rstage.append({"save": k, "t": st.t, "terms": reading})
        (out / f"{name}_rstage.json").write_text(json.dumps(rstage))
        return {"rstage_file": f"{name}_rstage{k}.npz"}

    _, rec = ctl.run(name, kernel, state, dt, out, with_cdot=False, t_end=T_END, on_save=on_save)
    rec["params"] = p
    (out / f"{name}.json").write_text(json.dumps(rec))


if __name__ == "__main__":
    if sys.argv[1] == "inputs":
        inputs(sys.argv[2])
    elif sys.argv[1] == "arm":
        arm(sys.argv[2], sys.argv[3])
    else:
        raise SystemExit(__doc__)
