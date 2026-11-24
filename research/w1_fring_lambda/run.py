"""dt-versus-dr origin trace (Jeremiah's go, Nov 24 04:43 UTC). Declarations: results/research/w1_fring_lambda/PRERUN.md.

Arms A (corrected W1 slice, #581) and B (background "none"), open edge, eight-node start, as in the f_ring smoke
(research/w1_fring_smoke/run.py, whose build/residuals are reused unchanged). Two (dr, dt) pairs, both at
lambda = dt/dr = 0.015 (half the smoke's cfl 0.03): (0.025, 3.75e-4) for 40 steps and (0.05, 7.5e-4) for 20 steps,
each to t = 0.015. Every step saves pointwise H and M at cells 0-5, plus max|H|, max|M| and their cell index and R.
No gating. Usage (repo root): OPENBLAS_NUM_THREADS=1 python3 research/w1_fring_lambda/run.py native_static.npz OUTDIR ARM
with ARM in A025 B025 A05 B05.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.argv = [sys.argv[0], sys.argv[1], sys.argv[2], "--lambda", sys.argv[3]]
src = Path("research/w1_fring_smoke/run.py").read_text().split("# ------------------------------------------------------------------ one arm")[0]
g = {"__file__": "research/w1_fring_smoke/run.py", "__name__": "lambda_trace"}
exec(compile(src, "run.py", "exec"), g)
arm, outdir = sys.argv[4], Path(sys.argv[2])
DT = {0.025: 3.75e-4, 0.05: 7.5e-4}
T_END = 0.015
dr = g["ARMS"][arm][0]
dt = DT[dr]
nsteps = int(round(T_END / dt))
k, adv = g["kernel"], g["advance_with_radiation_admissibility_retries"]
from engine.worldtube import areal_radius  # noqa: E402


def row(state, step):
    grid, geom, tot, raw, H, M = g["residuals"](state)
    R = areal_radius(grid, geom)
    iH, iM = int(np.argmax(np.abs(H))), int(np.argmax(np.abs(M)))
    return {"step": step, "t": float(state.t), "finite": g["finite"](state, H, M),
            "H_cells0_5": [float(x) for x in H[:6]], "M_cells0_5": [float(x) for x in M[:6]],
            "H_max": float(abs(H[iH])), "H_argcell": iH, "H_argR": float(R[iH]),
            "M_max": float(abs(M[iM])), "M_argcell": iM, "M_argR": float(R[iM])}


st = g["build"](arm)
rows, walls = [row(st, 0)], []
for n in range(1, nsteps + 1):
    t0 = time.perf_counter()
    st = adv(k, st, dt).state
    walls.append(time.perf_counter() - t0)
    rows.append(row(st, n))
    if not rows[-1]["finite"]:
        break
(outdir / f"arm_{arm}.json").write_text(json.dumps({"arm": arm, "dr": dr, "dt": dt, "lambda": dt / dr,
                                                    "nsteps": nsteps, "rows": rows, "step_wall": walls}, indent=0))
print(arm, dt, rows[-1]["t"], rows[-1]["H_cells0_5"][0], rows[-1]["H_max"], rows[-1]["H_argcell"])
