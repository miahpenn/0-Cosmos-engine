"""Timing test (no gating): ~10 steps per arm at dr 0.025 on the corrected slice.
Step time includes the kernel's own per-step diagnostics (kernel.step calls kernel.diagnostics; it cannot be switched
off without an engine change), so the bare step = step - diagnostics, with diagnostics timed separately on the same
state. The smoke's extra readings (measure) are timed separately too.
Usage: OPENBLAS_NUM_THREADS=1 python3 research/w1_fring_smoke/timing.py native_static.npz OUTDIR ARM
"""
import json
import sys
import time
from pathlib import Path

sys.argv = [sys.argv[0], sys.argv[1], sys.argv[2], "--timing", sys.argv[3]]
src = Path("research/w1_fring_smoke/run.py").read_text().split("# ------------------------------------------------------------------ one arm")[0]
g = {"__file__": "research/w1_fring_smoke/run.py", "__name__": "timing"}
exec(compile(src, "run.py", "exec"), g)
arm, outdir = sys.argv[4], Path(sys.argv[2])
NS = 10
st = g["build"](arm)
dt = g["CFL"] * float(st.grid.dr)
k = g["kernel"]
adv = g["advance_with_radiation_admissibility_retries"]
st = adv(k, st, dt).state  # warm-up step, not timed
tstep, tdiag, tmeas, tfin = [], [], [], []
for _ in range(NS):
    t0 = time.perf_counter(); st = adv(k, st, dt).state; tstep.append(time.perf_counter() - t0)
for _ in range(NS):
    t0 = time.perf_counter(); k.diagnostics(st); tdiag.append(time.perf_counter() - t0)
    t0 = time.perf_counter(); g["measure"](st); tmeas.append(time.perf_counter() - t0)
    t0 = time.perf_counter(); g["assert_finite_state"](st); tfin.append(time.perf_counter() - t0)
med = lambda x: sorted(x)[len(x) // 2]  # noqa: E731
res = {"arm": arm, "n_steps_timed": NS, "step_incl_kernel_diagnostics": med(tstep), "kernel_diagnostics": med(tdiag),
       "smoke_measure": med(tmeas), "assert_finite": med(tfin), "raw": {"step": tstep, "diag": tdiag, "measure": tmeas}}
res["bare_step"] = res["step_incl_kernel_diagnostics"] - res["kernel_diagnostics"]
(outdir / f"timing_{arm}.json").write_text(json.dumps(res, indent=1))
print(json.dumps({k_: v for k_, v in res.items() if k_ != "raw"}))
