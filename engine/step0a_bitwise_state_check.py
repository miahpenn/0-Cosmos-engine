#!/usr/bin/env python3
"""Step 0a: repeated full-state diagnostic-invariance comparison."""
from __future__ import annotations
import hashlib, json, os, subprocess, sys
from pathlib import Path

HARNESS = r"""
import hashlib, json, numpy as np
from engine.production_kernel import V55ProductionKernel

def checksum(state):
    h = hashlib.sha256()
    for label, obj in (("geometry", state.geometry),
                       ("scalars", state.scalars),
                       ("matter", state.matter)):
        h.update(label.encode())
        for name, value in vars(obj).items():
            if isinstance(value, np.ndarray):
                a = np.ascontiguousarray(value)
                h.update(name.encode())
                h.update(str(a.dtype).encode())
                h.update(str(a.shape).encode())
                h.update(a.tobytes(order="C"))
    for name in ("t", "tau", "e_folds"):
        h.update(name.encode())
        h.update(np.asarray([float(getattr(state, name))], dtype=np.float64).tobytes())
    return h.hexdigest()

def checkpoint(k, s, step):
    before = checksum(s)
    d = k.diagnostics(s, profiles=False)
    after = checksum(s)
    if before != after:
        raise RuntimeError(f"diagnostics mutated state at step {step}")
    return {"step": step, "before": before, "after": after, "diagnostics_called": True}

k = V55ProductionKernel()
s = k.initialize(
    resolution=160, r_max=80.0, amplitude=0.01, width=7.0,
    D_amplitude=1.0e-4, include_radiation=True,
)
checks = []
dt = 0.0075
for step in range(21):
    if step % 4 == 0:
        checks.append(checkpoint(k, s, step))
    if step < 20:
        s = k.step(s, dt)
if checks[-1]["step"] != 20:
    checks.append(checkpoint(k, s, 20))
print(json.dumps(checks))
"""

def run(repo: str):
    env = os.environ.copy()
    env["PYTHONPATH"] = repo
    p = subprocess.run(
        [sys.executable, "-c", HARNESS],
        cwd=repo, env=env, text=True, capture_output=True, check=True,
    )
    return json.loads(p.stdout.strip().splitlines()[-1])

def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: step0a_bitwise_state_check.py OLD_REPO NEW_REPO")
    old, new = map(lambda x: Path(x).resolve(), sys.argv[1:])
    a = run(str(old))
    b = run(str(new))
    if a != b:
        for oa, nb in zip(a, b):
            print(f"step {oa['step']}: old_before={oa['before']} new_before={nb['before']} "
                  f"old_after={oa['after']} new_after={nb['after']}")
        raise SystemExit("FAIL: diagnostic invariance or evolved-state checksum mismatch")
    for row in a:
        print(f"step {row['step']}: {row['before']} -> {row['after']}")
    print("PASS: diagnostics called and full-state checksums match at every sampled step")

if __name__ == "__main__":
    main()
