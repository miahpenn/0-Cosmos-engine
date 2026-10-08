#!/usr/bin/env python3
"""Step 0a: compare evolved numerical state checksums across kernel revisions.

This is diagnostic-only. It compares the parent revision with the candidate
revision at repeated intermediate steps. The diagnostic bookkeeping added in
Step 0a must not change any evolved numerical state.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np


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

k = V55ProductionKernel()
s = k.initialize(
    resolution=32, r_max=16.0, amplitude=0.01, width=7.0,
    D_amplitude=1.0e-4, include_radiation=True,
)
checks = []
for step in range(25):
    if step % 4 == 0:
        checks.append((step, checksum(s)))
    s = k.step(s, 0.01)
checks.append((25, checksum(s)))
print(json.dumps(checks))
"""

def run(repo: str):
    env = os.environ.copy()
    env["PYTHONPATH"] = repo
    p = subprocess.run(
        [sys.executable, "-c", HARNESS],
        cwd=repo,
        env=env,
        text=True,
        capture_output=True,
        check=True,
    )
    return json.loads(p.stdout.strip().splitlines()[-1])

def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: step0a_bitwise_state_check.py OLD_REPO NEW_REPO")
    old, new = map(str, map(Path, sys.argv[1:]))
    a = run(old)
    b = run(new)
    if a != b:
        for (sa, ca), (sb, cb) in zip(a, b):
            print(f"step {sa}: old={ca} new={cb} match={ca == cb}")
        raise SystemExit("FAIL: evolved-state checksum diverged")
    for step, digest in a:
        print(f"step {step}: {digest}")
    print("PASS: full numerical-state checksums match at every sampled step")

if __name__ == "__main__":
    main()
