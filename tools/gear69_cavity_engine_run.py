import sys, json; sys.path.insert(0, ".")
import numpy as np
from research import gear69_cavity_engine_harness_v1 as h
N = int(sys.argv[1]); T = float(sys.argv[2]); out = {}
for A in map(float, sys.argv[4:]):
    rows = h.run(N, A, T, every=0.1)
    out[str(A)] = rows.tolist()
    print(f"N={N} A={A}: rows={len(rows)} D_c end {rows[-1][2]:.3e} Gamma_tau end {rows[-1][3]:.8f}", flush=True)
json.dump(out, open(sys.argv[3], "w"))
