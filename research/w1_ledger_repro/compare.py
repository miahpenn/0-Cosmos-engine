"""Compare reproduced W1 ledgers with Aster's archived ledger.json files (no computation of physics here).

Usage: python3 research/w1_ledger_repro/compare.py X_JSON C_JSON  -> prints JSON: per-target values, bitwise equality,
relative difference, and whole-ledger leaf comparison (count of numeric leaves differing, max relative difference).
"""
import json
import sys

REF = {"x": "results/research/w1_fring_x_endpoint/ledger.json",
       "c": "results/research/w1_fring_conformal_ledger/ledger.json"}
TARGETS = {
    "nine-point H, nominal": ("x", ["runs", "1", "readback", "cells", "H_metric_9", 0]),
    "nine-point H, half": ("x", ["runs", "2", "readback", "cells", "H_metric_9", 0]),
    "nine-point H, quarter": ("x", ["runs", "4", "readback", "cells", "H_metric_9", 0]),
    "X writes, nominal-half": ("x", ["pairs", "1-2", "curvature", "9", "write", 0]),
    "X writes, half-quarter": ("x", ["pairs", "2-4", "curvature", "9", "write", 0]),
    "tracefree rate, nominal-half": ("c", ["pairs", "1-2", "curvature", "9", "contributions", "tracefree", "tangent", 0]),
    "tracefree rate, half-quarter": ("c", ["pairs", "2-4", "curvature", "9", "contributions", "tracefree", "tangent", 0]),
    "summed writes, nominal-half": ("c", ["pairs", "1-2", "curvature", "9", "combined_metric_write", 0]),
    "summed writes, half-quarter": ("c", ["pairs", "2-4", "curvature", "9", "combined_metric_write", 0]),
}


def get(o, path):
    for k in path:
        o = o[k]
    return o


def leaves(o, p=()):
    if isinstance(o, dict):
        for k, v in o.items():
            yield from leaves(v, p + (k,))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from leaves(v, p + (i,))
    else:
        yield p, o


def rel(a, b):
    return 0.0 if a == b else abs(a - b) / max(abs(a), abs(b))


new = {"x": json.load(open(sys.argv[1])), "c": json.load(open(sys.argv[2]))}
ref = {k: json.load(open(v)) for k, v in REF.items()}
out = {"targets": {}, "ledgers": {}}
for name, (k, path) in TARGETS.items():
    a, b = get(ref[k], path), get(new[k], path)
    out["targets"][name] = {"archived": a, "reproduced": b, "bitwise_equal": a == b, "rel_diff": rel(a, b)}
for k in ("x", "c"):
    R, N = dict(leaves(ref[k])), dict(leaves(new[k]))
    num = [p for p in R if isinstance(R[p], float) and isinstance(N.get(p), (int, float))]
    differ = [p for p in R if R[p] != N.get(p)]
    out["ledgers"][k] = {"leaves": len(R), "same_keys": set(R) == set(N), "leaves_differing": len(differ),
                         "first_differing": [list(map(str, p)) for p in differ[:5]],
                         "max_rel_diff_numeric": max((rel(R[p], N[p]) for p in num), default=0.0)}
print(json.dumps(out, indent=1))
