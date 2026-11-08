"""Snapshot of the V5.5 driver ledger (engine.bidirectional_production.run) at N=40 to T=6.

Run once on the old drivers (reference_pirk_unified stack) and once on the drivers rebuilt over
V55ProductionKernel; nothing else differs. Records the ledger at t = 0 and t = T, the Misner-Sharp
residual summary, and the initial COSMOS scalar values the run starts from.
Usage: python <this file> output.json
"""
import json, sys
import numpy as np
sys.path.insert(0, ".")
from engine import bidirectional_production as bp

out, meta = bp.run(40, 6.0, .06)
keys = [k for k in out if k != "t"]
snap = {"meta": {k: (float(v) if isinstance(v, (int, float, np.floating)) else v) for k, v in meta.items()},
        "t0": {k: float(out[k][0]) for k in keys}, "tT": {k: float(out[k][-1]) for k in keys}, "T": float(out["t"][-1])}
json.dump(snap, open(sys.argv[1], "w"), indent=1)
for k in sorted(keys):
    print(f"{k:12s} t0 {snap['t0'][k]: .6e}  tT {snap['tT'][k]: .6e}")
print({k: v for k, v in snap["meta"].items()})
