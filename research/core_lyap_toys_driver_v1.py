"""RESEARCH ONLY. #523 re-check with the toy's S (coordination log, Nov 19 2026). Usage: python research/core_lyap_toys_driver_v1.py baseline|claim|rest OUT
Runs run(n, rtol, k) of core2dof_lyap_v1 (old S) or core2dof_lyap_toyS_v1 (toy S) over phases k in range(20) in a process pool."""
import sys, math, numpy as np
from multiprocessing import Pool
sys.path.insert(0, "research")
import core2dof_lyap_v1 as old, core2dof_lyap_toyS_v1 as toy
def job(a):
    which, n, rtol, k = a
    m = old if which == "old" else toy
    return a, m.run(n, rtol, k)
def table(rs, label):
    keys = sorted({(a[1], a[2]) for a, _ in rs}); out = []
    for n, rt in keys:
        v = np.array([r for a, r in rs if a[1] == n and a[2] == rt]); kap2 = 1 - 1.5 * n
        out.append(f"{label} n {n:.4f} rtol {rt:g} kappa {math.sqrt(kap2) if kap2>0 else float('nan'):.4f} | lambda full T: mean {v[:,0].mean():.4f} sd {v[:,0].std():.4f} min {v[:,0].min():.4f} max {v[:,0].max():.4f} | 2nd half mean {v[:,1].mean():.4f} | n_phases {len(v)}")
    return out
if __name__ == "__main__":
    mode, outp = sys.argv[1], sys.argv[2]
    CLAIM = (0.3, 0.5, 0.6, 0.65, 2 / 3); REST = (0.1, 0.2, 0.4, 0.55, 0.62, 0.64, 0.7, 0.75, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)
    with Pool() as P:
        if mode == "baseline":
            rs = P.map(job, [("old", n, 1e-12, 0) for n in CLAIM]); lines = ["BASELINE old S, one phase (k=0), rtol 1e-12: lambda full-T, 2nd half"] + [f"  n {a[1]:.4f}: {r[0]:.4f} {r[1]:.4f}" for a, r in rs]
        elif mode == "oldrest":
            rs = P.map(job, [("old", n, 1e-12, k) for n in (0.65,) + REST for k in range(20)], chunksize=1); lines = table(rs, "OLD-S")
        else:
            ns = CLAIM if mode == "claim" else REST
            rs = P.map(job, [("toy", n, rt, k) for n in ns for rt in (1e-12, 1e-13) for k in range(20)], chunksize=1); lines = table(rs, "TOY-S")
    open(outp, "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
