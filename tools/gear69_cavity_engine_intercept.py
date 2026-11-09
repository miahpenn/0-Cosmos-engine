import json, sys, numpy as np


def _observed_order(h, e1, e2):
    """Observed convergence order from three spacings h0 > h1 > h2 (any ratios)."""
    from scipy.optimize import brentq
    if e2 == 0 or e1 / e2 <= 0:
        return float("nan")
    target = e1 / e2
    f = lambda p: (h[0] ** p - h[1] ** p) / (h[1] ** p - h[2] ** p) - target
    try:
        return brentq(f, 0.1, 12.0)
    except ValueError:
        return float("nan")

F = sys.argv[1]
def intercepts(N):
    d = json.load(open(f"{F}/cav{N}c.json")); lin = np.array(d["1e-05"])
    m = (lin[:, 0] >= 1.5) & (lin[:, 0] <= 3.0); g = float(np.mean(lin[m, 3] / (1 + 4.58 * lin[m, 2] ** 2)))
    out = {"gamma": g}
    for A in ("0.001", "0.002"):
        a = np.array(d[A]); x = a[:, 2] ** 2; y = (a[:, 3] - g) / g / x
        w = (x >= 5e-4) & (x <= 5e-3) & (a[:, 0] >= 1.5)
        out[A] = {deg: tuple(np.polyfit(x[w], y[w], deg)[::-1][:2]) for deg in (1, 2, 3)}
    return out
res = {N: intercepts(N) for N in sys.argv[2].split(",")}
for N, r in res.items():
    print(f"N={N}: gamma_num {r['gamma']:.10f} (exact 0.9493702942, diff {r['gamma'] - 0.9493702942:+.1e})")
    for A in ("0.001", "0.002"):
        print("   A=" + A + "  " + "  ".join(f"deg{d}: I={r[A][d][0]:.5f} ({r[A][d][0] - 4.58259:+.1e}) s={r[A][d][1]:.2f}" for d in (1, 2, 3)))
Ns = list(res)
if len(Ns) == 3:
    for A in ("0.001", "0.002"):
        I = [res[N][A][2][0] for N in Ns]; h = [40.0 / int(N) for N in Ns]
        e1, e2 = I[0] - I[1], I[1] - I[2]
        # Observed order for UNEQUAL refinement ratios: solve
        # (h0^p - h1^p) / (h1^p - h2^p) = e1 / e2. The closed form
        # log(e1/e2)/log(h0/h1) holds only when h0/h1 == h1/h2.
        p = _observed_order(h, e1, e2)
        # Richardson with ratio r = h1/h2 between the two finest
        r = h[1] / h[2]; Iinf = I[2] + (I[2] - I[1]) / (r ** p - 1)
        print(f"A={A} deg2: intercepts {[f'{x:.5f}' for x in I]}  observed order p={p:.2f}  Richardson I_inf={Iinf:.5f} (diff {Iinf - 4.58259:+.1e})")
