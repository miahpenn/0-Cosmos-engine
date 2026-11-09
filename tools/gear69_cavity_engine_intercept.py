import json, sys, numpy as np
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
        p = np.log(abs(e1 / e2)) / np.log(h[0] / h[1]) if e2 != 0 else float("nan")
        # Richardson with ratio r = h1/h2 between the two finest
        r = h[1] / h[2]; Iinf = I[2] + (I[2] - I[1]) / (r ** p - 1)
        print(f"A={A} deg2: intercepts {[f'{x:.5f}' for x in I]}  observed order p={p:.2f}  Richardson I_inf={Iinf:.5f} (diff {Iinf - 4.58259:+.1e})")
