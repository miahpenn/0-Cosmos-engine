"""RESEARCH ONLY. Floquet scan of the k=0 D mode above the band (design + amendment: coordination log, Nov 18 2026). Toy S: s''=-(s+0.5 s^3); D: d''=(1-1.5 s^2) d.
mu = amplitude growth per time of the Hill equation over the real S orbit (no friction). Descent uses the exact adiabatic law from the orbit, plus pure dust (-3Hx) and radiation (-2Hx) bounds."""
import math, numpy as np
from scipy.integrate import solve_ivp
def orbit(x):
    A = math.sqrt(x)
    def f(t, y):
        s, p, a, b, c, d, q = y            # Phi=[[a,b],[c,d]] columns evolve under D'' = (1-1.5 s^2) D
        w = 1 - 1.5 * s * s
        return [p, -(s + 0.5 * s ** 3), c, d, w * a, w * b, p * p]
    ev = lambda t, y: y[1]; ev.terminal = True; ev.direction = 1
    r = solve_ivp(f, [0, 10 / math.sqrt(1 + 0.5 * x) * 4], [A, 0, 1, 0, 0, 1, 0], method="DOP853", rtol=1e-12, atol=1e-14, events=ev)
    th = r.t_events[0][0]; y = r.y_events[0][0]; tr = y[2] + y[5]
    # Phi columns: (a,c) from IC (1,0) and (b,d) from IC (0,1): monodromy [[a,b],[c,d]] -> trace a+d
    mu = 0.0 if abs(tr) <= 2 else math.log((abs(tr) + math.sqrt(tr * tr - 4)) / 2) / th
    return dict(x=x, T=2 * th, mu=mu, muT=mu * 2 * th, K=y[6] / th, tr=tr)
GRID = (1.5, 2.0, 2.67, 3.0, 3.27, 4, 5, 5.33, 6, 8, 12, 16, 24, 32, 48, 64, 100, 157, 160)
if __name__ == "__main__":
    rows = {x: orbit(x) for x in GRID}
    print("x, mu (amp growth/time), mu*T (per S period), mu/sqrt(x), trace")
    for x in GRID:
        r = rows[x]; print(f"{x:7.2f} mu {r['mu']:.5g} muT {r['muT']:.4g} mu/sqrt(x) {r['mu']/math.sqrt(x):.4g} tr {r['tr']:.8f}")
    print("METHOD CHECK (inside band, expect ~0):", {x: round(rows[x]['mu'], 5) for x in (1.5, 2.0, 2.67)})
    above = [x for x in GRID if x > 2.67]; mus = [rows[x]['mu'] for x in above]
    print("(a) SIGN:", "POSITIVE-ALL" if min(mus) > 1e-6 else ("ZERO-somewhere" if any(abs(m) < 1e-8 for m in mus) else "MIXED"), "| min mu above band", f"{min(mus):.4g}", "at x", above[int(np.argmin(mus))])
    m = rows[3.27]['mu']; ok1 = abs(m - 0.35) <= 0.105; ok2 = abs(2 * m - 0.35) <= 0.105
    print(f"(b) SOURCE: mu(3.27)={m:.4f}, 2mu={2*m:.4f}; MATCH {'amplitude reading' if ok1 else ('energy-rate reading' if ok2 else 'NEITHER (GEAR-100 0.35 came from something else)')}")
    xs = np.geomspace(2.7, 160, 220); tab = [orbit(x) for x in xs]; ly = np.log(xs)
    mu_i = lambda y: np.interp(y, ly, [t['mu'] for t in tab]); K_i = lambda y: np.interp(y, ly, [t['K'] for t in tab])
    yy = np.linspace(math.log(157), math.log(8 / 3), 4000)
    print("(c) descent from x0=157 to 8/3, H=H_b*sqrt(x/x_b): law | t_total | N_gross | N_net(mu-1.5H) | frac time mu<1.5H | min mu | verdict")
    for Hb in (0.0057, 0.01, 0.02):
        out = []
        for law in ("exact", "dust", "radiation"):
            x = np.exp(yy); H = Hb * np.sqrt(x / (8 / 3)); A = np.sqrt(x)
            if law == "exact": rate = 3 * H * K_i(yy) / ((A + A ** 3 / 2)) / A * 2      # -dy/dt = 2 (dA/dt)/A
            elif law == "dust": rate = 3 * H
            else: rate = 2 * H
            dt = np.abs(np.diff(yy)) / (0.5 * (rate[1:] + rate[:-1])); mu = mu_i(yy); mm = 0.5 * (mu[1:] + mu[:-1]); hh = 0.5 * (H[1:] + H[:-1])
            Ng = float(np.sum(mm * dt)); Nn = float(np.sum((mm - 1.5 * hh) * dt)); fr = float(np.sum(dt[mm < 1.5 * hh]) / np.sum(dt))
            v = "SATURATES ABOVE" if (mu.min() > 0.05 and Nn >= 15) else ("DAMPED" if fr > 0.5 else "intermediate")
            out.append(f"  H_b {Hb} {law:9s} t {dt.sum():8.1f} N_gross {Ng:9.1f} N_net {Nn:9.1f} frac_damped {fr:.2f} min_mu {mu.min():.4g} -> {v}")
        print("\n".join(out))
    sel = (48, 64, 100, 157, 160); mm = [rows[x]['mu'] for x in sel]; mr = [rows[x]['mu'] / math.sqrt(x) for x in sel]
    sp = lambda v: max(v) / min(v) - 1
    print(f"(d) SCALING: mu spread over x=48..160 {sp(mm):.3f}; mu/sqrt(x) spread {sp(mr):.3f}; verdict:", "mu ~ const (my estimate)" if sp(mm) < 0.10 and sp(mr) >= 0.10 else ("mu ~ sqrt(x) (chat)" if sp(mr) < 0.10 and sp(mm) >= 0.10 else "no verdict"))
