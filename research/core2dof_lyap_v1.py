"""Gradient-free 2-DOF core of the S-D bag: does it reproduce lambda ~ 0.055?  (research; forced = no gradients, radiation, expansion)
V = 2 s^2 + s^4/8 + (3/4) s^2 d^2 - d^2/2 + d^4/8 ; kinetic 1/2 each.
s'' = -(4s + s^3/2 + 1.5 s d^2) ; d'' = d - d^3/2 - 1.5 s^2 d
Start: s = A, A^2 = 2n, d = 0, d' = 1e-3, advanced by phase_k * (pi/20) before measuring (20 phases; a choice).
Benettin: twin 1e-10 apart (random direction fixed seed), renormalise every TAU. Two integrator tolerances."""
import sys, math, numpy as np
from scipy.integrate import solve_ivp
T, TAU, D0, NPH = 2000.0, 10.0, 1e-10, 20
def f(t, y):
    s, d, ps, pd = y
    return [ps, pd, -(4*s + 0.5*s**3 + 1.5*s*d*d), d - 0.5*d**3 - 1.5*s*s*d]
def F2(t, Y): return np.concatenate([f(t, Y[:4]), f(t, Y[4:])])
def run(n, rtol, k):
    A = math.sqrt(2*n); y = np.array([A, 0, 0, 1e-3])
    if k:
        y = solve_ivp(f, (0, k*math.pi/NPH), y, method="DOP853", rtol=rtol, atol=1e-14).y[:, -1]
    rng = np.random.default_rng(1000 + k); u = rng.normal(size=4); u /= np.linalg.norm(u)
    Y = np.concatenate([y, y + D0*u]); acc = []; 
    for i in range(int(T/TAU)):
        Y = solve_ivp(F2, (0, TAU), Y, method="DOP853", rtol=rtol, atol=1e-14).y[:, -1]
        dv = Y[4:] - Y[:4]; dd = np.linalg.norm(dv); acc.append(math.log(dd/D0))
        Y[4:] = Y[:4] + dv*(D0/dd)
    a = np.array(acc); return a.sum()/T, a[len(a)//2:].sum()/(T/2)
if __name__ == "__main__":
    n = float(sys.argv[1]); rtol = float(sys.argv[2]); ks = range(NPH)
    r = np.array([run(n, rtol, k) for k in ks])
    kap2 = 1 - 1.5*n; kap = math.sqrt(kap2) if kap2 > 0 else float("nan")
    print(f"n {n:.4f} rtol {rtol:g} kappa {kap:.4f} | lambda full T: mean {r[:,0].mean():.4f} sd {r[:,0].std():.4f} min {r[:,0].min():.4f} max {r[:,0].max():.4f}"
          f" | 2nd half mean {r[:,1].mean():.4f} | lambda/kappa {r[:,0].mean()/kap:.3f}", flush=True)
