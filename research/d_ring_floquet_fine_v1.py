"""RESEARCH ONLY. GEAR-105's Floquet computation (chat script GEAR-105_floquet.py, same equations and units lambda = 2,
D* = 1) on a fine wavenumber grid, to compare with the field runs of research/d_ring_box_v1.py (HANDOFF Nov 14 (3),
Task A). GEAR-105's table scanned k^2 in steps of 0.01, which under-resolves the band at eps = 0.05 (band top 0.015).
Usage: python research/d_ring_floquet_fine_v1.py
"""
import numpy as np
from scipy.integrate import solve_ivp


def ring(eps):
    f = lambda t, y: [y[1], y[0] - y[0] ** 3]  # noqa: E731
    ev = lambda t, y: y[1]  # noqa: E731
    ev.direction = -1
    s = solve_ivp(f, (0, 50), [1 + eps, 0], rtol=1e-12, atol=1e-14, events=ev, dense_output=True)
    T = [t for t in s.t_events[0] if t > 1][0]
    return s.sol, T


def mu(sol, T, k2, kind):
    c = 1.0 if kind == "S" else -1.0

    def f(t, y):
        w2 = k2 + c + 3 * sol(t)[0] ** 2
        return [y[1], -w2 * y[0], y[3], -w2 * y[2]]
    r = solve_ivp(f, (0, T), [1, 0, 0, 1], rtol=1e-11, atol=1e-13)
    M = np.array([[r.y[0, -1], r.y[2, -1]], [r.y[1, -1], r.y[3, -1]]])
    return np.log(max(abs(np.linalg.eigvals(M)))) / T


if __name__ == "__main__":
    for eps in (0.05, 0.1, 0.2):
        sol, T = ring(eps)
        k2 = np.linspace(0, 8 * eps**2, 161)
        m = np.array([mu(sol, T, q, "D") for q in k2])
        i = int(m.argmax()); on = k2[m > 1e-6]
        mS = max(mu(sol, T, q, "S") for q in np.linspace(0, 16, 321))
        print(f"eps = {eps}: D band k^2 = {on.min():.5f}-{on.max():.5f} (k = {np.sqrt(on.min()):.4f}-{np.sqrt(on.max()):.4f}); peak rate {m[i]:.4e}"
              f" at k^2 = {k2[i]:.5f} (k = {np.sqrt(k2[i]):.4f}); rate/eps^2 = {m[i]/eps**2:.3f}; k^2/eps^2 = {k2[i]/eps**2:.2f}; "
              f"largest S rate over k^2 in [0, 16]: {mS:.1e}")
