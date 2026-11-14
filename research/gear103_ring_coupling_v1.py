"""RESEARCH ONLY. GEAR-103 check (HANDOFF to projects Nov 14 2026 (3), Task D).

Flat space, no gravity, S = 0. D's uniform ring under GEAR-91's coupling (research/gear91_coupling_v1.potential):
    V(D, phi) = -e^{-2 beta_D phi} D^2/2 + (lambda/8) D^4,      D'' = e^{-2 beta_D phi} D - (lambda/2) D^3,
with phi PRESCRIBED (FORCED, AUDIT §D): a slow uniform ramp phi(t) = phi_0 + v t. phi is not evolved, so nothing here
tests phi's own dynamics, only the two statements of GEAR-103 §2:
  (1) the ring's frequency follows sqrt2 e^{-beta_D phi}: per ring period, omega = 2 pi / T between successive maxima;
      least-squares slope of ln omega against phi (prediction: -beta_D);
  (2) the ring's term in phi's equation. phi'' = ... - 8 pi dV/dphi (the 8 pi of GEAR-91). Per period,
      G = <dV/dphi> - dV/dphi|_vacuum = beta_D e^{-2 beta_D phi} <D^2 - D*^2>,   D*^2 = 2 e^{-2 beta_D phi}/lambda,
      against the ring's energy density rho_ring = <D_t^2/2 + V(D) - V(D*)>. Prediction: G = -beta_D rho_ring, i.e. the
      source in phi's equation is -8 pi G = +8 pi beta_D rho_ring (Amendola's form for a mass m ~ e^{-beta phi}).
Uniform ring = one ODE; integrated with DOP853. Rulers: integrator tolerance (1e-11 / 1e-8), ramp speed v, ring size eps.
Negative controls: beta_D = +0.04 (signs flip) and beta_D = 0 (no dependence, G = 0).
Unpaid inputs: lambda = 5000 (scales out), beta_D = -0.04 (value anchored by observation; form GEAR-91), eps, phi_0 = 0,
v. G_reset not involved.
Usage: python research/gear103_ring_coupling_v1.py
"""
import math

import numpy as np
from scipy.integrate import solve_ivp

LAM = 5000.0


def run(beta, eps, v, phi_end=5.0, rtol=1e-11):
    T = phi_end / v
    eD = lambda t: math.exp(-2 * beta * v * t)  # noqa: E731
    D0 = math.sqrt(2 / LAM) * (1 + eps)

    def f(t, y):
        return [y[1], eD(t) * y[0] - 0.5 * LAM * y[0] ** 3]
    ev = lambda t, y: y[1]  # noqa: E731
    ev.direction = -1                                     # maxima of D
    sol = solve_ivp(f, (0, T), [D0, 0.0], method="DOP853", rtol=rtol, atol=rtol * 1e-4, events=ev, dense_output=True)
    tm = sol.t_events[0]
    rows = []
    for a, b in zip(tm[:-1], tm[1:]):
        ts = np.linspace(a, b, 400, endpoint=False)
        D, Dt = sol.sol(ts)
        e = np.exp(-2 * beta * v * ts)
        Ds2 = 2 * e / LAM
        Vv = -e * e / (2 * LAM)
        rho = np.mean(0.5 * Dt**2 - 0.5 * e * D**2 + LAM / 8 * D**4 - Vv)
        G = np.mean(beta * e * (D**2 - Ds2))
        amp = (D.max() - D.min()) / (2 * math.sqrt(Ds2.mean()))
        rows.append((0.5 * (a + b), v * 0.5 * (a + b), 2 * math.pi / (b - a), rho, G, amp))
    return np.array(rows)


if __name__ == "__main__":
    print("GEAR-103 check: uniform D ring, phi prescribed (ramp 0 -> 5), lambda = 5000. Predictions: d ln(omega)/d phi = -beta_D;")
    print("G/rho_ring = -beta_D (G = ring part of dV/dphi; the source in phi's equation is -8 pi G).")
    print("  beta_D   eps    v      rtol   | omega/sqrt2 at phi=0 / 2.5 / 5 | slope dln(omega)/dphi | G/rho_ring: start / mid / end | ring size start -> end | rho_ring end/start")
    for beta in (-0.04, 0.04, 0.0):
        for eps in (0.05, 0.10, 0.20):
            for v, rtol in ((0.005, 1e-11), (0.0025, 1e-11), (0.005, 1e-8)):
                r = run(beta, eps, v, rtol=rtol)
                om = r[:, 2] / math.sqrt(2)
                slope = np.polyfit(r[:, 1], np.log(r[:, 2]), 1)[0]
                q = r[:, 4] / r[:, 3]
                n = len(r)
                print(f"  {beta:+.2f}   {eps:.2f}  {v:.4f}  {rtol:.0e}  |  {om[0]:.5f} / {om[n//2]:.5f} / {om[-1]:.5f}   |   {slope:+.5f}"
                      f"          |  {q[0]:+.5f} / {q[n//2]:+.5f} / {q[-1]:+.5f}   |  {r[0,5]:.4f} -> {r[-1,5]:.4f}   |  {r[-1,3]/r[0,3]:.4f}")
