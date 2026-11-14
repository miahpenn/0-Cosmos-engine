"""RESEARCH ONLY. The uniform one-budget background as ordinary differential equations, for acceptance test T1 of the
lattice-cell wall (HANDOFF Nov 14 2026 (4)). Same content and couplings as the machine's local equations, in the
archive normalisation (8 pi G = 1 densities; H^2 = rho/3):
    phi' = Pi,   Pi' = -3 H Pi - V'(phi) + beta_dm rho_dm,   rho_dm' = -3 H rho_dm - beta_dm rho_dm Pi,
    rho_b' = -3 H rho_b,   rho_r' = -4 H rho_r,   H' = -(Pi^2 + rho_dm + rho_b + (4/3) rho_r)/2,
    V(phi) = V0 (e^{-LAMBDA phi} - D_POT),  plus the constant 8 pi V_D of D at rest in its vacuum.
H is integrated from its own equation (so the run passes through the turnaround H = 0); the Friedmann constraint
H^2 = rho/3 is returned as a check.
"""
import math

import numpy as np
from scipy.integrate import solve_ivp


def background(H0, phi0, Pi0, dm0, b0, r0, V0, LAMBDA, D_POT, VD8pi, beta_dm, T, n_out=2001):
    def f(t, y):
        phi, Pi, dm, b, rr, H, N = y
        Vp = -V0 * LAMBDA * math.exp(-LAMBDA * phi)
        return [Pi, -3 * H * Pi - Vp + beta_dm * dm, -3 * H * dm - beta_dm * dm * Pi, -3 * H * b, -4 * H * rr,
                -0.5 * (Pi * Pi + dm + b + 4.0 / 3.0 * rr), H]
    ts = np.linspace(0, T, n_out)
    s = solve_ivp(f, (0, T), [phi0, Pi0, dm0, b0, r0, H0, 0.0], method="DOP853", rtol=1e-13, atol=1e-16, t_eval=ts)
    phi, Pi, dm, b, rr, H, N = s.y
    rho = 0.5 * Pi**2 + V0 * (np.exp(-LAMBDA * phi) - D_POT) + VD8pi + dm + b + rr
    return dict(t=ts, phi=phi, Pi=Pi, dm=dm, b=b, r=rr, H=H, N=N, friedmann=H**2 - rho / 3)
