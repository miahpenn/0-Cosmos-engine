#!/usr/bin/env python3
"""Independent equation-level check of the GEAR-69 / GEAR-70B weak-field BVP.

Research-only: does not import or modify COSMOS Engine production code. Assumes
polar-areal metric, canonical real scalar, V=-mu^2 phi^2/2, and the displayed
GEAR-69 O(A^2) Einstein constraints and O(A^3) scalar equation. It checks a
regular-center solve_bvp against an independently integrated shooting result.
"""
from __future__ import annotations
import argparse, json, math
import numpy as np
from scipy.integrate import cumulative_simpson, solve_bvp, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.special import spherical_jn

GEAR69_TARGET_C = 1.443222
GEAR70B_REPORTED_C = 1.443200051347
GEAR70B_REPORTED_W0 = -72.160002567356

def verify(G: float = 1.0, R: float = 10.0, mu: float = 1.0,
           source_points: int = 20001, bvp_points: int = 1000) -> dict:
    if G <= 0 or R <= 0 or mu <= 0 or source_points < 1001 or bvp_points < 100:
        raise ValueError("Require G,R,mu>0, source_points>=1001, bvp_points>=100")
    k = math.pi / R
    gamma2 = mu * mu - k * k
    if gamma2 <= 0:
        raise ValueError("The fundamental mode is not unstable for these parameters")
    gamma = math.sqrt(gamma2)
    r = np.linspace(0.0, R, source_points)
    u = spherical_jn(0, k * r)
    up = -k * spherical_jn(1, k * r)
    upp = np.empty_like(r)
    upp[1:] = -2.0 * up[1:] / r[1:] - k * k * u[1:]
    upp[0] = -k * k / 3.0

    # O(A^2) density, mass, and metric constraints from GEAR-69.
    rho2 = 0.5 * (gamma2 * u * u + up * up - mu * mu * u * u)
    m2 = cumulative_simpson(4.0 * math.pi * G * r * r * rho2, x=r, initial=0.0)
    a2 = np.zeros_like(r)
    a2[1:] = m2[1:] / r[1:]
    alphap = np.zeros_like(r)
    alphap[1:] = (a2[1:] / r[1:] + 2.0 * math.pi * G * r[1:] *
                  (gamma2 * u[1:]**2 + up[1:]**2 + mu * mu * u[1:]**2))
    ialpha = cumulative_simpson(alphap, x=r, initial=0.0)
    alpha2 = -a2[-1] - (ialpha[-1] - ialpha)  # alpha2(R)=-a2(R)
    q, s = alpha2 - a2, alpha2 + a2
    qprime = CubicSpline(r, q).derivative()(r)
    div = np.empty_like(r)
    div[1:] = (2.0*q[1:]*up[1:]/r[1:] + qprime[1:]*up[1:] + q[1:]*upp[1:])
    div[0] = -k*k*q[0]  # regular r->0 limit of r^-2 d_r(r^2 q u')
    S3 = -div - mu*mu*s*u - 3.0*gamma2*q*u
    source = CubicSpline(r, S3)
    lam = mu*mu - 9.0*gamma2

    # Regular-center form of w'' + 2w'/r + lam*w = S3; w'(0)=0, w(R)=0.
    def rhs(x, y):
        curvature = source(x) - lam*y[0]
        x_safe = np.where(x == 0.0, 1.0, x)
        wpp = curvature - 2.0*y[1]/x_safe
        wpp = np.where(x == 0.0, curvature/3.0, wpp)
        return np.vstack((y[1], wpp))
    def bc(ya, yb):
        return np.array([ya[1], yb[0]])
    mesh = np.linspace(0.0, R, bvp_points)
    sol = solve_bvp(rhs, bc, mesh, np.zeros((2, len(mesh))), tol=1e-8, max_nodes=20000)
    if not sol.success:
        raise RuntimeError(f"independent solve_bvp failed: {sol.message}")
    w0_bvp = float(sol.sol(0.0)[0])
    C_bvp = -2.0*w0_bvp/(G*R*R)

    # Separate shooting cross-check, initialized from regular Taylor expansions.
    eps = min(1e-4, R*1e-5)
    S0 = float(S3[0])
    def hom_rhs(x, y): return [y[1], -lam*y[0] - 2.0*y[1]/x]
    def part_rhs(x, y): return [y[1], float(source(x)) - lam*y[0] - 2.0*y[1]/x]
    v = solve_ivp(hom_rhs, (eps, R), [1.0-lam*eps**2/6.0, -lam*eps/3.0],
                  rtol=2e-11, atol=2e-13)
    yp = solve_ivp(part_rhs, (eps, R), [S0*eps**2/6.0, S0*eps/3.0],
                   rtol=2e-11, atol=2e-13)
    if not v.success or not yp.success or abs(v.y[0, -1]) < 1e-12:
        raise RuntimeError("independent shooting integration failed")
    w0_shoot = float(-yp.y[0, -1]/v.y[0, -1])
    C_shoot = -2.0*w0_shoot/(G*R*R)
    return {
        "assumptions": "canonical real scalar; V=-mu^2 phi^2/2; polar-areal metric; stated exterior-clock boundary",
        "G": G, "R": R, "mu": mu, "source_points": source_points,
        "bvp_points": bvp_points, "k": k, "gamma": gamma,
        "m2_at_R": float(m2[-1]), "alpha2_at_R": float(alpha2[-1]),
        "solve_bvp_success": bool(sol.success), "solve_bvp_message": sol.message,
        "solve_bvp_adaptive_nodes": int(len(sol.x)),
        "w0_bvp": w0_bvp, "C_bvp": C_bvp,
        "w0_shooting": w0_shoot, "C_shooting": C_shoot,
        "abs_delta_C_bvp_vs_shooting": abs(C_bvp-C_shoot),
        "GEAR69_reported_target_C": GEAR69_TARGET_C,
        "delta_C_vs_GEAR69_target": C_bvp-GEAR69_TARGET_C,
        "relative_percent_vs_GEAR69_target": 100.0*(C_bvp-GEAR69_TARGET_C)/GEAR69_TARGET_C,
        "w0_implied_by_GEAR69_target": -GEAR69_TARGET_C*G*R*R/2.0,
        "GEAR70B_reported_w0": GEAR70B_REPORTED_W0,
        "GEAR70B_reported_C": GEAR70B_REPORTED_C,
        "delta_C_vs_GEAR70B_report": C_bvp-GEAR70B_REPORTED_C,
    }

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--G", type=float, default=1.0)
    parser.add_argument("--R", type=float, default=10.0)
    parser.add_argument("--mu", type=float, default=1.0)
    parser.add_argument("--source-points", type=int, default=20001)
    parser.add_argument("--bvp-points", type=int, default=1000)
    args=parser.parse_args()
    result=verify(args.G,args.R,args.mu,args.source_points,args.bvp_points)
    print(json.dumps(result, indent=2, sort_keys=True))
    if result["abs_delta_C_bvp_vs_shooting"] > 1e-6:
        raise SystemExit("FAIL: solve_bvp and shooting cross-check disagree")
if __name__ == "__main__": main()
