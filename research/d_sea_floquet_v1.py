#!/usr/bin/env python3
"""Research-only Floquet bridge: linear D perturbations about the prescribed sea.

This is Step 2 of the stone/Mathieu bridge.  It does not import the engine,
read stone exports, fit parameters, run a nonlinear stone, or change production.

Starting from (u=1)
    d_tt - lap(d) - d + 1/2 d^3 + 3/2 s^2 d = 0,
with s=A cos(omega_S t), x=A^2, and the adiabatic reference
    d0^2 = 2 - 3x/2,
write d=d0+delta.  For a Fourier mode k, the homogeneous perturbation obeys
    delta_tt + [k^2 + 2 - 3x/2 + (3x/4) cos(2 omega_S t)] delta = 0.
The linearization also has the separate additive term
    -(3/4) d0 x cos(2 omega_S t).
That forcing is reported as F below, but is deliberately excluded from the
homogeneous monodromy: parametric instability and additive response are
different questions.

With tau=omega_S t this is Mathieu's equation
    y'' + (a - 2q cos(2 tau)) y = 0,
    a=(k^2+2-3x/2)/omega_S^2, q=-3x/(8 omega_S^2).
The coefficient period is pi in tau, hence T=pi/omega_S in physical time.
The monodromy is integrated over one coefficient period.  For this
undamped second-order equation its determinant should be 1 (Wronskian
conservation); the script records the numerical deviation.

Scope: stability about a *prescribed* sea and periodic S, not a
self-consistent nonlinear stone.  d0 is an adiabatic reference, not an exact
periodic solution.  The averaged omega_S^2=4-15x/8 is an optional comparison
only; omega_S remains an independent input in the free-frequency sweep.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


def monodromy(x: float, k: float, omega_s: float, *, rtol: float, atol: float):
    """Return Mathieu parameters, Floquet multipliers, growth rate and det(M)."""
    if not (0.0 < x < 4.0 / 3.0):
        raise ValueError("x must lie in the open positive-sea range (0, 4/3)")
    if k < 0.0 or omega_s <= 0.0:
        raise ValueError("k must be nonnegative and omega_s positive")

    d0 = math.sqrt(2.0 - 1.5 * x)
    omega2 = omega_s * omega_s
    a = (k * k + 2.0 - 1.5 * x) / omega2
    q = -3.0 * x / (8.0 * omega2)
    forcing_F = -0.75 * d0 * x / omega2

    # Fundamental matrix columns are the solutions from (1,0) and (0,1).
    # Integrate in dimensionless tau over the coefficient period [0, pi].
    def rhs(tau, z):
        y11, y21, y12, y22 = z
        coeff = a - 2.0 * q * math.cos(2.0 * tau)
        return (y21, -coeff * y11, y22, -coeff * y12)

    sol = solve_ivp(
        rhs,
        (0.0, math.pi),
        (1.0, 0.0, 0.0, 1.0),
        method="DOP853",
        rtol=rtol,
        atol=atol,
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    y11, y21, y12, y22 = sol.y[:, -1]
    M = np.array([[y11, y12], [y21, y22]], dtype=float)
    det_m = float(np.linalg.det(M))
    multipliers = np.linalg.eigvals(M.astype(complex))
    max_mod = float(np.max(np.abs(multipliers)))
    period = math.pi / omega_s
    gamma = math.log(max_mod) / period
    return {
        "x": x,
        "k": k,
        "omega_S": omega_s,
        "d0": d0,
        "a": a,
        "q": q,
        "F_additive": forcing_F,
        "mu_plus_abs": float(np.max(np.abs(multipliers))),
        "mu_1_real": float(multipliers[0].real),
        "mu_1_imag": float(multipliers[0].imag),
        "mu_2_real": float(multipliers[1].real),
        "mu_2_imag": float(multipliers[1].imag),
        "gamma_F": gamma,
        "det_M": det_m,
        "det_M_minus_1": det_m - 1.0,
        "unstable": bool(max_mod > 1.0 + 1e-8),
        "omega_S_averaged": math.sqrt(4.0 - 15.0 * x / 8.0),
    }


def build_parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", default="results/d-sea-floquet/floquet.csv")
    p.add_argument("--x-min", type=float, default=0.05)
    p.add_argument("--x-max", type=float, default=1.30)
    p.add_argument("--nx", type=int, default=26)
    p.add_argument("--k-min", type=float, default=0.0)
    p.add_argument("--k-max", type=float, default=3.0)
    p.add_argument("--nk", type=int, default=61)
    p.add_argument(
        "--omega",
        type=float,
        nargs="+",
        default=[1.5, 1.75, 2.0, 2.25, 2.5],
        help="independent fixed omega_S values for the free-frequency sweep",
    )
    p.add_argument(
        "--include-averaged",
        action="store_true",
        help="also evaluate omega_S(x)=sqrt(4-15x/8), as a separate comparison",
    )
    p.add_argument("--rtol", type=float, default=2e-10)
    p.add_argument("--atol", type=float, default=2e-12)
    return p


def main():
    args = build_parser().parse_args()
    if args.nx < 1 or args.nk < 1:
        raise SystemExit("--nx and --nk must be positive")
    if not (0.0 < args.x_min <= args.x_max < 4.0 / 3.0):
        raise SystemExit("x bounds must lie inside (0, 4/3)")
    if args.k_min < 0.0 or args.k_max < args.k_min:
        raise SystemExit("invalid k bounds")
    if any(w <= 0.0 for w in args.omega):
        raise SystemExit("all --omega values must be positive")

    xs = np.linspace(args.x_min, args.x_max, args.nx)
    ks = np.linspace(args.k_min, args.k_max, args.nk)
    rows = []
    for x in xs:
        frequencies = [(w, "free") for w in args.omega]
        if args.include_averaged:
            frequencies.append((math.sqrt(4.0 - 15.0 * x / 8.0), "averaged"))
        for omega_s, source in frequencies:
            for k in ks:
                row = monodromy(float(x), float(k), float(omega_s),
                                rtol=args.rtol, atol=args.atol)
                row["omega_source"] = source
                rows.append(row)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    unstable = [r for r in rows if r["unstable"]]
    max_det_error = max(abs(r["det_M_minus_1"]) for r in rows)
    max_growth = max((r["gamma_F"] for r in rows), default=0.0)
    print(f"wrote {len(rows)} rows: {output}")
    print(f"unstable grid points: {len(unstable)}")
    print(f"maximum gamma_F: {max_growth:.9e}")
    print(f"maximum |det(M)-1|: {max_det_error:.3e}")
    if unstable:
        peak = max(unstable, key=lambda r: r["gamma_F"])
        print(
            "largest-growth point: "
            f"x={peak['x']:.6g}, k={peak['k']:.6g}, "
            f"omega_S={peak['omega_S']:.6g}, "
            f"a={peak['a']:.6g}, q={peak['q']:.6g}, "
            f"gamma_F={peak['gamma_F']:.9e}"
        )


if __name__ == "__main__":
    main()
