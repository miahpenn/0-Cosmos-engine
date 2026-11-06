"""High-precision reference for the Issue #4 split-derivative clean rerun.

This is a test-only mpmath transcription of the independent discrete CMC
fixed-point map. It does not import production or vendor code.
"""
from __future__ import annotations

import mpmath as mp


def _mp_vector(values):
    return [mp.mpf(float(value)) for value in values]


def _derivative(values, dr, *, parity=1, order=1):
    n = len(values)
    extended = [mp.mpf("0")] * (n + 4)
    extended[2:-2] = list(values)
    extended[1] = parity * values[0]
    extended[0] = parity * values[1]
    out = [mp.mpf("0")] * n
    for i in range(n - 2):
        s = extended[i:i + 5]
        if order == 1:
            out[i] = (s[0] - 8*s[1] + 8*s[3] - s[4]) / (12*dr)
        else:
            out[i] = (-s[0] + 16*s[1] - 30*s[2] + 16*s[3] - s[4]) / (12*dr**2)
    indices = list(range(n - 5, n))
    for i in (n - 2, n - 1):
        offsets = [mp.mpf(j - i) for j in indices]
        matrix = mp.matrix([[x**power for x in offsets] for power in range(5)])
        target = mp.matrix([mp.factorial(order) if power == order else 0
                            for power in range(5)])
        weights = mp.lu_solve(matrix, target)
        out[i] = sum(weights[j] * values[indices[j]] for j in range(5)) / dr**order
    return out


def _solve_map(state, alpha, *, dps=70):
    """One exact-precision application of the discrete incoming->output map."""
    with mp.workdps(dps):
        r, dr = state["r"], state["dr"]
        a, b, X = state["a"], state["b"], state["X"]
        Aa, K = state["Aa"], state["K"]
        rho, pr, pt = state["rho"], state["pr"], state["pt"]
        n = len(r)
        ap = _derivative(a, dr, parity=1, order=1)
        bp = _derivative(b, dr, parity=1, order=1)
        Xp = _derivative(X, dr, parity=1, order=1)
        inv = [X[i]**2 / a[i] for i in range(n)]
        Q = [mp.mpf("1.5")*Aa[i]**2 + K[i]**2/3
             + 4*mp.pi*(rho[i] + pr[i] + 2*pt[i]) for i in range(n)]
        c = [-mp.mpf("0.5")*ap[i]/a[i] + bp[i]/b[i]
             - Xp[i]/X[i] + 2/r[i] for i in range(n)]
        mass = [Q[i]/inv[i] for i in range(n)]
        operator = mp.matrix(n, n)
        operator[0, 0], operator[0, 1] = 1, -1
        for i in range(1, n - 1):
            ci = c[i]/(2*dr)
            operator[i, i-1] = 1/dr**2 - ci
            operator[i, i] = -2/dr**2 - mass[i]
            operator[i, i+1] = 1/dr**2 + ci
        operator[n-1, n-1] = 1

        alpha_d1 = _derivative(alpha, dr, parity=1, order=1)
        alpha_d2 = _derivative(alpha, dr, parity=1, order=2)
        h = [mp.mpf("0.5")*ap[i]/a[i] - bp[i]/b[i]
             + Xp[i]/X[i] - 2/r[i] for i in range(n)]
        k_rhs = [-inv[i]*(alpha_d2[i] - alpha_d1[i]*h[i])
                 + alpha[i]*Q[i] for i in range(n)]
        weights = [r[i]**2 * mp.sqrt(a[i])*b[i]/X[i]**3*dr
                   for i in range(n)]
        kbar = sum(weights[i]*k_rhs[i] for i in range(n))/sum(weights)
        rhs = mp.matrix(n, 1)
        for i in range(1, n - 1):
            rhs[i] = -kbar/inv[i]
        rhs[n-1] = 1
        result = mp.lu_solve(operator, rhs)
        return [result[i] for i in range(n)]


def exact_fixed_point(state, *, dps=70):
    """Return high-precision fixed-point root and Jacobian of the affine map."""
    with mp.workdps(dps):
        n = len(state["r"])
        zero = [mp.mpf("0")] * n
        f0 = _solve_map(state, zero, dps=dps)
        columns = []
        for j in range(n):
            basis = [mp.mpf("0")] * n
            basis[j] = 1
            fj = _solve_map(state, basis, dps=dps)
            columns.append([fj[i] - f0[i] for i in range(n)])
        jac = mp.matrix(n, n)
        for j in range(n):
            for i in range(n):
                jac[i, j] = columns[j][i]
        system = mp.eye(n) - jac
        root = mp.lu_solve(system, mp.matrix(f0))
        return [root[i] for i in range(n)], jac


def as_mp_state(base):
    """Convert the fixed-grid state arrays to exact binary-float mp inputs."""
    return {
        "r": _mp_vector(base["r"]),
        "dr": mp.mpf(float(base["dr"])),
        **{key: _mp_vector(base[key]) for key in
           ("a", "b", "X", "Aa", "K", "rho", "pr", "pt")},
    }
