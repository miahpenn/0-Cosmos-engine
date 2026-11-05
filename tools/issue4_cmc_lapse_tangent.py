#!/usr/bin/env python3
"""Independent discrete CMC lapse solve and directional-equation witness.

Pure research implementation. It does not import engine modules, production
CMC code, vendor modules, or any time integrator. Inputs are already reconstructed
base arrays and independently derived instantaneous directional arrays.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np


@dataclass(frozen=True)
class CMCDiscreteResult:
    alpha: np.ndarray
    kdot_cmc: float
    operator: np.ndarray
    rhs: np.ndarray
    residual: np.ndarray
    weights: np.ndarray
    proper_volume_total: float
    potential: np.ndarray
    mass_coefficient: np.ndarray
    radial_coefficient: np.ndarray


@dataclass(frozen=True)
class CMCDirectionalResult:
    alpha_dot: np.ndarray
    kdot_cmc_dot: float
    operator_dot: np.ndarray
    rhs_dot: np.ndarray
    differentiated_residual: np.ndarray
    weight_dot: np.ndarray
    proper_volume_total_dot: float
    potential_dot: np.ndarray
    mass_coefficient_dot: np.ndarray
    radial_coefficient_dot: np.ndarray


def _vec(name, value, n=None):
    x = np.asarray(value, dtype=np.float64)
    if x.ndim != 1 or (n is not None and x.shape != (n,)):
        raise ValueError(f"{name} must be a one-dimensional array of length {n}")
    if not np.all(np.isfinite(x)):
        raise FloatingPointError(f"{name} contains non-finite values")
    return x


def _validate_grid(r, dr):
    r = _vec("radial centres", r)
    if r.size < 3 or np.any(r <= 0.0):
        raise ValueError("CMC grid requires at least three positive cell centres")
    if not math.isfinite(float(dr)) or dr <= 0.0:
        raise ValueError("CMC spacing must be finite and positive")
    expected = (np.arange(r.size, dtype=np.float64) + 0.5) * float(dr)
    tol = 64.0 * math.ulp(max(float(np.max(np.abs(r))), float(expected[-1]), float(dr)))
    if np.any(np.diff(r) <= 0.0) or float(np.max(np.abs(r - expected))) > tol:
        raise ValueError("CMC centres must be uniform and cell-centred")
    return r


def _d1(d1, values, parity, n, label):
    out = _vec(label, d1(values, parity=parity), n)
    return out


def _measure(r, dr, a, b, X):
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        w = r*r*np.sqrt(a)*b/(X**3)*dr
        total = float(np.sum(w))
    if not np.all(np.isfinite(w)) or np.any(w <= 0.0) or not math.isfinite(total) or total <= 0.0:
        raise FloatingPointError("CMC proper-volume measure is invalid")
    return w, total


def _assemble_operator(r, dr, a, b, X, Aa, K, rho, pr, pt, d1):
    n = r.size
    ap = _d1(d1, a, 1, n, "D_r a")
    bp = _d1(d1, b, 1, n, "D_r b")
    Xp = _d1(d1, X, 1, n, "D_r X")
    c = -0.5*ap/a + bp/b - Xp/X + 2.0/r
    inv = X*X/a
    Q = 1.5*Aa*Aa + K*K/3.0 + 4.0*np.pi*(rho+pr+2.0*pt)
    m = Q/inv
    if not np.all(np.isfinite(c)) or not np.all(np.isfinite(m)):
        raise FloatingPointError("CMC operator coefficients are non-finite")
    h = float(dr)
    A = np.zeros((n, n), dtype=np.float64)
    A[0, 0], A[0, 1] = 1.0, -1.0
    for i in range(1, n-1):
        ci = c[i]/(2.0*h)
        A[i, i-1] = 1.0/h**2-ci
        A[i, i] = -2.0/h**2-m[i]
        A[i, i+1] = 1.0/h**2+ci
    A[-1, -1] = 1.0
    return A, c, m, inv, Q


def solve_cmc_discrete(*, r, dr, a, b, X, Aa, K, rho, pr, pt,
                       k_rhs, d1) -> CMCDiscreteResult:
    """Solve the pinned tridiagonal CMC equation from independently supplied arrays."""
    r = _validate_grid(r, dr)
    n = r.size
    fields = {name: _vec(name, value, n) for name, value in {
        "a": a, "b": b, "X": X, "Aa": Aa, "K": K,
        "rho": rho, "pr": pr, "pt": pt, "K RHS": k_rhs,
    }.items()}
    a,b,X,Aa,K,rho,pr,pt,k_rhs = (fields[k] for k in
        ("a","b","X","Aa","K","rho","pr","pt","K RHS"))
    if np.any(a <= 0.0) or np.any(b <= 0.0) or np.any(X <= 0.0):
        raise FloatingPointError("CMC metric coefficients must be positive")
    w, wsum = _measure(r, dr, a,b,X)
    kdot = float(np.sum(w*k_rhs)/wsum)
    A,c,m,inv,Q = _assemble_operator(r,dr,a,b,X,Aa,K,rho,pr,pt,d1)
    rhs = np.zeros(n, dtype=np.float64)
    rhs[1:-1] = -kdot/inv[1:-1]
    rhs[-1] = 1.0
    try:
        alpha = np.linalg.solve(A, rhs)
    except np.linalg.LinAlgError as exc:
        raise FloatingPointError("independent CMC discrete solve is singular") from exc
    residual = A@alpha-rhs
    if not np.all(np.isfinite(alpha)) or np.min(alpha) <= 0.0:
        raise FloatingPointError("independent CMC solve did not produce a positive finite lapse")
    return CMCDiscreteResult(alpha,kdot,A,rhs,residual,w,wsum,Q,m,c)


def differentiate_cmc_discrete(*, base: CMCDiscreteResult, r, dr,
        a,b,X,Aa,K,rho,pr,pt, k_rhs,
        a_dot,b_dot,X_dot,Aa_dot,K_dot,rho_dot,pr_dot,pt_dot,k_rhs_dot,
        d1) -> CMCDirectionalResult:
    """Differentiate the *discrete* CMC system, including its proper-volume projection.

    The directional K-RHS array is an independently derived instantaneous
    derivative of the pinned geometric K equation; it is never inferred from
    the lapse residual.
    """
    r = _validate_grid(r, dr); n = r.size
    names = ("a","b","X","Aa","K","rho","pr","pt","K RHS",
             "a_dot","b_dot","X_dot","Aa_dot","K_dot","rho_dot","pr_dot","pt_dot","K RHS dot")
    vals = (a,b,X,Aa,K,rho,pr,pt,k_rhs,
            a_dot,b_dot,X_dot,Aa_dot,K_dot,rho_dot,pr_dot,pt_dot,k_rhs_dot)
    f = {k:_vec(k,v,n) for k,v in zip(names,vals)}
    a,b,X,Aa,K,rho,pr,pt,kr = (f[k] for k in names[:9])
    ad,bd,Xd,Aad,Kd,rhod,prd,ptd,krd = (f[k] for k in names[9:])
    if np.any(a<=0.0) or np.any(b<=0.0) or np.any(X<=0.0):
        raise FloatingPointError("CMC metric coefficients must be positive")

    w, wsum = _measure(r,dr,a,b,X)
    wdot = w*(0.5*ad/a+bd/b-3.0*Xd/X)
    wsumdot = float(np.sum(wdot))
    kbar = float(np.sum(w*kr)/wsum)
    kbardot = float((np.sum(wdot*kr+w*krd)-kbar*wsumdot)/wsum)

    A,c,m,inv,Q = _assemble_operator(r,dr,a,b,X,Aa,K,rho,pr,pt,d1)
    ap = _d1(d1,a,1,n,"D_r a")
    bp = _d1(d1,b,1,n,"D_r b")
    Xp = _d1(d1,X,1,n,"D_r X")
    apd = _d1(d1,ad,1,n,"D_r a dot")
    bpd = _d1(d1,bd,1,n,"D_r b dot")
    Xpd = _d1(d1,Xd,1,n,"D_r X dot")
    cdot = (-0.5*(apd/a-ap*ad/(a*a)) +
            bpd/b-bp*bd/(b*b)-Xpd/X+Xp*Xd/(X*X))
    Qdot = 3.0*Aa*Aad + (2.0/3.0)*K*Kd + 4.0*np.pi*(rhod+prd+2.0*ptd)
    invdot = inv*(2.0*Xd/X-ad/a)
    mdot = Qdot/inv-m*invdot/inv

    Adot = np.zeros((n,n),dtype=np.float64)
    for i in range(1,n-1):
        Adot[i,i-1] = -cdot[i]/(2.0*dr)
        Adot[i,i] = -mdot[i]
        Adot[i,i+1] = cdot[i]/(2.0*dr)
    rhsdot = np.zeros(n,dtype=np.float64)
    rhsdot[1:-1] = -kbardot/inv[1:-1] + base.kdot_cmc*invdot[1:-1]/(inv[1:-1]**2)
    # Fixed regular-centre and outer-normalization rows have zero derivative.
    alpha_dot = np.linalg.solve(A, rhsdot-Adot@base.alpha)
    residual = A@alpha_dot+Adot@base.alpha-rhsdot
    if not np.all(np.isfinite(alpha_dot)):
        raise FloatingPointError("CMC lapse tangent is non-finite")
    return CMCDirectionalResult(alpha_dot,kbardot,Adot,rhsdot,residual,
        wdot,wsumdot,mdot,cdot)


def central_difference_witness(base_plus: CMCDiscreteResult,
                                base_minus: CMCDiscreteResult,
                                epsilon: float,
                                analytic_dot: np.ndarray) -> dict:
    """Compare independently solved +/- systems with the differentiated witness."""
    if not math.isfinite(epsilon) or epsilon <= 0.0:
        raise ValueError("epsilon must be finite and positive")
    fd=(base_plus.alpha-base_minus.alpha)/(2.0*epsilon)
    analytic=_vec("analytic lapse tangent",analytic_dot,fd.size)
    delta=fd-analytic
    return {"finite_difference":fd,"analytic":analytic,
            "difference":delta,"max_abs_difference":float(np.max(np.abs(delta)))}
