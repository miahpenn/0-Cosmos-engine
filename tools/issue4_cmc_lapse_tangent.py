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



def fourth_order_cell_derivative(values, dr: float, *, parity: int, order: int = 1) -> np.ndarray:
    """Pinned cell-centred fourth-order derivative with parity ghosts.

    Independent transcription of vendor commit
    d6052d605673ce9d82cdc99b0df79855fa2ea215, r0_operators.py,
    SphericalCellGrid._cell_derivative_fourth. No vendor module is imported.
    """
    if parity not in (-1, 1):
        raise ValueError("parity must be +1 (even) or -1 (odd)")
    if order not in (1, 2):
        raise ValueError("derivative order must be 1 or 2")
    if not math.isfinite(float(dr)) or dr <= 0.0:
        raise ValueError("grid spacing must be finite and positive")
    u = _vec("derivative input", values)
    n = u.size
    if n < 5:
        raise ValueError("fourth-order derivative requires at least five cells")
    result = np.empty_like(u)
    extended = np.empty(n + 4)
    extended[2:-2] = u
    extended[1] = parity * u[0]
    extended[0] = parity * u[1]
    for i in range(n - 2):
        stencil = extended[i:i+5]
        if order == 1:
            result[i] = (stencil[0]-8.0*stencil[1]+8.0*stencil[3]-stencil[4])/(12.0*dr)
        else:
            result[i] = (-stencil[0]+16.0*stencil[1]-30.0*stencil[2]+16.0*stencil[3]-stencil[4])/(12.0*dr**2)
    indices = np.arange(n-5,n)
    for i in (n-2,n-1):
        offsets = indices.astype(float)-float(i)
        matrix = np.vstack([offsets**power for power in range(5)])
        target = np.zeros(5)
        target[order] = float(math.factorial(order))
        weights = np.linalg.solve(matrix,target)
        result[i] = float(weights@u[indices])/dr**order
    return result

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
        wdot,wsumdot,Qdot,mdot,cdot)


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


@dataclass(frozen=True)
class StageAwareCMCTangentResult:
    """Stage-aware target projection and output-lapse tangent.

    alpha_in and alpha_in_dot construct the physical incoming K RHS. The
    returned base.alpha and tangent.alpha_dot are the distinct CMC output
    lapse and its direction.
    """
    base: CMCDiscreteResult
    tangent: CMCDirectionalResult
    k_rhs_in: np.ndarray
    k_rhs_in_dot: np.ndarray


def physical_k_rhs_and_direction(*, r, dr, a, b, X, alpha_in,
        Aa, K, rho, pr, pt, a_dot, b_dot, X_dot, alpha_in_dot,
        Aa_dot, K_dot, rho_dot, pr_dot, pt_dot, d1, d2):
    """Independently transcribe the pinned discrete L2 lapse Laplacian and its direction.

    No production/native/vendor imports. The physical K RHS is reconstructed
    from the pinned vacuum_bssn L2 lapse-Hessian stencil plus the pinned L3
    geometric and six-sector source terms. Its direction is derived from
    incoming-lapse, geometry, primary and matter directions. No k_rhs_dot
    argument exists.
    """
    r = _validate_grid(r, dr); n = r.size
    names=("a","b","X","alpha_in","Aa","K","rho","pr","pt",
           "a_dot","b_dot","X_dot","alpha_in_dot","Aa_dot","K_dot",
           "rho_dot","pr_dot","pt_dot")
    vals=(a,b,X,alpha_in,Aa,K,rho,pr,pt,a_dot,b_dot,X_dot,
          alpha_in_dot,Aa_dot,K_dot,rho_dot,pr_dot,pt_dot)
    f={k:_vec(k,v,n) for k,v in zip(names,vals)}
    (a,b,X,al,Aa,K,rho,pr,pt,ad,bd,Xd,ald,Aad,Kd,rhod,prd,ptd)=(
        f[k] for k in names)
    if np.any(a<=0.0) or np.any(b<=0.0) or np.any(X<=0.0) or np.any(al<=0.0):
        raise FloatingPointError("physical K-RHS metric and incoming lapse must be positive")
    D1a=lambda x,p: _vec("D1",d1(x,parity=p),n)
    D2a=lambda x,p: _vec("D2",d2(x,parity=p),n)
    with np.errstate(over="ignore",divide="ignore",invalid="ignore"):
        # Exact pinned vacuum_bssn.geometry_terms()["lap_alpha"] transcription:
        # inv=X^2/a; chip=-D1(X)/(2X);
        # lap=inv*(D2(alpha)-D1(alpha)*(D1(a)/(2a)-D1(b)/b-2chip-2/r)).
        # Keep derivatives on the original fields. Replacing this with
        # derivatives of A=a/X^2 and R=r*sqrt(b)/X is only continuum-equivalent;
        # discrete product/quotient rules do not commute with the stencil.
        inv=X**2/a
        inv_dot=inv*(2.0*Xd/X-ad/a)
        ap=D1a(a,1); bp=D1a(b,1); Xp=D1a(X,1)
        ap_dot=D1a(ad,1); bp_dot=D1a(bd,1); Xp_dot=D1a(Xd,1)
        dal=D1a(al,1); dald=D1a(ald,1)
        d2al=D2a(al,1); d2ald=D2a(ald,1)
        h=0.5*ap/a-bp/b+Xp/X-2.0/r
        h_dot=(0.5*(ap_dot/a-ap*ad/(a*a))
               -(bp_dot/b-bp*bd/(b*b))
               +(Xp_dot/X-Xp*Xd/(X*X)))
        lap=inv*(d2al-dal*h)
        lap_dot=(inv_dot*(d2al-dal*h)
                 +inv*(d2ald-dald*h-dal*h_dot))
        Q=1.5*Aa**2+K**2/3.0+4.0*np.pi*(rho+pr+2.0*pt)
        Qdot=3.0*Aa*Aad+(2.0/3.0)*K*Kd+4.0*np.pi*(rhod+prd+2.0*ptd)
        kr=-lap+al*Q
        krd=-lap_dot+ald*Q+al*Qdot
    if not np.all(np.isfinite(kr)) or not np.all(np.isfinite(krd)):
        raise FloatingPointError("physical incoming K RHS or its direction is non-finite")
    return kr,krd


def solve_stage_aware_cmc_tangent(*, r, dr, a,b,X,alpha_in,
        Aa,K,rho,pr,pt,a_dot,b_dot,X_dot,alpha_in_dot,
        Aa_dot,K_dot,rho_dot,pr_dot,pt_dot,d1,d2):
    """Build physical incoming K RHS, project it, solve output lapse and tangent."""
    kr,krd=physical_k_rhs_and_direction(
        r=r,dr=dr,a=a,b=b,X=X,alpha_in=alpha_in,Aa=Aa,K=K,
        rho=rho,pr=pr,pt=pt,a_dot=a_dot,b_dot=b_dot,X_dot=X_dot,
        alpha_in_dot=alpha_in_dot,Aa_dot=Aa_dot,K_dot=K_dot,
        rho_dot=rho_dot,pr_dot=pr_dot,pt_dot=pt_dot,d1=d1,d2=d2)
    base=solve_cmc_discrete(r=r,dr=dr,a=a,b=b,X=X,Aa=Aa,K=K,
        rho=rho,pr=pr,pt=pt,k_rhs=kr,d1=d1)
    tangent=differentiate_cmc_discrete(base=base,r=r,dr=dr,
        a=a,b=b,X=X,Aa=Aa,K=K,rho=rho,pr=pr,pt=pt,k_rhs=kr,
        a_dot=a_dot,b_dot=b_dot,X_dot=X_dot,Aa_dot=Aa_dot,K_dot=K_dot,
        rho_dot=rho_dot,pr_dot=pr_dot,pt_dot=pt_dot,k_rhs_dot=krd,d1=d1)
    return StageAwareCMCTangentResult(base,tangent,kr,krd)
