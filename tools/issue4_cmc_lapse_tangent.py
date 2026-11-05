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



@dataclass(frozen=True)
class CMCFixedPointResult:
    """Separately labelled self-consistent input/output lapse control."""
    base: CMCDiscreteResult
    alpha_in: np.ndarray
    linearization_evaluations: int
    fixed_point_residual: float


@dataclass(frozen=True)
class CMCFixedPointDirectionalResult:
    """Implicit directional equality for the self-consistent lapse control."""
    alpha_dot: np.ndarray
    tangent: CMCDirectionalResult
    implicit_residual: np.ndarray
    implicit_residual_max: float


def solve_cmc_lapse_fixed_point(*, r, dr, a,b,X,Aa,K,rho,pr,pt,
        initial_alpha=None, d1, d2, tolerance=2.0e-12):
    """Solve alpha_in = alpha_out by the exact affine lapse-map equation.

    At fixed geometry and matter the pinned physical K RHS and its proper-
    volume projection are linear in alpha_in, so the output map is affine.
    Its Jacobian is assembled from the exact differentiated equations, then
    the fixed-point linear system is solved directly. No iterative relaxation,
    fitting, or residual tuning is used. This is a separate analytic control,
    not accepted PIRK history.
    """
    r=_validate_grid(r,dr); n=r.size
    reference=np.ones(n,dtype=np.float64) if initial_alpha is None else _vec(
        "fixed-point reference lapse",initial_alpha,n).copy()
    if np.min(reference)<=0.0:
        raise ValueError("fixed-point reference lapse must be positive")
    zeros=np.zeros(n,dtype=np.float64)

    def map_base(alpha_in):
        kr,_=physical_k_rhs_and_direction(
            r=r,dr=dr,a=a,b=b,X=X,alpha_in=alpha_in,Aa=Aa,K=K,
            rho=rho,pr=pr,pt=pt,a_dot=zeros,b_dot=zeros,X_dot=zeros,
            alpha_in_dot=zeros,Aa_dot=zeros,K_dot=zeros,
            rho_dot=zeros,pr_dot=zeros,pt_dot=zeros,d1=d1,d2=d2)
        base=solve_cmc_discrete(
            r=r,dr=dr,a=a,b=b,X=X,Aa=Aa,K=K,rho=rho,pr=pr,pt=pt,
            k_rhs=kr,d1=d1)
        return base,kr

    reference_base,reference_kr=map_base(reference)
    jac=np.empty((n,n),dtype=np.float64)
    for j in range(n):
        basis=np.zeros(n,dtype=np.float64); basis[j]=1.0
        _,kr_dot=physical_k_rhs_and_direction(
            r=r,dr=dr,a=a,b=b,X=X,alpha_in=reference,Aa=Aa,K=K,
            rho=rho,pr=pr,pt=pt,a_dot=zeros,b_dot=zeros,X_dot=zeros,
            alpha_in_dot=basis,Aa_dot=zeros,K_dot=zeros,
            rho_dot=zeros,pr_dot=zeros,pt_dot=zeros,d1=d1,d2=d2)
        tangent=differentiate_cmc_discrete(
            base=reference_base,r=r,dr=dr,a=a,b=b,X=X,Aa=Aa,K=K,
            rho=rho,pr=pr,pt=pt,k_rhs=reference_kr,
            a_dot=zeros,b_dot=zeros,X_dot=zeros,Aa_dot=zeros,K_dot=zeros,
            rho_dot=zeros,pr_dot=zeros,pt_dot=zeros,k_rhs_dot=kr_dot,d1=d1)
        jac[:,j]=tangent.alpha_dot
    offset=reference_base.alpha-jac@reference
    system=np.eye(n)-jac
    try:
        alpha=np.linalg.solve(system,offset)
    except np.linalg.LinAlgError as exc:
        raise FloatingPointError("CMC lapse fixed-point affine system is singular") from exc
    if not np.all(np.isfinite(alpha)) or np.min(alpha)<=0.0:
        raise FloatingPointError("CMC lapse fixed point is not finite and positive")
    base,_=map_base(alpha)
    residual=float(np.max(np.abs(base.alpha-alpha)))
    if not math.isfinite(residual) or residual>tolerance:
        raise RuntimeError(
            f"CMC lapse affine fixed point failed residual gate: {residual:.6e}")
    return CMCFixedPointResult(base,alpha.copy(),n+1,residual)


def differentiate_cmc_lapse_fixed_point(*, fixed_point, r,dr,a,b,X,Aa,K,
        rho,pr,pt,a_dot,b_dot,X_dot,Aa_dot,K_dot,rho_dot,pr_dot,pt_dot,
        d1,d2):
    """Solve alpha_dot = D_state(F) + D_alpha(F)[alpha_dot].

    The state direction is supplied independently. The incoming-lapse tangent
    is not set to zero: the implicit linear system enforces equality of the
    input and output directional lapse. The Jacobian with respect to alpha_in
    is assembled by exact directional basis actions of the differentiated
    pinned equations, not by residual fitting.
    """
    r=_validate_grid(r,dr); n=r.size
    alpha=fixed_point.alpha_in
    zeros=np.zeros(n,dtype=np.float64)

    def output_direction(alpha_in_dot):
        kr,krd=physical_k_rhs_and_direction(
            r=r,dr=dr,a=a,b=b,X=X,alpha_in=alpha,Aa=Aa,K=K,
            rho=rho,pr=pr,pt=pt,a_dot=a_dot,b_dot=b_dot,X_dot=X_dot,
            alpha_in_dot=alpha_in_dot,Aa_dot=Aa_dot,K_dot=K_dot,
            rho_dot=rho_dot,pr_dot=pr_dot,pt_dot=pt_dot,d1=d1,d2=d2)
        return differentiate_cmc_discrete(
            base=fixed_point.base,r=r,dr=dr,a=a,b=b,X=X,Aa=Aa,K=K,
            rho=rho,pr=pr,pt=pt,k_rhs=kr,a_dot=a_dot,b_dot=b_dot,
            X_dot=X_dot,Aa_dot=Aa_dot,K_dot=K_dot,rho_dot=rho_dot,
            pr_dot=pr_dot,pt_dot=pt_dot,k_rhs_dot=krd,d1=d1)

    affine=output_direction(zeros).alpha_dot
    jac=np.empty((n,n),dtype=np.float64)
    for j in range(n):
        basis=np.zeros(n,dtype=np.float64); basis[j]=1.0
        jac[:,j]=output_direction(basis).alpha_dot-affine
    system=np.eye(n)-jac
    try:
        alpha_dot=np.linalg.solve(system,affine)
    except np.linalg.LinAlgError as exc:
        raise FloatingPointError("implicit CMC lapse tangent system is singular") from exc
    tangent=output_direction(alpha_dot)
    residual=alpha_dot-tangent.alpha_dot
    maximum=float(np.max(np.abs(residual)))
    if not np.all(np.isfinite(alpha_dot)) or not math.isfinite(maximum):
        raise FloatingPointError("implicit CMC lapse tangent is non-finite")
    return CMCFixedPointDirectionalResult(alpha_dot,tangent,residual,maximum)


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



SIX_SECTOR_NAMES = (
    "dark_matter", "baryons", "radiation", "S", "D", "phi",
)




def project_six_sector_stress_from_arrays(*, r, a, b, X, scalars, fluids, d1):
    """Independently reconstruct six named Eulerian stress projections.

    Inputs are frozen accepted-state arrays: scalar fields/momenta and the
    three Valencia conserved triples. This transcribes the pinned scalar and
    Valencia primitive/projection algebra without importing engine/vendor code.
    The returned per-sector arrays are rho/pr/pt/j in the locked normalization.
    """
    r=_vec("r",r); n=r.size
    if n<5 or np.any(r<=0.0) or np.any(np.diff(r)<=0.0):
        raise ValueError("six-sector projection requires increasing positive cell centres")
    a=_vec("a",a,n); b=_vec("b",b,n); X=_vec("X",X,n)
    if np.any(a<=0.0) or np.any(b<=0.0) or np.any(X<=0.0):
        raise FloatingPointError("six-sector projection requires positive metric fields")
    if not isinstance(scalars,dict) or set(scalars)!={"S","PS","D","PD","phi","Pi"}:
        raise ValueError("scalar inputs must contain exactly S, D, phi")
    if not isinstance(fluids,dict) or set(fluids)!={"dark_matter","baryons","radiation"}:
        raise ValueError("fluid inputs must contain exactly dark_matter, baryons, radiation")
    kappa=8.0*math.pi
    inv=X*X/a
    parts={}
    scalar_data=(
        ("S",0.5*np.asarray(scalars["S"],dtype=np.float64)**2,1.0,"PS"),
        ("D",-0.5*np.asarray(scalars["D"],dtype=np.float64)**2,1.0,"PD"),
        ("phi",8.242522415500654e-5*(np.exp(-np.asarray(scalars["phi"],dtype=np.float64))-0.10),kappa,"Pi"),
    )
    for name,potential,norm,momentum_key in scalar_data:
        field=_vec(f"{name} field",scalars[name],n)
        momentum=_vec(f"{name} momentum",scalars[momentum_key],n)
        field_r=_vec(f"{name} radial derivative",d1(field,parity=1),n)
        grad2=inv*field_r*field_r
        parts[name]={
            "rho":(0.5*momentum**2+0.5*grad2+potential)/norm,
            "pr":(0.5*momentum**2+0.5*grad2-potential)/norm,
            "pt":(0.5*momentum**2-0.5*grad2-potential)/norm,
            "j":(-momentum*field_r)/norm,
        }
    sqrt_gamma=np.sqrt(a)*b*r*r/(X**3)
    gamma_rr=a/(X*X); gamma_rr_inv=1.0/gamma_rr
    for name in ("dark_matter","baryons","radiation"):
        state=fluids[name]
        if not isinstance(state,dict) or set(state)!={"rest","energy_t","momentum_r"}:
            raise ValueError(f"{name} conserved state must contain rest/energy_t/momentum_r")
        rest=_vec(f"{name} rest",state["rest"],n)
        energy_t=_vec(f"{name} energy_t",state["energy_t"],n)
        momentum_r=_vec(f"{name} momentum_r",state["momentum_r"],n)
        E=energy_t/sqrt_gamma
        Sr=momentum_r/sqrt_gamma
        S2=gamma_rr_inv*Sr*Sr
        rho_fluid=np.empty(n); pressure=np.empty(n); v2=np.empty(n)
        for i in range(n):
            if name!="radiation":
                if rest[i]<=0.0:
                    raise ValueError(f"{name} rest-current must be positive at cell {i}")
                Dcons=rest[i]/sqrt_gamma[i]
                s2=S2[i]
                v2[i]=s2/(Dcons*Dcons+s2)
                W=1.0/math.sqrt(1.0-v2[i])
                rho_fluid[i]=Dcons/W
                pressure[i]=0.0
            else:
                if E[i]<0.0:
                    raise ValueError(f"radiation energy is negative at cell {i}")
                if S2[i]==0.0:
                    rho_fluid[i]=E[i]; pressure[i]=E[i]/3.0; v2[i]=0.0
                    continue
                if math.sqrt(S2[i])>=E[i]:
                    raise ValueError(f"radiation state is outside strict E>|S| cone at cell {i}")
                lo=0.0; hi=E[i]/3.0
                for _ in range(80):
                    mid=0.5*(lo+hi)
                    trial_v2=S2[i]/((E[i]+mid)**2)
                    correction=4.0*mid*trial_v2/(1.0-trial_v2)
                    residual=(3.0*mid-E[i])+correction
                    if abs(residual)<=1.0e-12*E[i]:
                        break
                    if residual>=0.0: hi=mid
                    else: lo=mid
                pressure[i]=mid; rho_fluid[i]=3.0*mid
                v2[i]=S2[i]/((E[i]+mid)**2)
        W=1.0/np.sqrt(1.0-v2)
        v_r=np.sign(Sr)*np.sqrt(v2/gamma_rr)
        h=rho_fluid+pressure
        parts[name]={
            "rho":(h*W*W-pressure)/kappa,
            "pr":(h*W*W*v2+pressure)/kappa,
            "pt":pressure/kappa,
            "j":(h*W*W*gamma_rr*v_r)/kappa,
        }
    for sector in parts:
        for key in ("rho","pr","pt","j"):
            parts[sector][key]=_vec(f"{sector} {key}",parts[sector][key],n)
    return parts

def sum_six_sector_projections(*, sector_projections, sector_projection_directions):
    """Validate and sum six named base and directional stress projections.

    Each base sector must provide rho/pr/pt/j; each direction must provide
    rho_dot/pr_dot/pt_dot/j_dot. No sector may be omitted or silently added.
    Arrays are expected to have already been independently projected under
    the locked six-sector convention.
    """
    expected=set(SIX_SECTOR_NAMES)
    if not isinstance(sector_projections,dict) or set(sector_projections)!=expected:
        raise ValueError("base sector projections must contain exactly the six locked sectors")
    if not isinstance(sector_projection_directions,dict) or set(sector_projection_directions)!=expected:
        raise ValueError("sector projection directions must contain exactly the six locked sectors")
    base_keys=("rho","pr","pt","j")
    dot_keys=("rho_dot","pr_dot","pt_dot","j_dot")
    base_total={k:None for k in base_keys}
    dot_total={k:None for k in dot_keys}
    shape=None
    for sector in SIX_SECTOR_NAMES:
        part=sector_projections[sector]
        direction=sector_projection_directions[sector]
        if not isinstance(part,dict) or set(part)!=set(base_keys):
            raise ValueError(f"{sector} base projection must contain rho/pr/pt/j")
        if not isinstance(direction,dict) or set(direction)!=set(dot_keys):
            raise ValueError(f"{sector} direction must contain rho_dot/pr_dot/pt_dot/j_dot")
        for key,value in part.items():
            arr=_vec(f"{sector} {key}",value)
            if shape is None: shape=arr.shape
            if arr.shape!=shape: raise ValueError("six-sector base projection shapes differ")
            base_total[key]=arr.copy() if base_total[key] is None else base_total[key]+arr
        for key,value in direction.items():
            arr=_vec(f"{sector} {key}",value)
            if arr.shape!=shape: raise ValueError("six-sector directional projection shapes differ")
            dot_total[key]=arr.copy() if dot_total[key] is None else dot_total[key]+arr
    return {"base":base_total,"direction":dot_total}



# Pinned instantaneous source-equation transcription. These routines are
# deliberately self-contained: no engine, production-kernel, adapter, or
# vendor import is permitted in this research module.
def _fluid_primitives_from_conserved(*, r, a,b,X,alpha,beta, state, radiation):
    r=_vec("fluid radii",r); n=r.size
    a=_vec("fluid a",a,n); b=_vec("fluid b",b,n); X=_vec("fluid X",X,n)
    alpha=_vec("fluid alpha",alpha,n); beta=_vec("fluid beta",beta,n)
    rest=_vec("fluid rest",state["rest"],n)
    energy=_vec("fluid energy_t",state["energy_t"],n)
    momentum=_vec("fluid momentum_r",state["momentum_r"],n)
    sqrtg=np.sqrt(a)*b*r*r/(X**3)
    grr=a/(X*X); invgrr=1.0/grr
    Sr=momentum/sqrtg
    q={"rho":np.empty(n),"pressure":np.zeros(n),"v_r":np.zeros(n),
       "v2":np.zeros(n),"W":np.ones(n),"E":np.empty(n),"Sr":Sr,
       "sqrtg":sqrtg,"grr":grr,"invgrr":invgrr}
    if not radiation:
        D=rest/sqrtg
        if np.any(D<=0.0): raise ValueError("dust rest density must be positive")
        s2=invgrr*Sr*Sr
        v2=s2/(D*D+s2)
        W=1.0/np.sqrt(1.0-v2)
        rho=D/W
        vr=np.sign(Sr)*np.sqrt(v2/grr)
        E=D*W
        q.update(rho=rho,v2=v2,W=W,v_r=vr,E=E)
        return q
    E=energy/sqrtg
    S2=invgrr*Sr*Sr
    p=np.empty(n); v2=np.empty(n)
    for i in range(n):
        if E[i]<0.0 or (S2[i]>0.0 and math.sqrt(S2[i])>=E[i]):
            raise ValueError(f"radiation state outside admissible cone at cell {i}")
        if S2[i]==0.0:
            p[i]=E[i]/3.0; v2[i]=0.0; continue
        lo=0.0; hi=E[i]/3.0
        for _ in range(100):
            mid=0.5*(lo+hi); vv=S2[i]/(E[i]+mid)**2
            residual=(3.0*mid-E[i])+4.0*mid*vv/(1.0-vv)
            if abs(residual)<=1.0e-12*E[i]: break
            if residual>=0.0: hi=mid
            else: lo=mid
        p[i]=mid; v2[i]=S2[i]/(E[i]+mid)**2
    W=1.0/np.sqrt(1.0-v2)
    vr=np.sign(Sr)*np.sqrt(v2/grr)
    q.update(rho=3.0*p,pressure=p,v2=v2,W=W,v_r=vr,E=E)
    return q


def _fluid_flux(metric, q, radiation):
    W=1.0/math.sqrt(1.0-metric["grr"]*q["v_r"]**2)
    v2=metric["grr"]*q["v_r"]**2
    h=q["rho"]+q["pressure"]
    E=h*W*W-q["pressure"]
    Sr=h*W*W*metric["grr"]*q["v_r"]
    Sup=Sr/metric["grr"]
    root=metric["sqrtg"]
    first=0.0 if radiation else root*q["rho"]*W*(metric["alpha"]*q["v_r"]-metric["beta"])
    return np.asarray([first,
        root*(metric["alpha"]*Sup-metric["beta"]*E),
        root*(metric["alpha"]*(h*W*W*v2+q["pressure"])-metric["beta"]*Sr)])


def _fluid_conserved_face(metric,q,radiation):
    W=1.0/math.sqrt(1.0-metric["grr"]*q["v_r"]**2)
    h=q["rho"]+q["pressure"]
    E=h*W*W-q["pressure"]
    Sr=h*W*W*metric["grr"]*q["v_r"]
    return np.asarray([0.0 if radiation else metric["sqrtg"]*q["rho"]*W,
                       metric["sqrtg"]*E,metric["sqrtg"]*Sr])


def _minmod_pair(x,y):
    if x*y<=0.0:return 0.0
    return math.copysign(min(abs(x),abs(y)),x)


def _reconstructed_fluid(q,i,side,grr):
    attrs={}
    for name in ("rho","pressure","v_r"):
        v=q[name]
        if i<=0: lo=hi=float(v[0])
        elif i>=len(v)-1: lo=hi=float(v[-1])
        else:
            slope=_minmod_pair(float(v[i]-v[i-1]),float(v[i+1]-v[i]))
            lo=float(v[i]-0.5*slope); hi=float(v[i]+0.5*slope)
        attrs[name]=lo if side=="left" else hi
    vmax=(1.0-1.0e-12)/math.sqrt(grr)
    attrs["rho"]=max(0.0,attrs["rho"])
    attrs["pressure"]=max(0.0,attrs["pressure"])
    attrs["v_r"]=max(-vmax,min(vmax,attrs["v_r"]))
    return attrs


def _hll_fluid_flux(metric,left,right,radiation):
    fl=_fluid_flux(metric,left,radiation); fr=_fluid_flux(metric,right,radiation)
    ul=_fluid_conserved_face(metric,left,radiation)
    ur=_fluid_conserved_face(metric,right,radiation)
    cs=1.0/math.sqrt(3.0) if radiation else 0.0
    def speed(q,sign):
        vhat=math.copysign(math.sqrt(metric["grr"]*q["v_r"]**2),q["v_r"])
        return metric["alpha"]*(q["v_r"]+sign*cs/math.sqrt(metric["grr"]))/(1.0+sign*vhat*cs)-metric["beta"]
    sm=min(0.0,speed(left,-1.0),speed(right,-1.0))
    sp=max(0.0,speed(left,1.0),speed(right,1.0))
    if sp<=0.0:return fr
    if sm>=0.0:return fl
    return (sp*fl-sm*fr+sm*sp*(ur-ul))/(sp-sm)


def _independent_fluid_rhs(*,r,a,b,X,alpha,beta,Aa,K,state,rho_phi_t,
                           dphi_r,d1,beta_dm=-0.04,radiation=False,dark_matter=False):
    r=_vec("r",r); n=r.size; dr=float(r[1]-r[0])
    a=_vec("a",a,n); b=_vec("b",b,n); X=_vec("X",X,n)
    alpha=_vec("alpha",alpha,n); beta=_vec("beta",beta,n)
    Aa=_vec("Aa",Aa,n); K=_vec("K",K,n)
    prim=_fluid_primitives_from_conserved(r=r,a=a,b=b,X=X,alpha=alpha,beta=beta,
                                         state=state,radiation=radiation)
    da=d1(a,parity=1); db=d1(b,parity=1); dX=d1(X,parity=1)
    dalpha=d1(alpha,parity=1); dbeta=d1(beta,parity=-1)
    grr=a/(X*X); gth=b*r*r/(X*X); sqrtg=np.sqrt(a)*b*r*r/(X**3)
    dgrr=grr*(da/a-2.0*dX/X)
    dgth=gth*(db/b+2.0/r-2.0*dX/X)
    face=np.zeros((n+1,3))
    for i in range(n-1):
        rf=0.5*(r[i]+r[i+1]); af=0.5*(a[i]+a[i+1])
        bf=0.5*(b[i]+b[i+1]); xf=0.5*(X[i]+X[i+1])
        alf=0.5*(alpha[i]+alpha[i+1]); betf=0.5*(beta[i]+beta[i+1])
        mf={"alpha":alf,"beta":betf,"grr":af/(xf*xf),
            "sqrtg":math.sqrt(af)*bf*rf*rf/(xf**3)}
        ql=_reconstructed_fluid(prim,i,"right",mf["grr"])
        qr=_reconstructed_fluid(prim,i+1,"left",mf["grr"])
        face[i+1]=_hll_fluid_flux(mf,ql,qr,radiation)
    if n>1:
        rf=r[-1]+0.5*dr
        mf={"alpha":alpha[-1],"beta":beta[-1],"grr":a[-1]/X[-1]**2,
            "sqrtg":math.sqrt(a[-1])*b[-1]*rf*rf/X[-1]**3}
    else:
        mf={"alpha":alpha[-1],"beta":beta[-1],"grr":grr[-1],"sqrtg":sqrtg[-1]}
    qout={k:float(prim[k][-1]) for k in ("rho","pressure","v_r")}
    face[-1]=_fluid_flux(mf,qout,radiation)
    W=prim["W"]; h=prim["rho"]+prim["pressure"]
    E=h*W*W-prim["pressure"]
    Sr=h*W*W*grr*prim["v_r"]
    Sup=Sr/grr
    Srr=(h*W*W*prim["v2"]+prim["pressure"])/grr
    Sth=prim["pressure"]/gth
    Krr=grr*(K/3.0+Aa); Kth=gth*(K/3.0-Aa/2.0)
    source_e=sqrtg*(alpha*(Krr*Srr+2.0*Kth*Sth)-Sup*dalpha)
    source_s=sqrtg*(-E*dalpha+Sr*dbeta+0.5*alpha*(Srr*dgrr+2.0*Sth*dgth))
    q_e=np.zeros(n); q_s=np.zeros(n); q_rest=np.zeros(n)
    if dark_matter:
        rho_raw=prim["rho"]
        q_t=beta_dm*rho_raw*rho_phi_t
        q_r=beta_dm*rho_raw*dphi_r
        ugrad=prim["W"]*((rho_phi_t-beta*dphi_r)/alpha+prim["v_r"]*dphi_r)
        q_rest=-alpha*sqrtg*beta_dm*rho_raw*ugrad
        q_e=-sqrtg*(q_t-beta*q_r)
        q_s=alpha*sqrtg*q_r
    rate=-(face[1:]-face[:-1])/dr
    if not radiation: rate[:,0]+=q_rest
    rate[:,1]+=source_e+q_e
    rate[:,2]+=source_s+q_s
    if radiation: rate[:,0]=0.0
    return {"rest":rate[:,0],"energy_t":rate[:,1],"momentum_r":rate[:,2]}


def independent_scalar_and_geometry_rates(*,r,dr,a,b,X,alpha,beta,Aa,K,Lambda,
        scalars,fluids,d1,d2,beta_dm=-0.04,lambda_m=2.0):
    """Derive the instantaneous scalar, fluid, metric, primary and Lambda rates.

    This is a direct transcription of the pinned scalar equations, Valencia
    finite-volume RHS, native CMC explicit metric equations, and pinned
    vacuum_bssn L2/L3 primary/Lambda equations. It calls no production RHS.
    """
    r=_validate_grid(r,dr); n=r.size
    a=_vec("a",a,n); b=_vec("b",b,n); X=_vec("X",X,n)
    alpha=_vec("alpha",alpha,n); beta=_vec("beta",beta,n)
    Aa=_vec("Aa",Aa,n); K=_vec("K",K,n); Lambda=_vec("Lambda",Lambda,n)
    if np.any(beta!=0.0): raise ValueError("locked CMC direction requires beta=0")
    inv=X*X/a
    ap=d1(a,parity=1); bp=d1(b,parity=1); Xp=d1(X,parity=1)
    alphap=d1(alpha,parity=1)
    sqrtg=np.sqrt(a)*b*r*r/(X**3)
    dm=fluids["dark_matter"]
    dm_D=_vec("DM rest",dm["rest"],n)/sqrtg
    dm_S=_vec("DM momentum",dm["momentum_r"],n)/sqrtg
    dm_s2=(X*X/a)*dm_S*dm_S
    dm_v2=dm_s2/(dm_D*dm_D+dm_s2)
    rho_dm=dm_D*np.sqrt(1.0-dm_v2)
    scalar_rates={}
    scalar_meta=(("S","PS",lambda f: f,lambda f:0.5*f*f,lambda f:f),
                 ("D","PD",lambda f:-f,lambda f:-0.5*f*f,lambda f:-f),
                 ("phi","Pi",lambda f:-8.242522415500654e-5*np.exp(-f),
                  lambda f:8.242522415500654e-5*(np.exp(-f)-0.10),
                  lambda f:-8.242522415500654e-5*np.exp(-f)))
    for field,momentum,Vprime,V,Vsecond in scalar_meta:
        f=_vec(field,scalars[field],n); p=_vec(momentum,scalars[momentum],n)
        fp=d1(f,parity=1); fpp=d2(f,parity=1)
        lap=inv*(fpp-fp*(0.5*ap/a-bp/b+Xp/X-2.0/r))
        ft=alpha*p
        pt=alpha*lap+alpha*K*p+inv*alphap*fp-alpha*Vprime(f)
        if field=="phi": pt=pt+alpha*beta_dm*rho_dm
        scalar_rates[field]=ft; scalar_rates[momentum]=pt
    dphi_t=scalar_rates["phi"]; dphi_r=d1(scalars["phi"],parity=1)
    fluid_rates={}
    for name in ("dark_matter","baryons","radiation"):
        fluid_rates[name]=_independent_fluid_rhs(
            r=r,a=a,b=b,X=X,alpha=alpha,beta=beta,Aa=Aa,K=K,
            state=fluids[name],rho_phi_t=dphi_t,dphi_r=dphi_r,d1=d1,
            beta_dm=beta_dm,radiation=(name=="radiation"),dark_matter=(name=="dark_matter"))
    # Reconstruct the six base projections before forming geometric sources.
    sectors=project_six_sector_stress_from_arrays(
        r=r,a=a,b=b,X=X,scalars=scalars,fluids=fluids,d1=d1)
    total={k:sum(sectors[s][k] for s in SIX_SECTOR_NAMES)
           for k in ("rho","pr","pt","j")}
    # Exact pinned vendor geometry_terms() transcription (zero KO, beta=0).
    app=d2(a,parity=1); bpp=d2(b,parity=1)
    chip=-0.5*Xp/X
    chipp=-0.5*(d2(X,parity=1)/X-(Xp/X)**2)
    Lp=d1(Lambda,parity=-1)
    metric_lambda=(1.0-a/b)/(r*r)
    radial_bracket=(0.5*app/a-a*Lp-0.75*(ap/a)**2+0.5*(bp/b)**2
        -0.5*Lambda*ap+ap/(r*b)+2.0*metric_lambda*(1.0+r*bp/b)
        +4.0*chipp-2.0*chip*(ap/a-bp/b-2.0/r))
    scalar_bracket=(0.5*app/a+bpp/b-a*Lp-(ap/a)**2+0.5*(bp/b)**2
        +2.0*(3.0-a/b)*bp/(r*b)+4.0*metric_lambda
        +8.0*(chipp+chip**2)-8.0*chip*(0.5*ap/a-bp/b-2.0/r))
    Rrr=-inv*radial_bracket; R=-inv*scalar_bracket
    alphapp=d2(alpha,parity=1)
    DrDr=inv*(alphapp-alphap*(0.5*ap/a+2.0*chip))
    lap=inv*(alphapp-alphap*(0.5*ap/a-bp/b-2.0*chip-2.0/r))
    ab=-0.5*Aa
    a_dot=-2.0*alpha*a*Aa
    b_dot=-2.0*alpha*b*ab
    X_dot=alpha*X*K/3.0
    Aa_dot=-(DrDr-lap/3.0)+alpha*(Rrr-R/3.0)+alpha*K*Aa-(16.0*math.pi/3.0)*alpha*(total["pr"]-total["pt"])
    K_dot=-lap+alpha*(Aa*Aa+2.0*ab*ab+K*K/3.0)+4.0*math.pi*alpha*(total["rho"]+total["pr"]+2.0*total["pt"])
    # Pinned Lambda L2 geometry plus adapter L3 matter source, with beta=0.
    Mgeom=(d1(Aa,parity=1)-(2.0/3.0)*d1(K,parity=1)
        +6.0*Aa*chip+(Aa-ab)*(2.0/r+bp/b))
    Aap=d1(Aa,parity=1)
    lambda_l2=(-2.0*(Aa*alphap+alpha*Aap)/a
        -4.0*alpha*(Aa-ab)/(r*b)+lambda_m*alpha*Mgeom/a)
    lambda_l3=2.0*alpha*Aa*Lambda-(8.0*math.pi*lambda_m)*alpha*total["j"]/a
    geometry_rates={"a":a_dot,"b":b_dot,"X":X_dot,
        "Aa":Aa_dot,"K":K_dot,"Lambda":lambda_l2+lambda_l3}
    return {"geometry":geometry_rates,"scalars":scalar_rates,
            "matter":fluid_rates,"sectors":sectors,"total":total,
            "diagnostics":{"M_geometry":Mgeom,"lap_alpha":lap,
                           "R":R,"Rr_r":Rrr,"rho_dm_raw":rho_dm}}


def six_sector_projection_directions(*,r,a,b,X,scalars,fluids,
        geometry_rates,scalar_rates,matter_rates,d1):
    """Analytic directional derivative of each locked sector projection.

    Includes metric dependence in scalar gradients and fluid volume/metric
    recovery. Dust and radiation are differentiated in their own conservative
    variables; no total-source residual allocation is used.
    """
    r=_vec("projection r",r); n=r.size
    a=_vec("projection a",a,n); b=_vec("projection b",b,n); X=_vec("projection X",X,n)
    ad=_vec("a_dot",geometry_rates["a"],n)
    bd=_vec("b_dot",geometry_rates["b"],n)
    Xd=_vec("X_dot",geometry_rates["X"],n)
    inv=X*X/a; invd=inv*(2.0*Xd/X-ad/a)
    out={}
    for name,field,momentum,norm,potential,potential_dot in (
        ("S","S","PS",1.0,0.5*np.asarray(scalars["S"])**2,
         np.asarray(scalars["S"])*scalar_rates["S"]),
        ("D","D","PD",1.0,-0.5*np.asarray(scalars["D"])**2,
         -np.asarray(scalars["D"])*scalar_rates["D"]),
        ("phi","phi","Pi",8.0*math.pi,
         8.242522415500654e-5*(np.exp(-np.asarray(scalars["phi"]))-0.10),
         -8.242522415500654e-5*np.exp(-np.asarray(scalars["phi"]))*scalar_rates["phi"]),
    ):
        f=_vec(f"{name} field",scalars[field],n); p=_vec(f"{name} momentum",scalars[momentum],n)
        fd=_vec(f"{name} field dot",scalar_rates[field],n); pd=_vec(f"{name} momentum dot",scalar_rates[momentum],n)
        fp=d1(f,parity=1); fpd=d1(fd,parity=1)
        grad=inv*fp*fp
        gradd=invd*fp*fp+2.0*inv*fp*fpd
        kin=0.5*p*p; kind=p*pd
        out[name]={
            "rho_dot":(kind+0.5*gradd+potential_dot)/norm,
            "pr_dot":(kind+0.5*gradd-potential_dot)/norm,
            "pt_dot":(kind-0.5*gradd-potential_dot)/norm,
            "j_dot":-(pd*fp+p*fpd)/norm,
        }
    sqrtg=np.sqrt(a)*b*r*r/(X**3)
    dlogg=0.5*ad/a+bd/b-3.0*Xd/X
    grr=a/(X*X); invgrr=1.0/grr
    dinvgrr=invgrr*(2.0*Xd/X-ad/a)
    kappa=8.0*math.pi
    for name in ("dark_matter","baryons","radiation"):
        U=fluids[name]; Ud=matter_rates[name]
        rest=_vec(f"{name} rest",U["rest"],n)
        energy=_vec(f"{name} energy",U["energy_t"],n)
        mom=_vec(f"{name} momentum",U["momentum_r"],n)
        restd=_vec(f"{name} rest dot",Ud["rest"],n)
        energyd=_vec(f"{name} energy dot",Ud["energy_t"],n)
        momd=_vec(f"{name} momentum dot",Ud["momentum_r"],n)
        S=mom/sqrtg; Sd=momd/sqrtg-S*dlogg
        S2=invgrr*S*S
        S2d=dinvgrr*S*S+2.0*invgrr*S*Sd
        if name!="radiation":
            D=rest/sqrtg; Dd=restd/sqrtg-D*dlogg
            den=D*D+S2; dend=2.0*D*Dd+S2d
            v2=S2/den
            v2d=(S2d*den-S2*dend)/(den*den)
            W=1.0/np.sqrt(1.0-v2)
            Wd=0.5*W**3*v2d
            E=D*W; Ed=Dd*W+D*Wd
            out[name]={
                "rho_dot":Ed/kappa,
                "pr_dot":(Ed*v2+E*v2d)/kappa,
                "pt_dot":np.zeros(n),
                "j_dot":Sd/kappa,
            }
        else:
            E=energy/sqrtg; Ed=energyd/sqrtg-E*dlogg
            p=np.empty(n); pd=np.empty(n); v2=np.empty(n); v2d=np.empty(n)
            for i in range(n):
                if S2[i]==0.0:
                    p[i]=E[i]/3.0
                    # At zero flux the first directional change of S^2 is zero.
                    pd[i]=Ed[i]/3.0
                    v2[i]=0.0; v2d[i]=0.0
                    continue
                lo=0.0; hi=E[i]/3.0
                for _ in range(100):
                    mid=0.5*(lo+hi); vv=S2[i]/(E[i]+mid)**2
                    F=(3.0*mid-E[i])+4.0*mid*vv/(1.0-vv)
                    if abs(F)<=1.0e-12*E[i]: break
                    if F>=0.0: hi=mid
                    else: lo=mid
                p[i]=mid; z=E[i]+p[i]; v2[i]=S2[i]/(z*z)
                vv=v2[i]
                Fp=3.0+4.0*vv/(1.0-vv)-8.0*p[i]*vv/(z*(1.0-vv)**2)
                FE=-1.0-8.0*p[i]*vv/(z*(1.0-vv)**2)
                FS=4.0*p[i]/(z*z*(1.0-vv)**2)
                pd[i]=-(FE*Ed[i]+FS*S2d[i])/Fp
                v2d[i]=S2d[i]/(z*z)-2.0*S2[i]*(Ed[i]+pd[i])/(z**3)
            out[name]={
                "rho_dot":Ed/kappa,
                "pr_dot":((Ed+pd)*v2+(E+p)*v2d+pd)/kappa,
                "pt_dot":pd/kappa,
                "j_dot":Sd/kappa,
            }
    for name in SIX_SECTOR_NAMES:
        for key in ("rho_dot","pr_dot","pt_dot","j_dot"):
            out[name][key]=_vec(f"{name} {key}",out[name][key],n)
    return out

def solve_stage_aware_cmc_tangent_six_sector(
        *, sector_projections, sector_projection_directions, **stage_args):
    """Stage-aware CMC tangent whose matter source is bound to six named sectors."""
    combined=sum_six_sector_projections(
        sector_projections=sector_projections,
        sector_projection_directions=sector_projection_directions)
    base=combined["base"]; dot=combined["direction"]
    return solve_stage_aware_cmc_tangent(
        **stage_args,rho=base["rho"],pr=base["pr"],pt=base["pt"],
        rho_dot=dot["rho_dot"],pr_dot=dot["pr_dot"],pt_dot=dot["pt_dot"])

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
