"""Opt-in diagnostic initializer: solve the production discrete H constraint."""
import math
import numpy as np
from .stress_energy import assemble_total_stress_energy
from . import v55_pirk_adapter as adapter

class ConvergenceError(RuntimeError):
    pass

def _residual(state,B,K_base,Aa_base):
    grid,G=state.grid,state.geometry
    r=np.asarray(grid.centers)
    _,vacuum,_=adapter.vendor_modules()
    A=1.0/np.sqrt(B)
    G.a[:]=A**(4.0/3.0); G.b[:]=A**(-2.0/3.0); G.X[:]=A**(-1.0/3.0)
    G.Lambda[:]=(np.asarray(grid.cell_derivative_fourth(G.a,parity=1))/(2*G.a**2)
      -np.asarray(grid.cell_derivative_fourth(G.b,parity=1))/(G.a*G.b)
      +2.0/r*(1.0/G.b-1.0/G.a))
    G.K[:]=K_base; G.Aa[:]=Aa_base
    raw=vacuum.constraints(grid,G)
    total=assemble_total_stress_energy(grid,G,state.scalars,state.matter)
    return np.asarray(raw["hamiltonian"])-16.0*math.pi*np.asarray(total.rho)

def discrete_consistent_state(kernel,resolution,r_max,amplitude=0.01,width=7.0,D_amplitude=1e-10,
 include_radiation=True,max_iter=12,tol=1e-13,max_cond=1e12,return_info=False):
    if isinstance(max_iter,bool) or not isinstance(max_iter,int) or max_iter<0: raise ValueError("invalid max_iter")
    if not np.isfinite(tol) or tol<=0 or not np.isfinite(max_cond) or max_cond<=1: raise ValueError("invalid solver tolerance/condition")
    state=kernel.initialize(resolution=resolution,r_max=r_max,amplitude=amplitude,width=width,D_amplitude=D_amplitude,include_radiation=include_radiation)
    K_base=np.asarray(state.geometry.K).copy(); Aa_base=np.asarray(state.geometry.Aa).copy()
    B=np.asarray(state.geometry.X,dtype=float)**6
    if B.ndim!=1 or not np.all(np.isfinite(B)) or np.any(B<=0): raise ConvergenceError("invalid initial B")
    def evaluate(candidate):
        value=np.asarray(_residual(state,candidate,K_base,Aa_base),dtype=float)
        if value.shape!=B.shape or not np.all(np.isfinite(value)): raise ConvergenceError("invalid residual")
        return value
    H=evaluate(B); hist=[float(np.max(np.abs(H)))]; max_condition=None; iterations=0
    for _ in range(max_iter):
        if hist[-1]<=tol: break
        J=np.zeros((B.size,B.size)); eps=1e-7
        for j in range(B.size):
            Bp=B.copy(); Bp[j]*=1+eps
            den=Bp[j]-B[j]
            if not np.isfinite(den) or den==0: raise ConvergenceError("invalid Jacobian step")
            J[:,j]=(evaluate(Bp)-H)/den
        condition=float(np.linalg.cond(J))
        max_condition=condition if max_condition is None else max(max_condition,condition)
        if not np.isfinite(condition) or condition>max_cond: raise ConvergenceError(f"ill-conditioned Jacobian: {condition}")
        delta=np.linalg.solve(J,-H)
        candidate=B+delta
        if not np.all(np.isfinite(candidate)) or np.any(candidate<=0): raise ConvergenceError("Newton update invalid")
        B=candidate; H=evaluate(B); hist.append(float(np.max(np.abs(H)))); iterations+=1
    if hist[-1]>tol: raise ConvergenceError(f"no convergence: {hist[-1]} > {tol}")
    evaluate(B)
    info={"iterations":iterations,"max_condition_number":max_condition,"final_max_residual":hist[-1],"tolerance":float(tol)}
    return (state,B,hist,info) if return_info else (state,B,hist)
