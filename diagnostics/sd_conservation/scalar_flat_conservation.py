"""Scalar-sector conservation in FIXED flat geometry (a=b=X=1, alpha=1, beta=0, K=0, Lambda=0).
Uses the production scalar RHS and the grid's own derivative operators.
Exact flat identity: d_t rho = (1/r^2) d_r (r^2 * P * dphi), rho=1/2(P^2+phi'^2)+V."""
import sys, copy, numpy as np
import pathlib as _p; sys.path.insert(0, str(_p.Path(__file__).resolve().parents[2]))
from engine.production_kernel import V55ProductionKernel
from engine.scalar_system import scalar_rhs_arrays, ScalarFields, cosmos_potential
from engine.discrete_consistent_initial import discrete_consistent_state
from engine import v55_pirk_adapter as ad
grid_ops, vacuum, _ = ad.vendor_modules()
k=V55ProductionKernel()
def flat_geometry(st):
    G=st.geometry; G.a[:]=1.0; G.b[:]=1.0; G.X[:]=1.0; G.alpha[:]=1.0; G.beta[:]=0.0
    G.K[:]=0.0; G.Aa[:]=0.0; G.Lambda[:]=0.0
    return G
def check(label, amp, D, phi0, pi0):
    st=k.initialize(resolution=40,r_max=40.0,amplitude=amp,width=7.0,D_amplitude=D,include_radiation=True)
    G=flat_geometry(st); g=st.grid; r=np.asarray(g.centers); n=len(r)
    z=np.zeros(n)
    sc=ScalarFields(S=np.asarray(st.scalars.S),PS=z.copy(),D=np.asarray(st.scalars.D),PD=z.copy(),
                    phi=np.full(n,phi0),Pi=np.full(n,pi0))
    rhs=scalar_rhs_arrays(g,G,sc,beta_dm=0.0,rho_dm=z.copy())
    dphi=lambda u: np.asarray(g.cell_derivative_fourth(u,parity=1))
    def rho_pair(f,p,V):
        d=dphi(f); return 0.5*(p*p+d*d)+V
    rhoS=lambda S,P: rho_pair(S,P,0.5*S*S)
    rhoD=lambda Dd,P: rho_pair(Dd,P,-0.5*Dd*Dd)
    res={}
    for name,(f,p,rr,V) in {"S":(sc.S,sc.PS,rhoS,None),"D":(sc.D,sc.PD,rhoD,None)}.items():
        eps=1e-7
        if name=="S":
            drho=(rhoS(sc.S+eps*rhs.S,sc.PS+eps*rhs.PS)-rhoS(sc.S,sc.PS))/eps
            flux=sc.PS*dphi(sc.S)
        else:
            drho=(rhoD(sc.D+eps*rhs.D,sc.PD+eps*rhs.PD)-rhoD(sc.D,sc.PD))/eps
            flux=sc.PD*dphi(sc.D)
        div=np.asarray(g.cell_derivative_fourth(r*r*flux,parity=-1))/(r*r)
        resid=drho-div
        res[name]=resid
        print(f"{label:<26} {name}: max|resid| cells0-4={np.max(np.abs(resid[:5])):.3e}  cells>=5={np.max(np.abs(resid[5:])):.3e}  "
              f"cell0={resid[0]:+.3e} cell1={resid[1]:+.3e} cell2={resid[2]:+.3e} cell3={resid[3]:+.3e}")
    return res
phi0=0.179055
pi0=0.0
check("flat, S amp 0.01, D 1e-10", 0.01, 1e-10, phi0, pi0)
check("flat, S amp 0.005, D 1e-10", 0.005, 1e-10, phi0, pi0)
check("flat, S amp 0.0 (background)", 0.0, 0.0, phi0, pi0)
