"""RESEARCH ONLY - imports no engine code. GEAR-87.
Independent nonlinear polar-areal evolution of the Gear 69/70B cavity
(G=1,R=10,mu=1): Dirichlet wall at fixed areal radius, growing-mode data to O(A^3).
Not engine code. Measures Gamma_tau = P(0)/phi(0) (a(0)=1) vs phi_c^2."""
import math, sys, numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import cumulative_simpson, solve_bvp
from scipy.special import spherical_jn
G,R,mu=1.0,10.0,1.0; k=math.pi/R; gam=math.sqrt(mu*mu-k*k)

def grid(N):
    h=R/N; return np.arange(N+1)*h, h

def metric(r,h,phi,P):
    """m, a, alpha from constraints. phi,P on vertex grid. returns a,alpha,dphi"""
    cs=CubicSpline(r,phi,bc_type=((1,0.0),(1,0.0)) if False else 'natural')
    dphi=cs.derivative()(r)
    return dphi
def solve_geom(r,h,phi,P,dphi):
    N=len(r)-1
    sp_phi=CubicSpline(r,phi); sp_d=CubicSpline(r,dphi); sp_P=CubicSpline(r,P)
    f=lambda rr,m: 4*math.pi*G*rr*rr*(0.5*sp_P(rr)**2*(1-2*G*m/rr)+0.5*(1-2*G*m/rr)*sp_d(rr)**2-0.5*mu*mu*sp_phi(rr)**2) if rr>0 else 0.0
    m=np.zeros(N+1)
    for j in range(N):
        r0=r[j]; mm=m[j]; r1=r[j]+h/2; r2=r[j+1]
        k1=f(r0,mm); k2=f(r1,mm+h/2*k1); k3=f(r1,mm+h/2*k2); k4=f(r2,mm+h*k3)
        m[j+1]=mm+h/6*(k1+2*k2+2*k3+k4)
    a=np.ones(N+1); a[1:]=1/np.sqrt(1-2*G*m[1:]/r[1:])
    # P here is a*Pi_N, Pi_N=P/a. ln alpha' = G m/(r(r-2Gm)) + 4 pi G r a^2 p_r, p_r=.5 Pi^2+.5 phi'^2/a^2+.5 mu^2 phi^2
    spm=CubicSpline(r,m)
    def g(rr,la):
        if rr<=0: return 0.0
        mm=spm(rr); a2=1/(1-2*G*mm/rr); Pi=sp_P(rr)/math.sqrt(a2)
        pr=0.5*Pi**2+0.5*sp_d(rr)**2/a2+0.5*mu*mu*sp_phi(rr)**2
        return G*mm/(rr*(rr-2*G*mm))+4*math.pi*G*rr*a2*pr
    la=np.zeros(N+1); la[N]=-math.log(a[N])
    for j in range(N,0,-1):
        r0=r[j]; r2=r[j-1]; r1=r0-h/2; hh=-h
        k1=g(r0,la[j]); k2=g(r1,la[j]+hh/2*k1); k3=g(r1,la[j]+hh/2*k2); k4=g(r2,la[j]+hh*k3)
        la[j-1]=la[j]+hh/6*(k1+2*k2+2*k3+k4)
    return m,a,np.exp(la)

def rhs(r,h,phi,P):
    dphi=metric(r,h,phi,P)
    m,a,al=solve_geom(r,h,phi,P,dphi)
    F=r*r*al*dphi/a
    dF=CubicSpline(r,F).derivative()(r)
    Pd=np.zeros_like(r); Pd[1:]=dF[1:]/r[1:]**2+a[1:]*al[1:]*mu*mu*phi[1:]
    Pd[0]=3*dF[0]/ (3*h)*0  # placeholder replaced below
    # centre: (1/r^2)F' -> limit via value at first cell extrapolated (F ~ r^3): use 3*F[1]/r1^3
    Pd[0]=3*F[1]/r[1]**3+a[0]*al[0]*mu*mu*phi[0]
    phid=al*P/a
    phid[-1]=0.0; Pd[-1]=0.0
    return phid,Pd,a,al

def w_profile(r):
    g2=gam**2
    rr=np.linspace(0,R,20001); u=spherical_jn(0,k*rr); up=-k*spherical_jn(1,k*rr)
    upp=np.empty_like(rr); upp[1:]=-2*up[1:]/rr[1:]-k*k*u[1:]; upp[0]=-k*k/3
    rho2=0.5*(g2*u*u+up*up-mu*mu*u*u)
    m2=cumulative_simpson(4*math.pi*G*rr*rr*rho2,x=rr,initial=0.0)
    a2=np.zeros_like(rr); a2[1:]=m2[1:]/rr[1:]
    ap=np.zeros_like(rr); ap[1:]=a2[1:]/rr[1:]+2*math.pi*G*rr[1:]*(g2*u[1:]**2+up[1:]**2+mu*mu*u[1:]**2)
    ia=cumulative_simpson(ap,x=rr,initial=0.0); al2=-a2[-1]-(ia[-1]-ia)
    q,s=al2-a2,al2+a2
    qp=CubicSpline(rr,q).derivative()(rr)
    div=np.empty_like(rr); div[1:]=2*q[1:]*up[1:]/rr[1:]+qp[1:]*up[1:]+q[1:]*upp[1:]; div[0]=-k*k*q[0]
    S3=-div-mu*mu*s*u-3*g2*q*u
    src=CubicSpline(rr,S3); lam=mu*mu-9*g2
    def f(x,y):
        c=src(x)-lam*y[0]; xs=np.where(x==0,1.0,x)
        return np.vstack((y[1],np.where(x==0,c/3,c-2*y[1]/xs)))
    mesh=np.linspace(0,R,1000)
    sol=solve_bvp(f,lambda ya,yb:np.array([ya[1],yb[0]]),mesh,np.zeros((2,1000)),tol=1e-8,max_nodes=20000)
    assert sol.success
    return sol.sol(r)[0], al2

def initial(r,h,A,use_w=True):
    u=spherical_jn(0,k*r); w,_=w_profile(r)
    if not use_w: w=0*w
    phi=A*u+A**3*w; phid=gam*A*u+3*gam*A**3*w
    P=phid.copy()
    for _ in range(4):
        dphi=metric(r,h,phi,P)
        m,a,al=solve_geom(r,h,phi,P,dphi)
        P=a*phid/al
    P[-1]=0.0
    return phi,P

def run(N,A,T,cfl=0.25,use_w=True,every=1):
    r,h=grid(N); phi,P=initial(r,h,A,use_w); dt=cfl*h
    out=[]; t=0.0
    n=int(round(T/dt))
    for s in range(n+1):
        if s%every==0: out.append((t,phi[0],P[0]))
        if s==n: break
        k1=rhs(r,h,phi,P); 
        k2=rhs(r,h,phi+dt/2*k1[0],P+dt/2*k1[1])
        k3=rhs(r,h,phi+dt/2*k2[0],P+dt/2*k2[1])
        k4=rhs(r,h,phi+dt*k3[0],P+dt*k3[1])
        phi=phi+dt/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]); P=P+dt/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]); t+=dt
    return np.array(out)

if __name__=="__main__":
    N=int(sys.argv[1]); T=float(sys.argv[2])
    for A in map(float,sys.argv[3:]):
        d=run(N,A,T)
        t,pc,Pc=d.T
        gt=Pc/pc   # Gamma_tau
        # report at a few times
        idx=[0,len(t)//4,len(t)//2,3*len(t)//4,len(t)-1]
        print("A",A,"N",N)
        for i in idx:
            x=pc[i]**2
            print("  t=%.3f phi_c=%.4e Gamma_tau=%.10f  (Gamma_tau/gam-1)/phi_c^2=%.4f"%(t[i],pc[i],gt[i],(gt[i]/gam-1)/x))

# Usage: python research/gear69_polar_areal_nonlinear_v1.py N T A1 A2 ...
# Reports Gamma_tau = P(0)/phi(0) (growth rate per centre proper time, a(0)=1) and
# (Gamma_tau/gamma - 1)/phi_c^2, to be compared with (2 w(0) - alpha2(0)) = 4.5826
# (G=1, R=10, mu=1; GEAR-86). Subtract the A->0 numerical rate for quantitative work
# (see GEAR-87 for the series procedure and its limits).
