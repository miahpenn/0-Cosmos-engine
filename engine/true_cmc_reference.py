import numpy as np, math, json, time
from scipy.linalg import solve_banded
import GEAR_V5_5_REFERENCE_PIRK_UNIFIED as q

# True stage-wise CMC slicing on the validated V5.5 matter/geometry system.
# Spatial shift is held at zero in this first isolated gauge gate. No physical
# equations are changed; only the lapse is obtained from the ADM CMC condition.

G=1.0

def matter_Q(g,s,fields):
    e,sr,sa,j=q.matter_projection(g,s,fields)
    A2=1.5*s.Aa*s.Aa
    return A2+s.K*s.K/3.0+4*np.pi*(e+sr+2*sa),e,sr,sa,j

def target_Kdot(g,s,fields,outer_frac=.20):
    l2=q.primary_l2(g,s); l3=q.primary_l3(g,s,fields)
    raw=l2['K']+l3['K']
    p0=int((1.0-outer_frac)*g.n)
    # Outer region supplies the dynamically determined cosmological CMC target.
    return float(np.mean(raw[p0:]))

def cmc_lapse(g,s,fields,outer_frac=.20):
    n=g.n; h=g.dr; r=g.centers
    ap=q.D(g,s.a,1); bp=q.D(g,s.b,1); Xp=q.D(g,s.X,1)
    inv=s.X**2/s.a
    c=-0.5*ap/s.a + bp/s.b - Xp/s.X + 2.0/r
    Q,_,_,_,_=matter_Q(g,s,fields)
    kd=target_Kdot(g,s,fields,outer_frac)
    # D^2 alpha + c alpha' - Q/inv alpha = -Kdot/inv
    A=np.zeros((3,n)); rhs=np.zeros(n)
    m=Q/inv; bvec=-kd/inv
    # regular center: alpha'(0)=0
    A[1,0]=1.0; A[0,1]=-1.0
    for i in range(1,n-1):
        cim=c[i]/(2*h)
        A[2,i-1]=1/h**2-cim
        A[1,i]=-2/h**2-m[i]
        A[0,i+1]=1/h**2+cim
        rhs[i]=bvec[i]
    # outer normalization fixes the remaining gauge scale
    A[1,-1]=1.0; rhs[-1]=1.0
    alpha=solve_banded((1,1),A,rhs,check_finite=False)
    return alpha,kd,float(np.min(alpha)),float(np.max(alpha))

def explicit_cmc(g,s,fields):
    out={}
    # beta=0 gauge: no artificial driver; the coordinate condition itself is the gauge.
    out['a']=-2*s.alpha*s.a*s.Aa+q.ko(g,s.a,1)
    Ab=-0.5*s.Aa
    out['b']=-2*s.alpha*s.b*Ab+q.ko(g,s.b,1)
    out['X']=s.X*(s.alpha*s.K)/3.0+q.ko(g,s.X,1)
    out['alpha']=np.zeros_like(s.alpha)
    out['beta']=np.zeros_like(s.beta)
    return out

def lambda_l2_cmc(g,s,lm=2.0):
    # beta=0 specialisation of the exact reference-metric Lambda L2 operator.
    Aa=s.Aa; Ab=-0.5*Aa; ap=q.D(g,s.alpha,1); Aap=q.D(g,Aa,1)
    return -2*(Aa*ap+s.alpha*Aap)/s.a -4*s.alpha*(Aa-Ab)/(g.centers*s.b)+lm*s.alpha*q.momentum_constraint(g,s)/s.a

def evolve_fields_stage(g,s,fields,dt):
    S,PS,Df,PD,phi,Pi,rdm,rb=fields
    q0=[]
    for f,p,Vp in ((S,PS,S),(Df,PD,-Df),(phi,Pi,q.dVc(phi))):
        q0.append(q.scalar_parts(g,s,f,p,Vp))
    # zero shift advection
    f1=[f+dt*z[0] for f,z in zip((S,Df,phi),q0)]
    rdm1=rdm+dt*(s.alpha*s.K*rdm)
    rb1=rb+dt*(s.alpha*s.K*rb)
    return [f1[0],f1[1],f1[2],None,None,None,rdm1,rb1], q0

def step_true_cmc(g,s,fields,dt,eta=2.0):
    S,PS,Df,PD,phi,Pi,rdm,rb=fields
    # Gauge at the beginning of the step.
    s=s.copy(); s.alpha=cmc_lapse(g,s,fields)[0]; s.beta.fill(0.0); s.B.fill(0.0)
    e0=explicit_cmc(g,s,fields)

    q0=[]
    for f,p,Vp in ((S,PS,S),(Df,PD,-Df),(phi,Pi,q.dVc(phi))):
        q0.append(q.scalar_parts(g,s,f,p,Vp))
    f1=[S+dt*q0[0][0], PS, Df+dt*q0[1][0], PD, phi+dt*q0[2][0], Pi, rdm+dt*(s.alpha*s.K*rdm), rb+dt*(s.alpha*s.K*rb)]

    # Explicit geometry stage, then solve CMC lapse for the actual stage.
    st=s.copy(); st.a=s.a+dt*e0['a']; st.b=s.b+dt*e0['b']; st.X=s.X+dt*e0['X']; st.beta.fill(0.0); st.B.fill(0.0)
    st.alpha=cmc_lapse(g,st,f1)[0]

    l20=q.primary_l2(g,s)
    l21=q.primary_l2(g,st)
    l30=q.primary_l3(g,s,fields)

    ns=s.copy(); ns.a=st.a; ns.b=st.b; ns.X=st.X; ns.alpha=st.alpha; ns.beta.fill(0.0); ns.B.fill(0.0)
    ns.Aa=s.Aa+dt*(0.5*l20['Aa']+0.5*l21['Aa']+l30['Aa'])
    ns.K=s.K+dt*(0.5*l20['K']+0.5*l21['K']+l30['K'])
    ns.Lambda=s.Lambda+dt*(0.5*lambda_l2_cmc(g,s)+0.5*lambda_l2_cmc(g,ns)+q.lambda_l3_matter(g,s,fields))
    ns.B.fill(0.0)

    # Stage momenta.
    ql20=[]
    for f,p,Vp in ((S,PS,S),(Df,PD,-Df),(phi,Pi,q.dVc(phi))): ql20.append(q.scalar_parts(g,s,f,p,Vp))
    ql21=[]
    for f,p,Vp in ((f1[0],PS,S),(f1[2],PD,-Df),(f1[4],Pi,q.dVc(f1[4]))): ql21.append(q.scalar_parts(g,ns,f,p,Vp))
    ps1=[PS+dt*(0.5*ql20[0][1]+0.5*ql21[0][1]+ql20[0][2]),
         PD+dt*(0.5*ql20[1][1]+0.5*ql21[1][1]+ql20[1][2]),
         Pi+dt*(0.5*ql20[2][1]+0.5*ql21[2][1]+ql20[2][2]+q.beta_dm*rdm)]
    fields_stage=[f1[0],ps1[0],f1[2],ps1[1],f1[4],ps1[2],f1[6],f1[7]]

    e1=explicit_cmc(g,ns,fields_stage)
    q1=[]
    for f,p,Vp in ((fields_stage[0],fields_stage[1],fields_stage[0]),(fields_stage[2],fields_stage[3],-fields_stage[2]),(fields_stage[4],fields_stage[5],q.dVc(fields_stage[4]))): q1.append(q.scalar_parts(g,ns,f,p,Vp))

    sf=s.copy(); sf.a=s.a+0.5*dt*(e0['a']+e1['a']); sf.b=s.b+0.5*dt*(e0['b']+e1['b']); sf.X=s.X+0.5*dt*(e0['X']+e1['X']); sf.beta.fill(0.0); sf.B.fill(0.0)
    # Use the stage lapse as predictor, then recompute the CMC lapse once the final U is formed.
    sf.alpha=ns.alpha.copy()

    ff=[S+0.5*dt*(q0[0][0]+q1[0][0]),None,Df+0.5*dt*(q0[1][0]+q1[1][0]),None,phi+0.5*dt*(q0[2][0]+q1[2][0]),None,
        rdm+0.5*dt*(s.alpha*s.K*rdm+ns.alpha*ns.K*fields_stage[6]),
        rb+0.5*dt*(s.alpha*s.K*rb+ns.alpha*ns.K*fields_stage[7])]

    oldv=ns.copy(); oldv.a=sf.a.copy(); oldv.b=sf.b.copy(); oldv.X=sf.X.copy(); oldv.alpha=sf.alpha.copy(); oldv.beta.fill(0.0)
    l2f=q.primary_l2(g,oldv); l31=q.primary_l3(g,ns,fields_stage)
    sf.Aa=s.Aa+0.5*dt*(l20['Aa']+l2f['Aa']+l30['Aa']+l31['Aa'])
    sf.K=s.K+0.5*dt*(l20['K']+l2f['K']+l30['K']+l31['K'])
    sf.Lambda=s.Lambda+0.5*dt*(lambda_l2_cmc(g,s)+lambda_l2_cmc(g,sf)+q.lambda_l3_matter(g,s,fields)+q.lambda_l3_matter(g,ns,fields_stage))
    sf.beta.fill(0.0); sf.B.fill(0.0)

    # Final matter momenta.
    ff[1]=PS+0.5*dt*(ql20[0][1]+q.scalar_parts(g,sf,ff[0],ps1[0],ff[0])[1]+ql20[0][2]+q1[0][2])
    ff[3]=PD+0.5*dt*(ql20[1][1]+q.scalar_parts(g,sf,ff[2],ps1[1],-ff[2])[1]+ql20[1][2]+q1[1][2])
    ff[5]=Pi+0.5*dt*(ql20[2][1]+q.scalar_parts(g,sf,ff[4],ps1[2],q.dVc(ff[4]))[1]+ql20[2][2]+q1[2][2]+q.beta_dm*(rdm+fields_stage[6]))

    # True end-of-step CMC solve: this is the gauge actually carried forward.
    sf.alpha=cmc_lapse(g,sf,tuple(ff))[0]
    sf.beta.fill(0.0); sf.B.fill(0.0)
    return sf,tuple(ff)

def run(N=40,T=24,cfl=.06):
    g=q.Grid(N,40.0); s,fields=q.make_initial(g,.01,7.); dt=cfl*g.dr; t=0.; steps=int(T/dt); hist=[]
    s.alpha=cmc_lapse(g,s,fields)[0]; s.beta.fill(0.0); s.B.fill(0.0)
    for n in range(steps+1):
        H,M,C,det=q.constraint(g,s,fields)
        if n%max(1,steps//12)==0 or n==steps:
            hist.append((t,float(np.min(s.alpha)),float(np.max(np.abs(H[2:]))),float(np.max(np.abs(M[2:]))),float(np.max(np.abs(C[2:]))),float(np.min(s.a)),float(np.min(s.b)),float(np.min(s.X))))
        if n==steps: break
        s,fields=step_true_cmc(g,s,fields,dt); t+=dt
        if not np.all(np.isfinite(s.alpha)) or np.min(s.alpha)<=0 or not np.all(np.isfinite(s.a)):
            return dict(ok=False,t=t,reason='invalid',hist=hist)
    return dict(ok=True,t=t,reason=None,hist=hist)

if __name__=='__main__':
    for N in (40,60):
        st=time.time(); o=run(N,24,.06); print('TRUE_CMC',N,json.dumps({'ok':o['ok'],'t':o['t'],'reason':o['reason'],'last':o['hist'][-1]})); print('sec',time.time()-st)