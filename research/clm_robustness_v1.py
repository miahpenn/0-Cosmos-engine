import sys; sys.path.insert(0, "research")
import numpy as np
import clm_solver_v1 as clm
nu=0.08; T=35.0; dt=5e-4; ic=lambda x:-np.sin(x)
def bis(fn,lo,hi,tol):
    while hi-lo>tol:
        m=0.5*(lo+hi); b=fn(m)
        if b: hi=m
        else: lo=m
    return round(lo,4),round(hi,4)
def study(N,blow):
    base=clm.run(N,nu,ic,1,T=T,dt=dt,gamma=0.95,rec=True,blow=blow)
    fr=np.array(base['frec']); inb=base['inb']>0
    def rep(s):
        def f(t):
            x=t/dt; i=min(int(x),len(fr)-2); a=x-i; v=np.zeros(N,complex); v[inb]=s*((1-a)*fr[i]+a*fr[i+1]); return v
        return clm.run(N,nu,ic,1,T=T,dt=dt,fext=f,blow=blow)['blew']
    cl=bis(lambda g:clm.run(N,nu,ic,1,T=T,dt=dt,gamma=g,blow=blow)['blew'],0.94,0.97,0.002)
    op=bis(rep,1.0,1.8,0.01)
    pk=max(m for t,m in base['out'])
    print(f"N={N} blow={blow:g}: gamma_0.95 peak {pk:.0f}; closed-loop gamma_crit {cl}; replay s_crit {op}",flush=True)
study(256,1e4)
study(512,1e3)
