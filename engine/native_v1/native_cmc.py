"""Native common weak lapse/trace operator, research-only, no timesteps.
Unchanged geometry assembly inherited. C1 piecewise cubic in y, physical flux
unknown at outer face; last-centre value is the sole lapse normalization.
"""
import numpy as np, scipy.linalg as la
from .geometry_core import Assembly, KAP, val
from scipy.special import roots_jacobi
from numpy.polynomial import polynomial as poly

class Lapse:
    def __init__(self,space,c,flux):self.space=space;self.coeff=c;self.flux=flux
    def __call__(self,y,order=0):return self.space.basis(y,order)@(self.coeff-self.space.one)+(1. if order==0 else 0.)
    def packed(self):return np.r_[self.coeff,self.flux]
    def minimum(self):
        v=[]
        for i,dy in enumerate(np.diff(self.space.ends)):
            c=self.coeff[2*i:2*i+4]; C=self.space.local(i); p=C@c
            t=[0.,1.]+[float(x.real) for x in poly.polyroots(poly.polyder(p)) if abs(x.imag)<1e-12 and 0<x.real<1]
            v.extend(poly.polyval(t,p))
        return float(min(v))

class Hermite:
    def __init__(self,A):
        self.ends=A.ends; self.n=A.n; self.dim=2*(A.n+1); widths=np.diff(self.ends);self.scale=np.r_[widths[0],np.minimum(widths[:-1],widths[1:]),widths[-1]]
        self.V=self.basis(A.y);self.Vy=self.basis(A.y,1);self.b=self.basis(A.cent[-1]);self.e=self.basis(A.Y);self.er=self.basis(A.Y,1);self.one=np.tile([1.,0.],A.n+1)
        self.sample_map=np.zeros((self.dim,A.n))
        # independent adapter from pinned centre samples; fixed local degree 3.
        for i,t in enumerate(self.ends):
            j=max(0,min(A.n-4,i-2));ids=np.arange(j,j+4);x=A.cent[ids]-t;sc=max(abs(x));V=np.array([(x/sc)**m for m in range(4)])
            self.sample_map[2*i,ids]=la.solve(V,[1.,0.,0.,0.]);self.sample_map[2*i+1,ids]=la.solve(V,[0.,self.scale[i]/sc,0.,0.])
    def local(self,i):
        d=self.ends[i+1]-self.ends[i];return np.array([[1,0,0,0],[0,d/self.scale[i],0,0],[-3,-2*d/self.scale[i],3,-d/self.scale[i+1]],[2,d/self.scale[i],-2,d/self.scale[i+1]]],float)
    def basis(self,y,order=0):
        yy=np.asarray(y);flat=yy.reshape(-1);out=np.zeros((len(flat),self.dim));ix=np.clip(np.searchsorted(self.ends,flat,side='right')-1,0,self.n-1)
        for i in np.unique(ix):
            mask=np.where(ix==i)[0];d=self.ends[i+1]-self.ends[i];t=(flat[mask]-self.ends[i])/d;C=poly.polyder(self.local(i),m=order,axis=0)/d**order
            out[np.ix_(mask,np.arange(2*i,2*i+4))]=poly.polyvander(t,len(C)-1)@C
        return out.reshape(yy.shape+(self.dim,))

class NativeAssembly(Assembly):
    def __init__(self,*a,**kw):
        super().__init__(*a,**kw);self.lapse=Hermite(self);self.pface=np.array([self.trace(v) for v in np.eye(self.np)])
        # Dedicated proper-volume quadrature: sqrt(y)/2 dy, not reweighting origin Jacobi rule by 1/y.
        yy=[];ww=[];parts=[]
        for c,mid in zip(self.cells,self.cent):
            qy=[];qw=[]
            for lo,hi in [(c['lo'],mid),(mid,c['hi'])]:
                t,w=roots_jacobi(len(c['y'])//2,0,.5 if lo==0 else 0);y=lo+(hi-lo)*(t+1)/2
                w=w*((hi-lo)/2)**1.5/2 if lo==0 else w*(hi-lo)/4*np.sqrt(y)
                qy.extend(y);qw.extend(w)
            qy=np.array(qy); sl=slice(len(yy),len(yy)+len(qy));yy.extend(qy);ww.extend(qw);parts.append((c,sl,qy))
        self.cy=np.array(yy);self.cw=np.array(ww);self.CP=np.zeros((len(yy),self.np))
        for c,sl,qy in parts:self.CP[sl,c['pi']]=val(c['cp'],2*(qy-c['lo'])/(c['hi']-c['lo'])-1)
        self.CV=self.lapse.basis(self.cy);self.CVy=self.lapse.basis(self.cy,1)
    def elliptic_system(self,k,ell,L0,p,K,incoming,source):
        B=self.lapse;y=self.cy;kv=self.CP@k;lv=self.CP@ell;pv=self.CP@p;L=L0+y*lv;J=np.exp(3*L-y*kv/2);s2=np.exp(-2*L+y*kv);rho,ani,jhat,tr=source(y);Q=1.5*(y*pv)**2+K*K/3+KAP/2*(rho+tr)
        mu=self.cw*J;grad=4*y*mu*s2; V=self.CV; Vy=self.CVy
        mass=V.T@(mu[:,None]*V); pot=V.T@(mu*Q);mat=Vy.T@(grad[:,None]*Vy)+V.T@((mu*Q)[:,None]*V);m=V.T@mu;volume=mu.sum()
        kR=self.pface@k;lR=self.pface@ell;LR=L0+self.Y*lR;JR=np.exp(3*LR-self.Y*kR/2);s2R=np.exp(-2*LR+self.Y*kR);flux_factor=2*self.Y**1.5*JR*s2R
        if len(incoming)==self.n:
            ic=B.sample_map@(incoming-incoming[-1])+incoming[-1]*B.one; fi=flux_factor*(B.er@(ic-incoming[-1]*B.one))
        elif len(incoming)==B.dim+1:ic=incoming[:-1];fi=incoming[-1]
        else:raise ValueError('incoming representation not identified')
        rate_matrix=(B.one@(mat@(ic-B.one)+pot-B.e*fi))/volume
        # constant test has exact zero derivative; evaluate equivalent rate separately as a check
        rate=(mu@(Q*(V@ic))-fi)/volume
        Lmat=np.zeros((B.dim+1,B.dim+1),dtype=np.result_type(mat,rate));Lmat[:-1,:-1]=mat;Lmat[:-1,-1]=-B.e;Lmat[-1,:-1]=B.b;rhs=np.r_[rate*m-pot,0.]
        return Lmat,rhs,rate,dict(mat=mat,mass=mass,m=m,volume=volume,mu=mu,Q=Q,s2=s2,flux_factor=flux_factor,ic=ic,fi=fi,rate_direct=rate_matrix,pot=pot)
    @staticmethod
    def scaled_solve(L,rhs):
        d=1/np.sqrt(np.diag(L[:-1,:-1]).real);c=1/np.linalg.norm(d*L[:-1,-1]);r=1/np.linalg.norm(L[-1,:-1]*d);cs=np.r_[d,c];rs=np.r_[d,r];return cs*la.solve(rs[:,None]*L*cs[None,:],rs*rhs)
    def cmc(self,k,ell,L0,p,K,incoming,source):
        B=self.lapse;L,f,rate,x=self.elliptic_system(k,ell,L0,p,K,incoming,source);sol=self.scaled_solve(L,f);al=Lapse(B,sol[:-1]+B.one,sol[-1]);minimum=al.minimum()
        if minimum<=0 or not np.isfinite(minimum):raise ValueError('nonpositive CMC lapse')
        weak=x['mat']@sol[:-1]+x['pot']-B.e*al.flux;local=la.solve(x['mass'],weak,assume_a='pos');res=local-rate*B.one
        # Strong cell interior operator: actual metric derivative, no reduction or weak rate reused.
        y=self.y;kv=self.P@k;lv=self.P@ell;pv=self.P@p;L0y=L0+y*lv;Js=np.exp(3*L0y-y*kv/2);ss=np.exp(-2*L0y+y*kv);rho,ani,jhat,tr=source(y);Qs=1.5*(y*pv)**2+K*K/3+KAP/2*(rho+tr);theta=kv+y*(self.Py@k);V=lv+y*(self.Py@ell);raw=-ss*(4*y*al(y,2)+(6+2*y*(theta+2*V))*al(y,1))+Qs*al(y);tau=rate-raw
        # Full weak projection uses the multiplier, while physical face derivative is separately measured.
        diag=dict(cmc_res=float(np.max(abs(L@sol-f))),weak_trace_norm=float(np.sqrt(res@x['mass']@res)),weak_trace_dual_res=float(np.max(abs(weak-rate*x['m']))),strong_trace_norm=float(np.sqrt((2*self.wh*Js/y)@(tau*tau))),alpha_min=minimum,alpha_face=float(al(self.Y)),last_centre=float(al(self.cent[-1])),face_flux=float(al.flux),physical_flux=float(x['flux_factor']*al(self.Y,1)),flux_defect=float(al.flux-x['flux_factor']*al(self.Y,1)),rate_formula_error=float(rate-x['rate_direct']),native_reprojection_error=float((B.one@weak)/x['volume']-rate),rate=float(rate))
        return al,rate,diag
