"""Factored scalar derivative and explicit wave port."""
import numpy as np, scipy.linalg as la
from numpy.polynomial.legendre import leggauss
from . import composed as C
KAP=C.KAP; mx=C.mx
class ScalarSpace:
 def __init__(self,y):
  self.y=y;self.n=len(y);self.E=np.diff(np.eye(self.n),axis=0)
  self.steps=[C.reconstruct(y,(np.arange(self.n)>i).astype(float)) for i in range(self.n-1)]
 def H(self,y,d=0):return np.stack([s(y,d) for s in self.steps],axis=-1)
 def B(self,y):return np.eye(self.n)[0]+self.H(y)@self.E
 def field(self,f):
  f=np.asarray(f);df=np.diff(f)
  def evaluate(y,d=0):return (f[0] if d==0 else 0.)+self.H(y,d)@df
  return evaluate
class StructuralRegistry(C.Registry):
 def __init__(self,A,S,raw,sc,mg):
  super().__init__(A,S,raw,sc,mg);self.space=ScalarSpace(A.cent);self.fields=[self.space.field(v) for v in sc]


def weak_port(A,S,raw,sc,mg,reg,al,scalar_boundary):
 V=reg.space;edges=np.unique(np.r_[A.ends,A.cent]);re=np.sqrt(edges);t,w=leggauss(24);r=((re[1:]+re[:-1])[:,None]/2+np.diff(re)[:,None]*t/2).ravel();wr=(np.diff(re)[:,None]*w/2).ravel();y=r*r
 shell=np.searchsorted(A.ends,y,side='right')-1
 shells=lambda z:np.bincount(shell,weights=wr*z,minlength=A.n)
 B=V.B(y);H=2*r[:,None]*V.H(y,1);J,s2,_,_,f,fy,_,g,_=reg.at(y);W=y*J;a=al(y);wv=wr*W;c=np.array([1,1,1/KAP])[:,None];p=f[1::2];fr=2*r*fy[::2];up=C.sh.up(f[::2]);L=a*(S['K']*p-up*np.array([1,1,KAP])[:,None]);L[2]+=a*g
 M=B.T@(wv[:,None]*B);Aa=B.T@((wv*a)[:,None]*B)
 def action(v,alpha):return V.E.T@(H.T@(wv*s2*alpha*(H@np.diff(v))))
 Fc=np.array(sc)[::2];Pc=np.array(sc)[1::2];R=np.sqrt(A.ends[-1]);e=V.B(R*R);JR,sR,*_=C.geom(A,S,np.array([R*R]));flux=al(R*R)*R**2*JR[0]*sR[0]*2*R*np.array([reg.fields[i](R*R,1) for i in (0,2,4)])
 physical_flux=flux.copy();flux=physical_flux if scalar_boundary is None else np.asarray(scalar_boundary,dtype=float);boundary_delta=flux-physical_flux
 rf=(Aa@Pc.T).T;rp=-np.array([action(v,a) for v in Fc])+(B.T@(wv*L).T).T+flux[:,None]*e
 ft=la.solve(M,rf.T,assume_a='pos').T;pt=la.solve(M,rp.T,assume_a='pos').T;sr=np.stack((ft,pt),axis=1).reshape(6,A.n)
 rec=[V.field(v) for v in sr];er=2*r*np.array([rec[i](y,1) for i in (0,2,4)])-(2*r*al(y,1)*p+a*2*r*fy[1::2]);work=[]
 for i in range(3):
  left=Pc[i]@M@pt[i]+Fc[i]@action(ft[i],1.);right=wr@(W*p[i]*L[i])+(Pc[i]@e)*flux[i]+wr@(W*s2*fr[i]*(2*r*al(y,1)*p[i]+er[i]));work.append(float(left-right))
 # Independent one-sided wave trace ledger on geometry interfaces.
 def trace(i,z):
  k=A.eval(S['k'],z,i=i);ell=A.eval(S['ell'],z,i=i);L0=S['L0']+z*ell;J=np.exp(3*L0-z*k/2);s2=np.exp(-2*L0+z*k)
  return al(z)*z*J*s2*2*np.sqrt(z)*np.array([reg.fields[j](z,1) for j in (0,2,4)])
 jumps=np.array([trace(i,z)-trace(i-1,z) for i,z in enumerate(A.ends[1:-1],1)])
 jump_load=(V.B(A.ends[1:-1]).T@jumps).T
 strong_load=(B.T@(wv*reg.scalar_rate(al,y)[1::2]).T).T
 wave_rem=rp-strong_load-jump_load-boundary_delta[:,None]*e
 X=a*g;xp=la.solve(M,B.T@(wv*X),assume_a='pos');exerr=B@xp-X;exE=W*p[2]*exerr/KAP;exM=-W*fr[2]*exerr/KAP
 ym=mg.q.ravel();Jm,_,_,_,fm,fym,_,gm,_=reg.at(ym);am=al(ym);dmE=am*ym*Jm*gm*fm[5]/KAP;dmM=-am*ym*Jm*gm*2*np.sqrt(ym)*fym[4]/KAP;XE=a*W*g*f[5]/KAP;XM=-a*W*g*fr[2]/KAP
 trE=shells(XE)-A.h*mg.average(dmE.reshape(mg.q.shape));trM=shells(XM)-A.h*mg.average(dmM.reshape(mg.q.shape))
 rng=np.random.default_rng(26011+A.n);vs=[np.ones(A.n),np.full(A.n,.123),A.cent,A.cent**2,*sc,rng.normal(size=A.n)];old=np.array([C.reconstruct(A.cent,v)(y,1) for v in vs]);new=np.array([V.field(v)(y,1) for v in vs]);cross=[float(la.norm(x-z)/max(la.norm(x),1.)) for x,z in zip(old,new)]
 one=np.ones(A.n); const=mx(action(one,a));assert const==0 and mx(V.field(one)(y,1))==0 and mx(V.field(one)(y,2))==0
 return sr,dict(constant_stiffness=const,constant_derivative=mx(V.field(one)(y,1)),constant_second_derivative=mx(V.field(one)(y,2)),derivative_crosswalk=max(cross),weak_field_residual=mx((M@ft.T).T-rf),weak_momentum_residual=mx((M@pt.T).T-rp),mass_min_eig=float(la.eigvalsh(M)[0]),paired_work=work,exchange_E_global=float(wr@exE),exchange_M_global=float(wr@exM),exchange_E_shell=shells(exE).tolist(),exchange_M_shell=shells(exM).tolist(),exchange_quadrature_E_shell=trE.tolist(),exchange_quadrature_M_shell=trM.tolist(),boundary_flux=flux.tolist(),physical_boundary_flux=physical_flux.tolist(),boundary_flux_delta=boundary_delta.tolist(),wave_interface_jumps=jumps.tolist(),wave_ibp_remainder=mx(wave_rem),weak_minus_strong_load=mx(rp-strong_load),scalar_rates=sr.tolist())
