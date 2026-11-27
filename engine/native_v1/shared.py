"""Research-only fixed mesh assembly and CMC tangent. No time step or coupled geometry RHS."""
import sys, pathlib, json, hashlib, platform
import numpy as np, scipy, scipy.linalg as la
from .candidate import Candidate
from .geometry_core import val,KAP
from .assembly import state
from .moments import MomentGrid, MomentSources, initial
from .matter_source import potential,fluid,V0,BETA_D,LAM,reconstruct
from .transport import scalar_rhs,geom

def mx(x):return float(np.max(np.abs(x)))
def relative(x,y):return float(la.norm(x-y)/max(la.norm(y),1e-300))
class Shared(Candidate):
 def __init__(self,n):
  super().__init__(n,True)
  # Rebuild all geometry-basis-dependent arrays in native_cmc.py and geometry_core.py.
  self.CP[:]=0
  for i,c in enumerate(self.cells):
   mask=(self.cy>=c['lo'])&(self.cy<c['hi']); t=2*(self.cy[mask]-c['lo'])/(c['hi']-c['lo'])-1
   self.CP[np.ix_(mask,c['pi'])]=val(c['cp'],t)
  self.pface=np.array([self.trace(v) for v in np.eye(self.np)])
  nf=n+1; self.S=self.N[:nf,:nf].copy();self.blocks=[]
  for c in self.cells:
   ids=c['ix'][2:];D=self.N[np.ix_(ids,ids)];E=self.N[:nf,ids];Di=la.cho_factor(D);self.S-=E@la.cho_solve(Di,E.T);self.blocks.append((ids,Di,E))
  self.sfscale=1/np.sqrt(np.diag(self.S));self.Schol=la.cho_factor(self.sfscale[:,None]*self.S*self.sfscale[None,:])
 def solveH(self,b):
  rhs=b[:self.n+1].copy();x=np.zeros(self.nh)
  for ids,Di,E in self.blocks:rhs-=E@la.cho_solve(Di,b[ids])
  x[:self.n+1]=self.sfscale*la.cho_solve(self.Schol,self.sfscale*rhs)
  for ids,Di,E in self.blocks:x[ids]=la.cho_solve(Di,b[ids]-E.T@x[:self.n+1])
  return x

class Source(MomentSources):
 def parts(self,y):
  J,s2,*_=geom(self.A,self.S,y);parts={}
  for name,fs in self.fluids.items():parts[name]=fluid('radiation' if name=='radiation' else 'dust',*[f(y) for f in fs],y,J,s2)[0]
  f=np.array([self.fields[i](y) for i in (0,2,4)]); P=np.array([self.fields[i](y) for i in (1,3,5)]);fy=np.array([self.fields[i](y,1) for i in (0,2,4)]); c=np.array([1,1,1/KAP])[:,None];a=np.sum(c*P**2,axis=0);b=4*s2*np.sum(c*fy**2,axis=0);u=potential(*f)
  parts['scalar']=(.5*(a+y*b)+u,b,-2*np.sum(c*P*fy,axis=0),1.5*a-.5*y*b-3*u)
  return parts
 def __call__(self,y):
  parts=self.parts(y);return tuple(sum(p[i] for p in parts.values()) for i in range(4))

def up(f):
 ed=np.exp(-2*BETA_D*f[2]);return np.array([f[0]+LAM/2*f[0]**3+1.5*LAM*f[0]*f[1]**2,-ed*f[1]+LAM/2*f[1]**3+1.5*LAM*f[1]*f[0]**2,(-V0*np.exp(-f[2])+KAP*BETA_D*ed*f[1]**2)/KAP])

def source_tangent(src,dr,ds,y,nu,lam):
 # Analytic tangent of rho+tr, same positive dust / strict radiation branches.
 A=src.A;J,s2,*_=geom(A,src.S,y);v=0.
 for name,fs in src.fluids.items():
  d,e,m=[f(y) for f in fs];dd,de,dm=[src.mg.reconstruction(z,j)(y) for z,j in zip(dr[name],(2,2,3))]
  D=d/J;E=e/J;u=m/J;Dt=dd/J-D*nu;Et=de/J-E*nu;ut=dm/J-u*nu
  q=s2*y*u*u;qt=2*lam*q+2*s2*y*u*ut
  if name=='radiation':
   er=E;ert=Et;P=(E**2-q)/(np.sqrt(4*E**2-3*q)+E);Pt=((E-P)*Et-qt/2)/(E+3*P)
  else:er=np.sqrt(D**2+q);ert=(D*Dt+qt/2)/er;P=Pt=np.zeros_like(y)
  den=er+P;v+=(ert+3*Pt+qt/den-q*(ert+Pt)/den**2)/KAP
 f=np.array([src.fields[i](y) for i in (0,2,4)]);P=np.array([src.fields[i](y) for i in (1,3,5)]);df=np.array([reconstruct(A.cent,ds[i])(y) for i in (0,2,4)]);dP=np.array([reconstruct(A.cent,ds[i])(y) for i in (1,3,5)]);c=np.array([1,1,1/KAP])[:,None]
 return v+2*(2*np.sum(c*P*dP,axis=0)-np.sum(up(f)*df,axis=0))

def run(n,mode):
 A=Shared(n);S,_=state(A,'compatible');mg=MomentGrid(A);raw,sc=initial(A,S,mg);src=Source(A,S,raw,sc,mg);B=A.lapse
 incoming=S['incoming'].copy();
 if mode=='packed':
  ic=B.sample_map@(incoming-incoming[-1])+incoming[-1]*B.one;incoming=np.r_[ic,.0013]
 S['incoming']=incoming
 before=hashlib.sha256(b''.join(np.asarray(x).tobytes() for x in list(S.values())+list(sc)+[a for vs in raw.values() for a in vs])).hexdigest()
 k,l,L0,p,K=[S[z] for z in ('k','ell','L0','p','K')];args=(k,l,L0,p,K,incoming,src)
 L,rhs,rate,x=A.elliptic_system(*args);sol=A.scaled_solve(L,rhs);al,_,diag=A.cmc(*args)
 rng=np.random.default_rng(12026+n);dk=A.projp(.002+.001*A.y);dl=A.projp(-.003+.002*A.y);dp=A.projp(.001-.0004*A.y);dL=.001;dK=.002
 dr={name:tuple(a*(.03+.01*A.cent) for a in vs) for name,vs in raw.items()};ds=[a*(.02-.01*A.cent)+1e-5*(i+1) for i,a in enumerate(sc)]
 di=.002*(A.cent-A.cent[-1]) if mode=='sampled' else np.r_[B.sample_map@(.002*(A.cent-A.cent[-1])),.00017]
 y=A.cy;dv=3*(dL+y*(A.CP@dl))-y*(A.CP@dk)/2;dm=-dL-y*(A.CP@dl)+y*(A.CP@dk)/2
 dQ=3*y*y*(A.CP@p)*(A.CP@dp)+2*K*dK/3+KAP/2*source_tangent(src,dr,ds,y,dv,dm)
 mu=x['mu'];Q=x['Q'];dmu=mu*dv;G=4*y*mu*x['s2'];dG=G*(dv+2*dm);dmQ=dmu*Q+mu*dQ;V=A.CV;Vy=A.CVy
 dmat=Vy.T@(dG[:,None]*Vy)+V.T@(dmQ[:,None]*V);dvol=dmu.sum();dvec=V.T@dmu;dpot=V.T@dmQ
 Y=A.Y;nvR=3*(dL+Y*(A.pface@dl))-Y*(A.pface@dk)/2;lmR=-dL-Y*(A.pface@dl)+Y*(A.pface@dk)/2
 if mode=='sampled':
  dic=B.sample_map@(di-di[-1])+di[-1]*B.one
  dfi=x['flux_factor']*((nvR+2*lmR)*(B.er@(x['ic']-incoming[-1]*B.one))+B.er@(dic-di[-1]*B.one))
 else:dic=di[:-1];dfi=di[-1]
 dRate=(dmu@(Q*(V@x['ic']))+mu@(dQ*(V@x['ic'])+Q*(V@dic))-dfi-rate*dvol)/x['volume']
 dLmat=np.zeros_like(L);dLmat[:-1,:-1]=dmat;drhs=np.r_[dRate*x['m']+rate*dvec-dpot,0.];dsol=A.scaled_solve(L,drhs-dLmat@sol)
 dphys=x['flux_factor']*((nvR+2*lmR)*(B.er@sol[:-1])+B.er@dsol[:-1]);ddefect=dsol[-1]-dphys
 def pert(t):
  st=dict(S,k=k+t*dk,ell=l+t*dl,L0=L0+t*dL,p=p+t*dp,K=K+t*dK,incoming=incoming+t*di);s=Source(A,st,{nm:tuple(a+t*b for a,b in zip(vs,dr[nm])) for nm,vs in raw.items()},[a+t*b for a,b in zip(sc,ds)],mg)
  z,b,r,xx=A.elliptic_system(*(st[j] for j in ('k','ell','L0','p','K','incoming')),s);solt=A.scaled_solve(z,b); phys=xx['flux_factor']*(B.er@solt[:-1]);return z,b,r,solt,xx['Q'],phys,solt[-1]-phys
 fd=[]
 for eps in (1e-3,5e-4,1e-4):
  a,b=pert(eps),pert(-eps);diff=[(u-v)/(2*eps) for u,v in zip(a,b)];fd.append(dict(eps=eps,matrix=relative(diff[0],dLmat),rhs=relative(diff[1],drhs),rate=abs(diff[2]-dRate),coeff=mx(diff[3]-dsol),Q=mx(diff[4]-dQ),physical_flux=abs(diff[5]-dphys),flux_defect=abs(diff[6]-ddefect)))
 # Scalar rate-transfer ledger: actual centre rates reconstructed, common solved lapse.
 q=mg.q.ravel();r=np.sqrt(q);J,s2,*_=geom(A,S,q);W=q*J;f=np.array([src.fields[i](q) for i in (0,2,4)]);P=np.array([src.fields[i](q) for i in (1,3,5)]);fr=2*r*np.array([src.fields[i](q,1) for i in (0,2,4)]);c=np.array([1,1,1/KAP])[:,None]
 strong=np.array(scalar_rhs(A,S,src,al,q));center=scalar_rhs(A,S,src,al,A.cent);rec=[reconstruct(A.cent,a) for a in center];hr=np.array([a(q) for a in rec]);hrf=2*r*np.array([rec[i](q,1) for i in (0,2,4)]);strf=2*r*(al(q,1)*P+al(q)*np.array([src.fields[i](q,1) for i in (1,3,5)]));ef=hr[::2]-strong[::2];ep=hr[1::2]-strong[1::2];er=hrf-strf
 DE=W*(np.sum(c*(P*ep+s2*fr*er),axis=0)+np.sum(up(f)*ef,axis=0));DM=-W*np.sum(c*(fr*ep+P*er),axis=0)
 def stress_rate(ft,Pt,frt):return W*(np.sum(c*(P*Pt+s2*fr*frt),axis=0)+np.sum(up(f)*ft,axis=0)), -W*np.sum(c*(fr*Pt+P*frt),axis=0)
 hs=stress_rate(hr[::2],hr[1::2],hrf);ss=stress_rate(strong[::2],strong[1::2],strf)
 rngv=rng.normal(size=A.nh);one=np.r_[np.ones(n+1),np.zeros(A.nh-n-1)];sparse=A.solveH(rngv);mono=la.solve(A.N,rngv,assume_a='pos');parts=src(q);orig=MomentSources(A,S,raw,sc,mg)(q)
 after=hashlib.sha256(b''.join(np.asarray(x).tobytes() for x in list(S.values())+list(sc)+[a for vs in raw.values() for a in vs])).hexdigest();assert before==after
 return dict(n=n,mode=mode,constant=mx(A.D@one),sample_crosswalk=mx(A.CP@k-np.array([A.eval(k,t) for t in y])),sparse_residual=relative(A.N@sparse,rngv),sparse_vs_monolithic=relative(sparse,mono),monolithic_residual=relative(A.N@mono,rngv),sparse_backward_error=float(la.norm(A.N@sparse-rngv)/(la.norm(A.N)*la.norm(sparse)+la.norm(rngv))),source_crosswalk=mx(np.array(parts)-np.array(orig)),alpha_min=diag['alpha_min'],cmc_res=diag['cmc_res'],cmc_strong=diag['strong_trace_norm'],cmc_weak=diag['weak_trace_dual_res'],cmc_flux_defect=diag['flux_defect'],cmc_rate=rate,tangent_rate=dRate,tangent_flux_multiplier=float(dsol[-1]),tangent_physical_flux=float(dphys),tangent_flux_defect=float(ddefect),tangent_equation=mx(L@dsol-(drhs-dLmat@sol)),fd=fd,scalar_E_closure=mx(hs[0]-ss[0]-DE),scalar_M_closure=mx(hs[1]-ss[1]-DM),scalar_E_max=mx(DE),scalar_M_max=mx(DM),scalar_E_shell=(A.h*mg.average(DE.reshape(mg.q.shape))).tolist(),scalar_M_shell=(A.h*mg.average(DM.reshape(mg.q.shape))).tolist(),immutable=before==after)
