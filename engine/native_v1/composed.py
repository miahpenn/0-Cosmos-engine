"""Static native research rate, restricted V1 witness; not a production backend."""
import hashlib
import numpy as np, scipy.linalg as la
from numpy.polynomial.legendre import legder,legval
from . import shared as sh, common as cm
from .geometry_core import Assembly,KAP
from . import geometry_core as gc
from .assembly import state
from .moments import MomentGrid,initial,MomentSources,matter_rhs
from .transport import geom,recover,beta_eff,scalar_rhs
from .matter_source import potential,reconstruct,V0,BETA_D,LAM
mx=sh.mx
class Integrated(sh.Shared):
 def __init__(self,n):super().__init__(n);cm.setup(self)
 def B(self,p,hat):
  i,d,*_=cm.split(self,p,hat);return i+d
 def joint(self,k,l,J,RV,kh,eh,alpha):
  a=-2*l/3+4*k/9;z=4*l-2*k/3
  ih,dh,il,*_=cm.split(self,a,-2*eh/3+4*kh/9);iz,dz,zil,zres,*_=cm.split(self,z,4*eh-2*kh/3)
  load=np.zeros(self.np)
  for i,c in enumerate(self.cells):
   n=len(c['ch']);dy=c['hi']-c['lo'];t=2*(c['y']-c['lo'])/dy-1
   v=cm.pad(il[i],n)+4*cm.hcoef(c,J)/3-8*cm.hcoef(c,RV)/3+cm.pad(cm.pcoef(c,z),n)
   d=legval(t,legder(v)*2/dy)+c['Hy']@dh[c['ix']]
   q=-legval(t,zres[i])+c['H']@dz[c['ix']];al=alpha(c['y'])
   load[c['pi']]=c['F'].T@(al*(c['wp']*d+c['wh']*q))
  return load
 rhs=gc.split_rhs
class Registry(sh.Source):
 def __init__(self,*args):super().__init__(*args);self.cache={}
 def at(self,y):
  y=np.asarray(y);key=(y.shape,y.tobytes())
  if key in self.cache:return self.cache[key]
  J,s2,V,th=geom(self.A,self.S,y);f=np.array([v(y) for v in self.fields]);fy=np.array([v(y,1) for v in self.fields]);flu={nm:recover(nm,v,y,J,s2) for nm,v in self.fluids.items()};g=beta_eff(f[4])*flu['dark_matter'][5]
  parts={}
  for nm,(D,E,Er,Sr,P,rho,v,lor) in flu.items():
   m=self.fluids[nm][2](y)/J;an=s2*m*m/(Er+P);parts[nm]=(Er/KAP,an/KAP,m/KAP,(3*P+y*an)/KAP)
  c=np.array([1,1,1/KAP])[:,None];ak=np.sum(c*f[1::2]**2,axis=0);b=4*s2*np.sum(c*fy[::2]**2,axis=0);u=potential(*f[::2]);parts['scalar']=(.5*(ak+y*b)+u,b,-2*np.sum(c*f[1::2]*fy[::2],axis=0),1.5*ak-.5*y*b-3*u)
  ans=(J,s2,V,th,f,fy,flu,g,parts);self.cache[key]=ans;return ans
 def parts(self,y):return self.at(y)[-1]
 def scalar_rate(self,al,y):
  J,s2,V,th,f,fy,flu,g,_=self.at(y);a=al(y);ay=al(y,1);up=sh.up(f[::2])*np.array([1,1,KAP])[:,None];up[2]-=g
  out=[]
  for i in range(3):
   lap=s2*(4*y*self.fields[2*i](y,2)+(6+2*y*(th+2*V))*fy[2*i]);out.extend([a*f[2*i+1],a*(lap+self.S['K']*f[2*i+1]-up[i])+4*y*s2*ay*fy[2*i]])
  return np.array(out)

def evaluate(A,S,raw,sc,mg,g_boundary,scalar_boundary):
 from .scalar import StructuralRegistry, weak_port
 n=A.n;reg=StructuralRegistry(A,S,raw,sc,mg)
 def digest():return hashlib.sha256(b''.join(np.asarray(x).tobytes() for x in list(S.values())+list(sc)+[a for v in raw.values() for a in v])).hexdigest()
 pre=digest();args=tuple(S[z] for z in ('k','ell','L0','p','K','incoming'));gauge=A.cmc(*args,reg);al,Kt,cd=gauge;G,gd=A.rhs(S,reg,g_boundary,prepared=gauge)
 # Same state and lapse, paired baseline (same split B for this comparison).
 G0,gd0=Assembly.rhs(A,S,reg,g_boundary,prepared=gauge)
 y=mg.q.ravel();r=np.sqrt(y);J,s2,V,th,f,fy,flu,g,_=reg.at(y);W=y*J;a=al(y);ay=al(y,1);ar=2*r*ay;sy=-2*V+th;ny=3*V-th/2;p=np.array([A.eval(S['p'],t) for t in y]);kr=S['K']/3+y*p;kt=S['K']/3-y*p/2
 avg=lambda z:mg.average(np.asarray(z).reshape(mg.q.shape));integ=lambda z:(A.h*avg(z)).tolist()
 yf=A.ends;Jf,s2f,_,_,ff,fyf,fluf,gf,_=reg.at(yf);af=al(yf)
 lv=-G['L0']-y*np.array([A.eval(G['ell'],t) for t in y])+y*np.array([A.eval(G['k'],t) for t in y])/2;nv=3*G['L0']+3*y*np.array([A.eval(G['ell'],t) for t in y])-y*np.array([A.eval(G['k'],t) for t in y])/2
 base,base_scalar,_,_=matter_rhs(A,S,reg,al,mg);mr={};rows={};exchange_dm=None
 for nm,(D,E,Er,Sr,P,rho,v,lor) in flu.items():
  d,e,m=[z(y) for z in reg.fluids[nm]];dy,ey,my=[z(y,1) for z in reg.fluids[nm]];u=m/J;uy=my/J-u*ny;Q=s2*y*u*u;Qy=s2*((1+y*sy)*u*u+2*y*u*uy)
  if nm=='radiation':Ery=ey/J-E*ny;Py=((E-P)*Ery-Qy/2)/(E+3*P)
  else:Dy=dy/J-D*ny;Ery=(D*Dy+Qy/2)/Er;Py=np.zeros_like(y)
  an=s2*u*u/(Er+P);any_=s2*(sy*u*u+2*u*uy)/(Er+P)-an*(Ery+Py)/(Er+P);pr=P+y*an
  Df,Ef,Erf,Srf,Pf,rhof,vf,lorf=fluf[nm];F=af*yf*Jf*np.array([Df*vf,s2f*Srf,Srf*vf+Pf]);ex=np.zeros((3,len(y)))
  if nm=='dark_matter':
   ex=np.array([-a*W*g*lor*(f[5]+v*2*r*fy[4]),-a*W*g*f[5],a*W*g*2*r*fy[4]]);exchange_dm=ex.copy()
  geoE=W*(a*(kr*pr+2*kt*P)-s2*Sr*ar);regM=-2*r*W*((Er+pr)*ay+a*(Py+2*an+y*any_+2*y*V*an));
  # Explicit one-sided physical traces and independent shared face flux.
  def oneface(i,t):
   k=A.eval(S['k'],t,i=i);ell=A.eval(S['ell'],t,i=i);L=S['L0']+t*ell;JJ=np.exp(3*L-t*k/2);ss=np.exp(-2*L+t*k);yy=np.array([t]);_,_,_,sr,pr0,_,vel0,_=recover(nm,reg.fluids[nm],yy,np.array([JJ]),np.array([ss]));return float(al(t)*t*JJ*(sr[0]*vel0[0]+pr0[0]))
  fl=np.array([oneface(i,yf[i]) for i in range(n)]);fright=np.array([oneface(i,yf[i+1]) for i in range(n)])
  ports=(F[2,:-1]-fl)-(F[2,1:]-fright)
  rates=np.array([-np.diff(F[0])/A.h+avg(ex[0]),-np.diff(F[1])/A.h+avg(geoE+ex[1]),avg(regM+ex[2])+ports/A.h]);mr[nm]=rates
  geoM=W*(-Er*ar+a*(pr*r*(2*V-th)+2*P*(1/r+2*r*V)))
  # Record quadrature IBP mismatch and physical momentum flux telescoping separately.
  ibp_m=A.h*avg(regM-geoM)+(fright-fl)
  tel=A.h*np.sum(rates-np.array([avg(ex[0]),avg(geoE+ex[1]),avg(geoM+ex[2])]),axis=1)+F[:,-1]-F[:,0]
  row=dict(rate_difference=mx(rates-np.array(base[nm])),center_flux=F[:,0].tolist(),outer_flux=F[:,-1].tolist(),momentum_flux_port=ports.tolist(),momentum_ibp_remainder=ibp_m.tolist(),flux_telescoping=tel.tolist(),telescoping_ledger_closure=abs(tel[2]-sum(ibp_m)),momentum_shell_difference=(A.h*(rates[2]-base[nm][2])).tolist())
  if nm!='radiation':
   er=J*Er;ery=J*(Ery+ny*Er);b=s2*m/er;by=s2*(sy*m+my)/er-b*ery/er;vel=a*b;vely=ay*b+a*by
   ds=-3*d*vel-2*y*(dy*vel+d*vely)+ex[0]/y;es=-3*er*vel-2*y*(ery*vel+er*vely)+(geoE+ex[1])/y;ms=(regM+ex[2])/r**3
   dh,eh,mh=[mg.reconstruction(z,j)(y) for z,j in zip(rates,(2,2,3))];w1=d/er;w2=s2*y*m/er;w3=J*pr;defect=eh-w1*dh-w2*mh-w3*lv;metric=w3*(a*kr-lv);rep=eh-es-w1*(dh-ds)-w2*(mh-ms)
   row.update(dust_ledger=mx(defect-metric-rep),dust_metric_shell=integ(y*metric),dust_rep_shell=integ(y*rep),dust_defect_shell=integ(y*defect))
  rows[nm]=row
 # Shared scalar rate, same cache and raw exchange. Compare original and reconstruct actual sample transfer.
 strong=reg.scalar_rate(al,y);sr,weakdiag=weak_port(A,S,raw,sc,mg,reg,al,scalar_boundary);rec=[reg.space.field(v) for v in sr];hr=np.array([z(y) for z in rec]);fr=2*r*fy[::2];P=f[1::2];c=np.array([1,1,1/KAP])[:,None];frt=ar*P+a*2*r*fy[1::2];hrfr=2*r*np.array([rec[i](y,1) for i in (0,2,4)]);ef=hr[::2]-strong[::2];ep=hr[1::2]-strong[1::2];er=hrfr-frt;up=sh.up(f[::2]);rho,an,jhat,tr=reg.parts(y)['scalar'];B=y*an;j=r*jhat
 DE=W*(np.sum(c*(P*ep+s2*fr*er),axis=0)+np.sum(up*ef,axis=0));DM=-W*np.sum(c*(fr*ep+P*er),axis=0)
 # Direct chain rule on represented and strong rates (difference ledger).
 def st(ft,pt,frate):return W*(np.sum(c*(P*pt+s2*fr*frate),axis=0)+np.sum(up*ft,axis=0)+lv*B+nv*rho),-W*np.sum(c*(fr*pt+P*frate),axis=0)+W*nv*j
 hE,hM=st(hr[::2],hr[1::2],hrfr);sE,sM=st(strong[::2],strong[1::2],frt)
 # Independent flux/source differentiation of scalar stress on smooth cell quadrature.
 pr=(tr+2*y*an)/3;pt=(tr-y*an)/3;fr_r=2*fy[::2]+4*y*np.array([reg.fields[i](y,2) for i in (0,2,4)]);Pr=2*r*fy[1::2];jr=-np.sum(c*(Pr*fr+P*fr_r),axis=0);Br=2*r*sy*B+2*s2*np.sum(c*fr*fr_r,axis=0);prr=np.sum(c*P*Pr,axis=0)+Br/2-np.sum(up*fr,axis=0);Wr=2*r*J+W*2*r*ny
 FE=(ar*W+a*Wr)*s2*j+a*W*s2*(2*r*sy*j+jr);FM=(ar*W+a*Wr)*pr+a*W*prr;GE=W*(a*(kr*pr+2*kt*pt)-s2*j*ar);GM=W*(-rho*ar+a*(pr*r*(2*V-th)+2*pt*(1/r+2*r*V)));XE=a*W*g*f[5]/KAP;XM=-a*W*g*fr[2]/KAP;ME=W*((nv+a*S['K'])*rho+(lv-a*kr)*B);MM=W*(nv+a*S['K'])*j
 # CMC constant test, physical flux jump, quadrature IBP separately on dedicated nodes.
 cy=A.cy;CJ,cs2,cV,cth,*_=reg.at(cy);Lmat,rhs,rate,x=A.elliptic_system(*args,reg);sol=np.r_[al.coeff-A.lapse.one,al.flux];mu=x['mu'];div=cs2*(4*cy*al(cy,2)+(6+2*cy*(cth+2*cV))*al(cy,1));T=rate-x['Q']*al(cy)+div
 Fphys=cd['physical_flux'];Fnum=al.flux;jumps=[]
 for i,t in enumerate(A.ends[1:-1],1):
  def traceflux(cell):
   k=A.eval(S['k'],t,i=cell);ell=A.eval(S['ell'],t,i=cell);L=S['L0']+t*ell;return 2*t**1.5*np.exp(L+t*k/2)*al(t,1)
  jumps.append(traceflux(i)-traceflux(i-1))
 ibp=float(mu@div-Fphys+sum(jumps));cw=float(A.lapse.one@(Lmat@sol-rhs)[:-1]);cmc=float(mu@T);cmc_cl=cmc-(Fphys-Fnum-sum(jumps)+ibp-cw)
 assert digest()==pre
 return dict(solved_lapse=al.packed().tolist(),K_rate=float(Kt),n=n,immutable=True,digest=pre,weak=weakdiag,matter_rates={k:v.tolist() for k,v in mr.items()},geometry_rates={k:np.asarray(G[k]).tolist() for k in ("p","k","ell","L0","K","J","RV","kh","eh","w")},geometry={k:gd[k] for k in ('principal_error','total_error','joint_load_difference','joint_work','joint_adjusted_energy_closure','zeta_origin_num','physical_flux_difference','H_m_norm','momentum_norm')},geometry_p_rate_diff=mx(G['p']-G0['p']),matter=rows,scalar=dict(weak_minus_sample_rate=mx(sr-np.array(base_scalar)),strong_rate_crosswalk=mx(strong-np.array(scalar_rhs(A,S,reg,al,y))),energy_transfer_closure=mx(hE-sE-DE),momentum_transfer_closure=mx(hM-sM-DM),energy_full_closure=mx(hE+FE-GE-XE-ME-DE),momentum_full_closure=mx(hM+FM-GM-XM-MM-DM),energy_metric_shell=integ(ME),momentum_metric_shell=integ(MM),energy_rep_shell=integ(DE),momentum_rep_shell=integ(DM)),exchange=dict(energy=mx(XE+exchange_dm[1]/KAP),momentum=mx(XM+exchange_dm[2]/KAP)),cmc=dict(strong_integral=cmc,physical_flux=Fphys,numerical_flux=Fnum,jumps=jumps,ibp=ibp,weak_constant=cw,closure=cmc_cl,alpha_min=cd['alpha_min']))
