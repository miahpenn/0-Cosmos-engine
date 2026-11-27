import sys,pathlib,json,numpy as np,scipy.linalg as la
from numpy.polynomial.legendre import legder,legmul,legval
# argv1 source, argv2 structural prototype directory
from .candidate import Candidate

from fractions import Fraction as F
from math import comb

def gram(n,i,d,a):
 l,r=F(i,n),F(i+1,n);dy=r*r-l*l;z=-(r*r+l*l)/dy;u=2/dy;w=F(2) if a==6 else F(1,2)
 m=[w*sum(F(comb(k,j))*z**(k-j)*u**j*(r**(a+2*j+1)-l**(a+2*j+1))/F(a+2*j+1) for j in range(k+1)) for k in range(2*d+1)]
 p=[[F(1)],[F(0),F(1)]]
 for j in range(1,d):
  q=[F(0)]+[v*F(2*j+1,j+1) for v in p[j]]
  for k,v in enumerate(p[j-1]):q[k]-=v*F(j,j+1)
  p.append(q)
 return np.array([[sum(v*t*m[j+k] for j,v in enumerate(b) if v for k,t in enumerate(c) if t) for c in p] for b in p],float)

def setup(A):
 for i,c in enumerate(A.cells):
  d=len(c['ch'])-1;gp=gram(A.n,i,d,6);gh=gram(A.n,i,d,4);cp=c['cp'];ch=c['ch'];cp=np.pad(cp,((0,d+1-len(cp)),(0,0)))
  c.update(Lh=ch.T@gh,Lp=cp.T@gp,Lq=cp.T@gh)
 return A

def pad(x,n):return np.pad(x,(0,n-len(x)))
def pcoef(c,f):return c['cp']@f[c['pi']]
def hcoef(c,f):
 a=f[c['ix']]; q=c['ch'][:,2:]@a[2:];q[:2]+=[(a[0]+a[1])/2,(a[1]-a[0])/2];return q

def strong(c,f):
 v=pcoef(c,f);dy=c['hi']-c['lo'];mid=(c['hi']+c['lo'])/2
 return 10*v+pad(legmul([8*mid/dy,4],legder(v)),len(v))

def split(A,f,hat):
 s=[strong(c,f) for c in A.cells]; pc=[pcoef(c,f) for c in A.cells]
 # Interpolation defined on every input: local low quadratic part, average its traces at shared faces.
 left=np.array([legval(-1,v[:3]) for v in s]);right=np.array([legval(1,v[:3]) for v in s]);faces=np.r_[left[0],right[:-1]+(left[1:]-right[:-1])/2,right[-1]]
 h=np.zeros(A.nh);h[:A.n+1]=faces; il=[];res=[]; load=np.zeros(A.nh); ports=np.zeros(A.nh)
 for i,c in enumerate(A.cells):
  a=s[i][2]; h[c['ix'][2]]=a/c['ch'][2,2]
  # use face-difference form for linear coeff, and known degree-2 coeff; remainder stays polynomial
  v=np.array([(faces[i]+faces[i+1])/2-a,(faces[i+1]-faces[i])/2,a]);il.append(v)
  r=s[i]-pad(v,len(s[i]));res.append(r)
  t=2*(c['y']-c['lo'])/(c['hi']-c['lo'])-1
  load[c['ix']]+=c['Lh']@pad(r,c['Lh'].shape[1])
 for i in range(1,A.n):ports[i]=A.ends[i]**2.5*(legval(-1,pc[i])-legval(1,pc[i-1]))
 ports[A.n]=hat-legval(1,pc[-1]); d=A.solveH(load+ports)
 return h,d,il,res,load,ports

def calc(A,l,k,j,rv,kh,eh):
 M,N,D=A.M,A.N,A.D;z=4*l-2*k/3;zh=4*eh-2*kh/3;a=-2*l/3+4*k/9;ah=-2*eh/3+4*kh/9
 ih,dh,il,_,_,_=split(A,a,ah);iz,dz,zil,zres,_,_=split(A,z,zh)
 hc=ih+dh+4*j/3-8*rv/3
 # Combine polynomial coefficients before differentiation; retain solved remainder as its own coordinate.
 load=np.zeros(A.np)
 for i,c in enumerate(A.cells):
  n=len(c['ch']); dy=c['hi']-c['lo']; t=2*(c['y']-c['lo'])/dy-1
  v=pad(il[i],n)+4*hcoef(c,j)/3-8*hcoef(c,rv)/3+pad(pcoef(c,z),n)
  deriv=legval(t,legder(v)*2/dy)
  load[c['pi']]=c['Lp']@pad(legder(v)*2/dy,n)-c['Lq']@pad(zres[i],n)
 pt=D@dh+la.solve(M,load+A.Q@dz,assume_a='pos')
 wt=-3*(A.e@hc+A.trace(z))
 # original paired moment formula (same operands and matrices)
 hc0=-2*A.B(l,eh)/3+4*A.B(k,kh)/9+4*j/3-8*rv/3; bz0=A.B(z,zh)
 pt0=D@hc0+la.solve(M,A.Q@bz0-10*A.W@z,assume_a='pos')
 wt0=-3*(A.e@hc0+A.trace(z))
 norm=lambda p:float(np.sqrt(p@M@p))
 bh=iz+dz; paired=-D.T@(M@z)+A.e*zh
 return dict(p_center=A.eval(pt,0),p_norm=norm(pt),w_rate=float(wt),paired_p_center=A.eval(pt0,0),paired_p_norm=norm(pt0),paired_w_rate=float(wt0),rate_diff=norm(pt-pt0),load_relative=float(la.norm(N@bh-paired)/max(1e-300,la.norm(paired))),origin_num=float(bh[0]-10*A.eval(z,0)),z_res_max=float(max(la.norm(v) for v in zres)),z_lift_norm=float(la.norm(dz)))
