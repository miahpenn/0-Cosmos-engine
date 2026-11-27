import sys,json, pathlib, numpy as np, scipy.linalg as la
from fractions import Fraction
from . import geometry_core as gc, native_cmc as nc
from numpy.polynomial.legendre import legmul
basis=json.load(open(pathlib.Path(__file__).with_name('basis.json'))); conv=lambda a:np.array([[float(Fraction(v)) for v in row] for row in a]); raw0=conv(basis['raw']); q0=conv(basis['quotient']);
class Candidate(nc.NativeAssembly):
 def __init__(self,n,exact=False):
  super().__init__(n)
  self.M[:]=0;self.N[:]=0;self.D[:]=0;self.Q=np.zeros((self.np,self.nh));self.W=np.zeros((self.np,self.np));self.Strong=np.zeros((self.nh,self.np)); mm=json.load(open(pathlib.Path(__file__).with_name("moments.json"))) if exact else None
  for i,c in enumerate(self.cells):
   pi,ix,sl=c['pi'],c['ix'],c['sl']; dy=c['hi']-c['lo']; t=2*(c['y']-c['lo'])/dy-1; wp,wh=c['wp'],c['wh']; raw=raw0 if i==0 else np.eye(3); m=raw.shape[1]; Q=q0 if i==0 else np.array([[-.25,0],[0,-.25]]) # integral P1/2=(t²-1)/4; P2/2=(t³-t)/4
   F0=gc.val(raw,t); L=la.cholesky(F0.T@(wp[:,None]*F0),lower=True); T=la.solve_triangular(L.T,np.eye(m)); cp=raw@T; F=gc.val(cp,t)
   B0=(1-t*t)[:,None]*gc.val(Q,t); Lb=la.cholesky(B0.T@(wh[:,None]*B0),lower=True); U=la.solve_triangular(Lb.T,np.eye(m-1)); qb=Q@U; H=np.column_stack([(1-t)/2,(1+t)/2,(1-t*t)[:,None]*gc.val(qb,t)])
   chb=np.column_stack([np.pad(legmul(qb[:,j],[2/3,0,-2/3]),(0,len(qb)+2-len(legmul(qb[:,j],[2/3,0,-2/3])))) for j in range(m-1)]); ch=np.zeros((len(chb),m+1));ch[:2,:2]=[[.5,.5],[-.5,.5]];ch[:,2:]=chb
   # derivative coordinates from exact primitive definition, no projected derivative solve
   A=np.zeros((m,m+1));A[0,0]=-1/dy; A[0,1]=1/dy; R=np.vstack([-raw[0,1:],np.eye(m-1)]); A[:,2:]=R@U/dy
   D=la.solve_triangular(T,A); D[:,0]=-D[:,1]; Hy=F@D
   M=F.T@(wp[:,None]*F);N=H.T@(wh[:,None]*H); Qm=F.T@(wh[:,None]*H); Wm=F.T@(wh[:,None]*F)
   if exact:
    j=mm[f"{n}/{i}"]; V=la.block_diag(np.eye(2),U); M=T.T@np.array(j["M"])@T; N=V.T@np.array(j["N"])@V; Qm=T.T@np.array(j["Q"])@V; Wm=T.T@np.array(j["W"])@T; self.Strong[np.ix_(ix,pi)]+=V.T@np.array(j["S"])@T
   self.Q[np.ix_(pi,ix)]=Qm;self.W[np.ix_(pi,pi)]=Wm;self.M[np.ix_(pi,pi)]=M;self.N[np.ix_(ix,ix)]+=N;self.D[np.ix_(pi,ix)]=D
   self.P[sl,pi]=F;self.H[sl,ix]=H;self.Py[sl,pi]=gc.val(cp,t,1,dy);self.Hy[sl,ix]=Hy;c.update(cp=cp,ch=ch,qb=qb,F=F,H=H,Hy=Hy,D=D)
  self.pface=np.array([self.trace(v) for v in np.eye(self.np)])
 def solveH(self,b):return la.solve(self.N,b,assume_a='pos')
 def eval(self,x,y,kind='p',der=0,i=None):
  i=min(np.searchsorted(self.ends,y,side='right')-1,self.n-1) if i is None else i;c=self.cells[i];dy=c['hi']-c['lo'];t=2*(y-c['lo'])/dy-1
  if kind=='H' and der==0:
   a=x[c['ix']];return float((1-t)/2*a[0]+(1+t)/2*a[1]+(1-t*t)*(gc.val(c['qb'],[t])[0]@a[2:]))
  if kind=='H' and der>=1:
   a=x[c['ix']]; v=c['D'][:,1]*(a[1]-a[0])+c['D'][:,2:]@a[2:];return float(gc.val(c['cp'],[t],der-1,dy)[0]@v)
  return super().eval(x,y,kind,der,i)
