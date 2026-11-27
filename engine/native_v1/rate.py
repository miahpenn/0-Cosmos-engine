"""Immutable rate interface, no temporal stages."""
import hashlib
import numpy as np
from . import composed as C, scalar as I
def snapshot(S,raw,sc):
 h=hashlib.sha256()
 for key,val in sorted(S.items()):
  a=np.asarray(val);h.update(key.encode());h.update(str(a.shape).encode());h.update(a.tobytes())
 for key,vs in sorted(raw.items()):
  h.update(key.encode())
  for v in vs:h.update(np.asarray(v).tobytes())
 for v in sc:h.update(np.asarray(v).tobytes())
 return h.hexdigest()
def bR(A,S):
 y=A.Y;L=S['L0']+y*(A.pface@S['ell']);k=A.pface@S['k'];return 2*y**1.5*np.exp(L+y*k/2)
def rate(A,S,raw,sc,z,g_boundary,scalar_boundary=None):
 pre=snapshot(S,raw,sc);z=np.asarray(z,dtype=float)
 if z.shape!=(A.lapse.dim+2,) or z[-1]!=1 or not np.all(np.isfinite(z)):raise ValueError('packed augmented gauge required')
 if not np.isfinite(g_boundary):raise ValueError('finite geometry boundary required')
 if scalar_boundary is not None and (np.shape(scalar_boundary)!=(3,) or not np.all(np.isfinite(scalar_boundary))):raise ValueError('three finite scalar wave fluxes required')
 st={k:np.array(v,copy=True) if np.ndim(v) else float(v) for k,v in S.items()};st['incoming']=np.r_[z[:-2],bR(A,st)*z[-2]]
 r={k:tuple(np.array(v,copy=True) for v in vs) for k,vs in raw.items()};f=[np.array(v,copy=True) for v in sc]
 d=C.evaluate(A,st,r,f,C.MomentGrid(A),float(g_boundary),scalar_boundary)
 l=np.array(d['solved_lapse']);d['outgoing_z']=np.r_[l[:-1],l[-1]/bR(A,st),1.].tolist();assert snapshot(S,raw,sc)==pre;d['external_digest']=pre;return d
def scalar_matrices(A,reg,al):
 V=reg.space;edges=np.sqrt(np.unique(np.r_[A.ends,A.cent]));t,w=I.leggauss(24);r=((edges[1:]+edges[:-1])[:,None]/2+np.diff(edges)[:,None]*t/2).ravel();wr=(np.diff(edges)[:,None]*w/2).ravel();y=r*r;J,s2,*_=reg.at(y);B=V.B(y);H=2*r[:,None]*V.H(y,1)@V.E;M=B.T@((wr*y*J)[:,None]*B);K=H.T@((wr*y*J*al(y)*s2)[:,None]*H);e=V.B(A.Y);d=2*np.sqrt(A.Y)*V.H(A.Y,1)@V.E;JR,sR,*_=C.geom(A,reg.S,np.array([A.Y]));factor=al(A.Y)*A.Y*JR[0]*sR[0];T=np.linalg.solve(np.linalg.cholesky(M),np.eye(A.n));port=T@np.outer(e,factor*d)@T.T;stiff=T@K@T.T
 return M,e,dict(boundary_mass_norm=float(np.linalg.norm(port,2)),stiffness_mass_norm=float(np.linalg.norm(stiff,2)),rank=int(np.linalg.matrix_rank(port)))
