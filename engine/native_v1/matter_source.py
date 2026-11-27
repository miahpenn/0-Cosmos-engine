"""Point-sample regular source adapter, not finite-volume transport.
Pinned physical source laws; raw densitized arrays remain authoritative.
"""
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.linalg import solve
from .assembly import KAP,val
V0=8.242522415500654e-5;BETA_D=-.04;PHI_TODAY=.179055;LAM=KAP/(2*.10*V0)*np.exp(-4*BETA_D*PHI_TODAY)
def potential(S,D,phi):return .5*S*S-.5*np.exp(-2*BETA_D*phi)*D*D+LAM/8*(S**4+6*S*S*D*D+D**4)+V0*np.exp(-phi)/KAP

def reconstruct(y,v):
    v=np.asarray(v);sl=[]
    for i,x in enumerate(y):
        j=max(0,min(len(y)-4,i-1));ids=np.arange(j,j+4);t=y[ids]-x;sc=max(abs(t));V=np.array([(t/sc)**m for m in range(4)]);sl.append(solve(V,[0.,1/sc,0.,0.])@(v[ids]-v[i]))
    return CubicHermiteSpline(y,v,sl,extrapolate=True)

def metric(A,S,y):
    y=np.asarray(y); k=np.array([A.eval(S['k'],x) for x in y]);ell=np.array([A.eval(S['ell'],x) for x in y]);L=S['L0']+y*ell;return np.exp(3*L-y*k/2),np.exp(-2*L+y*k)

def fluid(kind,d,e,m,y,J,s2):
    D=d/J; E=e/J; v=m/J; S2=s2*y*v*v
    if kind=='dust':
        if np.any(D<=0):raise ValueError('dust rest not positive')
        Er=np.sqrt(D*D+S2);P=np.zeros_like(E);defect=E-Er
    elif kind=='radiation':
        if np.any(E<0) or np.any((S2>=E*E)&(S2>0)):raise ValueError('radiation outside strict cone')
        if np.any(d!=0):raise ValueError('radiation rest nonzero')
        Er=E;P=np.zeros_like(E);active=E>0;P[active]=(E[active]**2-S2[active])/(np.sqrt(4*E[active]**2-3*S2[active])+E[active]);defect=np.zeros_like(E)
    else:raise ValueError('unknown species')
    den=Er+P;active=den>0;anis=np.zeros_like(E);anis[active]=s2[active]*v[active]**2/den[active];u=np.zeros_like(E);u[active]=s2[active]*v[active]/den[active];v2=y*u*u/s2
    if np.any(v2>=1):raise ValueError('noncausal recovery')
    return (Er/KAP,anis/KAP,v/KAP,(3*P+y*anis)/KAP),defect,v2

class Sources:
    def __init__(self,A,S,raw,scalar):
        self.A=A;self.S=S;self.raw=raw;self.scalar=scalar;y=A.cent
        self.fluids={name:[reconstruct(y,vs/y**power) for vs,power in zip(values,(1,1,1.5))] for name,values in raw.items()};self.fields=[reconstruct(y,v) for v in scalar]
    def __call__(self,y):return self.evaluate(y)[0]
    def evaluate(self,y):
        y=np.asarray(y);J,s2=metric(self.A,self.S,y);parts=[];defects={};vel={}
        for name,f in self.fluids.items():
            p,d,v=fluid('radiation' if name=='radiation' else 'dust',*[x(y) for x in f],y,J,s2);parts.append(p);defects[name]=d;vel[name]=v
        S,PS,D,PD,phi,Pi=[f(y) for f in self.fields];grad=np.array([self.fields[i](y,1) for i in (0,2,4)]);mom=np.array([PS,PD,Pi]);w=np.array([1,1,1/KAP])[:,None];kin=np.sum(w*mom*mom,axis=0);gg=np.sum(w*grad*grad,axis=0);V=potential(S,D,phi);an=4*s2*gg
        parts.append((.5*kin+.5*y*an+V,an,-2*np.sum(w*mom*grad,axis=0),1.5*kin-.5*y*an-3*V))
        return tuple(sum(p[i] for p in parts) for i in range(4)),defects,vel

def fixture(A,S):
    y=A.cent;J,s2=metric(A,S,y);raw={}
    # Independent manufactured raw physical sample arrays, no new model parameters.
    for name,d,m in [('dark_matter',.003+.001*y,.0002+.0001*y),('baryons',.001+.0002*y,-.0001+.00002*y)]:
        e=np.sqrt(d*d+s2*y*m*m);raw[name]=(y*d,y*e,y**1.5*m)
    e=.001+.0001*y;m=.00005+.00001*y;raw['radiation']=(np.zeros_like(y),y*e,y**1.5*m)
    scalar=[.0001+.00003*y,.00002+.00001*y,.0008+.00004*y,-.00003+.00001*y,.1+.002*y,.0003+.0001*y]
    return raw,scalar
