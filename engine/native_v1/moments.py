"""Research-only weighted moment reconstruction, separate from legacy point storage."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import solve
from scipy.interpolate import PPoly
from .native_cmc import NativeAssembly
from .assembly import state,KAP
from .matter_source import Sources,fixture,metric,reconstruct
from .transport import geom,recover,beta_eff,scalar_rhs
import json

class MomentGrid:
    def __init__(self,A,order=24):
        self.A=A; self.n=A.n; self.y=A.ends;self.h=A.h
        t,w=leggauss(order); r=np.sqrt(self.y);self.r=(r[:-1,None]+r[1:,None])/2+self.h*t/2;self.q=self.r**2;self.w=np.broadcast_to(w/2,self.q.shape)
        self.weights={j:np.diff(r**(j+1))/((j+1)*self.h) for j in (2,3)}
        self.map={j:self.build(j) for j in (2,3)}
    def average(self,values):return np.sum(self.w*values,axis=-1)
    def build(self,j):
        n=self.n;y=self.y;W=self.weights[j];q=self.q;rw=self.w*self.r**j/W[:,None];face=np.zeros((n+1,n));slope=face.copy()
        for i in range(n+1):
            start=max(0,min(n-4,i-2));ids=np.arange(start,start+4);sc=max(abs(y[ids[0]]-y[i]),abs(y[ids[-1]+1]-y[i]));v=(q[ids]-y[i])/sc
            H=np.stack([np.sum(rw[ids]*v**k,axis=1) for k in range(4)],axis=1)
            inv=solve(H,np.eye(4));face[i,ids]=inv[0];slope[i,ids]=inv[1]/sc
        bubble=np.array([0.,0.,1.,-2.,1.]); C=np.zeros((n,5,n))
        for i,dy in enumerate(np.diff(y)):
            a,b=face[i],face[i+1];c,d=dy*slope[i],dy*slope[i+1]
            C[i,:4]=np.stack([a,c,3*(b-a)-2*c-d,2*(a-b)+c+d]); t=(q[i]-y[i])/dy
            av=np.array([np.sum(rw[i]*t**k) for k in range(5)]);res=np.eye(n)[i]-av@C[i];C[i]+=bubble[:,None]*res/(av@bubble)
        return C
    def reconstruction(self,raw,j):
        coeff=self.map[j]@(raw/self.weights[j]);dy=np.diff(self.y)
        return PPoly((coeff/dy[:,None]**np.arange(5)).T[::-1],self.y)
    def moment(self,f,j):return self.average(self.r**j*f(self.q))

class MomentSources(Sources):
    def __init__(self,A,S,raw,scalar,mg):
        self.A=A;self.S=S;self.raw=raw;self.scalar=scalar;self.mg=mg
        self.fluids={name:[mg.reconstruction(v,j) for v,j in zip(vs,(2,2,3))] for name,vs in raw.items()}
        self.fields=[reconstruct(A.cent,v) for v in scalar]

def initial(A,S,mg):
    y=mg.q.ravel();J,s2=metric(A,S,y);raw={}
    for name,d,m in [('dark_matter',.003+.001*y,.0002+.0001*y),('baryons',.001+.0002*y,-.0001+.00002*y)]:
        e=np.sqrt(d*d+s2*y*m*m);raw[name]=tuple(mg.average((v*y**(j/2)).reshape(mg.q.shape)) for v,j in zip((d,e,m),(2,2,3)))
    e=.001+.0001*y;m=.00005+.00001*y;raw['radiation']=(np.zeros(A.n),mg.average((y*e).reshape(mg.q.shape)),mg.average((y**1.5*m).reshape(mg.q.shape)))
    return raw,fixture(A,S)[1]

def matter_rhs(A,S,src,al,mg):
    y=mg.q.ravel();r=np.sqrt(y);J,s2,V,theta=geom(A,S,y);a=al(y);ar=2*r*al(y,1);p=np.array([A.eval(S['p'],x) for x in y]);Aa=y*p;phi=src.fields[4](y);Pi=src.fields[5](y);pr=2*r*src.fields[4](y,1);be=beta_eff(phi)
    yf=A.ends;Jf,s2f,_,_=geom(A,S,yf);af=al(yf);out={};led={};ex=[]
    for name,f in src.fluids.items():
        D,E,Er,Sr,P,rho,v,W=recover(name,f,y,J,s2);Df,Ef,Erf,Srf,Pf,rhof,vf,Wf=recover(name,f,yf,Jf,s2f)
        flux=af*yf*Jf*np.array([Df*vf,s2f*Srf,Srf*vf+Pf]);root=y*J;stress_r=Sr*v+P
        geoE=root*(a*((S['K']/3+Aa)*stress_r+2*(S['K']/3-Aa/2)*P)-s2*Sr*ar)
        geoM=root*(-Er*ar+a*(stress_r*r*(2*V-theta)+2*P*(1/r+2*r*V)))
        exchange=np.zeros((3,len(y)))
        if name=='dark_matter':
            exchange=np.array([-a*root*be*rho*W*(Pi+v*pr),-a*root*be*rho*Pi,a*root*be*rho*pr]);ex=exchange
        source=np.array([np.zeros_like(y),geoE,geoM])+exchange;av=mg.average(source.reshape(3,A.n,-1));du=-np.diff(flux,axis=1)/A.h+av;out[name]=tuple(du)
        led[name]=dict(telescoping=(A.h*np.sum(du-av,axis=1)+flux[:,-1]-flux[:,0]).tolist(),center_flux=flux[:,0].tolist())
    rho=recover('dark_matter',src.fluids['dark_matter'],y,J,s2)[5];q=be*rho;sc=np.array([a*y*J*q*Pi,-a*y*J*q*pr]);net=A.h*np.sum(mg.average(((sc+ex[1:])/KAP).reshape(2,A.n,-1)),axis=1)
    return out,scalar_rhs(A,S,src,al,A.cent),led,net

def cons(A,S,src,y):
    J,s2=metric(A,S,y);rho,ani,jhat,tr=src(y);jr=np.sqrt(y)*jhat;pt=(tr-y*ani)/3;pr=pt+y*ani
    return y*J*np.array([rho,jr]),(J,s2,rho,jr,pr,pt)

def flux(A,S,src,al,r):
    y=r*r;_,(J,s2,rho,jr,pr,pt)=cons(A,S,src,y);return al(y)*y*J*np.array([s2*jr,pr])

def geo_source(A,S,src,al,y):
    _,(J,s2,rho,jr,pr,pt)=cons(A,S,src,y);r=np.sqrt(y);_,_,V,theta=geom(A,S,y);aa=y*np.array([A.eval(S['p'],x) for x in y]);ar=2*r*al(y,1);a=al(y);root=y*J
    return root*np.array([a*((S['K']/3+aa)*pr+2*(S['K']/3-aa/2)*pt)-s2*jr*ar,-rho*ar+a*(pr*r*(2*V-theta)+2*pt*(1/r+2*r*V))])

def compare(A,S,mg,raw,scalar):
    src=MomentSources(A,S,raw,scalar,mg);gauge=A.cmc(*(S[k] for k in ('k','ell','L0','p','K','incoming')),src);G,d=A.rhs(S,src,.13,prepared=gauge);al=G['alpha'];mr,sr,ledger,ex=matter_rhs(A,S,src,al,mg)
    y=A.cent;r=A.r;q=mg.q.ravel();J,s2=metric(A,S,q)
    dust={}
    for name in ('dark_matter','baryons'):
        f=src.fluids[name];D,E,M=[z(q) for z in f];dt,et,mt=[mg.reconstruction(z,j)(q) for z,j in zip(mr[name],(2,2,3))];logst=-G['L0']-q*np.array([A.eval(G['ell'],x) for x in q])+q/2*np.array([A.eval(G['k'],x) for x in q]);rec=np.sqrt(D**2+s2*q*M**2);rt=(D*dt+s2*q*M*mt+s2*q*logst*M*M)/rec;dust[name]={'max_stored_minus_recovered_moment':float(np.max(abs(mg.average((q*(E-rec)).reshape(mg.q.shape))))),'max_tangent_moment':float(np.max(abs(mg.average((q*(et-rt)).reshape(mg.q.shape)))))}
    def direction(e,yy):
        Q=dict(S)
        for k in ('k','ell','L0','p','K','J','RV','kh','eh','w'):Q[k]=S[k]+e*G[k]
        U={k:tuple(v+e*dv for v,dv in zip(raw[k],mr[k])) for k in raw};F=[v+e*dv for v,dv in zip(scalar,sr)];return cons(A,Q,MomentSources(A,Q,U,F,mg),yy)[0]
    result=[]
    avs=mg.average(geo_source(A,S,src,al,q).reshape(2,A.n,-1));fv=np.diff(flux(A,S,src,al,np.sqrt(A.ends)),axis=1)/A.h
    vol=mg.average((q*J).reshape(mg.q.shape));alpha_vol=mg.average((al(q)*q*J).reshape(mg.q.shape));Jc,_=metric(A,S,y)
    for eps,dr in ((1e-4,1e-5*A.h),(1e-5,5e-6*A.h)):
        td=(direction(eps,y)-direction(-eps,y))/(2*eps);fd=(flux(A,S,src,al,r+dr)-flux(A,S,src,al,r-dr))/(2*dr);res=(td+fd-geo_source(A,S,src,al,y))/(al(y)*y*Jc)
        qd=(direction(eps,q)-direction(-eps,q))/(2*eps);mres=(mg.average(qd.reshape(2,A.n,-1))+fv-avs)/alpha_vol
        result.append(dict(epsilon=eps,first_point=res[:,0].tolist(),max_point=np.max(abs(res),axis=1).tolist(),norm_point=np.sqrt(np.sum(y*Jc*A.h*res**2,axis=1)).tolist(),first_shell=mres[:,0].tolist(),max_shell=np.max(abs(mres),axis=1).tolist(),norm_shell=np.sqrt(np.sum(vol*A.h*mres**2,axis=1)).tolist()))
    roundtrip=max(np.max(abs(mg.moment(f,j)-v)) for name,fs in src.fluids.items() for f,j,v in zip(fs,(2,2,3),raw[name]))
    return dict(moment_roundtrip=float(roundtrip),dust=dust,exchange=ex.tolist(),telescoping=max(abs(x) for d0 in ledger.values() for x in d0['telescoping']),center_flux=max(abs(x) for d0 in ledger.values() for x in d0['center_flux']),geometry_energy=d['total_error'],conservation=result)

def main():
    rows=[]
    for n in (8,16,32,64):
        A=NativeAssembly(n);S,_=state(A,'compatible');mg=MomentGrid(A);raw,scalar=initial(A,S,mg)
        f=lambda y: 1+.2*y-.3*y*y+.4*y**3;df=lambda y:.2-.6*y+1.2*y*y;polychecks={}
        for j in (2,3):
            rec=mg.reconstruction(mg.moment(f,j),j);test=np.linspace(0,A.Y,501);polychecks[j]=dict(value=float(max(abs(rec(test)-f(test)))),derivative=float(max(abs(rec(test,1)-df(test)))),origin=float(rec(0)-1))
        fflux=np.sqrt(A.ends)**3;rate=-np.diff(fflux)/A.h;re=mg.reconstruction(rate,2);flux_witness=float(max(abs(re(np.linspace(0,A.Y,101))+3)))
        point,_=fixture(A,S);legacy_delta={k:[float(A.h*np.sum(a-b)) for a,b in zip(point[k],raw[k])] for k in raw};row=dict(n=n,polynomial=polychecks,cubic_flux_witness=flux_witness,legacy_minus_true_moment_sums=legacy_delta,native=compare(A,S,mg,raw,scalar),legacy_sum_preserving=compare(A,S,mg,point,scalar));rows.append(row);print(json.dumps(row),flush=True)
    open('results.json','w').write(json.dumps(rows,indent=2))
