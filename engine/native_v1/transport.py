"""Research-only zero-step conservative transport port, no evolution."""
import numpy as np
import scipy.linalg as la
from .native_cmc import NativeAssembly
from .assembly import state,KAP
from .matter_source import Sources,fixture,metric,V0,BETA_D,LAM

def beta_eff(phi):
    u=np.exp(-2*BETA_D*phi);u0=np.exp(-2*BETA_D*.179055)
    a=.93/np.sqrt(1+3*u0);b=.07/np.sqrt(u0)
    return BETA_D*(a*3*u/np.sqrt(1+3*u)+b*np.sqrt(u))/(a*np.sqrt(1+3*u)+b*np.sqrt(u))

def geom(A,S,y):
    k=np.array([A.eval(S['k'],x) for x in y]);ky=np.array([A.eval(S['k'],x,der=1) for x in y]);l=np.array([A.eval(S['ell'],x) for x in y]);ly=np.array([A.eval(S['ell'],x,der=1) for x in y]);L=S['L0']+y*l
    J=np.exp(3*L-y*k/2);s2=np.exp(-2*L+y*k);V=l+y*ly;theta=k+y*ky
    return J,s2,V,theta

def recover(name,f,y,J,s2):
    d,e,m=[v(y) for v in f];D=d/J;E=e/J;Sr=np.sqrt(y)*m/J;S2=s2*Sr**2
    if name=='radiation':
        if np.any(d!=0) or np.any(E<0) or np.any((S2>=E**2)&(S2>0)):raise ValueError('inadmissible radiation')
        P=np.zeros_like(E);mask=E>0;P[mask]=(E[mask]**2-S2[mask])/(np.sqrt(4*E[mask]**2-3*S2[mask])+E[mask]);Er=E;rho=3*P
    else:
        if np.any(D<=0):raise ValueError('inadmissible dust')
        Er=np.sqrt(D**2+S2);P=np.zeros_like(E);rho=D**2/Er
    v=np.zeros_like(E);mask=Er+P>0;v[mask]=s2[mask]*Sr[mask]/(Er+P)[mask];v2=v*v/s2
    if np.any(v2>=1):raise ValueError('noncausal velocity')
    return D,E,Er,Sr,P,rho,v,1/np.sqrt(1-v2)

def scalar_rhs(A,S,src,al,y):
    J,s2,V,theta=geom(A,S,y);a=al(y);ay=al(y,1);vals=[f(y) for f in src.fields];SS,PS,DD,PD,phi,Pi=vals
    ed=np.exp(-2*BETA_D*phi);vp=[SS+LAM/2*SS**3+1.5*LAM*SS*DD**2,-ed*DD+LAM/2*DD**3+1.5*LAM*DD*SS**2,-V0*np.exp(-phi)+KAP*BETA_D*ed*DD**2]
    rho_dm=recover('dark_matter',src.fluids['dark_matter'],y,J,s2)[5];vp[2]-=beta_eff(phi)*rho_dm
    out=[]
    for i in range(3):
        f=src.fields[2*i];mom=vals[2*i+1];fy=f(y,1);fyy=f(y,2);lap=s2*(4*y*fyy+(6+2*y*(theta+2*V))*fy)
        out.extend([a*mom,a*(lap+S['K']*mom-vp[i])+4*y*s2*ay*fy])
    return out

def rhs(A,S,src,al):
    y=A.cent;r=np.sqrt(y);J,s2,V,theta=geom(A,S,y);a=al(y);ar=2*r*al(y,1);p=np.array([A.eval(S['p'],x) for x in y]);Aa=y*p;phi=src.fields[4](y);Pi=src.fields[5](y);pr=2*r*src.fields[4](y,1);be=beta_eff(phi)
    yf=A.ends;rf=np.sqrt(yf);Jf,s2f,_,_=geom(A,S,yf);af=al(yf);out={};led={}
    for name,f in src.fluids.items():
        D,E,Er,Sr,P,rho,v,W=recover(name,f,y,J,s2);Df,Ef,Erf,Srf,Pf,rhof,vf,Wf=recover(name,f,yf,Jf,s2f)
        flux=af*yf*Jf*np.array([Df*vf,s2f*Srf,Srf*vf+Pf]);root=y*J;ani=Sr*v;stress_r=ani+P
        geoE=root*(a*((S['K']/3+Aa)*stress_r+2*(S['K']/3-Aa/2)*P)-s2*Sr*ar)
        # radial metric log derivative = 2r(2V-theta); angular log derivative = 2/r+4rV
        geoM=root*(-Er*ar+a*(stress_r*r*(2*V-theta)+2*P*(1/r+2*r*V)))
        exchange=np.zeros((3,A.n))
        if name=='dark_matter':exchange=np.array([-a*root*be*rho*W*(Pi+v*pr),-a*root*be*rho*Pi,a*root*be*rho*pr])
        source=np.array([np.zeros_like(y),geoE,geoM])+exchange;du=-np.diff(flux,axis=1)/A.h+source;out[name]=tuple(du)
        led[name]=dict(telescoping=(A.h*np.sum(du-source,axis=1)+flux[:,-1]-flux[:,0]).tolist(),center_flux=flux[:,0].tolist(),outer_flux=flux[:,-1].tolist(),exchange_sum=(A.h*exchange.sum(axis=1)).tolist())
    return out,scalar_rhs(A,S,src,al,y),led

def main():
    import json
    results=[]
    for n in (8,16,32,64):
        A=NativeAssembly(n);S,_=state(A,'compatible');raw,scalar=fixture(A,S);src=Sources(A,S,raw,scalar);gauge=A.cmc(*(S[k] for k in ('k','ell','L0','p','K','incoming')),src);G,d=A.rhs(S,src,.13,prepared=gauge);al=G['alpha'];mr,sr,ledger=rhs(A,S,src,al);y=A.cent;J,s2=metric(A,S,y)
        st=np.array([-G['L0']-x*A.eval(G['ell'],x)+x*A.eval(G['k'],x)/2 for x in y]);dust={}
        for name in ('dark_matter','baryons'):
            D,E,M=raw[name];Dt,Et,Mt=mr[name];Er=np.sqrt(D**2+s2*M**2);Ert=(D*Dt+s2*M*Mt+s2*st*M*M)/Er;dust[name]=dict(max_abs=float(max(abs(Et-Ert))),signed_sum=float(A.h*np.sum(Et-Ert)))
        # Independently compare scalar exchange to fluid energy and momentum source in Einstein units.
        rho=recover('dark_matter',src.fluids['dark_matter'],y,J,s2)[5];phi=src.fields[4](y);Pi=src.fields[5](y);phir=2*np.sqrt(y)*src.fields[4](y,1);q=beta_eff(phi)*rho;scE=al(y)*y*J*q*Pi/KAP;scM=-al(y)*y*J*q*phir/KAP;ex=np.array(ledger['dark_matter']['exchange_sum'])/KAP
        # Exact additive split at frozen gauge. Local explicit term alpha K p is projected with physical shear mass.
        kv=A.P@S['k'];lv=A.P@S['ell'];ss=np.exp(-2*(S['L0']+A.y*lv)+A.y*kv);Ms=A.P.T@((A.wp/ss)[:,None]*A.P);local=lambda p:la.solve(Ms,A.P.T@(A.wp*al(A.y)*S['K']*(A.P@p)/ss),assume_a='pos');stiff=G['p']-local(S['p']);S2=dict(S,p=S['p']+A.projp(.001+.002*A.y));G2,_=A.rhs(S2,src,.13,prepared=gauge);split=float(np.linalg.norm(G2['p']-local(S2['p'])-stiff))
        rec=dict(n=n,fluid=ledger,dust_energy_tangent=dust,scalar_exchange_net=[float(A.h*sum(scE)+ex[1]),float(A.h*sum(scM)+ex[2])],shear_partition_fixed_gauge_defect=split,geometry_energy_error=d['total_error'],scalar_rhs_max=float(np.max(abs(np.array(sr)))))
        results.append(rec);print(json.dumps(rec),flush=True)
    open('transport.json','w').write(json.dumps(results,indent=2))
