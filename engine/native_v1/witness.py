"""Homogeneous initial witness and reference RHS only; no evolution entry point."""
import numpy as np
from . import rate as R
C,I=R.C,R.I
from .matter_source import potential,BETA_D
from .transport import beta_eff
KAP=C.KAP
SC=np.array([.0001,.00002,.0008,-.00003,.1,.0003])
WITNESS_T=.001
STEPS=(.00025,.000125)
mx=lambda v:float(np.max(np.abs(v)))

def pressure_density(b):
    L,K,d,e,db,eb,er=b[:7];f=b[7::2];P=b[8::2];J=np.exp(3*L);c=np.array([1.,1.,1/KAP]);kin=c@(P*P);U=float(potential(*f))
    return (d+db+er)/(J*KAP)+kin/2+U,er/(3*J*KAP)+kin/2-U

def background_initial():
    b=np.r_[0.,0.,.003,.003,.001,.001,.001,SC];b[1]=-np.sqrt(3*KAP*pressure_density(b)[0]);return b

def background_rhs(t,b):
    L,K,d,e,db,eb,er=b[:7];f=b[7::2];P=b[8::2];rho,pressure=pressure_density(b);up=np.asarray(C.sh.up(f[:,None])).reshape(3);g=beta_eff(f[2])*d/np.exp(3*L);ex=-beta_eff(f[2])*d*P[2];p_dot=K*P-up*np.array([1,1,KAP]);p_dot[2]+=g
    return np.r_[-K/3,1.5*KAP*(rho+pressure),ex,ex,0.,0.,K*er/3,np.stack((P,p_dot),axis=1).ravel()]

def invariants(b):
    L,K,d,e,db,eb,er=b[:7];u=np.exp(-2*BETA_D*b[11]);u0=np.exp(-2*BETA_D*.179055);B=.93/np.sqrt(1+3*u0)*np.sqrt(1+3*u)+.07/np.sqrt(u0)*np.sqrt(u)
    return np.array([db,eb,er*np.exp(L),e-d,d/B,2*K*K/3-2*KAP*pressure_density(b)[0]])

def initial(n):
    A=C.Integrated(n);S,_=C.state(A,'flat');mg=C.MomentGrid(A);b=background_initial();S['K']=b[1];w=mg.weights[2];raw={'dark_matter':(b[2]*w,b[3]*w,np.zeros(n)),'baryons':(b[4]*w,b[5]*w,np.zeros(n)),'radiation':(np.zeros(n),b[6]*w,np.zeros(n))};sc=[np.full(n,v) for v in SC];z=np.r_[A.lapse.one,0.,1.];return A,(S,raw,sc),z,b

def refresh(A,x,z):
    S,raw,sc=x;mg=C.MomentGrid(A);reg=I.StructuralRegistry(A,S,raw,sc,mg);bR=R.bR(A,S);incoming=np.r_[z[:-2],bR*z[-2]];al,chi,diag=A.cmc(*(S[k] for k in ('k','ell','L0','p','K')),incoming,reg);return np.r_[al.coeff,al.flux/bR,1.],al,reg

def evaluate(A,x,z,b):
    S,raw,sc=x;before=R.snapshot(*x);zr,al,reg=refresh(A,x,z);J,s2,*_=C.geom(A,S,np.array([A.Y]));P=np.array([reg.fields[i](A.Y) for i in (1,3,5)]);flux=al(A.Y)*A.Y*J[0]*np.sqrt(s2[0])*(b[8::2]-P)
    # Use original incoming z, not refreshed z, for identical CMC solve in full rate call.
    res=R.rate(A,S,raw,sc,z,0.,flux);assert R.snapshot(*x)==before;assert np.array_equal(zr,np.array(res['outgoing_z']));return res,flux
