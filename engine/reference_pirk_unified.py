import numpy as np, math
G=1.0; KAPPA=8*np.pi; V0=8.242522415500654e-5; Dpot=0.10; beta_dm=-0.04
class Grid:
    def __init__(self,n,rmax=40.0):
        self.n=n; self.r_max=rmax; self.dr=rmax/n
        self.centers=(np.arange(n,dtype=float)+0.5)*self.dr
        self.faces=np.arange(n+1,dtype=float)*self.dr
        self.face_areas=self.faces**2
        f=self.faces; self.volumes=(f[1:]**3-f[:-1]**3)/3.0
    def _d4(self,u,parity,order=1):
        u=np.asarray(u,float); n=self.n; h=self.dr
        if u.shape!=(n,): raise ValueError
        x=np.empty(n+4); x[2:-2]=u; x[1]=parity*u[0]; x[0]=parity*u[1]
        out=np.empty(n)
        for i in range(n-2):
            s=x[i:i+5]
            if order==1: out[i]=(s[0]-8*s[1]+8*s[3]-s[4])/(12*h)
            else: out[i]=(-s[0]+16*s[1]-30*s[2]+16*s[3]-s[4])/(12*h*h)
        for i in (n-2,n-1):
            ids=np.arange(n-5,n); off=ids.astype(float)-float(i)
            mat=np.vstack([off**p for p in range(5)]); tar=np.zeros(5); tar[order]=math.factorial(order)
            w=np.linalg.solve(mat,tar); out[i]=(w@u[ids])/(h**order)
        return out
    def d1(self,u,p=1): return self._d4(u,p,1)
    def d2(self,u,p=1): return self._d4(u,p,2)
    def gradient_faces(self,u,inner=0.,outer=0.):
        g=np.empty(self.n+1); g[0]=inner;g[-1]=outer;g[1:-1]=(u[1:]-u[:-1])/self.dr; return g
    def divergence_faces(self,f): return (self.face_areas[1:]*f[1:]-self.face_areas[:-1]*f[:-1])/self.volumes

def Vc(phi): return V0*(np.exp(-phi)-Dpot)
def dVc(phi): return -V0*np.exp(-phi)
class State:
    names=('a','b','X','alpha','beta','Aa','K','Lambda','B')
    def __init__(self,n):
        for x in self.names: setattr(self,x,np.zeros(n))
        self.alpha[:]=1.; self.a[:]=1.; self.b[:]=1.; self.X[:]=1.
    def copy(self):
        s=State(len(self.a))
        for x in self.names: setattr(s,x,getattr(self,x).copy())
        return s

def D(g,u,p): return g.d1(u,p)
def D2(g,u,p): return g.d2(u,p)
def ricci_terms(g,s):
    r=g.centers;a=s.a;b=s.b;X=s.X
    ap=D(g,a,1);bp=D(g,b,1);app=D2(g,a,1);bpp=D2(g,b,1)
    chip=-0.5*D(g,X,1)/X
    chipp=-0.5*(D2(g,X,1)/X-(D(g,X,1)/X)**2)
    Lp=D(g,s.Lambda,-1); metric_lambda=(1-a/b)/r**2; inv=X**2/a
    rb=(0.5*app/a-a*Lp-0.75*(ap/a)**2+0.5*(bp/b)**2-0.5*s.Lambda*ap+ap/(r*b)+2*metric_lambda*(1+r*bp/b)+4*chipp-2*chip*(ap/a-bp/b-2/r))
    sb=(0.5*app/a+bpp/b-a*Lp-(ap/a)**2+0.5*(bp/b)**2+2*(3-a/b)*bp/(r*b)+4*metric_lambda+8*(chipp+chip**2)-8*chip*(0.5*ap/a-bp/b-2/r))
    Rrr=-inv*rb; R=-inv*sb
    alphap=D(g,s.alpha,1); alphapp=D2(g,s.alpha,1)
    DrDr=inv*(alphapp-alphap*(0.5*ap/a+2*chip))
    lap=inv*(alphapp-alphap*(0.5*ap/a-bp/b-2*chip-2/r))
    return R,Rrr,DrDr,lap

def ko(g,u,parity,epsilon=0.1):
    ext=np.empty(g.n+4); ext[2:-2]=u; ext[1]=parity*u[0]; ext[0]=parity*u[1]; ext[-2]=2*u[-1]-u[-2]; ext[-1]=3*u[-1]-2*u[-2]
    st=ext[:-4]-4*ext[1:-3]+6*ext[2:-2]-4*ext[3:-1]+ext[4:]
    out=-epsilon*st/(16*g.dr); out[-2:]=0.0; return out

def divbeta(g,s):
    r=g.centers; ap=D(g,s.a,1);bp=D(g,s.b,1); betap=D(g,s.beta,-1)
    return betap+s.beta*(0.5*ap/s.a+bp/s.b+2/r)

def momentum_constraint(g,s):
    r=g.centers; Aa=s.Aa;Ab=-0.5*Aa;chip=-0.5*D(g,s.X,1)/s.X
    return D(g,Aa,1)-(2/3)*D(g,s.K,1)+6*Aa*chip+(Aa-Ab)*(2/r+D(g,s.b,1)/s.b)

def lambda_l2(g,s,lm=2.0):
    r=g.centers; db=divbeta(g,s); Aa=s.Aa;Ab=-0.5*Aa; ap=D(g,s.alpha,1); Aap=D(g,Aa,1)
    return D2(g,s.beta,-1)/s.a+2*D(g,s.beta/r,1)/s.b+D(g,db,1)/(3*s.a)-2*(Aa*ap+s.alpha*Aap)/s.a-4*s.alpha*(Aa-Ab)/(r*s.b)+lm*s.alpha*momentum_constraint(g,s)/s.a

def lambda_l3(g,s):
    db=divbeta(g,s); return s.beta*D(g,s.Lambda,-1)-s.Lambda*D(g,s.beta,-1)+(2/3)*s.Lambda*db+2*s.alpha*s.Aa*s.Lambda+ko(g,s.Lambda,-1)

def matter_projection(g,s,fields):
    S,PS,Df,PD,phi,Pi,rdm,rb=fields
    e=np.zeros(g.n); sr=np.zeros(g.n); sa=np.zeros(g.n); j=np.zeros(g.n)
    invr=s.X**2/s.a
    for f,p,V,Vp,sc in ((S,PS,.5*S*S,S,1.),(Df,PD,-.5*Df*Df,-Df,1.),(phi,Pi,Vc(phi),dVc(phi),KAPPA)):
        fp=D(g,f,1); gr2=invr*fp*fp
        e += (.5*p*p+.5*gr2+V)/sc
        sr += (.5*p*p+.5*gr2-V)/sc
        sa += (.5*p*p-.5*gr2-V)/sc
        j += -p*fp/sc
    e += (rdm+rb)/KAPPA
    return e,sr,sa,j

def scalar_parts(g,s,f,p,Vp):
    fp=D(g,f,1); fpp=D2(g,f,1); invr=s.X**2/s.a
    lap=invr*(fpp-fp*(0.5*D(g,s.a,1)/s.a-D(g,s.b,1)/s.b+D(g,s.X,1)/s.X-2/g.centers))
    l2=s.alpha*lap+(s.X**-2)*D(g,s.alpha,1)*fp
    adv=s.beta*D(g,f,1)
    return s.alpha*p+adv, l2, s.alpha*s.K*p-s.alpha*Vp+s.beta*D(g,p,1)

def primary_l2(g,s):
    R,Rrr,DrDr,lap=ricci_terms(g,s)
    return {'Aa':-(DrDr-lap/3)+s.alpha*(Rrr-R/3),'K':-lap}

def primary_l3(g,s,fields):
    Aa=s.Aa;Ab=-0.5*Aa;e,sr,sa,j=matter_projection(g,s,fields)
    return {'Aa':s.beta*D(g,Aa,1)+s.alpha*s.K*Aa-(16*np.pi/3)*s.alpha*(sr-sa)+ko(g,Aa,1),
            'K':s.beta*D(g,s.K,1)+s.alpha*(Aa*Aa+2*Ab*Ab+s.K*s.K/3)+4*np.pi*s.alpha*(e+sr+2*sa)+ko(g,s.K,1)}

def lambda_l3_matter(g,s,fields,lm=2.0):
    e,sr,sa,j=matter_projection(g,s,fields)
    return lambda_l3(g,s)-8*np.pi*lm*s.alpha*j/s.a

def make_initial(g,Aamp=.01,width=7.):
    s=State(g.n); r=g.centers
    H0=math.sqrt((.5*0.179055**2+float(Vc(np.array([33.8983]))[0])+2.5857e-5+4.0306e-6)/3)
    S=Aamp*np.exp(-(r/width)**2); Df=1e-10*np.exp(-(r/width)**2); PS=np.zeros_like(r); PD=Df.copy()
    phi=np.full_like(r,33.8983); Pi=np.full_like(r,0.179055); rdm=np.full_like(r,2.5857e-5); rb=np.full_like(r,4.0306e-6)
    B=np.ones_like(r); Sp=-2*r/width**2*S; Dp=-2*r/width**2*Df
    for _ in range(6):
        Kr=-H0-4*np.pi*r*PD*Dp
        rho=.5*(PS*PS+B*Sp*Sp)+.5*S*S+.5*(PD*PD+B*Dp*Dp)-.5*Df*Df+(.5*Pi*Pi+Vc(phi))/KAPPA+(rdm+rb)/KAPPA
        src=1-8*np.pi*r*r*rho+2*r*r*Kr*(-H0)+r*r*H0*H0
        rr=np.r_[0.,r]; ss=np.r_[1.,src]; I=np.cumsum(0.5*(rr[1:]-rr[:-1])*(ss[1:]+ss[:-1])); B=I/r
    A=1/np.sqrt(B); Kt=np.full_like(r,-H0); Kr=Kt-4*np.pi*r*PD*Dp; K=Kr+2*Kt
    s.a=A**(4/3); s.b=A**(-2/3); s.X=A**(-1/3); s.alpha[:]=1.; s.K=K
    s.Aa=(Kr-K/3); s.beta[:]=0.; s.Lambda=D(g,s.a,1)/(2*s.a*s.a)-D(g,s.b,1)/(s.a*s.b)+2/r*(1/s.b-1/s.a); s.B=.75*s.Lambda
    return s,[S,PS,Df,PD,phi,Pi,rdm,rb]

def constraint(g,s,fields):
    R,Rrr,DrDr,lap=ricci_terms(g,s); Aa=s.Aa;Ab=-.5*Aa; H=R-(Aa*Aa+2*Ab*Ab)+(2/3)*s.K*s.K-16*np.pi*matter_projection(g,s,fields)[0]
    M=momentum_constraint(g,s)-8*np.pi*matter_projection(g,s,fields)[3]
    conn=s.Lambda-(D(g,s.a,1)/(2*s.a*s.a)-D(g,s.b,1)/(s.a*s.b)+2/g.centers*(1/s.b-1/s.a))
    det=s.a*s.b*s.b-1
    return H,M,conn,det