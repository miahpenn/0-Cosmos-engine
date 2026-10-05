import numpy as np, math, time, json, importlib.util

# V5.5 physics constants (archive locked)
G=1.0
eight_pi_G=8*np.pi*G
KAPPA=8*np.pi*G
V0=8.242522415500654e-5
Dpot=0.10
Y0=np.array([33.8983,0.179055,0.00341913,2.5857e-5,4.0306e-6,0.,0.])
NG=3

# ---------------- exact Engrenage spherical background ----------------
class Background:
    def __init__(self,r):
        N=len(r); self.r=r
        sv=np.zeros((N,3)); sv[:,0]=1; sv[:,1]=r; sv[:,2]=r
        self.scaling_vector=sv; self.inverse_scaling_vector=1/sv
        dsv=np.zeros((N,3,3)); dsv[:,1,0]=1; dsv[:,2,0]=1
        self.d1_scaling_vector=dsv
        disv=np.zeros_like(dsv); disv[:,1,0]=-1/r**2; disv[:,2,0]=-1/r**2
        self.d1_inverse_scaling_vector=disv
        d2isv=np.zeros((N,3,3,3)); d2isv[:,1,0,0]=2/r**3; d2isv[:,2,0,0]=2/r**3
        self.d2_inverse_scaling_vector=d2isv
        d2sv=np.zeros((N,3,3,3)); self.d2_scaling_vector=d2sv
        # at equator: only theta-theta component has nonzero second derivative
        d2sv[:,2,1,1]=-r
        sm=sv[:,:,None]*sv[:,None,:]; self.scaling_matrix=sm; self.inverse_scaling_matrix=1/sm
        d1sm=np.zeros((N,3,3,3))
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    d1sm[:,i,j,k]=sv[:,i]*dsv[:,j,k]+sv[:,j]*dsv[:,i,k]
        self.d1_scaling_matrix=d1sm
        d2sm=np.zeros((N,3,3,3,3))
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        d2sm[:,i,j,k,l]=(dsv[:,i,l]*dsv[:,j,k]+sv[:,i]*d2sv[:,j,k,l]+dsv[:,j,l]*dsv[:,i,k]+sv[:,j]*d2sv[:,i,k,l])
        self.d2_scaling_matrix=d2sm
        hg=np.zeros((N,3,3)); hg[:,0,0]=1;hg[:,1,1]=r*r;hg[:,2,2]=r*r;self.hat_gamma_LL=hg
        hc=np.zeros((N,3,3,3));
        hc[:,0,1,1]=-r;hc[:,0,2,2]=-r
        hc[:,1,0,1]=hc[:,1,1,0]=1/r
        hc[:,2,0,2]=hc[:,2,2,0]=1/r
        self.hat_christoffel=hc
        dhc=np.zeros((N,3,3,3,3))
        dhc[:,0,1,1,0]=-1;dhc[:,0,2,2,0]=-1
        dhc[:,1,2,2,1]=1
        dhc[:,1,0,1,0]=-1/r**2;dhc[:,1,1,0,0]=-1/r**2
        dhc[:,2,0,2,0]=-1/r**2;dhc[:,2,2,0,0]=-1/r**2
        dhc[:,2,1,2,1]=-1;dhc[:,2,2,1,1]=-1
        self.d1_hat_christoffel=dhc
        self.det_hat_gamma=r**4
        ddet=np.zeros((N,3));ddet[:,0]=4*r**3;self.d1_det_hat_gamma=ddet
        d2det=np.zeros((N,3,3));d2det[:,0,0]=12*r*r;self.d2_det_hat_gamma=d2det

# ---------------- fourth-order derivatives / parity ghosts ----------------
def derivative_matrices(n):
    I=np.eye(n)
    d=[None]*7
    d[0]=I.copy()
    d[1]=(1/12)*np.eye(n,k=-2)-2/3*np.eye(n,k=-1)+2/3*np.eye(n,k=1)-(1/12)*np.eye(n,k=2)
    d[2]=-(1/12)*np.eye(n,k=-2)+4/3*np.eye(n,k=-1)-2.5*I+4/3*np.eye(n,k=1)-(1/12)*np.eye(n,k=2)
    d[3]=-.5*np.eye(n,k=-2)+np.eye(n,k=-1)-np.eye(n,k=1)+.5*np.eye(n,k=2)
    d[4]=np.eye(n,k=-2)-4*np.eye(n,k=-1)+6*I-4*np.eye(n,k=1)+np.eye(n,k=2)
    d[5]=-.5*np.eye(n,k=-3)+2*np.eye(n,k=-2)-2.5*np.eye(n,k=-1)+2.5*np.eye(n,k=1)-2*np.eye(n,k=2)+.5*np.eye(n,k=3)
    d[6]=np.eye(n,k=-3)-6*np.eye(n,k=-2)+15*np.eye(n,k=-1)-20*I+15*np.eye(n,k=1)-6*np.eye(n,k=2)+np.eye(n,k=3)
    for k in range(1,7): d[k][:NG]=0; d[k][-NG:]=0
    return d

class Grid:
    def __init__(self,nphys=100,rmax=40.):
        self.nphys=nphys; self.N=nphys+2*NG
        allp=2*(self.N-NG)
        x=np.linspace(-rmax,rmax,allp)
        self.r=x[-self.N:]
        self.dr=x[1]-x[0]
        self.D=derivative_matrices(self.N)
        L=(-1/3)*np.eye(self.N,k=-3)+1.5*np.eye(self.N,k=-2)-3*np.eye(self.N,k=-1)+11/6*np.eye(self.N)
        R=(-11/6)*np.eye(self.N)+3*np.eye(self.N,k=1)-1.5*np.eye(self.N,k=2)+(1/3)*np.eye(self.N,k=3)
        L[:NG]=0; L[-NG:]=0; R[:NG]=0;R[-NG:]=0
        self.advec=np.stack([L,R])
    def d1(self,x): return x@self.D[1].T/self.dr
    def d2(self,x): return x@self.D[2].T/self.dr**2
    def d6(self,x): return x@self.D[6].T/(64*self.dr)
    def adv(self,x,direction):
        return x@self.advec[direction.astype(int),np.arange(self.N)].T/self.dr

PAR=np.array([1,1,1,1,1,1,1,1,-1,-1,-1,1, 1,1,1,1,1,1,1,1],int)

def fill_boundaries(U):
    # Exact inner parity convention for the three ghost points.
    U[:,0]=PAR*U[:,5]; U[:,1]=PAR*U[:,4]; U[:,2]=PAR*U[:,3]
    # Outer values are smooth/asymptotically constant over the short domain used here.
    U[:,-3]=U[:,-4]; U[:,-2]=U[:,-4]; U[:,-1]=U[:,-4]
    return U

# ---------------- BSSN tensor algebra ----------------
def bgamma(r,h,b): return h*b.scaling_matrix+b.hat_gamma_LL
def rgbg(r,h,b): return h+b.hat_gamma_LL*b.inverse_scaling_matrix
def bgu(r,h,b): return np.linalg.inv(bgamma(r,h,b))
def rbgU(r,h,b): return np.linalg.inv(rgbg(r,h,b))
def barA_LL(r,a,b): return b.scaling_matrix*a
def barA_UU(r,a,h,b):
    gu=bgu(r,h,b); return gu@barA_LL(r,a,b)@gu
def barA_sq(r,a,h,b):
    gu=rbgU(r,h,b); au=gu@a@gu; return np.einsum('xij,xij->x',a,au)

def hat_D_bar_gamma(h,d1h,b):
    out=d1h*b.scaling_matrix[:,:,:,None]+b.d1_scaling_matrix*h[:,:,:,None]
    eps=h*b.scaling_matrix
    out-=np.einsum('xlik,xlj->xijk',b.hat_christoffel,eps)+np.einsum('xljk,xil->xijk',b.hat_christoffel,eps)
    return out

def connections(h,d1h,b):
    bg=bgamma(b.r,h,b); gu=np.linalg.inv(bg); hd=hat_D_bar_gamma(h,d1h,b)
    Dull=.5*np.einsum('xil,xklj->xijk',gu,hd)+.5*np.einsum('xil,xjlk->xijk',gu,hd)-.5*np.einsum('xil,xjkl->xijk',gu,hd)
    Du=np.einsum('xjk,xijk->xi',gu,Dull); Dll=np.einsum('xil,xljk->xijk',bg,Dull)
    return Du,Dull,Dll

def hat_D2_bar_gamma(h,d1h,d2h,b):
    r=b.r; bg=bgamma(r,h,b);gu=np.linalg.inv(bg); hd=hat_D_bar_gamma(h,d1h,b);eps=h*b.scaling_matrix
    d1eps=d1h*b.scaling_matrix[:,:,:,None]+b.d1_scaling_matrix*h[:,:,:,None]
    dm=np.zeros((len(r),3,3,3,3))
    for k in range(3):
        for l in range(3): dm[:,:,:,k,l]=b.d1_scaling_matrix[:,:,:,k]*d1h[:,:,:,l]
    out=(np.einsum('xkl,xijkl->xij',gu,d2h)*b.scaling_matrix+np.einsum('xkl,xijkl->xij',gu,b.d2_scaling_matrix)*h+2*np.einsum('xkl,xijkl->xij',gu,dm)
      -np.einsum('xkl,xmlik,xmj->xij',gu,b.d1_hat_christoffel,eps)-np.einsum('xkl,xmljk,xim->xij',gu,b.d1_hat_christoffel,eps)
      -np.einsum('xkl,xmli,xmjk->xij',gu,b.hat_christoffel,d1eps)-np.einsum('xkl,xmlj,ximk->xij',gu,b.hat_christoffel,d1eps)
      -np.einsum('xkl,xijm,xmlk->xij',gu,hd,b.hat_christoffel)-np.einsum('xkl,xmjl,xmik->xij',gu,hd,b.hat_christoffel)-np.einsum('xkl,ximl,xmjk->xij',gu,hd,b.hat_christoffel))
    return out

def bar_ricci(h,d1h,d2h,L,d1L,b):
    r=b.r;bg=bgamma(r,h,b);gu=np.linalg.inv(bg);D=connections(h,d1h,b);Du,Dull,Dll=D
    Lam=b.inverse_scaling_vector*L
    hatDL=d1L*b.inverse_scaling_vector[:,:,None]+b.d1_inverse_scaling_vector*L[:,:,None]+np.einsum('xijk,xk->xij',b.hat_christoffel,Lam)
    hd2=hat_D2_bar_gamma(h,d1h,d2h,b)
    return (-.5*hd2+.5*np.einsum('xki,xkj->xij',bg,hatDL)+.5*np.einsum('xkj,xki->xij',bg,hatDL)+.5*np.einsum('xk,xijk->xij',Du,Dll)+.5*np.einsum('xk,xjik->xij',Du,Dll)
      +np.einsum('xkl,xmki,xjml->xij',gu,Dull,Dll)+np.einsum('xkl,xmkj,ximl->xij',gu,Dull,Dll)+np.einsum('xkl,xmik,xmjl->xij',gu,Dull,Dll))

# ---------------- V5.5 matter ----------------
def Vc(phi): return V0*(np.exp(-phi)-Dpot)
def dVc(phi): return -V0*np.exp(-phi)

def matter_sources(U,v,grid,b):
    N=len(b.r); em=type('EM',(),{})(); em.rho=np.zeros(N);em.Sij=np.zeros((N,3,3));em.Si=np.zeros((N,3));em.S=np.zeros(N)
    em4=np.exp(-4*v.phi); gu=bgu(b.r,v.h,b); bg=bgamma(b.r,v.h,b)
    S,PS,D,PD,phi,Pi,rdm,rb=U[12],U[13],U[14],U[15],U[16],U[17],U[18],U[19]
    fields=[(S,PS,.5*S*S,S,1.),(D,PD,-.5*D*D,-D,1.),(phi,Pi,Vc(phi),dVc(phi),KAPPA)]
    for f,p,V,Vp,sc in fields:
        fp=grid.d1(f)
        vec=np.zeros((N,3));vec[:,0]=fp
        grad2=em4*np.einsum('xij,xi,xj->x',gu,vec,vec)
        em.rho += (.5*p*p+.5*grad2+V)/sc
        em.Si += (-p[:,None]*vec)/sc
        Vt=-p*p+grad2
        sf=-(.5*Vt+V)/em4
        em.Sij += (sf[:,None,None]*bg+np.einsum('xi,xj->xij',vec,vec))/sc
        em.S += em4*np.einsum('xij,xij->x',gu,(sf[:,None,None]*bg+np.einsum('xi,xj->xij',vec,vec)))/sc
    em.rho += (rdm+rb)/KAPPA
    return em

class Vars:
    def __init__(self,N):
        self.N=N; self.h=np.zeros((N,3,3));self.a=np.zeros((N,3,3));self.L=np.zeros((N,3));self.shift=np.zeros((N,3));self.b=np.zeros((N,3));self.phi=np.zeros(N);self.K=np.zeros(N);self.lapse=np.ones(N)

def unpack(U):
    N=U.shape[1];v=Vars(N)
    v.phi=U[0].copy();v.h[:,0,0]=U[1];v.h[:,1,1]=U[2];v.h[:,2,2]=U[3];v.K=U[4].copy();v.a[:,0,0]=U[5];v.a[:,1,1]=U[6];v.a[:,2,2]=U[7];v.L[:,0]=U[8];v.shift[:,0]=U[9];v.b[:,0]=U[10];v.lapse=U[11].copy();return v

def derivatives(U,grid):
    N=grid.N;d1=D1=N
    a=type('D1',(),{})();d=type('D2',(),{})();ad=type('AD',(),{})()
    for name,shape in [('h',(N,3,3,3)),('a',(N,3,3,3)),('shift',(N,3,3)),('L',(N,3,3)),('phi',(N,3)),('K',(N,3)),('lapse',(N,3))]: setattr(a,name,np.zeros(shape))
    for name,shape in [('h',(N,3,3,3,3)),('shift',(N,3,3,3)),('phi',(N,3,3)),('lapse',(N,3,3))]: setattr(d,name,np.zeros(shape))
    for name,shape in [('h',(N,3,3,3)),('a',(N,3,3,3)),('L',(N,3,3)),('phi',(N,3)),('K',(N,3)),('lapse',(N,3))]: setattr(ad,name,np.zeros(shape))
    for Uidx,targ in [(1,a.h),(2,a.h),(3,a.h)]: pass
    a.h[:,0,0,0]=grid.d1(U[1]);a.h[:,1,1,0]=grid.d1(U[2]);a.h[:,2,2,0]=grid.d1(U[3]);a.a[:,0,0,0]=grid.d1(U[5]);a.a[:,1,1,0]=grid.d1(U[6]);a.a[:,2,2,0]=grid.d1(U[7]);a.shift[:,0,0]=grid.d1(U[9]);a.L[:,0,0]=grid.d1(U[8]);a.phi[:,0]=grid.d1(U[0]);a.K[:,0]=grid.d1(U[4]);a.lapse[:,0]=grid.d1(U[11])
    d.h[:,0,0,0,0]=grid.d2(U[1]);d.h[:,1,1,0,0]=grid.d2(U[2]);d.h[:,2,2,0,0]=grid.d2(U[3]);d.shift[:,0,0,0]=grid.d2(U[9]);d.phi[:,0,0]=grid.d2(U[0]);d.lapse[:,0,0]=grid.d2(U[11])
    return a,d,ad

def get_rhs(U,grid,b,return_geom=False):
    fill_boundaries(U)
    v=unpack(U);d1,d2,ad=derivatives(U,grid)
    # exact determinant projection from Engrenage
    det=np.linalg.det(bgamma(grid.r,v.h,b));rf=(det/b.det_hat_gamma if hasattr(b,'det_hat_gamma') else b.det_hat_gamma)**(-1/3)
    bg=bgamma(grid.r,v.h,b);v.h=(rf[:,None,None]*bg-b.hat_gamma_LL)*b.inverse_scaling_matrix
    d1,d2,ad=derivatives(U,grid)
    em=matter_sources(U,v,grid,b)
    # matter_sources expects d1 as radial arrays; pass wrapper
    class DerWrap: pass
    em=matter_sources(U,v,grid,b)
    r=grid.r; gu=bgu(r,v.h,b);bg=bgamma(r,v.h,b);D=connections(v.h,d1.h,b);ch=bar_chris_local(r,D,b)
    Shift=b.inverse_scaling_vector*v.shift
    d1Shift=b.d1_inverse_scaling_vector*v.shift[:,:,None]+d1.shift*b.inverse_scaling_vector[:,:,None]
    div=np.einsum('xii->x',d1Shift)+np.einsum('xiij,xj->x',ch,Shift)
    trA=np.einsum('xij,xij->x',gu,barA_LL(r,v.a,b)); A2=barA_sq(r,v.a,v.h,b); A_LL=barA_LL(r,v.a,b);A_UU=barA_UU(r,v.a,v.h,b)
    phid=-v.lapse*v.K/6+div/6
    hD=(2/3*(v.lapse*trA-div))[:,None,None]*rgbg(r,v.h,b)-2*v.lapse[:,None,None]*v.a
    # exact shift tensor term
    hshift=np.einsum('xjk,xki->xij',b.hat_gamma_LL,d1Shift)+np.einsum('xjk,xkil,xl->xij',b.hat_gamma_LL,b.hat_christoffel,Shift)
    hD += b.inverse_scaling_matrix*hshift+np.transpose(b.inverse_scaling_matrix*hshift,(0,2,1))
    barD2l=np.einsum('xij,xij->x',gu,d2.lapse)-np.einsum('xij,xkij,xk->x',gu,ch,d1.lapse)
    Kd=v.lapse*(v.K*v.K/3+A2+.5*eight_pi_G*(em.rho+em.S))-np.exp(-4*v.phi)*(barD2l+2*np.einsum('xij,xi,xj->x',gu,d1.lapse,d1.phi))
    Rbar=bar_ricci(v.h,d1.h,d2.h,v.L,d1.L,b); Aik=np.einsum('xkl,xik,xlj->xij',gu,A_LL,A_LL)
    TF=v.lapse[:,None,None]*(-2*d2.phi[:,None,None]*0) # placeholder
    dAdTF=v.lapse[:,None,None]*(-2*np.array([[[]]]) if False else 0)
    gradphi=np.zeros((grid.N,3));gradphi[:,0]=d1.phi[:,0];gradl=np.zeros((grid.N,3));gradl[:,0]=d1.lapse[:,0]
    hessphi=d2.phi.copy(); # embedded already
    hesslap=d2.lapse.copy()
    dAdTF=v.lapse[:,None,None]*(-2*np.array([np.diag([x,0,0]) for x in d2.phi[:,0,0]]) + 4*np.einsum('xi,xj->xij',gradphi,gradphi)+2*np.einsum('xkij,xk->xij',ch,gradphi)+Rbar-eight_pi_G*em.Sij)-hesslap+np.einsum('xkij,xk->xij',ch,gradl)+2*np.einsum('xi,xj->xij',gradphi,gradl)+2*np.einsum('xj,xi->xij',gradphi,gradl)
    tr=np.einsum('xij,xij->x',gu,dAdTF)[:,None,None]
    dadt=(b.inverse_scaling_matrix*dAdTF)
    rAik=b.inverse_scaling_matrix*Aik
    adt=-2/3*div[:,None,None]*v.a+v.lapse[:,None,None]*(-2*rAik+v.K[:,None,None]*v.a)+np.exp(-4*v.phi)[:,None,None]*(dadt-tr*rgbg(r,v.h,b)/3)
    # lambda; zero shift at first stage, full expression for later
    # For now exact zero-shift initial gauge; use the full shift terms only after shift nonzero is needed.
    dlamb=-2*np.einsum('xij,xj->xi',A_UU,d1.lapse)+12*v.lapse[:,None]*np.einsum('xij,xj->xi',A_UU,d1.phi)-4/3*v.lapse[:,None]*np.einsum('xij,xj->xi',gu,d1.K)-2*eight_pi_G*v.lapse[:,None]*np.einsum('xij,xj->xi',gu,em.Si)
    dlamb += 2*v.lapse[:,None]*np.einsum('xjk,xijk->xi',A_UU,D[1])
    dlamb *= b.scaling_vector
    out=np.zeros((12,grid.N));out[0]=phid;out[1]=hD[:,0,0];out[2]=hD[:,1,1];out[3]=hD[:,2,2];out[4]=Kd;out[5]=adt[:,0,0];out[6]=adt[:,1,1];out[7]=adt[:,2,2];out[8]=dlamb[:,0];out[9]=0.0;out[10]=0.0;out[11]=0.0
    return out,em,v,D,dualmask_local(grid.r,v,b),Rbar

def bar_chris_local(r,D,b): return b.hat_christoffel+D[1]
def dualmask_local(r,v,b): return np.zeros(len(r))

def initial(nphys=100,rmax=40,Aamp=.01,width=7):
    grid=Grid(nphys,rmax);r=grid.r; phys=slice(NG,-NG);rp=r[phys]
    H0=math.sqrt((.5*Y0[2]**2+Vc(Y0[1])+Y0[3]+Y0[4])/3)
    S=Aamp*np.exp(-(rp/width)**2);D=1e-10*np.exp(-(rp/width)**2);PS=np.zeros_like(rp);PD=D.copy();Phi=np.full_like(rp,Y0[1]);PC=np.full_like(rp,Y0[2]);rdm=np.full_like(rp,Y0[3]);rb=np.full_like(rp,Y0[4])
    # fourth-order-compatible initial gradients from analytic functions
    Sp=-2*rp/width**2*S; Dp=-2*rp/width**2*D
    B=np.ones_like(rp)
    for _ in range(3):
        Kr=np.full_like(rp,-H0)-4*np.pi*G*rp*PD*Dp
        rho=.5*(PS*PS+B*Sp*Sp)+.5*S*S+.5*(PD*PD+B*Dp*Dp)-.5*D*D+(.5*PC*PC+Vc(Phi))/KAPPA+(rdm+rb)/KAPPA
        src=1-8*np.pi*G*rp*rp*rho+2*rp*rp*Kr*(-H0)+rp*rp*H0*H0
        rr0=np.array([0.]);src0=np.array([1.0])
        Rall=np.concatenate([rr0,rp]); Sall=np.concatenate([src0,src])
        I=np.zeros_like(rp);I[:]=0
        vals=.5*(Rall[1:]-Rall[:-1])*(Sall[1:]+Sall[:-1]);I=np.cumsum(vals)
        B=I/rp
    A=1/np.sqrt(B);Kt=np.full_like(rp,-H0);Kr=Kt-4*np.pi*G*rp*PD*Dp;K=Kr+2*Kt;pc=np.log(A)/6
    hrr=A**(4/3)-1; htt=A**(-2/3)-1;hpp=htt.copy();arr=A**(4/3)*(Kr-K/3);att=A**(-2/3)*(Kt-K/3);app=att.copy()
    U=np.zeros((20,grid.N)); vals=[pc,hrr,htt,hpp,K,arr,att,app,np.zeros_like(rp),np.zeros_like(rp),np.zeros_like(rp),np.ones_like(rp),S,PS,D,PD,Phi,PC,rdm,rb]
    for i,val in enumerate(vals):U[i,NG:-NG]=val
    fill_boundaries(U);b=Background(r)
    v=unpack(U);dd,_,_=derivatives(U,grid);Delta=connections(v.h,dd.h,b);U[8]=Delta[0][:,0];fill_boundaries(U)
    return U,grid,b

def H_constraint(U,grid,b):
    v=unpack(U);d1,d2,_=derivatives(U,grid);em=matter_sources(U,v,grid,b);D=connections(v.h,d1.h,b);Rbar=bar_ricci(v.h,d1.h,d2.h,v.L,d1.L,b);gu=bgu(grid.r,v.h,b);ch=bar_chris_local(grid.r,D,b)
    barlap=np.einsum('xij,xij->x',gu,d2.phi)-np.einsum('xij,xkij,xk->x',gu,ch,d1.phi)
    grad2=np.einsum('xij,xi,xj->x',gu,d1.phi,d1.phi);R=np.exp(-4*v.phi)*(np.einsum('xij,xij->x',gu,Rbar)-8*barlap-8*grad2)
    H=R-np.exp(-8*v.phi)*barA_sq(grid.r,v.a,v.h,b)+2*v.K*v.K/3-16*np.pi*G*em.rho
    # lower-index momentum residual
    A_LL=barA_LL(grid.r,v.a,b);A_UU=barA_UU(grid.r,v.a,v.h,b)
    bg=bgamma(grid.r,v.h,b);ch=bar_chris_local(grid.r,D,b)
    # div_j A^j_r for diagonal tensor: d_r A^r_r + Gamma^j_jr A^r_r - Gamma^r_jr A^j_r
    # equivalent: d_r A^rr lowered? use mixed A^j_i, with i=r
    Amix=np.einsum('xjk,xki->xji',gu,A_LL)  # A^j_i
    q=Amix[:,0,0]; qt=Amix[:,1,1];qp=Amix[:,2,2]
    div=grid.d1(q)+(ch[:,0,0,0]+ch[:,1,1,0]+ch[:,2,2,0])*q-ch[:,0,0,0]*q-ch[:,0,1,1]*qt-ch[:,0,2,2]*qp
    M=div+6*q*d1.phi[:,0]-2/3*d1.K[:,0]-8*np.pi*G*np.exp(4*v.phi)*em.Si[:,0]
    return H,M,R,em

if __name__=='__main__':
  for n in (80,120,160,240):
    t=time.time();U,g,b=initial(n,40,.01,7);H,M,R,em=H_constraint(U,g,b);sl=slice(NG+4,-NG-4);print('INIT',n,'Hmax',np.max(abs(H[sl])),'Mmax',np.max(abs(M[sl])),'Rmax',np.max(abs(R[sl])),'rho',np.max(em.rho[sl]),'elapsed',time.time()-t)