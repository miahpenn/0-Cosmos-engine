"""Zero-step nonlinear geometry research operator, A8 rational centre.
No time stepping, no production imports. NumPy/SciPy; no physical tuning.
"""
import json, numpy as np, scipy.linalg as la
from scipy.special import roots_jacobi, eval_jacobi, hyp1f1, exprel
from scipy.interpolate import CubicSpline
from numpy.polynomial.legendre import legvander, legint, legder
from fractions import Fraction
from pathlib import Path
KAP=8*np.pi
Z=np.array([[float(Fraction(v)) for v in row] for row in json.load(open(Path(__file__).with_name('A8_certificate.json')))['center']])
t,w=roots_jacobi(48,0,0); L=legvander(t,24); J=np.array([eval_jacobi(k,0,2.5,t) for k in range(25)]).T
C=(L.T@(w[:,None]*np.column_stack([J[:,:5],J[:,5:]@Z])))*((2*np.arange(25)+1)/2)[:,None]
for j in range(5): C[j+1:,j]=0
C[:,0]=0; C[0,0]=1

def val(c,t,order=0,dy=2):
    cc=legder(c,m=order,axis=0)*(2/dy)**order if order else c
    return legvander(np.atleast_1d(t),len(cc)-1)@cc

def derivative(n,h,order):
    D=np.zeros((n,n)); weights=np.array([1,-8,0,8,-1] if order==1 else [-1,16,-30,16,-1])/12/h**order
    for i in range(n-2):
        for j,w in enumerate(weights):
            k=i+j-2; D[i,-k-1 if k<0 else k]+=w
    for i in (n-2,n-1):
        ids=np.arange(n-5,n); A=np.array([(ids-i)**k for k in range(5)],float); b=np.zeros(5); b[order]=1 if order==1 else 2; D[i,ids]=la.solve(A,b)/h**order
    return D

class Assembly:
    def __init__(self,n=8,order=36,R=1.):
        self.n=n; self.R=R; self.Y=R*R; self.h=R/n; self.ends=np.linspace(0,R,n+1)**2;self.r=(np.arange(n)+.5)*self.h; self.cent=self.r**2
        self.np=7+3*(n-1); self.nh=self.np+1; self.e=np.eye(self.nh)[n]; self.cells=[]; self.wp=[]; self.wh=[];self.ys=[]
        self.M=np.zeros((self.np,self.np)); self.N=np.zeros((self.nh,self.nh)); self.D=np.zeros((self.np,self.nh))
        p0=0; b0=n+1
        for i,(lo,hi) in enumerate(zip(self.ends[:-1],self.ends[1:])):
            dy=hi-lo; raw=C.copy() if i==0 else np.eye(3); m=raw.shape[1]; pi=np.arange(p0,p0+m); ix=np.r_[i,i+1,np.arange(b0,b0+m-1)];p0+=m;b0+=m-1
            # Split at spline knot, use y^1.5 Jacobi weight on origin segment, standard Gaussian elsewhere.
            yy=[];ww=[]
            for l,h in [(lo,self.cent[i]),(self.cent[i],hi)]:
                u,v=roots_jacobi(order,0,1.5 if l==0 else 0); y=l+(h-l)*(u+1)/2
                weight=v*((h-l)/2)**2.5 if l==0 else v*(h-l)/2*y**1.5
                yy.extend(y);ww.extend(weight/4)
            y=np.array(yy); wh=np.array(ww);wp=4*y*wh; t=2*(y-lo)/dy-1
            rawF=val(raw,t); Ch=la.cholesky(rawF.T@(wp[:,None]*rawF),lower=True); cp=raw@la.solve_triangular(Ch.T,np.eye(m))
            # Geometry shape: two face functions and primitive bubbles of original nonconstant shear.
            integ=legint(raw[:,1:],lbnd=-1,axis=0)/2; mean=raw[0,1:]; integ[:2]-=np.array([.5,.5])[:,None]*mean
            face=np.zeros((len(integ),2));face[:2]=[[.5,.5],[-.5,.5]]
            # bubble weighted scaling (retains support and trace)
            bv=val(integ,t); chb=la.cholesky(bv.T@(wh[:,None]*bv),lower=True);cb=integ@la.solve_triangular(chb.T,np.eye(m-1));ch=np.column_stack([face,cb]);
            # Force exact endpoint zeros for numerical bubble evaluations; polynomial formula unchanged inside.
            F=val(cp,t); H=val(ch,t); Hy=val(ch,t,1,dy); M=F.T@(wp[:,None]*F); N=H.T@(wh[:,None]*H)
            self.M[np.ix_(pi,pi)]=M;self.N[np.ix_(ix,ix)]+=N;self.D[np.ix_(pi,ix)]=la.solve(M,F.T@(wp[:,None]*Hy),assume_a='pos')
            self.cells.append(dict(pi=pi,ix=ix,cp=cp,ch=ch,lo=lo,hi=hi,y=y,wh=wh,wp=wp,F=F,H=H,Hy=Hy,sl=slice(len(self.ys),len(self.ys)+len(y))))
            self.ys.extend(y);self.wp.extend(wp);self.wh.extend(wh)
        self.y=np.array(self.ys);self.wp=np.array(self.wp);self.wh=np.array(self.wh)
        self.P=np.zeros((len(self.y),self.np));self.H=np.zeros((len(self.y),self.nh));self.Py=self.P.copy();self.Hy=self.H.copy()
        for c in self.cells:
            sl=c['sl']; t=2*(c['y']-c['lo'])/(c['hi']-c['lo'])-1; self.P[sl,c['pi']]=c['F'];self.H[sl,c['ix']]=c['H']; self.Py[sl,c['pi']]=val(c['cp'],t,1,c['hi']-c['lo']);self.Hy[sl,c['ix']]=c['Hy']
        # Exact algebraic Schur condensation, floating arithmetic. Face scaling for solve.
        nf=n+1;self.S=self.N[:nf,:nf].copy();self.blocks=[]
        for c in self.cells:
            ids=c['ix'][2:];D=self.N[np.ix_(ids,ids)]; E=self.N[:nf,ids]; Di=la.cho_factor(D);self.S-=E@la.cho_solve(Di,E.T);self.blocks.append((ids,Di,E))
        self.sfscale=1/np.sqrt(np.diag(self.S));self.Schol=la.cho_factor(self.sfscale[:,None]*self.S*self.sfscale[None,:]);self.D1=derivative(n,self.h,1);self.D2=derivative(n,self.h,2)
        assert np.max(np.abs(self.S[np.abs(np.arange(nf)[:,None]-np.arange(nf)[None,:])>1]),initial=0)==0
    def solveH(self,b):
        rhs=b[:self.n+1].copy();x=np.zeros(self.nh)
        for ids,Di,E in self.blocks:rhs-=E@la.cho_solve(Di,b[ids])
        x[:self.n+1]=self.sfscale*la.cho_solve(self.Schol,self.sfscale*rhs)
        for ids,Di,E in self.blocks:x[ids]=la.cho_solve(Di,b[ids]-E.T@x[:self.n+1])
        return x
    def projp(self,f):return la.solve(self.M,self.P.T@(self.wp*f),assume_a='pos')
    def projH(self,f):return self.solveH(self.H.T@(self.wh*f))
    def B(self,p,hat):return self.solveH(-self.D.T@(self.M@p)+self.e*self.R**5*hat)
    def eval(self,x,y,kind='p',der=0,i=None):
        i=min(np.searchsorted(self.ends,y,side='right')-1,self.n-1) if i is None else i;c=self.cells[i];t=2*(y-c['lo'])/(c['hi']-c['lo'])-1
        if kind=='H' and der==0 and y in (c['lo'],c['hi']): return x[i if y==c['lo'] else i+1]
        return float(val(c['cp' if kind=='p' else 'ch'],[t],der,c['hi']-c['lo'])[0]@x[c['pi' if kind=='p' else 'ix']])
    def trace(self,x,kind='p',der=0):return self.eval(x,self.Y,kind,der)
    def cmc(self,k,ell,L0,p,K,incoming,source):
        y=self.cent;r=self.r;kv=np.array([self.eval(k,t) for t in y]);lv=np.array([self.eval(ell,t) for t in y]);pv=np.array([self.eval(p,t) for t in y]);
        a=np.exp(-2*y*kv/3);b=np.exp(y*kv/3);X=np.exp(y*kv/6-L0-y*lv);inv=X*X/a
        # Same pinned CMC centre row equations and incoming lapse projection; algebraic stencil reproduction, not bitwise vendor RHS.
        D=self.D1; D2=self.D2; ap=D@(a-a[-1]);bp=D@(b-b[-1]);Xp=D@(X-X[-1]);cc=-.5*ap/a+bp/b-Xp/X+2/r
        rho,anis,jhat,trace=source(y);Q=1.5*(y*pv)**2+K*K/3+KAP/2*(rho+trace);wt=r*r*np.sqrt(a)*b/X**3*self.h
        raw=Q*incoming-inv*(D2@(incoming-incoming[-1])+cc*(D@(incoming-incoming[-1])));rate=wt@raw/wt.sum()
        L=np.zeros((self.n,self.n));rhs=-rate/inv;L[0,0]=-1/self.h**2-cc[0]/(2*self.h)-Q[0]/inv[0];L[0,1]=1/self.h**2+cc[0]/(2*self.h)
        for i in range(1,self.n-1):L[i,i-1]=1/self.h**2-cc[i]/(2*self.h);L[i,i]=-2/self.h**2-Q[i]/inv[i];L[i,i+1]=1/self.h**2+cc[i]/(2*self.h)
        L[-1,-1]=1;rhs[-1]=1;alpha=la.solve(L,rhs)
        if np.any(alpha<=0):raise ValueError('nonpositive CMC centre lapse')
        spline=CubicSpline(y,alpha,extrapolate=True);assert alpha[-1]==1 or abs(alpha[-1]-1)<1e-12
        # Incoming affine tangent at frozen state, rank-one denominator measured only.
        bvec=-1/inv;bvec[-1]=0;response=la.solve(L,bvec);direction=(y[-1]-y)/6
        d_raw=Q*direction-inv*(D2@direction+cc*(D@direction));d_rate=wt@d_raw/wt.sum();tan=response*d_rate
        denom=1-wt@(Q*response-inv*(D2@response+cc*(D@response)))/wt.sum()
        eps=1e-4;shifted=la.solve(L,rhs+eps*d_rate*bvec);tan_error=float(np.max(abs((shifted-alpha)/eps-tan)));native_raw=[];sp_in=CubicSpline(y,incoming)
        for i,t in enumerate(y):
            kval=kv[i];ky=self.eval(k,t,der=1);el=lv[i];ly=self.eval(ell,t,der=1);V=el+t*ly;dlog=-kval-t*ky;aa=sp_in(t);ay=sp_in(t,1);ayy=sp_in(t,2)
            native_raw.append(-inv[i]*(4*t*ayy+(6+2*t*(-dlog+2*V))*ay)+Q[i]*aa)
        return spline,rate,dict(cmc_res=float(np.max(abs(L@alpha-rhs))),alpha_min=float(min(alpha)),alpha_face=float(spline(self.Y)),denom=float(denom),incoming_tangent=float(np.max(abs(tan))),incoming_tangent_error=tan_error,native_rate_minus_pinned=float(wt@native_raw/wt.sum()-rate))
    def rhs(self,S,source,g=0.,prepared=None):
        p,k,ell,J,RV,L0,K,kh,eh,w,incoming=(S[t] for t in ('p','k','ell','J','RV','L0','K','kh','eh','w','incoming'));y=self.y;Y=self.Y;R=self.R
        pp=self.P@p;kk=self.P@k;ll=self.P@ell;L=L0+y*ll;s=np.exp(-L+y*kk/2);s2=s*s;X2=np.exp(y*kk/3-2*L);a=np.exp(-2*y*kk/3)
        if not np.all(np.isfinite(s2)) or np.any(s2<=0):raise ValueError('invalid metric')
        bk=self.B(k,kh);be=self.B(ell,eh);z=J+bk/3;Vh=be/4+RV;V=self.H@Vh-1.5*ll;zc=self.H@z
        Hc=-2*be/3+4*bk/9+4*J/3-8*RV/3;H=self.H@Hc;Hy=self.Hy@Hc
        zet=4*ell-2*k/3;zeth=4*eh-2*kh/3;Bzet=self.B(zet,zeth);numer=self.H@Bzet-10*(self.P@zet)
        # k+y*k_y via the same weak derivative, without division by y.
        kd=(self.H@bk-6*kk)/4
        ar,Kt,diag=self.cmc(k,ell,L0,p,K,incoming,source) if prepared is None else prepared;diag=dict(diag);al=ar(y);ay=ar(y,1);ayy=ar(y,2);al0=float(ar(0));alR=float(ar(Y));
        if np.any(al<=0) or al0<=0 or alR<=0:raise ValueError('nonpositive reconstructed lapse')
        beta=(al-al0)/y;betaR=(alR-al0)/Y # quotient of smooth known spline; evaluate stably near 0 below
        # cubic segment Taylor for first interval (including extrapolation to y=0)
        near=y<=self.cent[0];beta[near]=ar(0,1)+y[near]*ar(0,2)/2+y[near]**2*ar(0,3)/6
        rho,anis,jhat,stress=source(y);ap=self.projp(al*pp);betap=self.projp(beta)
        HR=self.e@Hc;kR=self.trace(k);lR=self.trace(ell);LR0=L0+Y*lR;sR=np.exp(-LR0+Y*kR/2);ph=(w-sR*HR)/(2*R)
        kt=3*ap;et=ap/2-K*betap/3;L0t=-K*al0/3;kht=3*alR*ph;eht=alR*ph/2-K*betaR/3
        Jsrc=12*al*y*pp*V-4*y*ay*pp-2*KAP*al*jhat;Jt=self.projH(Jsrc)
        bkt=self.B(kt,kht);bet=self.B(et,eht);zt=Jt+bkt/3;Vt=bet/4
        Ht=-2*bet/3+4*bkt/9+4*Jt/3;G=self.B(ap,alR*ph);Gsrc=2*K/9*self.B(betap,betaR)+4*Jt/3
        zett=4*et-2*kt/3;assert np.linalg.norm(zett+4*K*betap/3)<1e-9
        # curvature remainder after extracting s^2*(Hc_y + D_B zeta).
        kmq=.5*kk*kk*hyp1f1(1,3,-y*kk)
        f=al*(-2/3*(kmq-4*V*V+2*V*kd)+2/9*kd*zc)-2/3*(4*ayy-2*(4*V-kd)*ay)+al*K*pp/s2-2*KAP/3*al*anis/s2
        Ms=self.P.T@((self.wp/s2)[:,None]*self.P)
        lp=self.P.T@(self.wp*al*Hy);lz=self.P.T@(self.wh*al*numer);lf=self.P.T@(self.wp*f);pt=la.solve(Ms,lp+lz+lf,assume_a='pos')
        # actual reconstructed metric tangent: projection defects must not disappear in Mdot.
        logs_t=-L0t-y*(self.P@et)+y*(self.P@kt)/2;logs_cont=al*K/3+al*y*pp
        masswork=-np.dot(self.wp*pp*pp/s2,logs_t);Md=self.P.T@((-2*self.wp*logs_t/s2)[:,None]*self.P)
        # outer physical reduction-based lower terms. Numerical p trace distinct.
        pR=self.trace(p);VR=self.e@Vh-1.5*lR;zR=self.e@z;aR=np.exp(-2*Y*kR/3);X2R=np.exp(Y*kR/3-2*LR0);cR=zR/aR
        kdR=kR+Y*self.trace(k,der=1);db=R*kdR/3;BC=(np.exp(-2*LR0)-sR*sR*(1+2*Y*VR)**2)/Y;lapTF=Y*2/3*sR*sR*(4*ar(Y,2)-2*(4*VR-kdR)*ar(Y,1))
        rhoR,anisR,jhatR,_=[float(v[0]) for v in source(np.array([Y]))]
        lower=-2*R*aR*sR*Y*pR*cR+2*BC+4*X2R*(1/R+db)*R*cR+2*lapTF/alR-2*Y*K*pR+2*K*K/9-(Y*pR)**2/2-2*KAP*rhoR/3+4*KAP*Y*anisR/3
        atlog=-2*Y*self.trace(kt)/3;ctR=(self.e@zt-atlog*zR)/aR;stR=sR*(-L0t-Y*self.trace(et)+Y*self.trace(kt)/2)
        ph_t=alR/(2*Y)*(g-R*aR*sR/alR*ctR-lower)
        wt=2*R*ph_t+sR*(self.e@Ht)+stR*HR
        rem=alR/R*(g-lower+3*HR)+sR*(self.e@Ht)-aR*sR*ctR+stR*HR
        # independently assembled total energy rate and ledger
        bulk_rate=p@Ms@pt+Hc@self.N@Ht+.5*p@Md@p
        bnd=R**5*alR*ph*HR;srcp=p@lf;srcz=p@lz;srch=Hc@self.N@Gsrc
        total=bulk_rate+R**5/6*w*wt;pred=srcp+srcz+srch+masswork-R**4*alR*sR/2*HR*HR+R**5/6*w*rem
        # projected continuous test gives exact weak distribution identity for zeta work.
        T=self.projH(al*pp);Tv=self.H@T; cell=np.dot(self.wp*Tv,self.Py@zet);jumps=0.
        for i in range(1,self.n):
            t=self.ends[i];jumps+=t**2.5*T[i]*(self.eval(zet,t,i=i)-self.eval(zet,t,i=i-1))
        outer=Y**2.5*T[self.n]*(zeth-self.trace(zet));projection=10*np.dot(self.wh*(Tv-al*pp),self.P@zet)
        # independently evaluate strong continuum physical shear, not used to set the weak load.
        ky=self.Py@k;ks=kk+y*ky;Vs=ll+y*(self.Py@ell);Vsy=2*(self.Py@ell)+y*np.concatenate([val(c['cp'],2*(c['y']-c['lo'])/(c['hi']-c['lo'])-1,2,c['hi']-c['lo'])@ell[c['pi']] for c in self.cells]);
        # Compare physical and reduction geometry separately: reduced V derivative within each cell.
        Vy=self.Hy@Vh-1.5*(self.Py@ell);zy=self.Hy@z;zety=self.Py@zet
        fs=al*(-2/3*(kmq-4*V*V+2*V*ks)+2/9*ks*zc)-2/3*(4*ayy-2*(4*V-ks)*ay)+al*K*pp/s2-2*KAP/3*al*anis/s2
        strong=s2*(al*(Hy+zety)+fs)
        # Physical constraint comparator derives C independently from the actual metric,
        # rather than reusing the weak metric connection residual.
        kyy=np.concatenate([val(c['cp'],2*(c['y']-c['lo'])/(c['hi']-c['lo'])-1,2,c['hi']-c['lo'])@k[c['pi']] for c in self.cells]);zphys=self.H@J+(10*kk+4*y*ky)/3;zphysy=self.Hy@J+(14*ky+4*y*kyy)/3
        strong_physical=s2*(al*(-2/3*(ky+kmq+4*Vsy-4*Vs*Vs+2*Vs*ks)+4/3*zphysy+2/9*ks*zphys)-2/3*(4*ayy-2*(4*Vs-ks)*ay))+al*K*pp-2*KAP/3*al*anis
        # physical metric constraint witnesses using cell-strong reconstructed metric, independent connection.
        Bcurv=-np.exp(-2*L)*kk*exprel(y*kk)-s2*(4*Vs+4*y*Vs**2)
        ric=-2/3*s2*(ky+kmq+4*Vsy-4*Vs*Vs+2*Vs*ks);Rscalar=6*Bcurv+6*y*ric
        cy=(zphysy+2/3*ks*zphys)/a; c=zphys/a;Zphys=X2*(3*c+2*y*cy);Hm=Rscalar+2*K*K/3-1.5*(y*pp)**2-2*KAP*rho;M= np.sqrt(y)*(5*pp+2*y*(self.Py@p)+6*y*pp*Vs-KAP*jhat)
        # independently reconstructed physical connection differs from the weak reduction-based metric connection.
        zphys=self.H@J+(10*kk+4*y*ky)/3
        # physical characteristic at actual outer face, from independent radial constraint formulas (reduction comparator separate).
        Vphys=lR+Y*self.trace(ell,der=1);Vphysy=2*self.trace(ell,der=1)+Y*self.trace(ell,der=2);kyR=self.trace(k,der=1);kmqR=.5*kR*kR*hyp1f1(1,3,-Y*kR)
        Bphys=(np.exp(-2*LR0)-sR*sR*(1+2*Y*Vphys)**2)/Y;ricR=-2/3*sR*sR*(kyR+kmqR+4*Vphysy-4*Vphys**2+2*Vphys*kdR);zphysR=self.trace(J,'H')+(10*kR+4*Y*kyR)/3;zphysyR=self.trace(J,'H',1)+(14*kyR+4*Y*self.trace(k,der=2))/3;Zr=X2R/aR*(3*zphysR+2*Y*(zphysyR+2/3*kdR*zphysR))
        Hmr=6*Bphys+6*Y*ricR+2*K*K/3-1.5*(Y*pR)**2-2*KAP*rhoR;Mr=R*(5*pR+2*Y*self.trace(p,der=1)+6*Y*pR*Vphys-KAP*jhatR)
        Wphys=2*sR*Mr+(Hmr+4*Zr)/3
        l2=lambda u,weight=self.wp:float(np.sqrt(np.dot(weight,u*u)))
        # origin quotient checks: raw rational-centre residual, not divided by y.
        origin0=self.eval(Bzet,0,'H')-10*self.eval(zet,0);origin1=self.eval(Bzet,0,'H',1)-14*self.eval(zet,0,der=1)
        # Actual shared interface work: local geometry load residual yields face coefficients.
        face=[]
        for c0 in self.cells:
            ix=c0['ix'];sl=c0['sl'];res=c0['H'].T@(c0['wh']*(self.H[sl]@G))+c0['Hy'].T@(c0['wp']*(self.P[sl]@ap));face.append(res)
        cell_errors=[]
        for i,c0 in enumerate(self.cells):
            sl=c0['sl'];ix=c0['ix'];pwork=np.dot(c0['wp']*al[sl]*pp[sl],Hy[sl]);hwork=np.dot(c0['wh']*H[sl],self.H[sl]@G);cell_errors.append(pwork+hwork-Hc[i]*face[i][0]-Hc[i+1]*face[i][1])
        jump=max(abs(face[i][1]+face[i+1][0]) for i in range(self.n-1));bubble=max(np.max(abs(fa[2:])) for fa in face)
        diag.update(principal_cell_error=float(max(abs(np.array(cell_errors)))),J_projection_work_error=float(Hc@self.N@Jt-np.dot(self.wh*H,Jsrc)),J_projection_residual=float(np.linalg.norm(self.N@Jt-self.H.T@(self.wh*Jsrc))),lapse_product_projection_error=float(np.linalg.norm(self.M@ap-self.P.T@(self.wp*al*pp))),principal_error=float(p@lp+Hc@self.N@G-bnd),Hc_kinematic=float(np.linalg.norm(Ht-G-Gsrc)),total_error=float(total-pred),total_rate=float(total),bulk_rate=float(bulk_rate),zeta_work=float(srcz),curvature_lapse_matter_work=float(srcp),geometry_source_work=float(srch),mass_work=float(masswork),boundary_remainder_work=float(R**5/6*w*rem),boundary_negative_term=float(-R**4*alR*sR*HR*HR/2),zeta_jump_work=float(jumps),zeta_outer_work=float(outer),zeta_test_projection_work=float(projection),zeta_pair_error=float(srcz-cell-jumps-outer-projection),mass_identity=float(.5*p@Md@p-masswork),mass_directional_error=float(abs(np.dot(self.wp*pp*pp/s2,np.sinh(-2*logs_t*1e-5))/(2e-5)-masswork)),mass_continuum_difference=l2(logs_t-logs_cont),pt_norm=l2(self.P@pt),strong_shear_defect=l2(self.P@pt-strong),physical_shear_defect=l2(self.P@pt-strong_physical),reduction_defect=l2(V-Vs),connection_rep_defect=l2(zc-zphys),H_m_norm=l2(Hm,self.wh),H_T_norm=l2(Hm+Zphys,self.wh),H_pin_norm=l2(Hm+X2*(c+2*y*cy),self.wh),momentum_norm=l2(M,self.wh),connection_norm=l2(c*np.sqrt(y),self.wh),physical_W=float(Wphys),temporal_boundary_error=float(R*aR*sR/alR*ctR+2*Y/alR*ph_t+lower-g),trace_delta=float(ph-pR),trace_delta_t=float(ph_t-self.trace(pt)),physical_flux_difference=float(R**5*alR*(ph-pR)*HR),k_trace_delta_t=float(kht-self.trace(kt)),q_trace_delta_t=float(np.exp(-Y*kh)*kht-np.exp(-Y*kR)*self.trace(kt)),ell_trace_delta_t=float(eht-self.trace(et)),zeta_origin_num=float(origin0),zeta_origin_jet_error=float(origin1),face_jump_error=float(jump),bubble_balance_error=float(bubble),outer_load_error=float(face[-1][1]-R**5*alR*ph),centre_load=float(face[0][0]),sparse_res=float(np.linalg.norm(self.N@G-(-self.D.T@(self.M@ap)+self.e*R**5*alR*ph))),sparse_monolithic=float(np.linalg.norm(G-la.solve(self.N,-self.D.T@(self.M@ap)+self.e*R**5*alR*ph,assume_a='pos'))))
        return dict(p=pt,k=kt,ell=et,J=Jt,RV=np.zeros_like(RV),L0=L0t,K=Kt,kh=kht,eh=eht,w=wt,numerical_trace=ph,numerical_trace_t=ph_t,Hc=Hc,Hc_t=Ht,alpha=ar),diag

def no_source(y):return (np.zeros_like(y),)*4

def sources(y):return .002+.001*y,.0003+.0002*y,.0004+.0001*y,.001+.002*y

def state(A,name):
    y=A.y; z=np.zeros_like(y); flat=name=='flat'; broken=name in ('broken','offset_jump'); pert=name=='independent'; matter=name=='matter';
    k=z if flat else .04+.01*y+.003*y*y;ell=z if flat else .025+.006*y+.002*y*y;p=z if flat else .03*np.exp(y/5)
    if broken:p=np.where((y>=A.ends[1])&(y<A.ends[2]),.06,0)
    if name=='offset_jump':ell=ell+np.where((y>=A.ends[1])&(y<A.ends[2]),.005,0)
    J=z if flat else .005+.004*y;RV=np.where(y>=0,.003,0) if pert else z
    if pert:J=J+.01*y
    K=0. if flat else .08;L0=0. if flat else .01;kh=0. if flat else .053;eh=0. if flat else .033
    st={f:A.projp(v) for f,v in [('p',p),('k',k),('ell',ell)]};st.update(J=A.projH(J),RV=A.projH(RV),L0=L0,K=K,kh=kh,eh=eh,incoming=np.ones(A.n) if flat else 1+.02*(A.cent-A.cent[-1]))
    Hc=-2*A.B(st['ell'],eh)/3+4*A.B(st['k'],kh)/9+4*st['J']/3-8*st['RV']/3;sR=np.exp(-L0-A.Y*A.trace(st['ell'])+A.Y*A.trace(st['k'])/2)
    ph=0 if broken or flat else .03*np.exp(.2);st['w']=2*A.R*ph+sR*(A.e@Hc)
    return st,sources if matter else no_source

def main():
    results=[]
    for n in (8,16,32,64):
        A=Assembly(n)
        for name in (('flat','smooth','compatible','broken','offset_jump','independent','matter') if n==8 else ('smooth','compatible')):
            S,src=state(A,name);before={k:(v.copy() if isinstance(v,np.ndarray) else v) for k,v in S.items()};g=.13 if name=='matter' else json.load(open('boundary_datum.json'))['g'] if name=='compatible' else 0.
            rhs,d=A.rhs(S,src,g);assert all(np.array_equal(v,S[k]) for k,v in before.items());
            # Verification errors are numerical measurements, not production thresholds.
            assert all(np.isfinite(v) for v in d.values());results.append(dict(n=n,state=name,**d))
            print(n,name,json.dumps(d,sort_keys=True),flush=True)
        if n==8:
            # same state and dimension at higher quadrature order; rebuild rather than reuse data.
            A2=Assembly(n,48);S2,src2=state(A2,'smooth');_,d2=A2.rhs(S2,src2);_,d1=A.rhs(*state(A,'smooth'))
            print('quadrature comparison',json.dumps({k:d2[k]-d1[k] for k in d2},sort_keys=True),flush=True)
    Path('results.json').write_text(json.dumps(results,indent=2))

def split_rhs(self,S,source,g=0.,prepared=None):
    p,k,ell,J,RV,L0,K,kh,eh,w,incoming=(S[t] for t in ('p','k','ell','J','RV','L0','K','kh','eh','w','incoming'));y=self.y;Y=self.Y;R=self.R
    pp=self.P@p;kk=self.P@k;ll=self.P@ell;L=L0+y*ll;s=np.exp(-L+y*kk/2);s2=s*s;X2=np.exp(y*kk/3-2*L);a=np.exp(-2*y*kk/3)
    if not np.all(np.isfinite(s2)) or np.any(s2<=0):raise ValueError('invalid metric')
    bk=self.B(k,kh);be=self.B(ell,eh);z=J+bk/3;Vh=be/4+RV;V=self.H@Vh-1.5*ll;zc=self.H@z
    Hc=-2*be/3+4*bk/9+4*J/3-8*RV/3;H=self.H@Hc;Hy=self.Hy@Hc
    zet=4*ell-2*k/3;zeth=4*eh-2*kh/3;Bzet=self.B(zet,zeth);numer=self.H@Bzet-10*(self.P@zet)
    # k+y*k_y via the same weak derivative, without division by y.
    kd=(self.H@bk-6*kk)/4
    ar,Kt,diag=self.cmc(k,ell,L0,p,K,incoming,source) if prepared is None else prepared;diag=dict(diag);al=ar(y);ay=ar(y,1);ayy=ar(y,2);al0=float(ar(0));alR=float(ar(Y));
    if np.any(al<=0) or al0<=0 or alR<=0:raise ValueError('nonpositive reconstructed lapse')
    beta=(al-al0)/y;betaR=(alR-al0)/Y # quotient of smooth known spline; evaluate stably near 0 below
    # cubic segment Taylor for first interval (including extrapolation to y=0)
    near=y<=self.cent[0];beta[near]=ar(0,1)+y[near]*ar(0,2)/2+y[near]**2*ar(0,3)/6
    rho,anis,jhat,stress=source(y);ap=self.projp(al*pp);betap=self.projp(beta)
    HR=self.e@Hc;kR=self.trace(k);lR=self.trace(ell);LR0=L0+Y*lR;sR=np.exp(-LR0+Y*kR/2);ph=(w-sR*HR)/(2*R)
    kt=3*ap;et=ap/2-K*betap/3;L0t=-K*al0/3;kht=3*alR*ph;eht=alR*ph/2-K*betaR/3
    Jsrc=12*al*y*pp*V-4*y*ay*pp-2*KAP*al*jhat;Jt=self.projH(Jsrc)
    bkt=self.B(kt,kht);bet=self.B(et,eht);zt=Jt+bkt/3;Vt=bet/4
    Ht=-2*bet/3+4*bkt/9+4*Jt/3;G=self.B(ap,alR*ph);Gsrc=2*K/9*self.B(betap,betaR)+4*Jt/3
    zett=4*et-2*kt/3;assert np.linalg.norm(zett+4*K*betap/3)<1e-9
    # curvature remainder after extracting s^2*(Hc_y + D_B zeta).
    kmq=.5*kk*kk*hyp1f1(1,3,-y*kk)
    f=al*(-2/3*(kmq-4*V*V+2*V*kd)+2/9*kd*zc)-2/3*(4*ayy-2*(4*V-kd)*ay)+al*K*pp/s2-2*KAP/3*al*anis/s2
    Ms=self.P.T@((self.wp/s2)[:,None]*self.P)
    lp=self.P.T@(self.wp*al*Hy);lz=self.P.T@(self.wh*al*numer);lf=self.P.T@(self.wp*f);joint=self.joint(k,ell,J,RV,kh,eh,ar);pt=la.solve(Ms,joint+lf,assume_a='pos');joint_delta=joint-lp-lz
    # actual reconstructed metric tangent: projection defects must not disappear in Mdot.
    logs_t=-L0t-y*(self.P@et)+y*(self.P@kt)/2;logs_cont=al*K/3+al*y*pp
    masswork=-np.dot(self.wp*pp*pp/s2,logs_t);Md=self.P.T@((-2*self.wp*logs_t/s2)[:,None]*self.P)
    # outer physical reduction-based lower terms. Numerical p trace distinct.
    pR=self.trace(p);VR=self.e@Vh-1.5*lR;zR=self.e@z;aR=np.exp(-2*Y*kR/3);X2R=np.exp(Y*kR/3-2*LR0);cR=zR/aR
    kdR=kR+Y*self.trace(k,der=1);db=R*kdR/3;BC=(np.exp(-2*LR0)-sR*sR*(1+2*Y*VR)**2)/Y;lapTF=Y*2/3*sR*sR*(4*ar(Y,2)-2*(4*VR-kdR)*ar(Y,1))
    rhoR,anisR,jhatR,_=[float(v[0]) for v in source(np.array([Y]))]
    lower=-2*R*aR*sR*Y*pR*cR+2*BC+4*X2R*(1/R+db)*R*cR+2*lapTF/alR-2*Y*K*pR+2*K*K/9-(Y*pR)**2/2-2*KAP*rhoR/3+4*KAP*Y*anisR/3
    atlog=-2*Y*self.trace(kt)/3;ctR=(self.e@zt-atlog*zR)/aR;stR=sR*(-L0t-Y*self.trace(et)+Y*self.trace(kt)/2)
    ph_t=alR/(2*Y)*(g-R*aR*sR/alR*ctR-lower)
    wt=2*R*ph_t+sR*(self.e@Ht)+stR*HR
    rem=alR/R*(g-lower+3*HR)+sR*(self.e@Ht)-aR*sR*ctR+stR*HR
    # independently assembled total energy rate and ledger
    bulk_rate=p@Ms@pt+Hc@self.N@Ht+.5*p@Md@p
    bnd=R**5*alR*ph*HR;srcp=p@lf;srcz=p@lz;srch=Hc@self.N@Gsrc
    total=bulk_rate+R**5/6*w*wt;pred=srcp+srcz+srch+masswork-R**4*alR*sR/2*HR*HR+R**5/6*w*rem
    # projected continuous test gives exact weak distribution identity for zeta work.
    T=self.projH(al*pp);Tv=self.H@T; cell=np.dot(self.wp*Tv,self.Py@zet);jumps=0.
    for i in range(1,self.n):
        t=self.ends[i];jumps+=t**2.5*T[i]*(self.eval(zet,t,i=i)-self.eval(zet,t,i=i-1))
    outer=Y**2.5*T[self.n]*(zeth-self.trace(zet));projection=10*np.dot(self.wh*(Tv-al*pp),self.P@zet)
    # independently evaluate strong continuum physical shear, not used to set the weak load.
    ky=self.Py@k;ks=kk+y*ky;Vs=ll+y*(self.Py@ell);Vsy=2*(self.Py@ell)+y*np.concatenate([val(c['cp'],2*(c['y']-c['lo'])/(c['hi']-c['lo'])-1,2,c['hi']-c['lo'])@ell[c['pi']] for c in self.cells]);
    # Compare physical and reduction geometry separately: reduced V derivative within each cell.
    Vy=self.Hy@Vh-1.5*(self.Py@ell);zy=self.Hy@z;zety=self.Py@zet
    fs=al*(-2/3*(kmq-4*V*V+2*V*ks)+2/9*ks*zc)-2/3*(4*ayy-2*(4*V-ks)*ay)+al*K*pp/s2-2*KAP/3*al*anis/s2
    strong=s2*(al*(Hy+zety)+fs)
    # Physical constraint comparator derives C independently from the actual metric,
    # rather than reusing the weak metric connection residual.
    kyy=np.concatenate([val(c['cp'],2*(c['y']-c['lo'])/(c['hi']-c['lo'])-1,2,c['hi']-c['lo'])@k[c['pi']] for c in self.cells]);zphys=self.H@J+(10*kk+4*y*ky)/3;zphysy=self.Hy@J+(14*ky+4*y*kyy)/3
    strong_physical=s2*(al*(-2/3*(ky+kmq+4*Vsy-4*Vs*Vs+2*Vs*ks)+4/3*zphysy+2/9*ks*zphys)-2/3*(4*ayy-2*(4*Vs-ks)*ay))+al*K*pp-2*KAP/3*al*anis
    # physical metric constraint witnesses using cell-strong reconstructed metric, independent connection.
    Bcurv=-np.exp(-2*L)*kk*exprel(y*kk)-s2*(4*Vs+4*y*Vs**2)
    ric=-2/3*s2*(ky+kmq+4*Vsy-4*Vs*Vs+2*Vs*ks);Rscalar=6*Bcurv+6*y*ric
    cy=(zphysy+2/3*ks*zphys)/a; c=zphys/a;Zphys=X2*(3*c+2*y*cy);Hm=Rscalar+2*K*K/3-1.5*(y*pp)**2-2*KAP*rho;M= np.sqrt(y)*(5*pp+2*y*(self.Py@p)+6*y*pp*Vs-KAP*jhat)
    # independently reconstructed physical connection differs from the weak reduction-based metric connection.
    zphys=self.H@J+(10*kk+4*y*ky)/3
    # physical characteristic at actual outer face, from independent radial constraint formulas (reduction comparator separate).
    Vphys=lR+Y*self.trace(ell,der=1);Vphysy=2*self.trace(ell,der=1)+Y*self.trace(ell,der=2);kyR=self.trace(k,der=1);kmqR=.5*kR*kR*hyp1f1(1,3,-Y*kR)
    Bphys=(np.exp(-2*LR0)-sR*sR*(1+2*Y*Vphys)**2)/Y;ricR=-2/3*sR*sR*(kyR+kmqR+4*Vphysy-4*Vphys**2+2*Vphys*kdR);zphysR=self.trace(J,'H')+(10*kR+4*Y*kyR)/3;zphysyR=self.trace(J,'H',1)+(14*kyR+4*Y*self.trace(k,der=2))/3;Zr=X2R/aR*(3*zphysR+2*Y*(zphysyR+2/3*kdR*zphysR))
    Hmr=6*Bphys+6*Y*ricR+2*K*K/3-1.5*(Y*pR)**2-2*KAP*rhoR;Mr=R*(5*pR+2*Y*self.trace(p,der=1)+6*Y*pR*Vphys-KAP*jhatR)
    Wphys=2*sR*Mr+(Hmr+4*Zr)/3
    l2=lambda u,weight=self.wp:float(np.sqrt(np.dot(weight,u*u)))
    # origin quotient checks: raw rational-centre residual, not divided by y.
    origin0=self.eval(Bzet,0,'H')-10*self.eval(zet,0);origin1=self.eval(Bzet,0,'H',1)-14*self.eval(zet,0,der=1)
    # Actual shared interface work: local geometry load residual yields face coefficients.
    face=[]
    for c0 in self.cells:
        ix=c0['ix'];sl=c0['sl'];res=c0['H'].T@(c0['wh']*(self.H[sl]@G))+c0['Hy'].T@(c0['wp']*(self.P[sl]@ap));face.append(res)
    cell_errors=[]
    for i,c0 in enumerate(self.cells):
        sl=c0['sl'];ix=c0['ix'];pwork=np.dot(c0['wp']*al[sl]*pp[sl],Hy[sl]);hwork=np.dot(c0['wh']*H[sl],self.H[sl]@G);cell_errors.append(pwork+hwork-Hc[i]*face[i][0]-Hc[i+1]*face[i][1])
    jump=max(abs(face[i][1]+face[i+1][0]) for i in range(self.n-1));bubble=max(np.max(abs(fa[2:])) for fa in face)
    diag.update(principal_cell_error=float(max(abs(np.array(cell_errors)))),J_projection_work_error=float(Hc@self.N@Jt-np.dot(self.wh*H,Jsrc)),J_projection_residual=float(np.linalg.norm(self.N@Jt-self.H.T@(self.wh*Jsrc))),lapse_product_projection_error=float(np.linalg.norm(self.M@ap-self.P.T@(self.wp*al*pp))),principal_error=float(p@lp+Hc@self.N@G-bnd),Hc_kinematic=float(np.linalg.norm(Ht-G-Gsrc)),total_error=float(total-pred),total_rate=float(total),bulk_rate=float(bulk_rate),zeta_work=float(srcz),curvature_lapse_matter_work=float(srcp),geometry_source_work=float(srch),mass_work=float(masswork),boundary_remainder_work=float(R**5/6*w*rem),boundary_negative_term=float(-R**4*alR*sR*HR*HR/2),zeta_jump_work=float(jumps),zeta_outer_work=float(outer),zeta_test_projection_work=float(projection),zeta_pair_error=float(srcz-cell-jumps-outer-projection),mass_identity=float(.5*p@Md@p-masswork),mass_directional_error=float(abs(np.dot(self.wp*pp*pp/s2,np.sinh(-2*logs_t*1e-5))/(2e-5)-masswork)),mass_continuum_difference=l2(logs_t-logs_cont),pt_norm=l2(self.P@pt),strong_shear_defect=l2(self.P@pt-strong),physical_shear_defect=l2(self.P@pt-strong_physical),reduction_defect=l2(V-Vs),connection_rep_defect=l2(zc-zphys),H_m_norm=l2(Hm,self.wh),H_T_norm=l2(Hm+Zphys,self.wh),H_pin_norm=l2(Hm+X2*(c+2*y*cy),self.wh),momentum_norm=l2(M,self.wh),connection_norm=l2(c*np.sqrt(y),self.wh),physical_W=float(Wphys),temporal_boundary_error=float(R*aR*sR/alR*ctR+2*Y/alR*ph_t+lower-g),trace_delta=float(ph-pR),trace_delta_t=float(ph_t-self.trace(pt)),physical_flux_difference=float(R**5*alR*(ph-pR)*HR),k_trace_delta_t=float(kht-self.trace(kt)),q_trace_delta_t=float(np.exp(-Y*kh)*kht-np.exp(-Y*kR)*self.trace(kt)),ell_trace_delta_t=float(eht-self.trace(et)),zeta_origin_num=float(origin0),zeta_origin_jet_error=float(origin1),face_jump_error=float(jump),bubble_balance_error=float(bubble),outer_load_error=float(face[-1][1]-R**5*alR*ph),centre_load=float(face[0][0]),sparse_res=float(np.linalg.norm(self.N@G-(-self.D.T@(self.M@ap)+self.e*R**5*alR*ph))),sparse_monolithic=float(np.linalg.norm(G-la.solve(self.N,-self.D.T@(self.M@ap)+self.e*R**5*alR*ph,assume_a='pos'))))
    diag.update(joint_load_difference=float(np.linalg.norm(joint_delta)),joint_work=float(p@joint_delta),joint_adjusted_energy_closure=float(total-pred-p@joint_delta));return dict(p=pt,k=kt,ell=et,J=Jt,RV=np.zeros_like(RV),L0=L0t,K=Kt,kh=kht,eh=eht,w=wt,numerical_trace=ph,numerical_trace_t=ph_t,Hc=Hc,Hc_t=Ht,alpha=ar),diag
