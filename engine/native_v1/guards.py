"""Frozen release preflight and diagnostics, no time integrator."""
import pathlib, json, os, traceback
import numpy as np
from numpy.polynomial.legendre import leggauss
H_LIMIT=1e-8; M_LIMIT=1e-8; DEPART_LIMIT=1e-7; TICKS=8; DTICK=.000125

class Stop(Exception): pass

def require(ok,reason):
    if not ok: raise Stop(reason)

def clean(x):
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    if isinstance(x,np.ndarray):return [clean(v) for v in x] if x.ndim else clean(x.item())
    if isinstance(x,np.generic):return clean(x.item())
    if isinstance(x,float) and not np.isfinite(x):return str(x)
    return x

def write(p,x):
    # Preserve failed attempts too; atomic replace within a new per-run directory.
    tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(clean(x),indent=2,allow_nan=False)+'\n');os.replace(tmp,p)

def finite(x):
    if isinstance(x,dict):return all(finite(v) for v in x.values())
    if isinstance(x,(list,tuple)):return all(finite(v) for v in x)
    if isinstance(x,(float,int,np.number,np.ndarray)):return bool(np.all(np.isfinite(x)))
    return True

def maxabs(x):return float(np.max(np.abs(x)))
def weak_nodes(A):
    edges=np.sqrt(np.unique(np.r_[A.ends,A.cent]));t,w=leggauss(24)
    r=((edges[1:]+edges[:-1])[:,None]/2+np.diff(edges)[:,None]*t/2).ravel();wr=(np.diff(edges)[:,None]*w/2).ravel();return r*r,wr

def guard_physical(D,A,x,z,b):
    require(finite((x,z,b)),'nonfinite input')
    require(np.shape(z)==(A.lapse.dim+2,) and z[-1]==1,'invalid augmented packed gauge')
    S,raw,sc=x;mg=D.C.MomentGrid(A);reg=D.I.StructuralRegistry(A,S,raw,sc,mg)
    ys=np.unique(np.r_[A.y,A.cy,A.ends,A.cent,mg.q.ravel(),weak_nodes(A)[0]])
    require(np.all(np.isfinite(ys)),'nonfinite sample nodes')
    # Check raw branches before invoking source recovery. Include both one-sided metric face traces.
    minima={};dep={};inv={}
    def inspect(y,J,s2,label):
        require(finite((J,s2)) and np.all(J>0) and np.all(s2>0),'invalid metric '+label)
        for name,fs in reg.fluids.items():
            d,e,m=[v(y) for v in fs];require(finite((d,e,m)),'nonfinite reconstruction '+name)
            if name=='radiation':
                gap=e*e-s2*y*m*m;require(np.all(d==0) and np.all(e>0) and np.all(gap>0),'radiation branch/cone '+label)
                vals={'min_energy':float(e.min()),'min_cone_gap':float(gap.min())}
            else:
                er=np.sqrt(d*d+s2*y*m*m);require(np.all(d>0) and np.all(e>0) and np.all(er>0),'dust branch '+name+' '+label)
                vals={'min_mass':float(d.min()),'min_evolved_energy':float(e.min()),'min_recovered_energy':float(er.min())}
            for k,v in vals.items():minima[name+'_'+k]=min(v,minima.get(name+'_'+k,np.inf))
    J,s2,*_=D.C.geom(A,S,ys);inspect(ys,J,s2,'union')
    for i in range(A.n):
        y=A.ends[i:i+2];k=np.array([A.eval(S['k'],q,i=i) for q in y]);l=np.array([A.eval(S['ell'],q,i=i) for q in y]);L=S['L0']+y*l
        inspect(y,np.exp(3*L-y*k/2),np.exp(-2*L+y*k),'cell_face')
    k=np.array([A.eval(S['k'],q) for q in ys]);l=np.array([A.eval(S['ell'],q) for q in ys]);p=np.array([A.eval(S['p'],q) for q in ys]);L=S['L0']+ys*l;nu=3*L-ys*k/2;lam=-L+ys*k/2;L_nodes=L.copy()
    dep.update(K=abs(S['K']-b[1]),L0=abs(S['L0']-b[0]),nu=maxabs(nu-3*b[0]),lambda_=maxabs(lam+b[0]),k=maxabs(k),ell=maxabs(l),p=maxabs(p))
    for i in range(A.n):
        for q in A.ends[i:i+2]:
            k,l,p=(A.eval(S[name],q,i=i) for name in ('k','ell','p'));L=S['L0']+q*l
            for name,value in [('nu',3*L-q*k/2-3*b[0]),('lambda_',-L+q*k/2+b[0]),('k',k),('ell',l),('p',p)]:dep[name]=max(dep[name],abs(value))
    for i,field in enumerate(reg.fields):dep['scalar_'+str(i)]=maxabs(field(ys)-b[7+i])
    for name,target in [('dark_matter',(b[2],b[3],0.)),('baryons',(b[4],b[5],0.)),('radiation',(0.,b[6],0.))]:
        for j,field in enumerate(reg.fluids[name]):dep[name+'_'+str(j)]=maxabs(field(ys)-target[j])
    d,e,m=[v(ys) for v in reg.fluids['dark_matter']];db,eb,_=[v(ys) for v in reg.fluids['baryons']];rad=reg.fluids['radiation'][1](ys)
    phi=reg.fields[4](ys);u=np.exp(-2*D.BETA_D*phi);u0=np.exp(-2*D.BETA_D*.179055);B=.93*np.sqrt(1+3*u)/np.sqrt(1+3*u0)+.07*np.sqrt(u)/np.sqrt(u0)
    inv=dict(baryon_mass=db,baryon_energy=eb,radiation=rad*np.exp(L_nodes),dm_e_minus_d=e-d,dm_mass_over_B=d/B)
    # Initial invariants independent of reference time; archive pointwise signed changes on fixed nodes.
    inv0=D.invariants(D.background_initial())[:5];inv_delta={name:(v-inv0[i]).tolist() for i,(name,v) in enumerate(inv.items())}
    require(finite(dep) and finite(inv_delta),'nonfinite diagnostics')
    return dict(branch=minima,departures=dep,max_departure=max(dep.values()),nodes=ys.tolist(),invariant_delta=inv_delta)

def guard_rates(r,pre):
    require(finite(r),'nonfinite rate or ledger')
    require(r['immutable'],'native input mutation')
    require(r['cmc']['alpha_min']>0,'nonpositive lapse')
    require(r['geometry']['H_m_norm']<=H_LIMIT,'Hamiltonian norm limit')
    require(r['geometry']['momentum_norm']<=M_LIMIT,'momentum norm limit')
    require(pre['max_departure']<=DEPART_LIMIT,'homogeneous departure limit')

def attempt(D,A,x,z,t,b,path,kind):
    # Archive input prior to any check; record failed attempted state even when no rate exists.
    z=np.asarray(z); b=np.asarray(b); before=D.R.snapshot(*x);zb=z.tobytes();bb=b.tobytes()
    row=dict(kind=kind,t=t,state=x,input_digest=before,incoming_z=z,reference=b,geometry_boundary=0.,status='attempted')
    write(path,row)
    try:
        pre=guard_physical(D,A,x,z,b);row['preflight']=pre;write(path,row)
        require(pre['max_departure']<=DEPART_LIMIT,'homogeneous departure limit')
        S,raw,sc=x;zr,al,reg=D.refresh(A,x,z);J,s2,*_=D.C.geom(A,S,np.array([A.Y]));P=np.array([reg.fields[i](A.Y) for i in (1,3,5)]);flux=al(A.Y)*A.Y*J[0]*np.sqrt(s2[0])*(b[8::2]-P)
        row.update(scalar_flux=flux,reference_incoming=b[8::2],interior_scalar_P=P,pre_solved_z=zr,pre_solved_physical_lapse_flux=float(D.R.bR(A,S)*al(A.Y,1)),pre_solved_numerical_lapse_flux=float(al.flux));write(path,row)
        require(finite((zr,flux)) and al.minimum()>0,'invalid pre-solved lapse or boundary')
        zz=zr.tobytes(); zrr,alrr,_=D.refresh(A,x,zr)
        require(finite(zrr) and alrr.minimum()>0,'invalid fixed-state second refresh')
        require(zr.tobytes()==zz,'fixed-state refresh mutated input gauge')
        row.update(fixed_state_refreshed_z=zrr,fixed_state_refresh_signed=zrr-zr,fixed_state_refresh_max=maxabs(zrr-zr));write(path,row)
        r=D.R.rate(A,S,raw,sc,z,0.,flux)
        row.update(rate=r,actual_radiation_outer_flux=r['matter']['radiation']['outer_flux']);write(path,row)
        row['same_incoming_cmc_max']=maxabs(zr-np.array(r['outgoing_z']));write(path,row)
        require(np.array_equal(zr,np.array(r['outgoing_z'])),'pre-solve versus rate lapse mismatch')
        require(D.R.snapshot(*x)==before and z.tobytes()==zb and b.tobytes()==bb,'input mutation')
        guard_rates(r,pre);row['status']='valid';write(path,row);return r,np.array(r['outgoing_z'])
    except Exception as exc:
        row.update(status='failed',error=type(exc).__name__+': '+str(exc),traceback=traceback.format_exc());write(path,row);raise

def energy_data(D,A,x,r):
    S,raw,sc=x;mg=D.C.MomentGrid(A);reg=D.I.StructuralRegistry(A,S,raw,sc,mg);G=r['geometry_rates'];y=A.y;J,s2,*_=D.C.geom(A,S,y);p=np.asarray(S['p']);k=np.asarray(S['k']);l=np.asarray(S['ell']);lv=-G['L0']-y*(A.P@np.array(G['ell']))+y/2*(A.P@np.array(G['k']));Ms=A.P.T@((A.wp/s2)[:,None]*A.P);Md=A.P.T@((-2*A.wp/s2*lv)[:,None]*A.P)
    hc=lambda v:-2*A.B(np.array(v['ell']),v['eh'])/3+4*A.B(np.array(v['k']),v['kh'])/9+4*np.array(v['J'])/3-8*np.array(v['RV'])/3
    blocks={'geometry_p':(p,Ms,np.asarray(G['p']),Md),'geometry_Hc':(hc(S),A.N,hc(G),np.zeros_like(A.N)),'geometry_w':(np.atleast_1d(S['w']),np.array([[A.R**5/6]]),np.atleast_1d(G['w']),np.zeros((1,1)))}
    yy,wr=weak_nodes(A);J,s2,*_=reg.at(yy);B=reg.space.B(yy);H=2*np.sqrt(yy)[:,None]*reg.space.H(yy,1)@reg.space.E;W=yy*J
    kdot=np.array([A.eval(np.array(G['k']),q) for q in yy]);ldot=np.array([A.eval(np.array(G['ell']),q) for q in yy]);nv=3*G['L0']+3*yy*ldot-yy*kdot/2;lv=-G['L0']-yy*ldot+yy*kdot/2
    M=B.T@((wr*W)[:,None]*B);K=H.T@((wr*W*s2)[:,None]*H);Mdot=B.T@((wr*W*nv)[:,None]*B);Kdot=H.T@((wr*W*s2*(nv+2*lv))[:,None]*H)
    sr=np.asarray(r['weak']['scalar_rates']);c=(1.,1.,1/D.KAP)
    for i in range(3):
        blocks[f'scalar_{i}_P']=(sc[2*i+1],c[i]*M,sr[2*i+1],c[i]*Mdot);blocks[f'scalar_{i}_f']=(sc[2*i],c[i]*K,sr[2*i],c[i]*Kdot)
    f=np.array([reg.fields[i](yy) for i in (0,2,4)]);ft=np.array([reg.space.field(sr[i])(yy) for i in (0,2,4)]);U=D.potential(*f);Ud=np.sum(D.C.sh.up(f)*ft,axis=0)
    potential=(W,U,W*nv,Ud,wr)
    return blocks,potential

def dust_data(D,A,x,r):
    S,raw,sc=x;mg=D.C.MomentGrid(A);y=mg.q.ravel();J,s2,*_=D.C.geom(A,S,y);G=r['geometry_rates'];lv=-G['L0']-y*np.array([A.eval(np.array(G['ell']),q) for q in y])+y/2*np.array([A.eval(np.array(G['k']),q) for q in y]);out={}
    for name in ('dark_matter','baryons'):
        vals=[mg.reconstruction(v,j)(y) for v,j in zip(raw[name],(2,2,3))];rates=[mg.reconstruction(np.array(v),j)(y) for v,j in zip(r['matter_rates'][name],(2,2,3))];out[name]=(np.array(vals),y*s2,np.array(rates),lv)
    return out,mg

def step_ledger(D,A,x0,r0,xm,rm,xp,rp,h):
    a,pa=energy_data(D,A,x0,r0);m,pm=energy_data(D,A,xm,rm);b,pb=energy_data(D,A,xp,rp);out={}
    for key,(u0,M0,_,_) in a.items():
        um,Mm,rate,Mdot=m[key];u1,M1,_,_=b[key];du=u1-u0;ub=(u1+u0)/2;Mb=(M1+M0)/2;dM=M1-M0
        actual=(u1@M1@u1-u0@M0@u0)/2;stage=h*(um@Mm@rate);mixed=h*((ub@Mb-um@Mm)@rate);mass=ub@dM@ub/2;cross=du@dM@du/8;increment=ub@Mb@(du-h*rate);stage_mass=h*(um@Mdot@um)/2
        out[key]=dict(delta=actual,stage_work=stage,mixed_work=mixed,mass_endpoint=mass,coupled_increment=cross,update_roundoff_work=increment,stage_mass_work=stage_mass,metric_quadrature_defect=mass+cross-stage_mass,closure=actual-stage-mixed-mass-cross-increment)
    W0,U0,*_=pa;W1,U1,*_=pb;Wm,Um,Wdot,Udot,wr=pm;actual=wr@(W1*U1-W0*U0);stage=h*(wr@(Wdot*Um+Wm*Udot));p1=wr@((W1+W0)/2*(U1-U0));p2=wr@((U1+U0)/2*(W1-W0));out['scalar_potential']=dict(delta=actual,stage_work=stage,temporal_defect=actual-stage,potential_work=p1,volume_work=p2,closure=actual-p1-p2)
    d0,mg=dust_data(D,A,x0,r0);dm,_=dust_data(D,A,xm,rm);d1,_=dust_data(D,A,xp,rp);y=mg.q.ravel();shell=lambda v:(A.h*mg.average((y*v).reshape(mg.q.shape))).tolist();dust={}
    for name,(v0,q0,_,_) in d0.items():
        v1,q1,_,_=d1[name];vm,qm,rt,lv=dm[name];a0,e0,m0=v0;a1,e1,m1=v1;am,em,mm=vm;da,de,dmom=v1-v0;ab=(a1+a0)/2;mb=(m1+m0)/2;qb=(q1+q0)/2;er0=np.sqrt(a0*a0+q0*m0*m0);er1=np.sqrt(a1*a1+q1*m1*m1);erm=np.sqrt(am*am+qm*mm*mm)
        delta_rec=(2*ab*da+2*qb*mb*dmom+(q1-q0)*(mb*mb+dmom*dmom/4))/(er1+er0);actual=(e1-er1)-(e0-er0);stable=de-delta_rec;stage=h*(rt[1]-(am*rt[0]+qm*mm*rt[2]+qm*mm*mm*lv)/erm)
        dust[name]=dict(endpoint_defect=shell(actual),stable_endpoint_defect=shell(stable),endpoint_identity_closure=shell(actual-stable),stage_defect=shell(stage),temporal_defect=shell(stable-stage))
    return dict(quadratic=out,dust=dust)
