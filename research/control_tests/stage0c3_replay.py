"""Authorized Stage 0c.3, flat coarse OFF only. Research readback; no engine edits.
Original snapshot recovered from Stage 0; exact per-step observer vs plain
comparison. No projection of evolved state. Archived projected geometry map
is a reference, not a closed full-engine tangent. Stage increments telescope
but are not independent counterfactual causal contributions.
"""
import argparse, copy, json, pathlib, hashlib
import numpy as np
from scipy.linalg import eig, svd, svdvals, null_space, qr
from research.control_tests import ctl, observer, stage0
from engine import lattice_wall as wall
from engine import production_kernel as pk

SOURCE="cbf4e6bd759271b2ca15295db8f13e570e8fa3aa"
VENDOR="d6052d605673ce9d82cdc99b0df79855fa2ea215"
G=("a","b","X","Aa","K","Lambda","alpha")
DT=.000375; DR=.05; N=160; STEPS=5600; SZ=7*N
SCALE=np.repeat([1.,1.,1.,DR,DR,DR,1.],N)
LABEL="stage0_flat_off_1a_dr05_dt0_r8"
def read(p): return json.loads(p.read_text())
def write(p,x): p.write_text(json.dumps(x,indent=2,allow_nan=False))
def norm(x): return float(np.linalg.norm(x))
def mx(x): return float(np.max(np.abs(x)))
def ratio(x,y): return float(x/y) if y else None
def cpair(x): return [float(x.real),float(x.imag)]
def rms(x): return float(np.sqrt(np.mean(np.abs(x)**2)))
def pack(g): return np.concatenate([getattr(g,k) for k in G])
def snap_stage(st,g,sc,m):
    q=pk.ProductionState(grid=st.grid,geometry=g,scalars=sc,matter=m,
                         t=st.t,tau=st.tau,e_folds=st.e_folds)
    return stage0.snapshot(q)
def compare_runtime(p,q):
    keys=("head","submodules","python","numpy","scipy","threads")
    return {k:p[k]==q[k] for k in keys}
def load(ref,spec,out):
    old=read(ref/"off"/(LABEL+".json"))
    sp=read(spec/"off"/"report.json")
    assert old["stage0"]["provenance"]["head"]==SOURCE
    assert sp["provenance"]["head"]==SOURCE and VENDOR in sp["provenance"]["submodules"]
    assert (old["n"],old["dr"],old["dt"],old["steps_done"])==(N,DR,DT,STEPS)
    raw=np.load(ref/"preflight/flat_dr0.05_r8_off.npz",allow_pickle=False)
    raw={k:raw[k] for k in raw.files}
    assert stage0.hashes(raw)==old["stage0"]["initial_hashes"]==sp["initial_hashes"]
    k,s=stage0.build(DR,8.,"flat",False)
    before=stage0.snapshot(s); rebuilt=stage0.hashes(before)
    mismatch=[f for f in raw if rebuilt[f]!=old["stage0"]["initial_hashes"][f]]
    # Restore ONLY initial archived arrays, not evolved state. Grid/metadata
    # must match constructor. No charged dust/EM in original control.
    assert np.array_equal(before["r"],raw["r"])
    assert s.electromagnetic is None and s.matter.charged_dust is None
    assert np.array_equal(before["clocks"],raw["clocks"])
    for f in ctl.VAC.ALL_FIELDS: setattr(s.geometry,f,raw["g_"+f].copy())
    for f in ctl.FIELDS_S: getattr(s.scalars,f)[:]=raw["s_"+f]
    for spc in ("dark_matter","baryons","radiation"):
        for f in ("rest","energy_t","momentum_r"):
            getattr(getattr(s.matter,spc),f)[:]=raw["m_"+spc+"_"+f]
    assert stage0.hashes(stage0.snapshot(s))==old["stage0"]["initial_hashes"]
    p=stage0.provenance(); assert p["head"]==SOURCE and not p["tracked_status"]
    checks=compare_runtime(p,old["stage0"]["provenance"]); assert all(checks.values()),checks
    assert float(ctl.ss.V0)==old["stage0"]["V0"] and float(ctl.ss.SD_QUARTIC_LAMBDA)==old["stage0"]["SD_QUARTIC_LAMBDA"]
    arrays=np.load(spec/"off/spectral_arrays.npz",allow_pickle=False)
    assert np.array_equal(pack(s.geometry)*SCALE,arrays["baseline_input"])
    report={"scientific_admission":False,"source":SOURCE,"vendor":VENDOR,"provenance":p,
      "archive_provenance":old["stage0"]["provenance"],"pinned_runtime_checks":checks,
      "runtime_text_equal":p["numpy_runtime"]==old["stage0"]["provenance"]["numpy_runtime"],
      "reconstructed_initial_hash_mismatches_before_archive_restore":mismatch,
      "archived_initial_state_verified":True,"scope":"flat OFF dr=.05 R=8 dt=.000375 t=2.1; frozen geometry block; no physical energy attribution"}
    write(out/"report.json",report)
    return k,s,old,arrays,report
def flat_constraints(s):
    H,M=ctl.constraints(s); raw=ctl.VAC.constraints(s.grid,s.geometry)
    return np.stack([H,M,raw["connection"],raw["determinant"]])
def pert(s,z,h):
    q=copy.deepcopy(s)
    for j,f in enumerate(G): setattr(q.geometry,f,getattr(s.geometry,f)+h*z[j*N:(j+1)*N]/SCALE[j*N])
    return q
def tangent(s):
    # Exact flat vacuum linearization of pinned vacuum_bssn constraints.
    # Frozen matter background is O(1e-150); separately check coupled
    # constraint finite differences in geometry directions, no matter columns.
    grid=s.grid; r=grid.centers; eye=np.eye(N); zero=np.zeros_like(eye)
    D=np.column_stack([grid.cell_derivative_fourth(v,parity=1) for v in eye])
    O=np.column_stack([grid.cell_derivative_fourth(v,parity=-1) for v in eye])
    # Use exact vendor derivative helper, not a new stencil.
    D2=np.column_stack([ctl.VAC._d2(grid,v,1) for v in eye])
    R=np.diag(1/r); R2=R@R
    H=np.hstack([-.5*D2+4*R2,-D2-4*R@D-4*R2,4*D2+8*R@D,zero,zero,O,zero])
    M=np.hstack([zero,zero,zero,D+3*R,-2*D/3,zero,zero])
    C=np.hstack([-.5*D-2*R,D+2*R,zero,zero,zero,eye,zero])
    det=np.hstack([eye,2*eye,zero,zero,zero,zero,zero])
    return {k:v/SCALE[None,:] for k,v in (("hamiltonian",H),("momentum",M),("connection",C),("determinant",det))}
def mode_data(P):
    A=(P-np.eye(SZ))/DT
    lam,L,V=eig(A,left=True,right=True)
    mu=1+DT*lam
    indices=[int(np.argmin(np.where(np.angle(mu)>0,abs(np.angle(mu)/DT-w),np.inf))) for w in (35.89750577213832,44.8718822151729,20.94021170041402)]
    V=V[:,indices]; L=L[:,indices]; mu=mu[indices]
    overlap=np.sum(L.conj()*V,axis=0); assert np.all(abs(overlap)>0)
    W=L/np.conj(overlap)[None,:]
    return mu,V,W
def restriction(k,s,arr,report,out):
    P=arr["step_1"]; P0=arr["step_0"]; A=arr["lapse_1"]
    J=tangent(s); E=np.zeros((N,SZ)); E[:,-N:]=np.eye(N)
    J["incoming_lapse_fixed_point"]=E-A
    J["incoming_lapse_outer_normalization"]=E[-1:,:]
    # A fixed-point-compatible incoming lapse is an implicit tangent
    # condition, not an inverse of I-Aalpha; null direction retained.
    # Ab=-Aa/2 is algebraic in this engine so trace constraint is identically zero.
    blocks={"incoming_lapse":np.vstack([J["incoming_lapse_fixed_point"],J["incoming_lapse_outer_normalization"]]),
       "determinant":J["determinant"],"connection":J["connection"],
       "hamiltonian":J["hamiltonian"],"momentum":J["momentum"],
       "hamiltonian_momentum":np.vstack([J["hamiltonian"],J["momentum"]]),
       "all_geometry_constraints":np.vstack([J[x] for x in ("determinant","connection","hamiltonian","momentum")]),
       "trace_algebraic_identity":np.zeros((1,SZ))}
    # Row equilibration for numerical null-space rank only, not dynamics.
    # Explicit tolerance sensitivity, not an admission gate.
    def basis(M,tol):
        z=np.linalg.norm(M,axis=1); M=M[z>0]/z[z>0,None]
        if not len(M): return np.eye(SZ)
        return null_space(M,rcond=tol)
    # Smooth analytic even/odd functions of y, with correct centre orders.
    r=s.grid.centers; x=r/r[-1]; y=2*x*x-1
    def smooth(degree):
        T=np.polynomial.chebyshev.chebvander(y,degree); cols=[]
        for f in ("isotropic_metric","anisotropy","X","Aa","K","Lambda","alpha"):
            for v in T.T:
                z=np.zeros((7,N))
                if f=="isotropic_metric": z[0]=v; z[1]=v
                elif f=="anisotropy": z[0]=x*x*v; z[1]=-.5*x*x*v
                elif f=="Aa": z[3]=x*x*v
                elif f=="Lambda": z[5]=x*v
                elif f=="alpha": z[6]=(1-x*x)*v
                else: z[G.index(f)]=v
                cols.append(z.ravel()*SCALE)
        return qr(np.array(cols).T,mode="economic")[0]
    bases={}
    for key,M in blocks.items():
        for tol in (1e-8,1e-6):
            Z=basis(M,tol); bases[f"{key}_rcond{tol:g}"]=Z
    for d in (4,8): bases[f"smooth_regular_degree{d}"]=smooth(d)
    # Combined smooth + determinant and connection, not necessarily invariant.
    Q=bases["smooth_regular_degree8"]; B=np.vstack([J["determinant"],J["connection"]])
    bases["smooth_and_det_connection"]=Q@basis_small(B@Q)
    def comp(v):
        w=v.reshape(7,N)
        return {"field_fraction":{f:float(np.sum(abs(w[j])**2)) for j,f in enumerate(G)},
          "centre3_fraction":float(np.sum(abs(w[:,:3])**2)),
          "outer3_fraction":float(np.sum(abs(w[:,-3:])**2)),
          "constraint_tangent_l2":{f:norm(M@v) for f,M in J.items()}}
    summaries={}; saved={"r":r,"scale":SCALE,**{"constraint_"+f:M for f,M in J.items()}}
    rng=np.random.default_rng(613); witnesses={}
    # Check exact analytic flat constraints vs actual production evaluator
    # at deterministic directions, not claim exact full coupled tangent.
    for tag,z in [("random",rng.normal(size=SZ)),("smooth",smooth(4)[:,2])]:
        z=z/np.linalg.norm(z); vals=[]
        for h in (1e-4,3e-5):
            diff=(flat_constraints(pert(s,z,h))-flat_constraints(pert(s,z,-h)))/(2*h)
            vals.append({f:{"absolute_l2":norm(diff[i]-J[f]@z),"relative_l2":ratio(norm(diff[i]-J[f]@z),norm(diff[i]))}
               for i,f in enumerate(("hamiltonian","momentum","connection","determinant"))})
        witnesses[tag]=vals
    for m in (1,80,800,1600,3200,STEPS):
        Q=np.linalg.matrix_power(P,m); U,ss,Vh=svd(Q,full_matrices=False); u=U[:,0]; v=Vh[0].conj()
        row={"t":m*DT,"gain":float(ss[0]),"input":comp(v),"output":comp(u),
             "restricted":{},"second_singular_value":float(ss[1])}
        for name,Z in bases.items():
            row["restricted"][name]={"dimension":Z.shape[1],"gain":float(svdvals(Q@Z)[0]) if Z.shape[1] else 0.0,
              "dominant_input_projection":norm(Z.T@v)}
        summaries[str(m)]=row
        saved[f"left_{m}"]=u; saved[f"right_{m}"]=v
        if m in (1,STEPS):
            row["other_epsilon_gain"]=float(svdvals(np.linalg.matrix_power(P0,m))[0])
    # Immediate one-step response along the leading right singular vector,
    # actual observer operation order, +/- perturbations. Also independently
    # verify flat quadratic lapse null shape at last-centre normalization.
    tests=[]; v=saved["right_1"]
    def trace(st):
        rec=observer.Recorder(); q=copy.deepcopy(st); new,log=rec.recorded_step(q,DT)
        seq=[("start",pack(st.geometry)*SCALE)]
        for i,e in enumerate(log):
            g=e[1]; seq.append((f"{i}_{e[0]}_input",pack(g)*SCALE))
            if e[0]=="solve":
                h=g.copy(); h.alpha=e[4].copy(); seq.append((f"{i}_solve_output",pack(h)*SCALE))
        seq.append(("accepted",pack(new.geometry)*SCALE)); return seq
    for h in (1e-4,3e-5):
        plus=trace(pert(s,v,h)); minus=trace(pert(s,v,-h)); assert [x[0] for x in plus]==[x[0] for x in minus]
        rows=[]
        for (name,p),(_,m) in zip(plus,minus):
            z=(p-m)/(2*h)
            rows.append({"stage":name,"scaled_gain":norm(z),"field_l2":{f:norm(z[j*N:(j+1)*N]) for j,f in enumerate(G)}})
        tests.append({"epsilon":h,"stage_gains":rows,"accepted_tangent_difference_l2":norm((plus[-1][1]-minus[-1][1])/(2*h)-P@v)})
    quadratic=[]; shape=(r[-1]**2-r**2)/6
    baseg=s.geometry.copy()
    a0,kd0=k._solve_lapse_for_state(s.grid,baseg,s.scalars,s.matter,s.electromagnetic)
    for h in (1e-4,3e-5,1e-6,1e-8,1e-10,1e-12):
        g=baseg.copy(); g.alpha=baseg.alpha+h*shape
        # actual incoming K RHS projection, band solve, then production L2.
        before=ctl.VAC.primary_l2_rhs(s.grid,g)
        alpha,kd=k._solve_lapse_for_state(s.grid,g,s.scalars,s.matter,s.electromagnetic)
        g.alpha=alpha
        after=k._geometry_stage_terms(s.grid,g,s.scalars,s.matter,s.electromagnetic)["primary_l2"]
        g0=baseg.copy(); g0.alpha=a0
        a_base=k._geometry_stage_terms(s.grid,g0,s.scalars,s.matter,s.electromagnetic)["primary_l2"]
        # The continuum zero is the shear RHS, not a tolerance waiver.
        quadratic.append({"epsilon":h,"kdot_response_over_epsilon":float((kd-kd0)/h),
          "solved_shape_error_over_epsilon":mx((alpha-a0)/h-shape),
          "incoming_Aa_rhs_absmax":mx(before["Aa"]),
          "outgoing_Aa_rhs_delta_absmax":mx(after["Aa"]-a_base["Aa"]),
          "outgoing_regular_shear_delta_absmax":mx((after["Aa"]-a_base["Aa"])/r**2),
          "outgoing_Aa_rhs_delta_cell0":float(after["Aa"][0]-a_base["Aa"][0])})
        saved[f"quadratic_alpha_{h:g}"]=alpha
        saved[f"quadratic_Aa_delta_{h:g}"]=after["Aa"]-a_base["Aa"]
    report.update({"status":"completed_restriction_diagnostic","gain":summaries,
       "tangent_checks":witnesses,"immediate_stage_direction":tests,"quadratic_lapse_null":quadratic,
       "restriction_scope":"geometry-only frozen flat vacuum tangent, matter held in input; no invariance claim; rank tolerance sensitivity retained; trace algebraic Ab=-Aa/2",
       "physical_admissibility_proven":False})
    np.savez_compressed(out/"singular_vectors.npz",**saved)
    write(out/"report.json",report)
    print(json.dumps({"gain":summaries,"quadratic_lapse_null":quadratic,"immediate_stage_direction":tests,"tangent_checks":witnesses},allow_nan=False),flush=True)
def basis_small(M):
    z=np.linalg.norm(M,axis=1); M=M[z>0]/z[z>0,None]
    return null_space(M,rcond=1e-8) if len(M) else np.eye(M.shape[1])
def replay(k,s,old,arr,report,out):
    P=arr["step_1"]; base=pack(s.geometry)
    b=arr["baseline_step"]-arr["baseline_input"]
    mu,V,W=mode_data(P)
    plain=copy.deepcopy(s); rec_state=copy.deepcopy(s); rec=observer.Recorder()
    report["status"]="replaying"; write(out/"report.json",report)
    # Store unscaled accepted arrays and exact per-step stage inputs/outputs.
    # Chunked outputs prevent RAM accumulation. No extra diagnostics feed back.
    chunks=[]; stage_chunks=[]; terms_chunks=[]; constraints=[]; hashes=[]; checks=[]
    states=[]; stage_geom=[]; archived_centre_mismatch=0; archived_save_mismatch=[]
    stage_labels=None; initial_hash=stage0.hashes(stage0.snapshot(s))
    archived_saves={r["step"]:r["state_hashes"] for r in old["saves"]}
    stage_times=[]
    def flush(i):
        if not chunks: return
        def stack_dicts(seq):
            return {k:np.stack([d[k] for d in seq]) for k in seq[0]}
        np.savez_compressed(out/f"state_{i:05d}.npz",**stack_dicts(chunks),constraints=np.stack(constraints))
        if stage_chunks:
            np.savez_compressed(out/f"stages_{i:05d}.npz",
              **stack_dicts(stage_chunks),**stack_dicts(terms_chunks))
        chunks.clear(); stage_chunks.clear(); terms_chunks.clear(); constraints.clear()
    for i in range(STEPS+1):
        q=stage0.snapshot(rec_state); h=stage0.hashes(q)
        hp=stage0.hashes(stage0.snapshot(plain))
        assert h==hp,("observer mismatch",i,[f for f in h if h[f]!=hp[f]])
        hashes.append(h); chunks.append(q); constraints.append(flat_constraints(rec_state))
        states.append((pack(rec_state.geometry)-base)*SCALE)
        oldrow=old["centre_trace"][i]; g=rec_state.geometry
        now=[(g.a[:3]-g.b[:3]).tolist(),g.Aa[:3].tolist(),g.X[:3].tolist(),g.Lambda[:3].tolist()]
        if now!=[oldrow[f] for f in ("a_minus_b","Aa","X","Lambda")]: archived_centre_mismatch+=1
        if i in archived_saves and h!=archived_saves[i]:
            archived_save_mismatch.append({"step":i,"fields":[f for f in h if h[f]!=archived_saves[i][f]]})
        if i==STEPS:
            flush(i); break
        new,log=rec.recorded_step(rec_state,DT)
        plain=k.step(plain,DT)
        sd={}; td={}; seq=[pack(rec_state.geometry)]
        names=["start"]
        for j,e in enumerate(log):
            name=f"{j:02d}_{e[0]}"; snap=snap_stage(rec_state,e[1],e[2],e[3])
            for f,v in snap.items():
                if f not in ("r","clocks"): sd[name+"_"+f]=v
            seq.append(pack(e[1])); names.append(name+"_input")
            if e[0]=="solve":
                td[name+"_alpha"]=e[4]; td[name+"_kdot"]=np.array(e[5])
                g=e[1].copy(); g.alpha=e[4].copy()
                seq.append(pack(g)); names.append(name+"_output")
            elif e[0]=="terms":
                for f,v in e[4].items():
                    if isinstance(v,dict):
                        for l,w in v.items(): td[name+"_"+f+"_"+l]=w
                    else: td[name+"_"+f]=v
        # Reconstruct exact final L2 operands following production copy order;
        # not claimed as an independent intercepted call.
        st,_=observer.stages(rec_state,new,log); fc=st["final_l2_call"][0]
        for f in ctl.VAC.ALL_FIELDS: sd["reconstructed_final_l2_g_"+f]=np.array(getattr(fc,f))
        final_l2=ctl.VAC.primary_l2_rhs(wall.bind(rec_state.grid,fc),fc)
        for f,v in final_l2.items(): td["reconstructed_final_l2_"+f]=np.asarray(v)
        final_primary=fc.copy()
        final_primary.Aa=log[8][1].Aa.copy(); final_primary.K=log[8][1].K.copy()
        for f in ctl.VAC.ALL_FIELDS: sd["reconstructed_final_lambda_l2_g_"+f]=np.array(getattr(final_primary,f))
        td["reconstructed_final_lambda_l2"]=ctl.VAC.lambda_l2_rhs(wall.bind(rec_state.grid,final_primary),final_primary,lambda_m=pk.LAMBDA_M)
        seq.append(pack(new.geometry)); names.append("accepted")
        if stage_labels is None: stage_labels=names
        assert names==stage_labels
        stage_chunks.append(sd); terms_chunks.append(td)
        # Differences taken in original physical coordinates before scaling.
        stage_geom.append(np.diff(np.array(seq),axis=0)*SCALE)
        rec_state=new
        if (i+1)%100==0:
            flush(i); print("replay steps",i+1,"archived centre differences",archived_centre_mismatch,flush=True)
    X=np.array(states); inc=np.array(stage_geom)
    # Ledger: subtract before scaling for accepted increment; avoid adding
    # nearly equal O(1) states. Reference b uses the archived scaled baseline.
    I=np.diff(X,axis=0); A=P-np.eye(SZ)
    predicted=X[:-1]@A.T
    residual=I-predicted-b
    z=X@W.conj(); f=residual@W.conj(); bj=b@W.conj()
    force=f+bj
    intrinsic=(abs(mu)**2-1)*abs(z[:-1])**2
    cross=2*np.real(np.conj(z[:-1]*mu)*force); squared=abs(force)**2
    closure=abs(z[1:])**2-abs(z[:-1])**2-intrinsic-cross-squared
    # Ordered stage changes sum to actual change; relative to first-step
    # stage changes, with explicit interaction remainder. Not counterfactual
    # independent source attribution. Report signs, not physical energy.
    modal_inc=inc@W.conj(); telescope=inc.sum(axis=1)-I
    stage_departures=modal_inc-modal_inc[:1]
    stage_interaction=f-(stage_departures.sum(axis=1)-predicted@W.conj())
    def coherence(v,j):
        demod=v*np.exp(-1j*np.angle(mu[j])*np.arange(len(v)))
        return ratio(abs(np.sum(demod)),np.sum(abs(demod)))
    def contribution(v,j):
        # exact forced recurrence from zero initial amplitude.
        out=np.zeros(len(v)+1,complex)
        for i,x in enumerate(v): out[i+1]=mu[j]*out[i]+x
        return out
    modes=[]
    for j in range(len(mu)):
        modes.append({"mu":cpair(mu[j]),"omega":float(np.angle(mu[j])/DT),
          "growth":float(np.log(abs(mu[j]))/DT),"modal_absmax":mx(z[:,j]),
          "constant_drift_abs":float(abs(bj[j])),"remainder_rms":rms(f[:,j]),
          "remainder_demodulated_coherence":coherence(f[:,j],j),
          "total_forcing_demodulated_coherence":coherence(force[:,j],j),
          "signed_amplitude_budget_sums":{"intrinsic":float(intrinsic[:,j].sum()),"cross":float(cross[:,j].sum()),"excitation_square":float(squared[:,j].sum()),"closure_absmax":mx(closure[:,j])},
          "drift_response_absmax":mx(contribution(np.full(STEPS,bj[j]),j)),
          "remainder_response_absmax":mx(contribution(f[:,j],j)),
          "ordered_stage_departure_response_absmax":{name:mx(contribution(stage_departures[:,i,j],j)) for i,name in enumerate(stage_labels[1:])}})
    # Full affine model vs exact replay. Decompose remainder response by
    # actual signed vector forcing (not treating all omissions as excitatory).
    affine=np.zeros_like(X); forced=np.zeros_like(X)
    for i in range(STEPS):
        affine[i+1]=P@affine[i]+b; forced[i+1]=P@forced[i]+residual[i]
    C=np.zeros((12,SZ))
    for j in range(3):
        C[j,j]=1; C[j,N+j]=-1; C[3+j,3*N+j]=1/DR
        C[6+j,2*N+j]=1; C[9+j,5*N+j]=1/DR
    obs=X@C.T; pred=affine@C.T; rem=forced@C.T
    obs_report={}
    for j,field in enumerate(("a_minus_b","Aa","X","Lambda")):
        sl=slice(3*j,3*j+3); o=obs[:,sl]; p=pred[:,sl]; q=rem[:,sl]
        obs_report[field]={"observed_rms":rms(o),"affine_rms":rms(p),"remainder_response_rms":rms(q),
          "signed_cross_mean_square":float(2*np.mean(p*q)),"closure_rms":rms(o-p-q),
          "affine_residual_over_observed":ratio(rms(o-p),rms(o))}
    # Alignment with actual departures and residuals, unrestricted geometry
    # norm; do not label these inputs physically admissible.
    alignment={}; singular={}
    for m in (1,800,STEPS):
        U,ss,Vh=svd(np.linalg.matrix_power(P,m),full_matrices=False)
        v=Vh[0]; u=U[:,0]
        dnorm=np.linalg.norm(X,axis=1); fnorm=np.linalg.norm(residual,axis=1)
        ax=np.divide(abs(X@v),dnorm,out=np.zeros(len(X)),where=dnorm>0)
        af=np.divide(abs(residual@v),fnorm,out=np.zeros(len(residual)),where=fnorm>0)
        alignment[str(m)]={"gain":float(ss[0]),"departure_cosine_max":mx(ax),"remainder_cosine_max":mx(af),
          "departure_cosine_rms":rms(ax),"remainder_cosine_rms":rms(af)}
        singular[f"right_{m}"]=v; singular[f"left_{m}"]=u
        singular[f"departure_alignment_{m}"]=ax; singular[f"remainder_alignment_{m}"]=af
    np.savez_compressed(out/"modal_ledger.npz",X=X,increment=I,baseline=b,predicted_increment=predicted,
      remainder=residual,modal=z,modal_remainder=f,modal_constant=bj,modal_force=force,
      intrinsic=intrinsic,cross=cross,excitation_squared=squared,closure=closure,
      ordered_stage_increments=inc,ordered_stage_modal=modal_inc,ordered_stage_departures=stage_departures,
      stage_interaction_remainder=stage_interaction,affine=affine,forced=forced,**singular)
    write(out/"state_hashes.json",hashes)
    report.update({"status":"completed_replay_diagnostic","steps":STEPS,
      "observer_off_on_digest_equal_every_step":True,"archived_centre_mismatch_count":archived_centre_mismatch,
      "archived_save_hash_mismatches":archived_save_mismatch,
      "exact_archived_trajectory_reproduced":not archived_centre_mismatch and not archived_save_mismatch,
      "stage_labels":stage_labels,"stage_increment_telescope_absmax":mx(telescope),
      "stage_interaction_remainder_absmax":mx(stage_interaction),
      "frozen_modal_budget":modes,"observable_comparison":obs_report,"singular_alignment":alignment,
      "initial_drift_vs_first_actual_increment_inf":mx(I[0]-b),
      "baseline_plus_remainder_closure_inf":mx(X-affine-forced),
      "stage_scope":"actual recorded operands, reconstructed final L2; signed ordered differences not independent source decomposition",
      "native_attribution":"held"})
    write(out/"report.json",report)
    print(json.dumps({k:report[k] for k in ("status","steps","observer_off_on_digest_equal_every_step","archived_centre_mismatch_count","archived_save_hash_mismatches","exact_archived_trajectory_reproduced","stage_increment_telescope_absmax","stage_interaction_remainder_absmax","frozen_modal_budget","observable_comparison","singular_alignment")},allow_nan=False),flush=True)
def main():
    p=argparse.ArgumentParser(); p.add_argument("task",choices=["restrict","replay"])
    p.add_argument("stage0",type=pathlib.Path); p.add_argument("stage0b",type=pathlib.Path); p.add_argument("out",type=pathlib.Path)
    a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    k,s,old,arr,report=load(a.stage0,a.stage0b,a.out)
    if a.task=="restrict": restriction(k,s,arr,report,a.out)
    else: replay(k,s,old,arr,report,a.out)
if __name__=="__main__": main()
