"""Stage 0c.1/.2 ONLY. Archive-only diagnostics; no engine imports, probes,
new trajectories or acceptance gates. Frozen projected geometry map, not a
closed full-engine tangent. The affine linear response is not production replay.
"""
import argparse, hashlib, json, pathlib
import numpy as np
from scipy.linalg import eig, schur, svdvals
from scipy.signal import detrend

SOURCE="cbf4e6bd759271b2ca15295db8f13e570e8fa3aa"
VENDOR="d6052d605673ce9d82cdc99b0df79855fa2ea215"
G=("a","b","X","Aa","K","Lambda","alpha")
def digest(v):
    return hashlib.sha256(str(v.dtype).encode()+str(v.shape).encode()+np.ascontiguousarray(v).tobytes()).hexdigest()
def read(p): return json.loads(p.read_text())
def write(p,x): p.write_text(json.dumps(x,indent=2,allow_nan=False))
def mx(v): return float(np.max(np.abs(v)))
def norm(v): return float(np.linalg.norm(v))
def ratio(a,b): return float(a/b) if b else None
def cpair(z): return [float(z.real),float(z.imag)]
def rms(v): return float(np.sqrt(np.mean(np.abs(v)**2)))
def provenance(p):
    assert p["head"]==SOURCE and not p["tracked_status"]
    assert VENDOR in p["submodules"]
def main():
    pa=argparse.ArgumentParser()
    pa.add_argument("mode",choices=["flat","native"])
    pa.add_argument("switch",choices=["off","on"])
    pa.add_argument("stage0",type=pathlib.Path)
    pa.add_argument("stage0b",type=pathlib.Path)
    pa.add_argument("out",type=pathlib.Path)
    a=pa.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    label=f"stage0_{a.mode}_{a.switch}_1a_dr05_dt0_r8"
    old=read(a.stage0/a.switch/(label+".json"))
    spec=read(a.stage0b/a.switch/"report.json")
    pre=read(a.stage0/"preflight"/"stage0_preflight.json")
    provenance(old["stage0"]["provenance"]); provenance(spec["provenance"]); provenance(pre["provenance"])
    assert spec["archive_provenance"]==old["stage0"]["provenance"]
    assert spec["status"]=="completed_diagnostic"
    n=old["n"]; dr=old["dr"]; dt=old["dt"]; N=7*n
    assert (n,dr,dt,old["r_max"])==(160,.05,.000375,8.)
    assert tuple(spec["coordinates"])==G
    scale=np.repeat(spec["scale_by_field"],n)
    assert np.array_equal(scale,np.repeat([1.,1.,1.,dr,dr,dr,1.],n))
    target=old["stage0"]["initial_hashes"]
    raw=np.load(a.stage0/"preflight"/f"{a.mode}_dr0.05_r8_{a.switch}.npz",allow_pickle=False)
    raw_hash={k:digest(raw[k]) for k in raw.files}
    assert raw_hash==target
    assert pre["pairs"][f"{a.mode}_dr0.05_r8"]["arms"][int(a.switch=="on")]["hashes"]==target
    snap=np.load(a.stage0/a.switch/(label+".npz"),allow_pickle=False)
    # Check every archived geometry and scalar profile against its exact state digest.
    saves=old["saves"]; steps=np.array([r["step"] for r in saves],dtype=int)
    assert steps[0]==0 and steps[-1]==old["steps_done"] and np.all(np.diff(steps)==800)
    for i,rec in enumerate(saves):
        for k in target:
            if k.startswith(("g_","s_")):
                assert digest(snap[f"t{i}_{k}"])==rec["state_hashes"][k],(i,k)
    for k in target:
        if k.startswith(("g_","s_")): assert digest(snap["t0_"+k])==target[k]
    arr=np.load(a.stage0b/a.switch/"spectral_arrays.npz",allow_pickle=False)
    base=arr["baseline_input"]; b=arr["baseline_step"]-base
    raw_scaled=np.concatenate([raw["g_"+k] for k in G])*scale
    mismatch=[k for k in target if target[k]!=spec["initial_hashes"][k]]
    assert mismatch==spec["initial_hash_mismatches"]
    exact=not mismatch and np.array_equal(raw_scaled,base)
    result={"scientific_admission":False,"scope":"archive-only frozen projected geometry; no engine replay or physical energy attribution",
        "mode":a.mode,"switch":a.switch,"coordinates":G,"scale":scale[::n].tolist(),
        "stage0_source":SOURCE,"stage0_vendor":VENDOR,"archived_full_initial_state_hashes_verified":True,
        "saved_geometry_scalar_hashes_verified":True,"stage0b_initial_hash_mismatches":mismatch,
        "stage0b_scaled_geometry_diff_by_field":{k:mx(raw_scaled[j*n:(j+1)*n]-base[j*n:(j+1)*n]) for j,k in enumerate(G)},
        "trajectory_projection_eligible":bool(exact),
        "native_mismatch_policy":"archive initial arrays recovered exactly, but no new tangent at those operands; unmatched trajectory projection held",
        "profile_sampling_dt":.3,"profile_angular_nyquist":float(np.pi/.3),
        "baseline_drift_inf":mx(b),"baseline_drift_l2_scaled":norm(b),
        "shared_diagnostics":{"metric_anisotropy_shear":"modal observation",
          "connection_residual":"not evaluable from archived spectra without operator crosswalk; Lambda trace is not connection residual",
          "momentum_residual":"archived M profiles only","CMC_lapse":"stored alpha and solved-lapse tangent",
          "matter_conservation":"no every-step conserved state archive; no attribution",
          "centre_derivatives":"no derivative modal attribution without crosswalk"}}
    write(a.out/"report.json",result)
    # 3 centre cells each: physical a-b,Aa,X,Lambda. Scale-free relative
    # comparisons are reported by field; the state norm is scaled Euclidean.
    C=np.zeros((12,N))
    for j in range(3):
        C[j,j]=1; C[j,n+j]=-1
        C[3+j,3*n+j]=1/dr
        C[6+j,2*n+j]=1
        C[9+j,5*n+j]=1/dr
    tr=old["centre_trace"]; ts=np.array([x["t"] for x in tr])
    obs=np.array([sum([v[f] for f in ("a_minus_b","Aa","X","Lambda")],[]) for v in tr])
    assert len(ts)==steps[-1]+1 and np.array_equal(np.array([v["step"] for v in tr]),np.arange(len(ts)))
    assert np.allclose(ts,np.arange(len(ts))*dt,rtol=0,atol=1e-14)
    # subtract observable in its original coordinates, before scaling
    obs_delta=obs-obs[0]
    obs_peak=spec["trace_peaks"]
    targets=[p["omega"] for p in obs_peak if p["field"]=="a-b"][:3]
    binw=2*np.pi/(len(ts)*dt)
    result["predeclared_selection"]="first three archived a-b spectral peaks, nearest positive-frequency mode; Schur cluster includes all modes within half FFT bin, both signs; no fitted amplitudes"
    result["fft_bin_width"]=float(binw); result["target_omega"]=targets
    out_arrays={"time":ts,"observed_centre_departures":obs_delta,"sample_steps":steps,"C":C,"baseline_drift":b}
    P=arr["step_1"]; P0=arr["step_0"]
    # Eigenvalues near one are better represented via (P-I)/dt.
    A=(P-np.eye(N))/dt; A0=(P0-np.eye(N))/dt
    vals,L,V=eig(A,left=True,right=True)
    vals0=eig(A0,left=False,right=False)
    mu=1+dt*vals; mu0=1+dt*vals0
    positive=np.flatnonzero((np.abs(mu)>1e-8)&(np.angle(mu)>0))
    lam=np.full(N,np.nan+1j*np.nan,dtype=complex)
    good=np.abs(mu)>1e-8; lam[good]=np.log(mu[good])/dt
    overlap=np.sum(L.conj()*V,axis=0)
    cond=np.divide(1.,np.abs(overlap),out=np.full(N,np.inf),where=np.abs(overlap)>0)
    # Frozen map nonnormality, not intermodal transfer. Exact dense SVD on
    # selected times in the declared coordinate norm, and observable row norm.
    power_rows=[]; powers={}
    for m in [1,80,800,1600,3200,int(steps[-1])]:
        Q=np.linalg.matrix_power(P,m); powers[m]=Q
        sv=svdvals(Q)
        row={"steps":m,"t":m*dt,"sigma_max_scaled":float(sv[0]),
             "spectral_radius_power":float(np.max(np.abs(mu))**m),
             "centre_observable_row_norms":np.linalg.norm(C@Q,axis=1).tolist()}
        power_rows.append(row)
    result["finite_time_amplification"]=power_rows
    result["map_max_real_effective"]=float(np.max(lam[good].real))
    result["map_zero_like_multipliers"]=int(np.sum(~good))
    result["epsilon_map_relative_difference"]=ratio(norm(P-P0),norm(P))
    result["geometry_to_non_geometry_output_norm"]=norm(arr["other_1"])
    result["full_map_closure_proven"]=False
    out_arrays.update({"map_eigenvalues":mu,"left_vectors_unit":L,"right_vectors_unit":V,"left_right_overlap":overlap,"eigenvalue_condition":cond})
    # All eight full geometry saves, no interpolation and no modal frequencies
    # inferred from these sub-Nyquist samples.
    X=np.array([np.concatenate([snap[f"t{i}_g_{k}"]-snap["t0_g_"+k] for k in G])*scale for i in range(len(steps))])
    out_arrays["archived_full_geometry_departures_at_saves"]=X
    # Affine propagation from archived baseline with no free coefficients.
    # For unmatched native arms only frozen self-response is reported.
    linear=np.empty((len(ts),12)); u=np.zeros(N)
    save_linear=[]
    for m in range(len(ts)):
        linear[m]=C@u
        if m%800==0: save_linear.append(u.copy())
        if m+1<len(ts): u=P@u+b
    out_arrays["affine_baseline_centre_response"]=linear
    save_linear=np.array(save_linear)
    out_arrays["affine_baseline_full_at_saves"]=save_linear
    if exact:
        result["affine_observable_comparison"]={f:{"observed_rms":rms(obs_delta[:,3*j:3*j+3]),
          "predicted_rms":rms(linear[:,3*j:3*j+3]),
          "unexplained_rms":rms(obs_delta[:,3*j:3*j+3]-linear[:,3*j:3*j+3]),
          "relative_unexplained_rms":ratio(rms(obs_delta[:,3*j:3*j+3]-linear[:,3*j:3*j+3]),rms(obs_delta[:,3*j:3*j+3]))}
          for j,f in enumerate(("a_minus_b","Aa","X","Lambda"))}
        result["affine_full_snapshot_residual_l2"]=[norm(x-y) for x,y in zip(X,save_linear)]
    # Select candidates without fitting amplitudes or treating small bin
    # differences as measured frequency precision. Keep left/right sensitivity.
    candidates=[]; seen=set()
    for targetw in targets:
        j=int(positive[np.argmin(abs(lam[positive].imag-targetw))])
        if j in seen: continue
        seen.add(j)
        v=V[:,j]; l=L[:,j]
        entry={"index":j,"target_omega":targetw,"omega":float(lam[j].imag),"growth":float(lam[j].real),
           "mu":cpair(mu[j]),"within_half_bin":bool(abs(lam[j].imag-targetw)<=binw/2),
           "left_right_overlap_abs":float(abs(overlap[j])),"condition":float(cond[j]) if np.isfinite(cond[j]) else None,
           "nearest_other_multiplier_distance":float(np.min(np.abs(np.delete(mu-mu[j],j)))),
           "nearest_coarse_epsilon_multiplier_distance":float(np.min(abs(mu0-mu[j]))),
           "right_residual":ratio(norm(P@v-mu[j]*v),norm(P)*norm(v)),
           "left_residual":ratio(norm(l.conj()@P-mu[j]*l.conj()),norm(P)*norm(l)),
           "observable_eigenvector": [cpair(x) for x in C@v]}
        if overlap[j]==0:
            entry["modal_projection"]="undefined"; candidates.append(entry); continue
        w=l/np.conj(overlap[j]); bj=w.conj()@b
        entry["biorthogonal_check"]=cpair(w.conj()@v)
        entry["projected_constant_drift"]=cpair(bj)
        z=np.empty(len(ts),complex); z[0]=0
        for m in range(1,len(ts)): z[m]=mu[j]*z[m-1]+bj
        # Coherence in demodulated phase, constant forcing only, not actual excitation.
        phase_terms=bj*np.exp(-1j*np.angle(mu[j])*np.arange(len(ts)-1))
        entry["constant_drift_demodulated_coherence"]=ratio(abs(np.sum(phase_terms)),np.sum(abs(phase_terms)))
        entry["constant_drift_modal_abs_max"]=mx(z)
        entry["constant_drift_modal_abs_final"]=float(abs(z[-1]))
        entry["constant_drift_centre_conjugate_pair_max"]=np.max(abs(2*np.real(z[:,None]*(C@v)[None,:])),axis=0).tolist()
        # Exact block residual over 800 steps relative to frozen reference,
        # not per-step forcing nor its phase coherence. Powers computed directly.
        if exact:
            observed=X@w.conj()
            affine=z[steps]
            block=observed[1:]-mu[j]**np.diff(steps)*observed[:-1]
            block_res=block-(affine[1:]-mu[j]**np.diff(steps)*affine[:-1])
            entry["archived_modal_at_saves"]=[cpair(x) for x in observed]
            entry["affine_modal_at_saves"]=[cpair(x) for x in affine]
            entry["unexplained_block_excitation"]=[cpair(x) for x in block_res]
            entry["signed_block_amplitude_budget"]=[]
            for z0,z1,d,f in zip(observed[:-1],observed[1:],np.diff(steps),block):
                propagated=mu[j]**d*z0
                terms=[(abs(mu[j]**d)**2-1)*abs(z0)**2,2*np.real(np.conj(propagated)*f),abs(f)**2]
                entry["signed_block_amplitude_budget"].append({"intrinsic":float(terms[0]),"cross":float(terms[1]),"excitation_square":float(terms[2]),"closure":float(abs(z1)**2-abs(z0)**2-sum(terms))})
        out_arrays[f"candidate_{j}_constant_drift_modal"]=z
        candidates.append(entry)
    result["candidate_modes"]=candidates
    result["actual_excitation_phase_coherence"]="not measurable from 0.3-spaced full-state archive; every-step centre signals are not full modal coordinates"
    # Schur leading invariant subspace selected by frequency windows,
    # orthogonal projection is not generally an autonomous modal coordinate.
    def select(z):
        m=1+dt*z
        return bool(abs(m)>1e-8 and any(abs(abs(np.angle(m)/dt)-t)<=binw/2 for t in targets))
    T,U,count=schur(A,output="complex",sort=select)
    Q=U[:,:count]; Tm=np.eye(count)+dt*T[:count,:count]
    sub={"dimension":int(count),"right_invariance_relative":ratio(norm(P@Q-Q@Tm),norm(P)*norm(Q)),
         "complement_to_subspace_coupling_norm":norm(Q.conj().T@P@U[:,count:]),
         "note":"orthogonal Schur coordinates can receive complement contributions; no physical energy or causal attribution"}
    if count:
        out_arrays["schur_basis"]=Q; out_arrays["schur_block"]=Tm
        if exact:
            projected=X@Q.conj()
            sub["orthogonal_subspace_norm_at_saves"]=np.linalg.norm(projected,axis=1).tolist()
            sub["orthogonal_subspace_observable_residual_rms_at_saves"]=rms(X@C.T-np.real(projected@Q.T)@C.T)
            sub["complement_state_norm_at_saves"]=np.linalg.norm(X-projected@Q.T,axis=1).tolist()
    result["ordered_schur_cluster"]=sub
    # Common observations explicitly retain signed constraint and lapse values
    # at saved times, without alleging a causal ledger.
    result["saved_constraint_and_lapse_readings"]=[{"t":float(r["t"]),"H":r["H"],"M":r["M"],
        "alpha_stored_cell0":r["cell0"]["alpha"],"alpha_solved_cell0":r["alpha_rhs_cell0"]} for r in saves]
    result["status"]="completed_archive_diagnostic"
    np.savez_compressed(a.out/"diagnostics.npz",**out_arrays)
    write(a.out/"report.json",result)
    print(json.dumps({k:result[k] for k in ("status","mode","switch","stage0b_initial_hash_mismatches","trajectory_projection_eligible","baseline_drift_inf","map_max_real_effective","finite_time_amplification","affine_observable_comparison","candidate_modes","ordered_schur_cluster") if k in result},indent=2,allow_nan=False),flush=True)
if __name__=="__main__": main()
