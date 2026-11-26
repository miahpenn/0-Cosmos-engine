"""Stage 0b pilot: frozen initial geometry-block spectra, not a full coupled admission.
No changes to engine, boundary, CMC or PIRK2. Separate +/- one-step probes.
Coordinate scaling: a,b,X,dr*Aa,dr*K,dr*Lambda,stored alpha.
Semi-discrete geometry RHS is conditional on incoming stored alpha. Actual step
map retains alpha and includes live matter/scalar stage response, but its
geometry-to-geometry block alone is not the full matter-coupled map.
"""
import argparse, copy, json, pathlib, time
import numpy as np
from scipy.linalg import eig, svdvals
from scipy.signal import detrend, find_peaks
from research.control_tests import ctl, stage0, test1

G=tuple(ctl.FIELDS_G)+("alpha",)
EPS=(1e-4,3e-5)
def write_json(path,obj):
    path.write_text(json.dumps(obj,indent=2,allow_nan=False))
def maxabs(x):
    return float(np.max(np.abs(x)))
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("arm",choices=[k for k,v in test1.ARMS.items() if v[0]=="1a"])
    ap.add_argument("mode",choices=["flat","native"])
    ap.add_argument("switch",choices=["off","on"])
    ap.add_argument("archive",type=pathlib.Path)
    ap.add_argument("out",type=pathlib.Path)
    a=ap.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    _,dr,dt,rm=test1.ARMS[a.arm]
    k,s=stage0.build(dr,rm,a.mode,a.switch=="on")
    n=s.grid.n; scale=np.repeat([1.,1.,1.,dr,dr,dr,1.],n); size=len(scale)
    assert s.electromagnetic is None and s.matter.charged_dust is None
    label=f"stage0_{a.mode}_{a.switch}_{a.arm}"
    archived=json.loads((a.archive/a.switch/(label+".json")).read_text())
    initial=stage0.hashes(stage0.snapshot(s))
    expected=archived["stage0"]["initial_hashes"]
    assert initial.keys()==expected.keys()
    differences=[f for f in initial if initial[f]!=expected[f]]
    # Never claim exact trace correspondence if recreated initial hashes differ.
    prov=stage0.provenance()
    report={"status":"assembling","scientific_admission":False,"provenance":prov,
      "archive_provenance":archived["stage0"]["provenance"],"initial_hashes":initial,
      "initial_hash_mismatches":differences,"dr":dr,"dt":dt,"r_max":rm,
      "scope":"initial-slice conditional geometry RHS and projected geometry block of actual accepted PIRK2 map; not full coupled spectrum",
      "coordinates":list(G),"scale_by_field":scale[::n].tolist(),"eps":list(EPS)}
    write_json(a.out/"report.json",report)
    def pack(st):
        return np.concatenate([getattr(st.geometry,f) for f in G])*scale
    x=pack(s)
    def state(z):
        q=copy.deepcopy(s)
        for i,f in enumerate(G):
            setattr(q.geometry,f,(z[i*n:(i+1)*n]/scale[i*n]).copy())
        return q
    # Input reconstruction itself must be bitwise exact: don't perturb untouched
    # coordinates by scale roundtrips. The scaled displacement is applied only
    # to the selected original array element below.
    def perturb(j,h):
        q=copy.deepcopy(s); f=G[j//n]
        getattr(q.geometry,f)[j%n]+=h/scale[j]
        return q
    dust_scales={sp:maxabs(getattr(s.matter,sp).rest) for sp in ("dark_matter","baryons")}
    def other(st):
        chunks=[getattr(st.scalars,f) for f in ctl.FIELDS_S]
        for sp in ("dark_matter","baryons","radiation"):
            den=dust_scales.get(sp,1.)
            chunks.extend(getattr(getattr(st.matter,sp),f)/den for f in ("rest","energy_t","momentum_r"))
        return np.concatenate(chunks)
    def sample(q):
        g0,kdot,R,_=ctl.semi(k,q)
        rhs=np.concatenate([R[f] for f in G[:-1]])*scale[:-n]
        solved=np.asarray(g0.alpha)
        new=k.step(q,dt)
        return rhs,solved,pack(new),other(new)
    base=sample(copy.deepcopy(s))
    arrays={"baseline_input":x,"baseline_rhs":base[0],"baseline_solved_lapse":base[1],"baseline_step":base[2],"r":np.asarray(s.grid.centers)}
    # Separate baseline drift from the perturbation response.
    report["baseline_step_drift_inf"]=maxabs(base[2]-x)
    report["baseline_rhs_inf"]=maxabs(base[0])
    mats=[]
    for i,h in enumerate(EPS):
        J=np.empty((6*n,size)); A=np.empty((n,size)); P=np.empty((size,size)); O=np.empty((len(base[3]),size))
        for j in range(size):
            plus=sample(perturb(j,h)); minus=sample(perturb(j,-h))
            for mat,p,m in zip((J,A,P,O),plus,minus):
                mat[:,j]=(p-m)/(2*h)
            if j%128==0:
                print(f"{label} eps={h} column {j}/{size}",flush=True)
        for name,M in zip(("rhs","lapse","step","other"),(J,A,P,O)):
            assert np.all(np.isfinite(M)),name
            arrays[f"{name}_{i}"]=M
        mats.append((J,A,P,O))
        np.savez_compressed(a.out/f"assembly_eps{i}.npz",J=J,A=A,P=P,O=O)
    assert stage0.hashes(stage0.snapshot(s))==initial,"base state changed"
    def compare(a,b):
        return {"relative_frobenius":float(np.linalg.norm(a-b)/max(np.linalg.norm(b),np.finfo(float).tiny)),"max_abs":maxabs(a-b)}
    report["epsilon_comparison"]={name:compare(v,w) for name,v,w in zip(("rhs","lapse","step","other"),mats[0],mats[1])}
    J,A,P,O=mats[-1]
    report["non_geometry_output_response"]={"coordinate_convention":"scalar physical; dust divided by initial species maximum rest; radiation physical","frobenius":float(np.linalg.norm(O)),"max_abs":maxabs(O),"not_a_closure_test":True}
    # Observed peaks: retain multiple bins, not just dominant peak. Same window
    # as Stage 0, in cell 0. Phase is Aa relative to a-b at each shared bin.
    tr=archived["centre_trace"]; t=np.array([v["t"] for v in tr])
    obs=np.array([[v[f][0] for v in tr] for f in ("a_minus_b","Aa")])
    w=2*np.pi*np.fft.rfftfreq(len(t),dt)
    ft=np.fft.rfft(detrend(obs,axis=1)*np.hanning(len(t))[None,:],axis=1)
    peaks=[]
    for z in range(2):
        inds=find_peaks(abs(ft[z]))[0]
        inds=sorted(inds,key=lambda j:abs(ft[z,j]),reverse=True)[:8]
        for j in inds:
            peaks.append({"field":("a-b","Aa")[z],"omega":float(w[j]),"omega_dr":float(w[j]*dr),
              "amplitude":float(abs(ft[z,j])),"phase_Aa_minus_a_b":float(np.angle(ft[1,j]*np.conj(ft[0,j])))})
    report["trace_peaks"]=peaks
    # Eigendecomposition in scaled coordinates; no mode selection by target.
    # Map eigenvalues are taken from (P-I)/dt for conditioning; reconstruct mu.
    # Boundary/gauge projections can yield mu approximately zero, for which
    # log(mu)/dt is not a physically meaningful continuum growth rate.
    def spectral(name,M,is_map=False):
        vals,V=eig(M,check_finite=True)
        arrays[name+"_eigvals"]=vals; arrays[name+"_eigvecs"]=V
        mu=1+dt*vals if is_map else vals
        valid=np.abs(mu)>1e-8 if is_map else np.ones(len(vals),dtype=bool)
        lam=np.full(len(vals),complex(np.nan,np.nan))
        if is_map: lam[valid]=np.log(mu[valid])/dt
        else: lam=vals
        aa=V[3*n,:]/dr; ab=V[0,:]-V[n,:]
        frac=np.sum(np.abs(V.reshape((-1,n,len(vals)))[:,:3,:])**2,axis=(0,1))/np.sum(abs(V)**2,axis=0)
        def row(j):
            v=V[:,j]
            res=np.linalg.norm(M@v-vals[j]*v)/(max(np.linalg.norm(M,ord=np.inf),abs(vals[j]),np.finfo(float).tiny)*np.linalg.norm(v))
            return {"index":int(j),"growth":float(lam[j].real),"omega":float(abs(lam[j].imag)),"omega_dr":float(abs(lam[j].imag)*dr),
              "centre3_scaled_norm_fraction":float(frac[j]),"a_minus_b_cell0_abs":float(abs(ab[j])),"Aa_cell0_abs":float(abs(aa[j])),
              "phase_Aa_minus_a_b":float(np.angle(aa[j]*np.conj(ab[j]))* (1 if lam[j].imag>=0 else -1)),
              "relative_eigen_residual":float(res),"map_mu_abs":float(abs(mu[j])) if is_map else None}
        ii=np.flatnonzero(valid)
        positive=ii[lam[ii].imag>0]
        near=[]
        for peak in peaks:
            jj=sorted(positive,key=lambda j:abs(abs(lam[j].imag)-peak["omega"]))[:3]
            near.append({"observed":peak,"nearest": [row(j) for j in jj]})
        return {"dimension":len(vals),"near_zero_map_multipliers":int(np.sum(~valid)),
           "max_real":float(np.max(lam[ii].real)),"max_abs_omega":float(np.max(abs(lam[ii].imag))),
           "top_growth":[row(j) for j in ii[np.argsort(lam[ii].real)[-8:][::-1]]],
           "top_centre_localization":[row(j) for j in ii[np.argsort(frac[ii])[-8:][::-1]]],
           "frequency_neighbours_not_identifications":near}
    report["conditional_rhs_fixed_incoming_lapse"]=spectral("rhs_conditional",J[:,:6*n])
    report["actual_step_projected_geometry"]=spectral("step_geometry",(P-np.eye(size))/dt,True)
    # A fixed-point elimination is only a separate hypothetical diagnostic,
    # since production uses its incoming stored lapse, not repeated solves.
    F=np.eye(n)-A[:,6*n:]; sv=svdvals(F)
    report["lapse_fixed_point_matrix"]={"sigma_min":float(sv[-1]),"sigma_max":float(sv[0]),
      "condition":float(sv[0]/sv[-1]) if sv[-1]>0 else None,"not_eliminated":True}
    # Repeat spectral summary at larger epsilon to expose unstable numerical identifications.
    report["coarser_epsilon_conditional_rhs"]=spectral("rhs_conditional_eps0",mats[0][0][:,:6*n])
    report["coarser_epsilon_actual_step"]=spectral("step_geometry_eps0",(mats[0][2]-np.eye(size))/dt,True)
    np.savez_compressed(a.out/"spectral_arrays.npz",**arrays)
    report["status"]="completed_diagnostic"; write_json(a.out/"report.json",report)
    print(json.dumps({k:report[k] for k in ("status","initial_hash_mismatches","epsilon_comparison","baseline_step_drift_inf","non_geometry_output_response","lapse_fixed_point_matrix")},indent=2),flush=True)
if __name__=="__main__": main()
