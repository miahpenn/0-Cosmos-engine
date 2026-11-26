"""Stage 0c.4a: single-host paired runtime diagnostic, no accepted step.
Run only under the exact Stage 0c.3 dependencies. No engine change.
"""
import contextlib, hashlib, io, json, os, pathlib, platform, subprocess, sys

MASK = "X86_V4,AVX512_ICL,AVX512_SPR"
SOURCE = "cbf4e6bd759271b2ca15295db8f13e570e8fa3aa"
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
def write(p, x):
    p.write_text(json.dumps(x, indent=2, allow_nan=False))
def sha(p):
    with open(p, "rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()

def child(ref, out, arm):
    import numpy as np
    import scipy
    import scipy.linalg
    from threadpoolctl import threadpool_info
    from engine import cmc_gauge as cmc, production_kernel as pk
    from research.control_tests import ctl, stage0
    assert platform.python_version() == "3.12.16"
    assert np.__version__ == "2.5.5" and scipy.__version__ == "1.18.1"
    p = stage0.provenance()
    assert p["head"] == SOURCE and VENDOR in p["submodules"] and not p["tracked_status"]
    assert all(os.environ.get(k) == "1" for k in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS"))
    features = {k: bool(v) for k,v in np._core._multiarray_umath.__cpu_features__.items()}
    runtime = io.StringIO()
    with contextlib.redirect_stdout(runtime):
        np.show_runtime()
    pools = threadpool_info()
    # A preflight record is retained even if the host cannot run the approved mask comparison.
    write(out/"preflight.json", {"features":features,"runtime":runtime.getvalue(),"pools":pools,
        "provenance":p,"mask":os.environ.get("NPY_DISABLE_CPU_FEATURES"),
        "cpu":subprocess.check_output(["lscpu"],text=True)})
    assert features.get("X86_V3"), "X86_V3 not available"
    if arm == "default":
        assert features.get("X86_V4"), "AVX512-capable default host required"
    else:
        assert not any(features.get(k,False) for k in MASK.split(",")), "mask ineffective"
    with np.load(ref/"preflight/flat_dr0.05_r8_off.npz",allow_pickle=False) as z:
        raw={k:z[k].copy() for k in z.files}
    old=json.loads((ref/"off/stage0_flat_off_1a_dr05_dt0_r8.json").read_text())
    expected=old["stage0"]["initial_hashes"]
    assert stage0.hashes(raw)==expected
    assert old["stage0"]["provenance"]["head"]==SOURCE
    k,s=stage0.build(.05,8.,"flat",False)
    before=stage0.snapshot(s)
    assert np.array_equal(before["r"],raw["r"]) and np.array_equal(before["clocks"],raw["clocks"])
    assert s.electromagnetic is None and s.matter.charged_dust is None
    for f in ctl.VAC.ALL_FIELDS: setattr(s.geometry,f,raw["g_"+f].copy())
    for f in ctl.FIELDS_S: getattr(s.scalars,f)[:]=raw["s_"+f]
    for species in ("dark_matter","baryons","radiation"):
        for f in ("rest","energy_t","momentum_r"):
            getattr(getattr(s.matter,species),f)[:]=raw["m_"+species+"_"+f]
    assert stage0.hashes(stage0.snapshot(s))==expected
    assert ctl.ss.V0==old["stage0"]["V0"] and ctl.ss.SD_QUARTIC_LAMBDA==old["stage0"]["SD_QUARTIC_LAMBDA"]
    events=[]; arrays={}
    def snapshot(g,sc,m):
        return stage0.snapshot(pk.ProductionState(grid=s.grid,geometry=g,scalars=sc,matter=m,
            t=s.t,tau=s.tau,e_folds=s.e_folds))
    def geo(grid,g):
        return {"r":np.asarray(grid.centers),**{f:np.asarray(getattr(g,f)) for f in ctl.VAC.ALL_FIELDS}}
    def flatten(x,prefix=""):
        if isinstance(x,dict):
            for key,v in x.items(): yield from flatten(v,prefix+str(key)+"/")
        elif isinstance(x,(tuple,list)):
            for i,v in enumerate(x): yield from flatten(v,prefix+str(i)+"/")
        else:
            a=np.asarray(x)
            if a.dtype.kind in "biufc": yield prefix.rstrip("/"),a.copy()
    def record(name,inputs,outputs):
        i=len(events); vals=dict(flatten(outputs)); assert vals
        assert all(np.all(np.isfinite(v)) for v in vals.values())
        keys={f:f"e{i:03d}_{j:03d}" for j,f in enumerate(vals)}
        arrays.update({keys[f]:v for f,v in vals.items()})
        events.append({"event":name,"input":stage0.hashes({f:np.asarray(v) for f,v in inputs.items()}),
            "output":stage0.hashes(vals),"array_keys":keys})
    # Instrument only call boundaries, using original functions and unchanged return values.
    original=cmc._proper_volume_weights
    def weights(grid,g):
        v=original_weights(grid,g); record("CMC volume weights",geo(grid,g),v); return v
    original_weights=original; cmc._proper_volume_weights=weights
    original_project=cmc._project_kdot_from_weights
    def project(grid,g,sc,m,w,total_weight,**kw):
        inputs=snapshot(g,sc,m); inputs.update(weights=np.asarray(w),total_weight=np.asarray(total_weight))
        v=original_project(grid,g,sc,m,w,total_weight,**kw)
        record("CMC projected Kdot",inputs,v); return v
    cmc._project_kdot_from_weights=project
    original_band=cmc.solve_banded
    def band(lu,ab,b,**kw):
        inputs={"band":np.asarray(ab).copy(),"rhs":np.asarray(b).copy(),"bandwidth":np.asarray(lu)}
        v=original_band(lu,ab,b,**kw); record("CMC banded solve",inputs,v); return v
    cmc.solve_banded=band
    _,vac,_=pk.adapter.vendor_modules()
    original_l2=vac.primary_l2_rhs
    def l2(grid,g,*a,**kw):
        inputs=geo(grid,g); v=original_l2(grid,g,*a,**kw)
        record("primary L2 K/Aa",inputs,{f:v[f] for f in ("K","Aa")}); return v
    vac.primary_l2_rhs=l2
    original_l3=pk.adapter.primary_l3_with_matter
    def l3(grid,g,sc,m,*a,**kw):
        inputs=snapshot(g,sc,m); v=original_l3(grid,g,sc,m,*a,**kw)
        record("primary L3 K/Aa",inputs,{f:v[f] for f in ("K","Aa")}); return v
    pk.adapter.primary_l3_with_matter=l3
    original_solve=k._solve_lapse_for_state
    def solve(grid,g,sc,m,em):
        inputs=snapshot(g,sc,m); v=original_solve(grid,g,sc,m,em)
        record("stage lapse refresh",inputs,v); return v
    k._solve_lapse_for_state=solve
    original_terms=k._geometry_stage_terms
    terms=[]
    def stage_terms(grid,g,sc,m,em):
        inputs=snapshot(g,sc,m); v=original_terms(grid,g,sc,m,em)
        record("stage geometry RHS",inputs,v)
        terms.append((inputs,v)); return v
    k._geometry_stage_terms=stage_terms
    class StopAtFirstPredictor(Exception): pass
    def stop(grid,g,sc,m,em):
        assert len(terms)==2
        start,first=terms[0]; pred,second=terms[1]
        inputs={f"base_{f}":start["g_"+f] for f in ("K","Aa")}
        outputs={}
        for f in ("K","Aa"):
            for name,v in (("L2_start",first["primary_l2"][f]),("L2_pred",second["primary_l2"][f]),("L3_start",first["primary_l3"][f])):
                inputs[f"{name}_{f}"]=v
            inputs["dt"]=np.asarray(.000375)
            outputs[f]=getattr(g,f).copy()
            exact=(start["g_"+f]+.000375*(.5*first["primary_l2"][f]+.5*second["primary_l2"][f]+first["primary_l3"][f]))
            assert np.array_equal(exact,outputs[f])
        record("first primary predictor K/Aa writes (pre-boundary)",inputs,outputs)
        raise StopAtFirstPredictor
    k._apply_outer_boundary_for_state=stop
    try:
        k.step(s,.000375)
        raise AssertionError("unexpected accepted step")
    except StopAtFirstPredictor:
        pass
    finally:
        assert stage0.hashes(stage0.snapshot(s))==expected
    pools=threadpool_info()
    for q in pools:
        q["binary_sha256"]=sha(q["filepath"])
    assert pools and all(q["num_threads"]==1 for q in pools if q["user_api"]=="blas")
    np.savez_compressed(out/"events.npz",**arrays)
    write(out/"report.json",{"arm":arm,"features":features,"pools":pools,"events":events,
        "initial_hashes":expected,"immutable":True,"accepted_steps":0,
        "numpy":np.__version__,"scipy":scipy.__version__,"python":platform.python_version(),
        "numpy_file":np.__file__,"scipy_file":scipy.__file__})

def parent(ref,out):
    # No imports of NumPy/SciPy before spawning fresh children.
    def run(arm,repeat):
        d=out/f"{arm}-{repeat}"; d.mkdir()
        env=os.environ.copy()
        env.pop("NPY_DISABLE_CPU_FEATURES",None)
        if arm=="masked": env["NPY_DISABLE_CPU_FEATURES"]=MASK
        env["OPENBLAS_VERBOSE"]="2"
        with (d/"stdout.txt").open("w") as f:
            p=subprocess.run([sys.executable,__file__,"child",str(ref),str(d),arm],env=env,stdout=f,stderr=subprocess.STDOUT)
        if p.returncode:
            write(out/"summary.json",{"status":"blocked_preflight_or_probe","arm":arm,"returncode":p.returncode,"accepted_steps":0})
            raise SystemExit(p.returncode)
        return json.loads((d/"report.json").read_text())
    a,b=run("default",1),run("masked",1)
    def match(a,b):
        return a["pools"]==b["pools"] and a["initial_hashes"]==b["initial_hashes"]
    def first(a,b):
        if len(a["events"])!=len(b["events"]): return {"status":"confounded_event_count"}
        for i,(x,y) in enumerate(zip(a["events"],b["events"])):
            if x["event"]!=y["event"] or x["input"]!=y["input"]:
                return {"status":"confounded_unobserved_operand_divergence","index":i,"event":x["event"]}
            if x["output"]!=y["output"]:
                return {"status":"first_matched_operand_output_divergence","index":i,"event":x["event"],
                        "fields":[f for f in x["output"] if x["output"][f]!=y["output"][f]]}
        return {"status":"no_difference_detected_in_sampled_operators"}
    if not match(a,b):
        result={"status":"confounded_blas_or_initial_operands"}
    else:
        result=first(a,b)
        if result["status"]=="first_matched_operand_output_divergence":
            c,d=run("default",2),run("masked",2)
            if not(match(a,c) and match(b,d) and first(c,d)==result and
                   c["events"]==a["events"] and d["events"]==b["events"]):
                result["status"]="divergence_not_reproduced_or_runtime_changed"
            else:
                import numpy as np
                i=result["index"]; x=a["events"][i]; y=b["events"][i]; differences={}
                with np.load(out/"default-1/events.npz",allow_pickle=False) as aa,np.load(out/"masked-1/events.npz",allow_pickle=False) as bb:
                    for f in result["fields"]:
                        av=aa[x["array_keys"][f]]; bv=bb[y["array_keys"][f]]
                        delta=bv-av
                        differences[f]={"signed_min":float(np.min(delta)),"signed_max":float(np.max(delta)),
                            "nonzero":int(np.count_nonzero(delta)),"shape":list(delta.shape),
                            "signed_delta":delta.tolist()}
                result["differences_masked_minus_default"]=differences
                result["repeat_verified"]=True
    result.update({"accepted_steps":0,"scientific_admission":False,"scope":"one host; first predictor only; no accepted evolution"})
    write(out/"summary.json",result)
    print(json.dumps(result,allow_nan=False),flush=True)

if __name__=="__main__":
    mode,ref,out,*rest=sys.argv[1:]
    (child if mode=="child" else parent)(pathlib.Path(ref),pathlib.Path(out),*rest)
