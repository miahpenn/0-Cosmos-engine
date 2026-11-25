"""Read signed W1 momentum budgets through unmodified production steps."""
import argparse
import io
import json
from pathlib import Path
import subprocess

import numpy as np

from engine.adaptive_step import advance_with_radiation_admissibility_retries
from engine.production_kernel import V55ProductionKernel
from engine.run_production import assert_finite_state
from engine.scalar_system import ScalarFields
from engine.stress_energy import assemble_total_stress_energy
from engine.v55_initial import SuppliedSlice
from research.w1_fring_smoke import origin_trace as pins
from research.w1_fring_smoke.frozen_centre_derivative import shell_rms
from research.w1_fring_smoke.momentum_precision import ARMS, read_git
from research.w1_fring_smoke.prepare_spatial import DR, sha
from research.w1_fring_smoke.spatial_replay import FOLDER

PARTS = ("Aa_prime", "Aa_over_r_term", "K_gradient", "metric_coupling", "matter_scalar", "matter_fluid")
STAGES = ("initial_lapse", "explicit_predictor", "primary_after_lapse")
LOCAL_CELLS = 24


class Observe(V55ProductionKernel):
    def __init__(self, vacuum):
        super().__init__()
        self.vacuum = vacuum
        self.active = False
        self.geometry_calls = 0
        self.boundary_calls = 0
        self.stages = {}

    def record(self, name, grid, g, scalars, matter):
        r = grid.centers
        d = lambda f: self.vacuum._d(grid,f,1)
        t = assemble_total_stress_energy(grid,g,scalars,matter)
        aa, ab = g.Aa, -.5*g.Aa
        geometry = self.vacuum.momentum_constraint(grid,g)
        terms = {
            "Aa_prime":d(aa), "Aa_over_r_term":(aa-ab)*2/r,
            "K_gradient":-2/3*d(g.K),
            "metric_coupling":6*aa*(-.5*d(g.X)/g.X)+(aa-ab)*d(g.b)/g.b,
            "matter_scalar":-8*np.pi*t.scalar_j, "matter_fluid":-8*np.pi*t.fluid_j,
        }
        terms["M"] = geometry-8*np.pi*t.j
        scale = np.maximum(1e-30,sum(abs(terms[key]) for key in PARTS))
        if np.any(abs(terms["M"]-sum(terms[k] for k in PARTS)) > 32*np.finfo(float).eps*scale):
            raise ValueError("momentum terms do not close")
        terms["K"] = g.K
        terms["j"] = t.j
        self.stages[name] = {key:value[:LOCAL_CELLS].copy() for key,value in terms.items()}

    def _geometry_stage_terms(self, grid, geometry, scalars, matter, electromagnetic):
        if self.active:
            self.record(STAGES[self.geometry_calls],grid,geometry,scalars,matter)
            self.geometry_calls += 1
        return super()._geometry_stage_terms(grid,geometry,scalars,matter,electromagnetic)

    def _apply_outer_boundary_for_state(self, grid, geometry, scalars, matter, electromagnetic):
        if self.active:
            name = ("primary","corrector")[self.boundary_calls]
            self.record(name+"_before_boundary",grid,geometry,scalars,matter)
        super()._apply_outer_boundary_for_state(grid,geometry,scalars,matter,electromagnetic)
        if self.active:
            self.record(name+"_after_boundary",grid,geometry,scalars,matter)
            self.boundary_calls += 1


def local_rms(r, values):
    if r[1]-r[0] <= .02500000000000001:
        return shell_rms(r,values)["rms"]
    mask = (r>=.05) & (r<.25)
    return float(np.sqrt(np.average(values[mask]**2,weights=r[mask]**2)))


def summarize(r, stages):
    result = {}
    reference = stages["initial_lapse"]
    for name, row in stages.items():
        terms = (*PARTS,"M")
        delta = {key:row[key]-reference[key] for key in terms}
        result[name] = {
            "rows":{key:v.tolist() for key,v in row.items()},
            "delta_from_refreshed_start":{key:v.tolist() for key,v in delta.items()},
            "band_rms":{key:local_rms(r,row[key]) for key in terms},
            "band_delta_rms":{key:local_rms(r,delta[key]) for key in terms},
        }
    return result


def provenance():
    parsed = json.loads((FOLDER/"input.json").read_text())
    vendor = subprocess.check_output(["git","-C","vendor/bb-palatini-unified-r0","rev-parse","HEAD"],text=True).strip()
    subprocess.run(["git","diff","--exit-code",pins.SOURCE,"--","engine","vendor"],check=True)
    if vendor != pins.VENDOR or sha(pins.INPUT) != pins.DIGEST:
        raise ValueError("unexpected production or fixture pins")
    return parsed


def advance(kernel,state,dt,watch):
    kernel.active = watch
    kernel.geometry_calls = 0
    kernel.boundary_calls = 0
    kernel.stages = {}
    if watch:
        kernel.record("accepted_start",state.grid,state.geometry,state.scalars,state.matter)
    start = state.t
    update = advance_with_radiation_admissibility_retries(kernel,state,dt)
    if update.retry_count:
        raise ValueError("stage replay retried")
    end = update.state
    if not watch:
        return end,None
    if kernel.geometry_calls!=3 or kernel.boundary_calls!=2:
        raise ValueError("unexpected stage order")
    assert_finite_state(end)
    kernel.record("accepted_end",end.grid,end.geometry,end.scalars,end.matter)
    r = end.grid.centers[:LOCAL_CELLS]
    rows = summarize(r,kernel.stages)
    return end,{"t_start":start,"t_end":end.t,"stages":rows}


def fine_replay(dr, go, vacuum, input_meta):
    dt = .000375
    n = round(60/dr)
    grid = go.SphericalCellGrid(n,60.)
    path = FOLDER/"static_0p0125.npz" if dr==DR else pins.INPUT
    digest = input_meta["input_sha256"] if dr==DR else pins.DIGEST
    prefix = "" if dr==DR else "dr0.025_"
    if sha(path)!=digest:
        raise ValueError("input digest does not match")
    with np.load(path,allow_pickle=False) as data:
        if dr==DR:
            np.testing.assert_array_equal(data["r"],grid.centers)
        fields = [data[prefix+k] for k in ("S","PS","D","PD")]
    gc = SuppliedSlice(ScalarFields(*fields,np.full(n,.179055),np.full(n,.00341913)),
                       grid.centers.copy(),input_meta["H0"],input_meta["edge_tol"],f"W1 static {digest}")
    kernel = Observe(vacuum)
    state = kernel.initialize(resolution=n,r_max=60.,outer_edge="open_edge",
                              eight_node_areal=True,initial_structure="supplied_slice",supplied_slice=gc)
    archived = json.loads((FOLDER/"comparison.json").read_text())["runs"][str(dr)][str(dt)]["rows"]
    reports = []
    for step in range(1,41):
        state,sampled = advance(kernel,state,dt,step in (1,20,40))
        if sampled is not None:
            reports.append(sampled)
        if step in (20,40):
            observed = sampled["stages"]["accepted_end"]["rows"]["M"][0]
            if observed != archived[step//10]["cell_0"]["M"]:
                raise ValueError("accepted centre M does not reproduce archived same-grid readback")
            print(f"dr={dr} step {step}/40 matches the archive",flush=True)
    if not np.isclose(state.t,.015,atol=2e-17,rtol=0):
        raise ValueError("unexpected endpoint")
    return {"dr":dr,"dt":dt,"steps":40,"r":grid.centers[:LOCAL_CELLS].tolist(),
            "input_sha256":digest,"observed_steps":reports}


def archived_replay(arm,vacuum):
    dt = .0015 if arm=="A1" else .00075
    kernel = Observe(vacuum)
    state = kernel.initialize(resolution=1200,r_max=60.,outer_edge="open_edge",
                              eight_node_areal=True,initial_structure="none")
    rec = json.loads(read_git(ARMS/f"{arm}.json"))
    path = ARMS/f"{arm}-t0.015.npz"
    content = read_git(path)
    with np.load(io.BytesIO(content),allow_pickle=False) as data:
        for key in vacuum.ALL_FIELDS:
            np.copyto(getattr(state.geometry,key),data[f"geometry_{key}"])
        state.scalars = ScalarFields(*(data[f"scalars_{k}"] for k in ("S","PS","D","PD","phi","Pi")))
        for species in ("dark_matter","baryons","radiation"):
            field = getattr(state.matter,species)
            for key in ("rest","energy_t","momentum_r"):
                np.copyto(getattr(field,key),data[f"matter_{species}_{key}"])
        state.t,state.tau,state.e_folds = map(float,data["clock"])
    steps = round(.0015/dt)
    reports = []
    for i in range(steps):
        state,row = advance(kernel,state,dt,True)
        reports.append(row)
        start_step = rec["saves"]["0.015"]["step"]
        for key,offset in (("accepted_start",i),("accepted_end",i+1)):
            np.testing.assert_array_equal(row["stages"][key]["rows"]["M"][:6],
                                          rec["per_step"][start_step+offset]["M"])
    import hashlib
    return {"dr":.05,"dt":dt,"steps":steps,"r":state.grid.centers[:LOCAL_CELLS].tolist(),
            "state_sha256":hashlib.sha256(content).hexdigest(),"observed_steps":reports}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case",choices=("0.025","0.0125","A1","A2"))
    parser.add_argument("output",type=Path)
    args = parser.parse_args()
    meta = provenance()
    go,vacuum,_ = pins.ad.vendor_modules()
    if vacuum.KO_EPSILON != 0:
        raise ValueError("expected KO off")
    report = (archived_replay(args.case,vacuum) if args.case in ("A1","A2")
              else fine_replay(float(args.case),go,vacuum,meta))
    report.update(engine=pins.SOURCE,vendor=pins.VENDOR,script_sha256=sha(Path(__file__)))
    args.output.write_text(json.dumps(report,allow_nan=False,separators=(",",":"))+"\n")
