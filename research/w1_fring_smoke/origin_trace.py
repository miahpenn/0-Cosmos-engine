"""Short, read-only H/M localization on the corrected W1 smoke inputs."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

from engine import scalar_system as ss, v55_initial as vi, v55_pirk_adapter as ad
from engine.adaptive_step import advance_with_radiation_admissibility_retries
from engine.production_kernel import V55ProductionKernel
from engine.run_production import assert_finite_state
from engine.stress_energy import assemble_total_stress_energy

SOURCE = "05c2a09aae5fc1224f110c055c0c5b4b3acf8b21"
VENDOR = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
INPUT = Path("results/research/w1_static_slice/native_static.npz")
DIGEST = "339652a2a7e1379b19812893e43bafcf9ab0211422f2b8685c951107b863a295"
T_END = .015
SAMPLES = [0, 1, 2, 3, 4, 5, -3, -2, -1]


def snapshot(grid, geometry, scalars, matter, vacuum):
    stress = assemble_total_stress_energy(grid, geometry, scalars, matter)
    raw = vacuum.constraints(grid, geometry)
    Ricci = vacuum.geometry_terms(grid, geometry)["R"]
    a, b, X, Aa, K = (getattr(geometry, name) for name in ("a", "b", "X", "Aa", "K"))
    r = grid.centers
    M_Aa = grid.cell_derivative_fourth(Aa, parity=1) + 3*Aa/r
    M_K = -(2/3)*grid.cell_derivative_fourth(K, parity=1)
    H = raw["hamiltonian"]-16*np.pi*stress.rho
    M = raw["momentum"]-8*np.pi*stress.j
    signed = dict(H=H, M=M, H_geometry=raw["hamiltonian"], Ricci=Ricci,
                  H_source=16*np.pi*stress.rho, H_kinematic=raw["hamiltonian"]-Ricci,
                  M_geometry=raw["momentum"], M_Aa=M_Aa, M_K=M_K,
                  M_other_geometry=raw["momentum"]-M_Aa-M_K,
                  M_source=8*np.pi*stress.j, M_scalar=8*np.pi*stress.scalar_j,
                  M_fluid=8*np.pi*stress.fluid_j, Aa=Aa, Aa_over_r2=Aa/r**2,
                  anisotropy_over_r2=(a/b-1)/r**2, C=raw["connection"], alpha=geometry.alpha)
    return dict(cells={k: v[SAMPLES].tolist() for k,v in signed.items()},
                maxima={k: {"cell": (i:=int(np.argmax(abs(v)))), "signed":float(v[i])}
                        for k,v in (("H", H),("M", M))})


class Observer(V55ProductionKernel):
    def __init__(self, vacuum):
        super().__init__()
        self.vacuum = vacuum
        self.enabled = False
        self.term_calls = 0
        self.boundary_calls = 0
        self.stages = []

    def record(self, stage, grid, geometry, scalars, matter):
        self.stages.append(dict(stage=stage, **snapshot(grid, geometry, scalars, matter, self.vacuum)))

    def _geometry_stage_terms(self, grid, geometry, scalars, matter, electromagnetic):
        if self.enabled:
            name = ("after_initial_lapse_refresh", "explicit_predictor", "primary_predictor_after_lapse")
            self.record(name[self.term_calls], grid, geometry, scalars, matter)
            self.term_calls += 1
        return super()._geometry_stage_terms(grid, geometry, scalars, matter, electromagnetic)

    def _apply_outer_boundary_for_state(self, grid, geometry, scalars, matter, electromagnetic):
        if self.enabled:
            name = ("primary_predictor", "corrected")[self.boundary_calls]
            self.record(name+"_before_boundary", grid, geometry, scalars, matter)
        super()._apply_outer_boundary_for_state(grid, geometry, scalars, matter, electromagnetic)
        if self.enabled:
            self.record(name+"_after_boundary", grid, geometry, scalars, matter)
            self.boundary_calls += 1


def run(dr, arm, half):
    parent = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    vendor = subprocess.check_output(["git", "-C", "vendor/bb-palatini-unified-r0", "rev-parse", "HEAD"], text=True).strip()
    if parent != SOURCE or vendor != VENDOR or hashlib.sha256(INPUT.read_bytes()).hexdigest() != DIGEST:
        raise ValueError("unexpected smoke source, vendor or static input")
    go,vac,_ = ad.vendor_modules()
    if vac.KO_EPSILON != 0:
        raise ValueError("expected KO off")
    n = round(60/dr)
    grid = go.SphericalCellGrid(n,60.)
    opts = dict(resolution=n,r_max=60.,outer_edge="open_edge",eight_node_areal=True)
    kernel = Observer(vac)
    if arm == "A":
        key = f"dr{dr:g}"
        labels = json.loads(Path("results/research/w1_c9_static/input_checks/report.json").read_text())
        with np.load(INPUT,allow_pickle=False) as z:
            values = [np.array(z[key+"_"+k],dtype=float) for k in ("S","PS","D","PD")]
        values.extend((np.full(n,.179055),np.full(n,.00341913)))
        supplied = vi.SuppliedSlice(ss.ScalarFields(*values),np.asarray(grid.centers).copy(),
                                   .0056689681106181795,{.025:1.943e-7,.05:1.949e-7}[dr],
                                   labels["grids"][key]["G0.1_label"])
        state = kernel.initialize(**opts,initial_structure="supplied_slice",supplied_slice=supplied)
    else:
        state = kernel.initialize(**opts,initial_structure="none")

    def accepted():
        assert_finite_state(state)
        return dict(t=state.t, **snapshot(state.grid,state.geometry,state.scalars,state.matter,vac))

    rows = [accepted()]
    dt = .03*dr/(2 if half else 1)
    attempts = 0
    for step in range(round(T_END/dt)):
        kernel.enabled = step == 0
        advanced = advance_with_radiation_admissibility_retries(kernel,state,dt)
        if advanced.retry_count:
            raise RuntimeError("retry would break the matched times")
        attempts += len(advanced.attempts)
        state = advanced.state
        rows.append(accepted())
    if kernel.term_calls!=3 or kernel.boundary_calls!=2 or not np.isclose(state.t,T_END,atol=1e-17,rtol=0):
        raise RuntimeError("unexpected stage count or end time")
    return dict(source=parent,vendor=vendor,input_sha256=DIGEST, script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                dr=dr,arm=arm,dt=dt,attempts=attempts,sample_cells=SAMPLES,
                sample_r=grid.centers[SAMPLES].tolist(),rows=rows,first_step_stages=kernel.stages)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("dr",type=float,choices=(.025,.05))
    p.add_argument("arm",choices=("A","B"))
    p.add_argument("--half",action="store_true")
    args = p.parse_args()
    print(json.dumps(run(args.dr,args.arm,args.half),indent=2,allow_nan=False))
