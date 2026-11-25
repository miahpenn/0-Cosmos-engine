"""Replay the fine and reference W1 grids at matched timesteps."""
import argparse
import json
import math
from pathlib import Path
import subprocess

import numpy as np

from engine import scalar_system as ss, v55_initial as vi
from engine.adaptive_step import advance_with_radiation_admissibility_retries
from engine.production_kernel import V55ProductionKernel
from engine.run_production import assert_finite_state
from engine.stress_energy import assemble_total_stress_energy
from research.w1_fring_smoke import areal_trace, origin_trace as trace
from research.w1_fring_smoke.prepare_spatial import DR, sha

FOLDER = Path("results/research/w1_fring_spatial")
SAMPLE_INTERVAL = .00375
WIDE = np.longdouble


class MetricWrites(V55ProductionKernel):
    def __init__(self):
        super().__init__()
        self.rates = []

    def _geometry_stage_terms(self, grid, geometry, scalars, matter, electromagnetic):
        result = super()._geometry_stage_terms(grid, geometry, scalars, matter, electromagnetic)
        self.rates.append((result["explicit"], geometry.alpha.astype(WIDE)*geometry.Aa))
        return result


def profiles(state, vacuum, tracefree, writes):
    grid, metric = state.grid, state.geometry
    stress = assemble_total_stress_energy(grid, metric, state.scalars, state.matter)
    raw = vacuum.constraints(grid, metric)
    R = vacuum.geometry_terms(grid, metric)["R"]
    H = raw["hamiltonian"] - 16*np.pi*stress.rho
    M = raw["momentum"] - 8*np.pi*stress.j
    terms = {"H": H, "M": M, "tracefree_rate": tracefree}
    terms.update({key: getattr(state.scalars, key) for key in ("S","D","PS","PD")})
    a, b, X = (getattr(metric, key).astype(WIDE) for key in ("a","b","X"))

    def direct(aa, bb, xx, width):
        x2_offset = (xx-1)*(xx+1)
        q = (x2_offset - (aa-1))/aa
        u = ((bb-1)-x2_offset)/(xx*(np.sqrt(bb)+xx))
        return areal_trace.areal_scalar(grid, q, u, width)

    # Hold the recorded stage rates fixed and remove all cumulative metric writes.
    # A separate evolved storage experiment would also change the stage rates.
    for width in areal_trace.WIDTHS:
        curvature = direct(a, b, X, width)
        terms[f"H_direct_{width}"] = curvature + (raw["hamiltonian"]-R) - 16*np.pi*stress.rho
        without_writes = direct(a-writes["a"], b-writes["b"], X-writes["X"], width)
        terms[f"curvature_writes_{width}"] = curvature-without_writes
        reduced_a = a * np.exp(2*tracefree)
        reduced_b = b * np.exp(-tracefree)
        terms[f"curvature_tracefree_{width}"] = curvature-direct(reduced_a,reduced_b,X,width)
    terms["R_areal"] = grid.centers*np.sqrt(b)/X
    return {key: np.asarray(value, dtype=float) for key, value in terms.items()}


def run(dr, dt, output, zero_only):
    vendor_path = "vendor/bb-palatini-unified-r0"
    vendor = subprocess.check_output(["git","-C",vendor_path,"rev-parse","HEAD"],text=True).strip()
    subprocess.run(["git","diff","--exit-code",trace.SOURCE,"--","engine","vendor"],check=True)
    if vendor != trace.VENDOR or sha(trace.INPUT) != trace.DIGEST:
        raise ValueError("unexpected source pins")
    go, vacuum, _ = trace.ad.vendor_modules()
    if vacuum.KO_EPSILON != 0 or np.finfo(WIDE).nmant < 63:
        raise ValueError("expected KO off and at least 64 significand bits")
    source = json.loads((FOLDER/"input.json").read_text())
    if dr == DR:
        input_path, prefix = FOLDER/"static_0p0125.npz", ""
        digest = source["input_sha256"]
    else:
        input_path, prefix = trace.INPUT, "dr0.025_"
        digest = trace.DIGEST
    if sha(input_path) != digest:
        raise ValueError("W1 input digest mismatch")
    n = round(60/dr)
    r = go.SphericalCellGrid(n, 60.).centers
    with np.load(input_path, allow_pickle=False) as data:
        if not prefix:
            np.testing.assert_array_equal(r,data["r"])
        values = [np.array(data[prefix+k],dtype=float) for k in ("S","PS","D","PD")]
    supplied = vi.SuppliedSlice(ss.ScalarFields(*values,np.full(n,.179055),np.full(n,.00341913)),
                               r.copy(),source["H0"],source["edge_tol"],
                               f"W1 static phase, native grid build {dr}, sha256 {digest}")
    kernel = MetricWrites()
    state = kernel.initialize(resolution=n,r_max=60.,outer_edge="open_edge",
                              eight_node_areal=True,initial_structure="supplied_slice",supplied_slice=supplied)
    writes = {key: np.zeros(n,dtype=WIDE) for key in ("a","b","X")}
    tracefree = np.zeros(n,dtype=WIDE)
    initial = profiles(state,vacuum,tracefree,writes)
    print(json.dumps({"dr":dr,"dt":dt,"initial_H_max":float(np.max(abs(initial["H"]))),
                      "initial_M_max":float(np.max(abs(initial["M"]))),
                      "initial_centre_H":float(initial["H"][0])}),flush=True)
    if zero_only:
        return
    rows = [initial]
    times = [state.t]
    every = round(SAMPLE_INTERVAL/dt)
    steps = 4*every
    if not math.isclose(every*dt,SAMPLE_INTERVAL,abs_tol=1e-17,rel_tol=0):
        raise ValueError("timestep must divide the sampling interval")
    for i in range(1,steps+1):
        kernel.rates.clear()
        old = state.geometry
        advanced = advance_with_radiation_admissibility_retries(kernel,state,dt)
        if advanced.retry_count or len(kernel.rates) != 3:
            raise RuntimeError("unexpected retry or stage count")
        state = advanced.state
        start, _, primary = kernel.rates
        tracefree += (WIDE(dt)/2)*(start[1]+primary[1])
        for key in writes:
            initial_field = getattr(old,key)
            final_field = getattr(state.geometry,key)
            increment = (dt/2)*(start[0][key]+primary[0][key])
            np.testing.assert_array_equal(final_field,initial_field+increment)
            writes[key] += (final_field.astype(WIDE)-initial_field.astype(WIDE))-increment.astype(WIDE)
        if i % every == 0:
            assert_finite_state(state)
            rows.append(profiles(state,vacuum,tracefree,writes))
            times.append(state.t)
            print(f"{dr:g} {dt:g} step {i}/{steps} t {state.t:.8f}",flush=True)
    np.testing.assert_allclose(times,np.arange(5)*SAMPLE_INTERVAL,atol=steps*np.spacing(.015),rtol=0)
    meta = {"dr":dr,"dt":dt,"steps":steps,"source":trace.SOURCE,"vendor":vendor,
            "input_sha256":digest,"script_sha256":sha(Path(__file__))}
    np.savez_compressed(output,metadata=json.dumps(meta),t=times,r=r,
                        **{key:np.stack([row[key] for row in rows]) for key in initial})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dr",type=float,choices=(.025,DR))
    parser.add_argument("dt",type=float,choices=(.000375,.0001875))
    parser.add_argument("output",type=Path)
    parser.add_argument("--zero-only",action="store_true")
    args = parser.parse_args()
    run(args.dr,args.dt,args.output,args.zero_only)
