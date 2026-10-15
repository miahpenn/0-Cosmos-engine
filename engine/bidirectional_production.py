"""GEAR V5.5 bidirectional production coupling.

There is deliberately NO added interface source term.
All local/COSMOS feedback is carried by the single metric + total Tmunu.
This module records the derived exchange ledger only.
"""
import argparse
from pathlib import Path

import numpy as np

from .cosmology_observables import effective_hubble
from . import true_cmc_reference as cmc
from . import reference_pirk_unified as q

def ledger_observables(g,s,fields,Rcoord=10.0):
    e,pr,pt,j=q.matter_projection(g,s,fields)
    r=g.centers; gamma_inv=s.X**2/s.a
    R=r*np.sqrt(s.b)/s.X; Rr=q.D(g,R,1)
    Rt=-s.alpha*R*(s.K/3.-s.Aa/2.)
    chi=gamma_inv*Rr**2-(Rt/s.alpha)**2
    M=.5*R*(1.-chi)
    k=int(np.argmin(abs(r-Rcoord)))
    # Exact spherical unified-first-law current at the coordinate worldtube.
    flux=4*np.pi*R[k]**2*(-s.alpha[k]*gamma_inv[k]*j[k]*Rr[k])
    work=4*np.pi*R[k]**2*(-pr[k]*Rt[k])
    outer=r>=.8*g.r_max
    inner=r<=r[k]
    return {
        'H_eff': effective_hubble(s, g.volumes),
        'tau_rate':float(s.alpha[0]),
        'R_sigma':float(R[k]),
        'M_MS':float(M[k]),
        'chi_sigma':float(chi[k]),
        'flux_T':float(flux),
        'work_pR':float(work),
        'rho_outer':float(np.mean(e[outer])),
        'p_outer':float(np.mean(pr[outer])),
        'j_outer':float(np.mean(j[outer])),
        'rho_inner':float(np.mean(e[inner])),
        'phi_outer':float(np.mean(fields[4][outer])),
        'Pi_outer':float(np.mean(fields[5][outer])),
    }

def run(N=40,T=24,cfl=.06,Rcoord=10.):
    g=q.Grid(N,40.); s,fields=q.make_initial(g,.01,7.)
    s.alpha=cmc.cmc_lapse(g,s,fields)[0]; s.beta.fill(0.); s.B.fill(0.)
    dt=T/round(T/(cfl*g.dr)); nsteps=round(T/dt); dt=T/nsteps; t=0.; hist=[]
    for n in range(nsteps+1):
        z=ledger_observables(g,s,fields,Rcoord); z['t']=t; hist.append(z)
        if n==nsteps: break
        s,fields=cmc.step_true_cmc(g,s,fields,dt); t+=dt
    out={k:np.asarray([z[k] for z in hist]) for k in hist[0]}
    dM=np.gradient(out['M_MS'],out['t'],edge_order=2)
    rhs=out['flux_T']+out['work_pR']
    mask=(out['t']>=2)&(out['t']<=out['t'][-1]-2)
    meta={
      'N':N,'T':float(out['t'][-1]),'dt':dt,
      'max_Misner_Sharp_residual':float(np.max(abs(dM[mask]-rhs[mask]))),
      'rms_Misner_Sharp_residual':float(np.sqrt(np.mean((dM[mask]-rhs[mask])**2))),
      'H_final':float(out['H_eff'][-1]),
      'tau_final':float(np.trapezoid(out['tau_rate'],out['t'])),
      'no_added_source_term':True,
    }
    return out,meta

def save_ledger_archive(out, N, output_dir):
    """Save a ledger using its recorded final time in the artifact filename.

    A filename is part of the evidence provenance: never label an arbitrary
    time span as T24 just because that is the CLI's current default.
    """
    if "t" not in out:
        raise ValueError("ledger output must include its sampled 't' history")
    times = np.asarray(out["t"])
    if (
        times.ndim != 1
        or times.size == 0
        or times.dtype.kind not in "iuf"
        or not np.all(np.isfinite(times))
    ):
        raise ValueError("ledger output 't' must be a non-empty finite numeric vector")
    final_time = float(times[-1])
    if final_time < 0.0:
        raise ValueError("ledger final time must be nonnegative")
    time_tag = f"{final_time:g}".replace("-", "m").replace(".", "p")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    target = output_dir / f"GEAR_V5_5_BIDIRECTIONAL_PRODUCTION_LEDGER_N{N}_T{time_tag}.npz"
    np.savez_compressed(target, **out)
    return target


if __name__ == "__main__":
    import json

    parser = argparse.ArgumentParser(
        description="Run the bounded GEAR V5.5 bidirectional research ledger."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("artifacts/gear_v55"),
        help="directory for NPZ outputs (default: ./artifacts/gear_v55)",
    )
    args = parser.parse_args()

    for N in (40, 60):
        out, meta = run(N)
        output_path = save_ledger_archive(out, N, args.output_dir)
        print(json.dumps({**meta, "output_path": str(output_path)}, indent=2))