"""Resolution audit for exact production constraint witnesses at strong D."""
from __future__ import annotations
import json
from pathlib import Path
from .production_kernel import V55ProductionKernel

def run_case(N):
    k=V55ProductionKernel()
    s=k.initialize(resolution=N,r_max=float(N),D_amplitude=1e-4,include_radiation=True)
    dt_nom=.03*s.grid.dr
    while s.t < 24.0-1e-14:
        s=k.step(s,min(dt_nom,24.0-s.t))
    o=s.history[-1]
    return {
        "resolution":N,"t":s.t,"tau":s.tau,
        "hamiltonian_max":o["hamiltonian_max"],
        "hamiltonian_normalized_max":o["hamiltonian_normalized_max"],
        "hamiltonian_l2_inner":o["hamiltonian_l2_inner"],
        "hamiltonian_l2_outer":o["hamiltonian_l2_outer"],
        "momentum_max":o["momentum_max"],
        "momentum_normalized_max":o["momentum_normalized_max"],
        "momentum_l2_inner":o["momentum_l2_inner"],
        "momentum_l2_outer":o["momentum_l2_outer"],
        "cmc_residual_outer_max":o["cmc_residual_outer_max"],
        "determinant_constraint_max":o["determinant_constraint_max"],
        "determinant_min":o["determinant_min"],
        "trapping_min":o["trapping_min"],
        "trapped_roots_count":len(o["trapped_roots"]),
        "lapse_min":o["lapse_min"],
        "lapse_min_r":o["lapse_min_r"],
    }

def main(output="runs/constraint-resolution/report.json"):
    rows=[run_case(160),run_case(320)]
    out={"status":"completed","D_amplitude":1e-4,"final_time":24.0,
         "diagnostic_only":True,"cases":rows}
    Path(output).parent.mkdir(parents=True,exist_ok=True)
    Path(output).write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))
if __name__=="__main__": main()
