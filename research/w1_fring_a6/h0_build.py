"""W1 starting-state H wiggle, part A (engine venv): rebuild the t = 0 states the archived spatial replays started from,
capture the areal loop's last source, and store the near-centre slices for the 60-digit part (h0_mp.py).

Zero steps: kernel.initialize only (the supplied-slice initial-data build and its lapse solve), exactly as
research/w1_fring_smoke/spatial_replay.py builds its start. The stored t = 0 H of the archived replays is checked
bitwise against the rebuilt H before anything is kept.

Usage: python -m research.w1_fring_a6.h0_build INPUT_DIR RUNS_DIR OUT.npz
  INPUT_DIR: holds static_0p0125.npz + input.json (fine) and native_static.npz (dr 0.025, 0.05)
  RUNS_DIR:  holds sp/<dr>-0.000375.npz (archived replays; t = 0 H for the bitwise check)
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from engine import scalar_system as ss, v55_initial as vi, v55_pirk_adapter as ad
from engine.production_kernel import V55ProductionKernel
from engine.stress_energy import assemble_total_stress_energy

N0 = {0.05: 24, 0.025: 34, 0.0125: 48}   # cells kept: the band (r <= 0.25) plus stencil and integral reach


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main(inp, runs, out):
    inp, runs = Path(inp), Path(runs)
    meta = json.loads((inp / "input.json").read_text())
    _, vac, _ = ad.vendor_modules()
    saved, report = {}, {"inputs": {}}
    for dr in (0.05, 0.025, 0.0125):
        n = round(60 / dr)
        if dr == 0.0125:
            path, prefix = inp / "static_0p0125.npz", ""
            assert sha(path) == meta["input_sha256"]
        else:
            path, prefix = inp / "native_static.npz", f"dr{dr}_"
        report["inputs"][str(dr)] = {"file": path.name, "sha256": sha(path)}
        with np.load(path, allow_pickle=False) as z:
            vals = [np.array(z[prefix + k], dtype=float) for k in ("S", "PS", "D", "PD")]
        srcs = []
        orig = vi._eight_node_origin_integral

        def cap(grid, rr, h, f):
            srcs.append(np.array(f, dtype=float))
            return orig(grid, rr, h, f)

        vi._eight_node_origin_integral = cap
        try:
            go = ad.vendor_modules()[0]
            r = go.SphericalCellGrid(n, 60.).centers
            sl = vi.SuppliedSlice(ss.ScalarFields(*vals, np.full(n, .179055), np.full(n, .00341913)), r.copy(),
                                  meta["H0"], meta["edge_tol"], f"W1 static phase, dr {dr}")
            state = V55ProductionKernel().initialize(resolution=n, r_max=60., outer_edge="open_edge",
                                                     eight_node_areal=True, initial_structure="supplied_slice",
                                                     supplied_slice=sl)
        finally:
            vi._eight_node_origin_integral = orig
        g = state.geometry
        st = assemble_total_stress_energy(state.grid, g, state.scalars, state.matter)
        H = vac.constraints(state.grid, g)["hamiltonian"] - 16.0 * np.pi * st.rho
        chk = None
        f = runs / "sp" / f"{dr}-0.000375.npz"
        if f.exists():
            with np.load(f) as z:
                chk = bool(np.array_equal(z["H"][0], H))
                report.setdefault("stored_H_bitwise", {})[str(dr)] = chk
                if not chk:
                    report.setdefault("stored_H_maxdiff", {})[str(dr)] = float(np.max(np.abs(z["H"][0] - H)))
        k = N0[dr]
        tag = f"dr{dr}_"
        for name in ("a", "b", "X", "Aa", "K", "Lambda"):
            saved[tag + name] = np.asarray(getattr(g, name))[:k]
        saved[tag + "rho"] = np.asarray(st.rho)[:k]
        saved[tag + "H"] = H[:k]
        saved[tag + "src_last"] = srcs[-1][:k]
        saved[tag + "src_prev"] = srcs[-2][:k]
        saved[tag + "r"] = np.asarray(r)[:k]
        report.setdefault("areal_loop_calls", {})[str(dr)] = len(srcs)
    np.savez_compressed(out, report=json.dumps(report), **saved)
    print(json.dumps(report))


if __name__ == "__main__":
    main(*sys.argv[1:4])
