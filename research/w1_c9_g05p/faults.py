"""Fault tests for the C9 W1 gate script (research/w1_c9_g05p/run.py), stored native.npz only, no simulation.

Each case writes a modified copy of native.npz and runs the gate script itself on it, then records which of
G0.5'a-e pass. Reversed tail (AMENDMENT 3, 9.5; Aster's case): the field's rest-frame momentum P_rest = P - H0 r f_r
is negated beyond r ~ 12 (smooth switch), the producer's chart term kept. D and S are reversed separately.
Usage: python3 research/w1_c9_g05p/faults.py native.npz members.json OUTDIR
"""
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, os.getcwd())
import numpy as np  # noqa: E402

from engine import v55_initial as vi, v55_pirk_adapter as ad  # noqa: E402

npz, members, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
out.mkdir(parents=True, exist_ok=True)
go, vac, _ = ad.vendor_modules()
z = dict(np.load(npz))
mod = {"baseline": dict(z), "S tail reversed": dict(z), "D tail reversed": dict(z)}
for dr, n in ((0.025, 2400), (0.05, 1200)):
    key = "dr%g_" % dr
    grid = go.SphericalCellGrid(n, 60.0)
    r = np.asarray(grid.centers)
    H0 = float(vi.build_initial_data(grid, vac.flat_state, initial_structure="none", eight_node_areal=True).H0)
    sw = 0.5 * (1 + np.tanh((r - 12.0) / 2))
    for f, p, case in (("S", "PS", "S tail reversed"), ("D", "PD", "D tail reversed")):
        F, P = np.asarray(z[key + f], dtype=np.float64), np.asarray(z[key + p], dtype=np.float64)
        rest = P - H0 * r * grid.cell_derivative_fourth(F, parity=1)
        mod[case][key + p] = P - 2 * sw * rest
rep = {}
with tempfile.TemporaryDirectory() as td:
    for case, arrs in mod.items():
        d = Path(td) / case.replace(" ", "_")
        d.mkdir()
        np.savez(d / "native.npz", **arrs)
        subprocess.run([sys.executable, "research/w1_c9_g05p/run.py", str(d / "native.npz"), str(members), str(d)],
                       check=True, capture_output=True)
        g = json.loads((d / "report.json").read_text())["grids"]
        rep[case] = {k: {"entry": v["entry"]["engine"][:60], "entry_reading": v["entry_reading"],
                         **{c: v["G0.5p"][c + "_PASS"] for c in "abcde"}, "e_values": v["G0.5p"]["e"]}
                     for k, v in g.items()}
(out / "faults.json").write_text(json.dumps(rep, indent=1, default=float))
for case, g in rep.items():
    for k, v in g.items():
        print(f"{case:16s} {k:7s} a={v['a']} b={v['b']} c={v['c']} d={v['d']} e={v['e']} entry={v['entry_reading']} "
              f"FD/P={v['e_values']['F_win_D/P_tail']:.4f} FS/P={v['e_values']['F_win_S/P_S_tail']:.4f}")
