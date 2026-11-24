"""Input checks and fault cases for the corrected W1 slice (W1 derivation AMENDMENTS 3 and 4). No simulation.

Runs the C9 gate script research/w1_c9_g05p/run.py with --static on: the corrected slice; the admitted
(uncorrected) slice; the build with the v^2 P_st piece removed; the build with the wrong-sign phase; and the
corrected slice with its D or S tail reversed (static-frame momentum P_st = P - H0 r f_r negated beyond r ~ 12).
Expected (A3 5, A4 3): corrected passes everything; every fault fails 'f; reversed tails also fail 'e.
Usage: python3 research/w1_static_slice/checks.py BUILD_DIR ADMITTED_NPZ MEMBERS_JSON OUTDIR
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

bdir, admitted, members, out = map(Path, sys.argv[1:5])
out.mkdir(parents=True, exist_ok=True)
build = bdir / "build.json"
go, vac, _ = ad.vendor_modules()
z = dict(np.load(bdir / "static" / "native.npz"))
cases = {"corrected": bdir / "static" / "native.npz", "uncorrected (admitted)": admitted,
         "v^2 piece missing": bdir / "fault_no_v2" / "native.npz",
         "wrong-sign phase": bdir / "fault_wrong_sign" / "native.npz"}
rev = {"D tail reversed": dict(z), "S tail reversed": dict(z)}
for dr, n in ((0.025, 2400), (0.05, 1200)):
    key = "dr%g_" % dr
    grid = go.SphericalCellGrid(n, 60.0)
    r = np.asarray(grid.centers)
    H0 = float(vi.build_initial_data(grid, vac.flat_state, initial_structure="none", eight_node_areal=True).H0)
    sw = 0.5 * (1 + np.tanh((r - 12.0) / 2))
    for f, p, case in (("S", "PS", "S tail reversed"), ("D", "PD", "D tail reversed")):
        F, P = np.asarray(z[key + f], dtype=np.float64), np.asarray(z[key + p], dtype=np.float64)
        pst = P - H0 * r * grid.cell_derivative_fourth(F, parity=1)
        rev[case][key + p] = P - 2 * sw * pst
rep = {}
with tempfile.TemporaryDirectory() as td:
    for case, arrs in rev.items():
        d = Path(td) / case.replace(" ", "_")
        d.mkdir()
        np.savez(d / "native.npz", **arrs)
        cases[case] = d / "native.npz"
    for case, npz in cases.items():
        d = out / case.split(" (")[0].replace(" ", "_").replace("^", "")
        subprocess.run([sys.executable, "research/w1_c9_g05p/run.py", str(npz), str(members), str(d), "--static", str(build)],
                       check=True, capture_output=True)
        R = json.loads((d / "report.json").read_text())
        rep[case] = {"outcome": R["outcome"], "grids": {k: {
            "G0.1-4": all(v[g] for g in ("G0.1_PASS", "G0.2_PASS", "G0.3_PASS", "G0.4_PASS")),
            "entry": v["entry_reading"], **{c: v["G0.5p"][c + "_PASS"] for c in "abcdef"},
            "P_tail_rel": v["G0.5p"]["P_tail_rel"], "m_D": v["G0.5p"]["D"]["m"], "m_S": v["G0.5p"]["S"]["m"],
            "PD_last5": v["G0.5p"]["d"]["PD_last5"], "PD_bound": v["G0.5p"]["d"]["PD_bound"],
            "e_D": v["G0.5p"]["e"]["F_win_D/P_tail"], "e_S": v["G0.5p"]["e"]["F_win_S/P_S_tail"],
            "f_D": v["G0.5p"]["f"]["F_st_D/P_tail"], "f_S": v["G0.5p"]["f"]["F_st_S/P_S_tail"]}
            for k, v in R["grids"].items()}}
(out / "checks.json").write_text(json.dumps(rep, indent=1, default=float))
for case, c in rep.items():
    for k, v in c["grids"].items():
        print(f"{case:24s} {k:7s} G0.1-4={v['G0.1-4']} entry={v['entry']} " + " ".join(f"{x}={v[x]}" for x in "abcdef")
              + f" Prel={v['P_tail_rel']:+.4f} mD={v['m_D']:.1e} mS={v['m_S']:.1e} e={v['e_D']:.3f}/{v['e_S']:.3f}"
              f" f={v['f_D']:.4f}/{v['f_S']:.4f}")
