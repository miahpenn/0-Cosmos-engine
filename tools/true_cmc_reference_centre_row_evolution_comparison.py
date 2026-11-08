"""Old vs repaired centre row of engine/true_cmc_reference.cmc_lapse over real step_true_cmc steps.

Run from a checkout carrying the repaired row. The old row (alpha[0] - alpha[1] = 0) is
recovered by loading engine/true_cmc_reference.py from git ref OLD_REF and swapping its
cmc_lapse in; nothing else differs. Reference initial data q.make_initial(g, 0.01, 7.0),
r_max = 40, dt = 0.06 dr, to t = 0.60, N = 40, 80, 160. Records H, M, C from
q.constraint at cells 0-2 and in 5 <= r <= 10, and the lapse at cells 0-1.
Usage: python <this file> OLD_REF output.json
"""
import importlib.util, json, subprocess, sys, tempfile
import numpy as np
sys.path.insert(0, ".")
from engine import reference_pirk_unified as q
from engine import true_cmc_reference as tcr

old_src = subprocess.check_output(["git", "show", f"{sys.argv[1]}:engine/true_cmc_reference.py"], text=True)
with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as f:
    f.write(old_src)
spec = importlib.util.spec_from_file_location("tcr_old", f.name); old = importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
NEW = tcr.cmc_lapse

def measure(g, s, fields):
    H, M, C, _ = q.constraint(g, s, fields); r = g.centers; ann = (r >= 5) & (r <= 10)
    out = {f"{k}_{w}": float(np.max(np.abs(v[sl]))) for k, v in (("H", H), ("M", M), ("C", C)) for w, sl in (("c02", slice(0, 3)), ("ann", ann))}
    out["alpha0"], out["alpha1_minus_alpha0"] = float(s.alpha[0]), float(s.alpha[1] - s.alpha[0])
    return out

def run(n, lapse):
    tcr.cmc_lapse = lapse
    try:
        g = q.Grid(n, 40.0); s, fields = q.make_initial(g, .01, 7.); dt = .06 * g.dr
        s.alpha = lapse(g, s, fields)[0]; s.beta.fill(0.0); s.B.fill(0.0)
        for _ in range(int(round(0.6 / dt))):
            s, fields = tcr.step_true_cmc(g, s, fields, dt)
        return measure(g, s, fields), s
    finally:
        tcr.cmc_lapse = NEW

res = {}
for n in (40, 80, 160):
    (mo, so), (mn, sn) = run(n, old.cmc_lapse), run(n, NEW)
    far = q.Grid(n, 40.0).centers >= 5
    res[n] = {"old": mo, "new": mn, "max_field_diff_r_ge_5": {k: float(np.max(np.abs(getattr(so, k) - getattr(sn, k))[far])) for k in ("a", "b", "X", "K", "Aa", "alpha")}}
    for lab, m in (("old", mo), ("new", mn)):
        print(f"N={n:>3} {lab}: M c0-2 {m['M_c02']:.2e} ann {m['M_ann']:.1e} | H c0-2 {m['H_c02']:.2e} ann {m['H_ann']:.1e} | C c0-2 {m['C_c02']:.2e} | a1-a0 {m['alpha1_minus_alpha0']:.2e}", flush=True)
    print(f"N={n:>3} r>=5 diffs: " + " ".join(f"{k} {v:.1e}" for k, v in res[n]["max_field_diff_r_ge_5"].items()), flush=True)
json.dump(res, open(sys.argv[2], "w"), indent=1)
