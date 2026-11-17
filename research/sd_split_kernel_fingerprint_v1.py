"""RESEARCH ONLY. Kernel fingerprints (state hashes after a few steps) for the split-follows-phi PR: run on main as
`python research/sd_split_kernel_fingerprint_v1.py main` and on the branch as `... ref` (the labelled references:
imported_dust_reference(), hand_floor_reference()); equal hashes = bit-identical."""
import hashlib, sys, contextlib
sys.path.insert(0, ".")
import numpy as np
import engine.scalar_system as ss
from engine.production_kernel import V55ProductionKernel
ref = sys.argv[1] == "ref"
def h(st):
    m = hashlib.sha256()
    for k in ("a","b","X","alpha","Aa","K","Lambda"): m.update(np.asarray(getattr(st.geometry,k)).tobytes())
    for k in ("S","PS","D","PD","phi","Pi"): m.update(np.asarray(getattr(st.scalars,k)).tobytes())
    for sp in ("dark_matter","baryons","radiation"):
        for f in ("rest","energy_t","momentum_r"): m.update(np.asarray(getattr(getattr(st.matter,sp),f)).tobytes())
    return m.hexdigest()[:16]
def run(label, ctx, kw, n):
    k = V55ProductionKernel()
    with ctx:
        st = k.initialize(**kw); dt = 0.06*float(st.grid.dr)
        for _ in range(n): st = k.step(st, dt)
    print(f"{label:45s} {h(st)}")
dust = ss.imported_dust_reference() if ref else contextlib.nullcontext()
run("background wall N20 x40 (dust reference)", dust, dict(resolution=20, r_max=40.0), 40)
dust = ss.imported_dust_reference() if ref else contextlib.nullcontext()
run("rescaled pocket wall N24 x6 (dust reference)", dust, dict(resolution=24, r_max=24.0, initial_structure="rescaled_pocket", pocket_x=2.0), 6)
run("hand-floor reference open edge N20 x10", ss.hand_floor_reference(), dict(resolution=20, r_max=40.0), 10)
