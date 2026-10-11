"""One-step comparison of R=80 and R=160 production initial states on shared cells.

Measures where the first step's field changes differ, to test whether the centre-first
divergence is already present after one step from the initial data. Measurement only.
"""
import sys
import pathlib as _p

import numpy as np

sys.path.insert(0, str(_p.Path(__file__).resolve().parents[2]))
from engine.production_kernel import V55ProductionKernel

k = V55ProductionKernel()


def init(r_max, N):
    return k.initialize(
        resolution=N, r_max=r_max, amplitude=0.01, width=7.0,
        D_amplitude=1e-4, include_radiation=True
    )


s80 = init(80.0, 160)
s160 = init(160.0, 320)
n = 160
dr = s80.grid.dr
dt = 0.0075 * dr


def fields(st):
    G = st.geometry
    S = st.scalars
    return {
        "alpha": np.asarray(G.alpha)[:n],
        "a": np.asarray(G.a)[:n],
        "b": np.asarray(G.b)[:n],
        "X": np.asarray(G.X)[:n],
        "S": np.asarray(S.S)[:n],
        "PS": np.asarray(S.PS)[:n],
        "D": np.asarray(S.D)[:n],
        "PD": np.asarray(S.PD)[:n],
        "phi": np.asarray(S.phi)[:n],
        "Pi": np.asarray(S.Pi)[:n],
    }


f80, f160 = fields(s80), fields(s160)
print("initial max difference on shared cells:")
for key in f80:
    d = np.abs(f80[key] - f160[key])
    print(f"  {key:6s} max={d.max():.3e} at cell {int(np.argmax(d))}")

a80 = k.step(s80, dt)
a160 = k.step(s160, dt)
g80, g160 = fields(a80), fields(a160)
print(f"after one step (dt={dt:.5f}), difference between domains:")
for key in g80:
    d = np.abs(g80[key] - g160[key])
    print(f"  {key:6s} centre cells0-2 max={d[:3].max():.3e}   all max={d.max():.3e} at cell {int(np.argmax(d))}")
print("change over the step, centre cell 0 (R=80, R=160):")
for key in ("S", "D", "phi", "a", "alpha"):
    print(f"  {key:6s} dR80={g80[key][0]-f80[key][0]:+.4e}  dR160={g160[key][0]-f160[key][0]:+.4e}")
