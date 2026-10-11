"""Initial-state reference check: production lapse at r=79.75 and at the centre for three
domain sizes at fixed dr=0.5. Measurement only; no trajectory."""
import sys
import pathlib as _p

import numpy as np

sys.path.insert(0, str(_p.Path(__file__).resolve().parents[2]))
from engine.production_kernel import V55ProductionKernel

k = V55ProductionKernel()
rows = {}
for rmax, N in ((80.0, 160), (160.0, 320), (240.0, 480), (320.0, 640)):
    st = k.initialize(
        resolution=N, r_max=rmax, amplitude=0.01, width=7.0,
        D_amplitude=1e-4, include_radiation=True
    )
    r = np.asarray(st.grid.centers)
    al = np.asarray(st.geometry.alpha)
    j = int(np.argmin(np.abs(r - 79.75)))
    rows[rmax] = (al, r, j)
    print(
        f"rmax={rmax:6.1f} N={N}: alpha(r=79.75)={al[j]:.9f} "
        f"(cell {j})  alpha(centre)={al[0]:.9f}  "
        f"alpha(outer cell)={al[-1]:.9f}",
        flush=True,
    )
a80, r80, j80 = rows[80.0]
a160, r160, j160 = rows[160.0]
a240, r240, j240 = rows[240.0]
a320, r320, j320 = rows[320.0]
print(
    f"difference at r=79.75: R80-R160 = {a80[j80]-a160[j160]:+.6e};  "
    f"R160-R240 = {a160[j160]-a240[j240]:+.6e};  "
    f"R240-R320 = {a240[j240]-a320[j320]:+.6e}"
)
print(
    f"centre: R80-R160 = {a80[0]-a160[0]:+.6e};  "
    f"R160-R240 = {a160[0]-a240[0]:+.6e};  "
    f"R240-R320 = {a240[0]-a320[0]:+.6e}"
)
