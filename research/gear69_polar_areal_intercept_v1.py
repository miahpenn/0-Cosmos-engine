"""RESEARCH ONLY - imports no engine code. GEAR-87 / AMENDMENT-3-GEAR-87.
Dense small-amplitude run of gear69_polar_areal_nonlinear_v1 and the phi_c^2 -> 0
intercept of (Gamma_tau - Gamma_num)/gamma/phi_c^2, to compare with 2 w(0) - alpha2(0)
= 4.58260 (G=1, R=10, mu=1). Usage: python research/gear69_polar_areal_intercept_v1.py N"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gear69_polar_areal_nonlinear_v1 as pe

N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
ref = pe.run(N, 1e-5, 3.0, every=4)
gnum = np.median(ref[:, 2] / ref[:, 1])
t, pc, Pc = pe.run(N, 0.001, 5.2, every=2).T
x = pc ** 2
y = (Pc / pc - gnum) / pe.gam / x
sel = (x >= 5e-4) & (x <= 5e-3)
for deg in (1, 2, 3):
    c = np.polyfit(x[sel], y[sel], deg)[::-1]
    print("N", N, "deg", deg, "intercept %.5f slope %.3f diff from 4.58260 %.2e"
          % (c[0], c[1], c[0] - 4.5826))
