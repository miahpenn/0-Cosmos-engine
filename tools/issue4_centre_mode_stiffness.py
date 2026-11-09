"""Issue #4: rate at which the Aa RHS (L2) responds to an irregular (a-b) component at the origin.
Perturb Aa by eps*g(r), g a narrow even bump with g(0) != 0 (width 3h), and by eps*r^2-weighted
regular bump; compare -d(L2 Aa)/d(Aa) at cell 0 and its scaling with h."""
import sys
import numpy as np
sys.path.insert(0, ".")
from engine.production_kernel import V55ProductionKernel
for N in (60, 120, 240, 480):
    k = V55ProductionKernel(); s = k.initialize(resolution=N, r_max=40.0)
    r = np.asarray(s.grid.centers); h = float(s.grid.dr)
    def l2(geo):
        g = geo.copy(); g.alpha = s.geometry.alpha
        return k._geometry_stage_terms(s.grid, g, s.scalars, s.matter, None)["primary_l2"]["Aa"]
    base = l2(s.geometry); eps = 1e-8
    out = []
    for name, bump in (("irregular g(0)=1", np.exp(-(r / (3 * h)) ** 2)),
                       ("regular r^2 bump", (r / (3 * h)) ** 2 * np.exp(-(r / (3 * h)) ** 2))):
        # traceless metric perturbation keeping a*b^2 fixed: da = 2 eps bump, db = -eps bump
        g = s.geometry.copy(); g.a = g.a * (1 + 2 * eps * bump); g.b = g.b * (1 - eps * bump)
        resp = (l2(g) - base) / eps
        out.append(f"{name}: response cell0 {resp[0]:+.3e} (x h^2 {resp[0] * h * h:+.3f})")
    print(f"N={N} h={h:.4f}: " + " | ".join(out), flush=True)
