"""Synthetic checks of the 60-digit transcription (python3 -m research.control_tests.check_mp60).

Uses only fields built here in exact arithmetic; no engine call and no stored run data, so no gate is read.
"""
import mpmath as mp

from research.control_tests import mp60 as M


def field(g, f):
    return [M.E(f(r.v)) for r in g.r]


def maxabs(xs, skip_edge=2):
    xs = xs[:len(xs) - skip_edge] if skip_edge else xs
    return max(abs(x.v) for x in xs if x is not None)


def main():
    g = M.Grid(40, 0.05)
    # 1. stencils exact on even quartics, everywhere (edge rows are exact 5-point Lagrange weights)
    p = lambda r: 1 + 3 * r**2 - 2 * r**4
    for order, exact in ((1, lambda r: 6 * r - 8 * r**3), (2, lambda r: 6 - 24 * r**2)):
        d = g.d(field(g, p), 1, order)
        err = max(abs(d[i].v - exact(g.r[i].v)) for i in range(g.n))
        assert err < mp.mpf(10) ** -45, (order, err)
    # 2. solve operator and Delta_4 agree exactly on an even quadratic alpha (both exact there)
    L = mp.mpf(30)
    w = lambda r: r**2 / (4 * L**2)
    s = {"a": field(g, lambda r: 1), "b": field(g, lambda r: 1), "X": field(g, lambda r: 1 - w(r)),
         "Lambda": field(g, lambda r: 0), "alpha": field(g, lambda r: mp.mpf("0.5") + mp.mpf("0.2") * r**2),
         "Aa": field(g, lambda r: 0), "K": field(g, lambda r: -3 / L)}
    d24 = M.delta2_minus_delta4(g, s)
    assert maxabs(d24) < mp.mpf(10) ** -45, maxabs(d24)
    # 3. 1b slice (X = 1 - r^2/4L^2, K = -3/L, a = b = 1): vacuum H and M vanish in exact arithmetic
    H = M.hamiltonian_vac(g, s)
    assert maxabs(H) < mp.mpf(10) ** -45, maxabs(H)
    mt = M.momentum_terms(g, s)
    assert maxabs(mt["M_vac"]) < mp.mpf(10) ** -45
    # 4. the rounding magnitude is positive and the depth is counted
    assert all(x.a > 0 and x.n > 10 for x in H[:-2])
    # 5. a quartic alpha gives the 3-point/2-point truncation, (dr^2/12) alpha'''' + (dr^2/6) c alpha''' to leading order
    s4 = dict(s, alpha=field(g, lambda r: 1 + r**4))
    d24 = M.delta2_minus_delta4(g, s4)
    i = 10
    r = g.r[i].v; h = g.dr.v
    inv = (1 - w(r))**2
    c = 2 * (r / (4 * L**2)) / (1 - w(r)) + 2 / r   # -X'/X + 2/r
    lead = inv * (h**2 / 12 * 24 + h**2 / 6 * c * 24 * r)
    assert abs(d24[i].v / lead - 1) < mp.mpf("1e-6"), d24[i].v / lead
    # 6. explicit X RHS on the 1b slice is alpha K X / 3; Kdot projection of a uniform l3 (alpha const, Aa = 0,
    #    vacuum) is alpha K^2/3 within its own bound
    sc = dict(s, alpha=field(g, lambda r: mp.mpf("0.7")))
    ex = M.explicit_rhs(g, sc)
    assert max(abs(x.v - mp.mpf("0.7") * (-3 / L) * X.v / 3) for x, X in zip(ex["X"], sc["X"])) < mp.mpf(10) ** -45
    m0 = {k: field(g, lambda r: 0) for k in ("rho", "pr", "pt", "j")}
    kv, kb = M.kdot_projection(g, sc, m0)
    assert abs(kv - mp.mpf("0.7") * (3 / L) ** 2 / 3) < mp.mpf(10) ** -45 and kb > 0
    # 7. Lambda l3 vanishes with Aa = 0 and j = 0; K' vanishes for constant K
    assert maxabs(M.lambda_l3(g, sc, m0), 0) == 0
    assert maxabs(M.k_prime(g, sc), 0) < mp.mpf(10) ** -45
    print("mp60 synthetic checks pass")


if __name__ == "__main__":
    main()
