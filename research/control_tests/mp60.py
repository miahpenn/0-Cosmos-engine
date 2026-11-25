"""60-digit transcription of the engine formulas used by the control tests (Amendment 3 section 5, A6 B, A9 section 7).

Runs under plain python3 with mpmath (no engine import). Inputs are the binary64 arrays a research run stores,
converted exactly to mpf; dr is the binary64 dr. Every operation is done in mp in the engine's own order, and each
value carries (v, a, n): v the 60-digit value, a the same expression evaluated with every term made non-negative
(the summand magnitude), n the floating-point operation depth. The rounding bound of the engine's binary64
evaluation is then gamma_n * a  (A6 B5 / A13 section 1: |stored - 60 digits| <= gamma_m Sum_k |s_k|).

Not widened (inputs): the lapse solve (alpha is the engine's), the matter totals rho, pr, pt, j (binary64 inputs
for Tests 1a/1b, where dust is 1e-150), and anything stored before the call.

Edge cells: the engine's last two rows use one-sided weights from a binary64 5x5 solve; here they are solved
exactly. Those cells are readings only (A6 B0, A9 section 1) and are returned but flagged.
"""
from fractions import Fraction
from math import factorial as _fact

import mpmath as mp
import numpy as np

mp.mp.dps = 60
U = mp.mpf(2) ** -53  # unit roundoff


def gamma(n):
    return n * U / (1 - n * U)


class E:
    __slots__ = ("v", "a", "n")

    def __init__(self, v, a=None, n=0):
        self.v = mp.mpf(v)
        self.a = abs(self.v) if a is None else a
        self.n = n

    @staticmethod
    def c(x):
        return x if isinstance(x, E) else E(x)

    def __add__(s, o):
        o = E.c(o); return E(s.v + o.v, s.a + o.a, max(s.n, o.n) + 1)

    __radd__ = __add__

    def __sub__(s, o):
        o = E.c(o); return E(s.v - o.v, s.a + o.a, max(s.n, o.n) + 1)

    def __rsub__(s, o):
        return E.c(o) - s

    def __mul__(s, o):
        o = E.c(o); return E(s.v * o.v, s.a * o.a, s.n + o.n + 1)

    __rmul__ = __mul__

    def __truediv__(s, o):
        o = E.c(o); return E(s.v / o.v, s.a / abs(o.v), s.n + o.n + 1)

    def __rtruediv__(s, o):
        return E.c(o) / s

    def __neg__(s):
        return E(-s.v, s.a, s.n)

    def __pow__(s, k):
        assert k in (2, 3)  # numpy x**2 is x*x; x**3 goes through pow (within one rounding of x*x*x)
        return s * s if k == 2 else s * s * s


def arr(x):
    """binary64 array -> list of exact E inputs."""
    return [E(mp.mpf(float(v))) for v in np.asarray(x, dtype=float)]


def const(x, n):
    return [E(x)] * n


def vals(xs):
    return np.array([float(x.v) for x in xs])


def bound(xs):
    """gamma_n * a per cell: the forward-error bound of the binary64 evaluation."""
    return np.array([float(gamma(x.n) * x.a) for x in xs])


def mpvals(xs):
    return [x.v for x in xs]


PI = E(mp.pi, None, 1)  # exact pi; the engine rounds it once


class Grid:
    def __init__(self, n, dr):
        self.n = n
        self.dr = E(mp.mpf(float(dr)))
        self.r = [E((mp.mpf(i) + mp.mpf(1) / 2) * self.dr.v) for i in range(n)]  # exact (i+1/2) dr
        self.dr2 = self.dr * self.dr
        self._edge = {}

    def edge_weights(self, i, order):
        key = (i, order)
        if key not in self._edge:
            offs = [Fraction(j - i) for j in range(self.n - 5, self.n)]
            M = [[o ** p for o in offs] + [Fraction(1 if p == order else 0) * _fact(order)] for p in range(5)]
            for c in range(5):  # exact Gauss-Jordan in rationals
                piv = next(rw for rw in range(c, 5) if M[rw][c] != 0)
                M[c], M[piv] = M[piv], M[c]
                M[c] = [x / M[c][c] for x in M[c]]
                for rw in range(5):
                    if rw != c and M[rw][c] != 0:
                        f = M[rw][c]
                        M[rw] = [x - f * y for x, y in zip(M[rw], M[c])]
            self._edge[key] = [mp.mpf(M[k][5].numerator) / M[k][5].denominator for k in range(5)]
        return self._edge[key]

    def d(self, u, parity, order):
        """cell_derivative_fourth / cell_second_derivative_fourth, engine operation order."""
        n = self.n
        ext = [parity * u[1], parity * u[0]] + list(u) + [None, None]
        out = [None] * n
        for i in range(n - 2):
            s = ext[i:i + 5]
            if order == 1:
                out[i] = (s[0] - 8.0 * s[1] + 8.0 * s[3] - s[4]) / (12.0 * self.dr)
            else:
                out[i] = (-s[0] + 16.0 * s[1] - 30.0 * s[2] + 16.0 * s[3] - s[4]) / (12.0 * self.dr2)
        for i in (n - 2, n - 1):
            w = self.edge_weights(i, order)
            acc = E(0)
            for k, j in enumerate(range(n - 5, n)):
                acc = acc + E(w[k]) * u[j]
            out[i] = acc / (self.dr if order == 1 else self.dr2)
        return out


def _ew(f, *cols):
    return [f(*xs) for xs in zip(*cols)]


def geometry_terms(g, s):
    """vendor vacuum_bssn.geometry_terms (R, Rr_r, DrDr_alpha, lap_alpha); s = dict of E lists."""
    a, b, X, Lam, al, r = s["a"], s["b"], s["X"], s["Lambda"], s["alpha"], g.r
    ap, bp = g.d(a, 1, 1), g.d(b, 1, 1)
    app, bpp = g.d(a, 1, 2), g.d(b, 1, 2)
    Xp, Xpp = g.d(X, 1, 1), g.d(X, 1, 2)
    chip = _ew(lambda xp, x: -0.5 * xp / x, Xp, X)
    chipp = _ew(lambda xpp, xp, x: -0.5 * (xpp / x - (xp / x) ** 2), Xpp, Xp, X)
    Lp = g.d(Lam, -1, 1)
    ml = _ew(lambda a_, b_, r_: (1.0 - a_ / b_) / r_ ** 2, a, b, r)
    inv = _ew(lambda x, a_: x ** 2 / a_, X, a)
    rb = _ew(lambda app_, a_, Lp_, ap_, bp_, b_, L_, r_, ml_, cpp, cp:
             0.5 * app_ / a_ - a_ * Lp_ - 0.75 * (ap_ / a_) ** 2 + 0.5 * (bp_ / b_) ** 2 - 0.5 * L_ * ap_
             + ap_ / (r_ * b_) + 2.0 * ml_ * (1.0 + r_ * bp_ / b_) + 4.0 * cpp - 2.0 * cp * (ap_ / a_ - bp_ / b_ - 2.0 / r_),
             app, a, Lp, ap, bp, b, Lam, r, ml, chipp, chip)
    sb = _ew(lambda app_, a_, bpp_, b_, Lp_, ap_, bp_, r_, ml_, cpp, cp:
             0.5 * app_ / a_ + bpp_ / b_ - a_ * Lp_ - (ap_ / a_) ** 2 + 0.5 * (bp_ / b_) ** 2
             + 2.0 * (3.0 - a_ / b_) * bp_ / (r_ * b_) + 4.0 * ml_ + 8.0 * (cpp + cp ** 2)
             - 8.0 * cp * (0.5 * ap_ / a_ - bp_ / b_ - 2.0 / r_),
             app, a, bpp, b, Lp, ap, bp, r, ml, chipp, chip)
    alp, alpp = g.d(al, 1, 1), g.d(al, 1, 2)
    out = {
        "R": _ew(lambda i, x: -i * x, inv, sb),
        "Rr_r": _ew(lambda i, x: -i * x, inv, rb),
        "DrDr_alpha": _ew(lambda i, pp, p, ap_, a_, cp: i * (pp - p * (0.5 * ap_ / a_ + 2.0 * cp)),
                          inv, alpp, alp, ap, a, chip),
        "lap_alpha": _ew(lambda i, pp, p, ap_, a_, bp_, b_, cp, r_: i * (pp - p * (0.5 * ap_ / a_ - bp_ / b_ - 2.0 * cp - 2.0 / r_)),
                         inv, alpp, alp, ap, a, bp, b, chip, r),
        "chip": chip, "ap": ap, "bp": bp, "inv": inv,
    }
    return out


def momentum_terms(g, s, gt=None):
    """The A4 R-stage split of the vacuum momentum: D_r Aa, 3Aa/r, -(2/3) D_r K, remaining metric couplings.
    (Aa - Ab)(2/r + b'/b) with Ab = -Aa/2 is split as 3Aa/r + 1.5 Aa b'/b, in the engine's grouping order for M."""
    gt = gt or geometry_terms(g, s)
    Aa, K, b, r = s["Aa"], s["K"], s["b"], g.r
    Aap, Kp = g.d(Aa, 1, 1), g.d(K, 1, 1)
    out = {
        "DrAa": Aap,
        "3Aa_r": _ew(lambda A, r_: 1.5 * A * (2.0 / r_), Aa, r),
        "DrK": _ew(lambda kp: -(E(2) / 3.0) * kp, Kp),
        "metric": _ew(lambda A, cp, bp_, b_: 6.0 * A * cp + 1.5 * A * (bp_ / b_), Aa, gt["chip"], gt["bp"], b),
    }
    # the engine's own grouping, for the independence check against vendor momentum_constraint
    out["M_vac"] = _ew(lambda dA, A, kp, cp, r_, bp_, b_: dA - (E(2) / 3.0) * kp + 6.0 * A * cp + (A - (-0.5 * A)) * (2.0 / r_ + bp_ / b_),
                       Aap, Aa, Kp, gt["chip"], r, gt["bp"], b)
    return out


def hamiltonian_vac(g, s, gt=None):
    gt = gt or geometry_terms(g, s)
    return _ew(lambda R, A, K: R - (A ** 2 + 2.0 * (-0.5 * A) ** 2) + (E(2) / 3.0) * K ** 2, gt["R"], s["Aa"], s["K"])


def rhs_primary(g, s, m, gt=None):
    """primary_l2 / primary_l3 (KO off, beta 0) with the adapter's matter terms; m = dict of E lists rho, pr, pt."""
    gt = gt or geometry_terms(g, s)
    al, Aa, K = s["alpha"], s["Aa"], s["K"]
    l2 = {"Aa": _ew(lambda D, L, al_, Rr, R: -(D - L / 3.0) + al_ * (Rr - R / 3.0),
                    gt["DrDr_alpha"], gt["lap_alpha"], al, gt["Rr_r"], gt["R"]),
          "K": _ew(lambda L: -L, gt["lap_alpha"])}
    l3 = {"Aa": _ew(lambda al_, k, A, pr, pt: al_ * k * A - (16.0 * PI / 3.0) * al_ * (pr - pt), al, K, Aa, m["pr"], m["pt"]),
          "K": _ew(lambda al_, A, k, rho, pr, pt: al_ * (A ** 2 + 2.0 * (-0.5 * A) ** 2 + k ** 2 / 3.0)
                   + 4.0 * PI * al_ * (rho + pr + 2.0 * pt), al, Aa, K, m["rho"], m["pr"], m["pt"])}
    return l2, l3


def lambda_l2(g, s, lambda_m=2.0, gt=None):
    """vendor lambda_l2_rhs with beta = 0 (shift terms vanish identically)."""
    gt = gt or geometry_terms(g, s)
    M = momentum_terms(g, s, gt)["M_vac"]
    alp = g.d(s["alpha"], 1, 1)
    Aap = g.d(s["Aa"], 1, 1)
    return _ew(lambda A, ap_, al_, Ap, a_, r_, b_, M_: -2.0 * (A * ap_ + al_ * Ap) / a_ - 4.0 * al_ * (A - (-0.5 * A)) / (r_ * b_)
               + lambda_m * al_ * M_ / a_, s["Aa"], alp, s["alpha"], Aap, s["a"], g.r, s["b"], M)


def solve_operator(g, s):
    """The CMC solve's own operator on alpha (cmc_gauge: 3-point alpha'', 2-point alpha', even ghost at cell 0),
    inv * (alpha'' + c alpha'), with c from the 4th-order d1 of a, b, X. Last cell: None (alpha = 1 row)."""
    a, b, X, al = s["a"], s["b"], s["X"], s["alpha"]
    ap, bp, Xp = g.d(a, 1, 1), g.d(b, 1, 1), g.d(X, 1, 1)
    c = _ew(lambda ap_, a_, bp_, b_, xp, x, r_: -0.5 * ap_ / a_ + bp_ / b_ - xp / x + 2.0 / r_, ap, a, bp, b, Xp, X, g.r)
    n = g.n
    out = [None] * n
    for i in range(n - 1):
        lo = al[i - 1] if i > 0 else al[0]
        hi = al[i + 1]
        inv = X[i] ** 2 / a[i]
        out[i] = inv * ((lo - 2.0 * al[i] + hi) / g.dr2 + c[i] * (hi - lo) / (2.0 * g.dr))
    return out


def delta2_minus_delta4(g, s, gt=None):
    gt = gt or geometry_terms(g, s)
    L2 = solve_operator(g, s)
    return [None if x is None else x - y for x, y in zip(L2, gt["lap_alpha"])]


def state_from(npz, prefix, n, extra=()):
    """Read one stored state (binary64) as exact E lists."""
    keys = ("a", "b", "X", "Aa", "K", "Lambda", "alpha") + tuple(extra)
    return {k: arr(npz[prefix + k]) for k in keys}


def explicit_rhs(g, s):
    """native CMC explicit a, b, X RHS with beta = 0 (the beta and div_beta terms are exact zeros in binary64)."""
    a, b, X, al, Aa, K = s["a"], s["b"], s["X"], s["alpha"], s["Aa"], s["K"]
    return {"a": _ew(lambda al_, a_, A: -(2.0 * al_ * a_ * A), al, a, Aa),
            "b": _ew(lambda al_, b_, A: -(2.0 * al_ * b_ * (-0.5 * A)), al, b, Aa),
            "X": _ew(lambda x, al_, k: x * (al_ * k) / 3.0, X, al, K)}


def lambda_l3(g, s, m, lambda_m=2.0):
    """adapter lambda_l3_with_matter, beta = 0, KO off: 2 alpha Aa Lambda - 8 pi lambda_m alpha j / a."""
    return _ew(lambda al_, A, L, j, a_: 2.0 * al_ * A * L - 8.0 * PI * lambda_m * al_ * j / a_,
               s["alpha"], s["Aa"], s["Lambda"], m["j"], s["a"])


def proper_weights(g, s):
    """cmc_gauge._proper_volume_weights: r^2 sqrt(a) b / X^3, times dr."""
    return _ew(lambda r_, a_, b_, x: (r_ ** 2 * E(mp.sqrt(a_.v), a_.a ** 0.5, a_.n + 1) * b_ / x ** 3) * g.dr,
               g.r, s["a"], s["b"], s["X"])


def kdot_projection(g, s, m):
    """cmc_gauge._project_kdot_from_weights on the slice passed to the solve (its alpha is the solve's input).
    Returns (value, bound): bound = gamma_N Sum |w_i K_i| / Sum w with N = depth + n (summation)."""
    l2, l3 = rhs_primary(g, s, m)
    w = proper_weights(g, s)
    terms = [wi * (x + y) for wi, x, y in zip(w, l2["K"], l3["K"])]
    num = sum(t.v for t in terms)
    den = sum(x.v for x in w)
    depth = max(t.n for t in terms) + g.n
    bnd = gamma(depth) * sum(t.a for t in terms) / den
    return num / den, bnd


def k_prime(g, s):
    """-(2/3) D_r K, the K term of the momentum constraint."""
    return _ew(lambda kp: -(E(2) / 3.0) * kp, g.d(s["K"], 1, 1))
