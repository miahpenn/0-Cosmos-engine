"""Research-only shared helper: extended-precision (mpmath) affine-root map chain for the synthetic
Issue #4 fixture.

Shared by two independent reference generators:
  * tests/data/issue4_affine_root/generate_reference.py       (affine root)
  * tests/data/issue4_implicit_derivative/generate_derivative_reference.py (exact-root derivative)

Not a production module. Not imported by any CI test: CI tests read committed reference data.
The port is validated against the float64 tool map by tests/data/issue4_research/validate_port.py.
Module state: call load_historical_fixture() before using the functions that need the fixture.
"""
import hashlib
import importlib.util
import pathlib

import mpmath as mp

HELPER_SHA256 = hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
T = None      # historical fixture module, set by load_historical_fixture()
DR = 0.5
mp.mp.dps = 30


def load_historical_fixture(repo_root):
    """Load the historical test module that defines the synthetic _stage_state fixture."""
    global T
    path = pathlib.Path(repo_root) / "tests" / "test_issue4_cmc_lapse_tangent.py"
    spec = importlib.util.spec_from_file_location("issue4_historical", path)
    T = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(T)
    return T


def set_dps(dps):
    mp.mp.dps = int(dps)


FIELDS = ("a", "b", "X", "Aa", "K", "rho", "pr", "pt")
PI = mp.pi

def cell_deriv(u, parity, order):
    """Fourth-order cell-centred derivative with parity ghosts (port of vendor transcription)."""
    n = len(u); dr = mp.mpf(DR)
    ext = [None] * (n + 4)
    for i in range(n):
        ext[i + 2] = u[i]
    ext[1] = parity * u[0]; ext[0] = parity * u[1]
    res = [None] * n
    for i in range(n - 2):
        s = ext[i:i + 5]
        if order == 1:
            res[i] = (s[0] - 8 * s[1] + 8 * s[3] - s[4]) / (12 * dr)
        else:
            res[i] = (-s[0] + 16 * s[1] - 30 * s[2] + 16 * s[3] - s[4]) / (12 * dr * dr)
    idx = list(range(n - 5, n))
    for i in (n - 2, n - 1):
        M = mp.matrix(5, 5)
        for col, j in enumerate(idx):
            for row in range(5):
                M[row, col] = mp.mpf(j - i) ** row
        tgt = mp.matrix(5, 1); tgt[order] = mp.factorial(order)
        w = mp.lu_solve(M, tgt)
        res[i] = sum(w[k] * u[idx[k]] for k in range(5)) / dr ** order
    return res

def map_output(state, alpha_in):
    """Extended-precision port of _fixed_point_map_output (returns output alpha list)."""
    a, b, X, Aa, K, rho, pr, pt = (state[k] for k in FIELDS)
    r = state["r"]; n = len(r)
    ap = cell_deriv(a, 1, 1); bp = cell_deriv(b, 1, 1); Xp = cell_deriv(X, 1, 1)
    dal = cell_deriv(alpha_in, 1, 1); d2al = cell_deriv(alpha_in, 1, 2)
    inv = [X[i] ** 2 / a[i] for i in range(n)]
    Q = [1.5 * Aa[i] ** 2 + K[i] ** 2 / 3 + 4 * PI * (rho[i] + pr[i] + 2 * pt[i]) for i in range(n)]
    hh = [0.5 * ap[i] / a[i] - bp[i] / b[i] + Xp[i] / X[i] - 2 / r[i] for i in range(n)]
    lap = [inv[i] * (d2al[i] - dal[i] * hh[i]) for i in range(n)]
    kr = [-lap[i] + alpha_in[i] * Q[i] for i in range(n)]
    # solve_cmc_discrete
    c = [-hh[i] for i in range(n)]
    m = [Q[i] / inv[i] for i in range(n)]
    w = [r[i] * r[i] * mp.sqrt(a[i]) * b[i] / X[i] ** 3 * mp.mpf(DR) for i in range(n)]
    wsum = sum(w)
    kdot = sum(w[i] * kr[i] for i in range(n)) / wsum
    h = mp.mpf(DR)
    A = mp.matrix(n, n); rhs = mp.matrix(n, 1)
    A[0, 0] = 1; A[0, 1] = -1
    for i in range(1, n - 1):
        ci = c[i] / (2 * h)
        A[i, i - 1] = 1 / h ** 2 - ci
        A[i, i] = -2 / h ** 2 - m[i]
        A[i, i + 1] = 1 / h ** 2 + ci
        rhs[i] = -kdot / inv[i]
    A[n - 1, n - 1] = 1; rhs[n - 1] = 1
    alpha = mp.lu_solve(A, rhs)
    out = [alpha[i] for i in range(n)]
    if min(out) <= 0:
        raise FloatingPointError("mp CMC solve did not produce a positive lapse")
    return out

def to_mp(x):
    return [mp.mpf(float(v)) for v in x]

def base_state(delta=None, scale=0.0):
    base, rates = T._stage_state()
    st = {"r": to_mp(base["r"])}
    for k in FIELDS:
        vals = to_mp(base[k])
        if delta is not None:
            d = to_mp(rates[k + "_dot"]) if delta == "rates" else None
            vals = [vals[i] + mp.mpf(scale) * d[i] for i in range(len(vals))]
        st[k] = vals
    return st, base, rates

def shifted_state(st, delta_vals):
    out = dict(st)
    for k in FIELDS:
        out[k] = [st[k][i] + delta_vals[k][i] for i in range(len(st[k]))]
    return out

def affine_root(shifted, step):
    n = len(shifted["r"])
    ref = [mp.mpf(1)] * n
    fref = map_output(shifted, ref)
    Jm = mp.matrix(n, n)
    h = mp.mpf(step)
    for j in range(n):
        plus = list(ref); minus = list(ref)
        plus[j] = ref[j] + h; minus[j] = ref[j] - h
        pp = map_output(shifted, plus); pm = map_output(shifted, minus)
        for i in range(n):
            Jm[i, j] = (pp[i] - pm[i]) / (2 * h)
    system = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            system[i, j] = (1 if i == j else 0) - Jm[i, j]
    rhs = mp.matrix(n, 1)
    for i in range(n):
        rhs[i] = fref[i] - sum(Jm[i, j] * ref[j] for j in range(n))
    root = mp.lu_solve(system, rhs)
    # infinity-norm condition number of (I - J), computed in extended precision
    inv = mp.inverse(system)
    def inf_norm(M):
        return max(sum(abs(M[i, j]) for j in range(n)) for i in range(n))
    cond = inf_norm(system) * inf_norm(inv)
    return [root[i] for i in range(n)], cond

def compare(root_mp, solver64):
    d = [root_mp[i] - mp.mpf(float(solver64[i])) for i in range(len(solver64))]
    return dict(cell0_signed=float(d[0]), cell0_abs=float(abs(d[0])),
                max_abs=float(max(abs(x) for x in d)))

