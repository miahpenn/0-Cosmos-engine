"""Generate the exact-root directional-derivative reference for the synthetic Issue #4 fixture.

Research/reference data only. Run from the repository root:
    python tests/data/issue4_implicit_derivative/generate_derivative_reference.py <dps> <out.json>
For epsilon in {2e-5, 1e-5}: central difference of the extended-precision affine root along the
fixture's rates for the eight fields (a, b, X, Aa, K, rho, pr, pt). Consumes the shared research helper.
"""
import json
import pathlib
import platform
import subprocess
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tests" / "data" / "issue4_research"))
import mp_fixture_helper as P  # noqa: E402


# Correctly rounded exp for the fixture. _stage_state() (and _state(), which it calls) use one
# transcendental function, np.exp, through the fixture module's global `np`. NumPy's float64 exp can
# differ in the last bit between SIMD code paths (AVX-512 vs not). While the fixture's own _stage_state()
# runs, the fixture module's `np` is replaced by a proxy whose exp evaluates each element in mpmath at
# 40 digits and rounds once to float64; every other attribute is NumPy's. The original `np` is restored
# after each call. Only this generator process is affected; the test fixture file is unchanged.
class _CorrectlyRoundedExpNumpy:
    def __getattr__(self, name):
        return getattr(np, name)

    @staticmethod
    def exp(x):
        import mpmath as mp
        arr = np.asarray(x, dtype=np.float64)
        with mp.workdps(40):
            out = [float(mp.exp(mp.mpf(float(v)))) for v in arr.ravel()]
        return np.array(out, dtype=np.float64).reshape(arr.shape)


def _use_correctly_rounded_fixture_exp():
    """Wrap the loaded fixture's _stage_state so np.exp inside it is correctly rounded."""
    original = P.T._stage_state

    def stage_state():
        saved = P.T.np
        P.T.np = _CorrectlyRoundedExpNumpy()
        try:
            return original()
        finally:
            P.T.np = saved

    P.T._stage_state = stage_state

EPSILONS = (2.0e-5, 1.0e-5)


def derivative(eps):
    import mpmath as mp
    st0 = P.base_state()[0]
    _, rates = P.T._stage_state()
    n = len(st0["r"]); e = mp.mpf(eps)
    pd = {k: [e * mp.mpf(float(rates[k + "_dot"][i])) for i in range(n)] for k in P.FIELDS}
    md = {k: [-x for x in pd[k]] for k in P.FIELDS}
    p, _ = P.affine_root(P.shifted_state(st0, pd), 1e-5)
    m, _ = P.affine_root(P.shifted_state(st0, md), 1e-5)
    return [(p[i] - m[i]) / (2 * e) for i in range(n)]


def main(dps, out):
    import mpmath as mp
    P.load_historical_fixture(REPO)
    _use_correctly_rounded_fixture_exp()
    P.set_dps(dps)
    data = dict(
        description="Exact-root central-difference directional derivative of the affine lapse root "
                    "for the synthetic _stage_state fixture (reference only)",
        dps=int(dps),
        epsilons=[f"{e:g}" for e in EPSILONS],
        direction_fields=list(P.FIELDS),
        derivative={f"{e:g}": [mp.nstr(x, 30) for x in derivative(e)] for e in EPSILONS},
        provenance=dict(
            python=platform.python_version(), numpy=np.__version__, mpmath=mp.__version__,
            repo_commit=subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                                       cwd=str(REPO)).stdout.strip(),
            helper_sha256=P.HELPER_SHA256,
            map_chain="shared extended-precision port; affine_root step 1e-5"))
    json.dump(data, open(out, "w"), indent=1)
    print("wrote", out, "mpmath", mp.__version__)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 30,
         sys.argv[2] if len(sys.argv) > 2 else str(HERE / "issue4_implicit_derivative_reference.json"))
