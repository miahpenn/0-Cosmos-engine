"""Generate the extended-precision affine-root reference for the synthetic Issue #4 fixture.

Research/reference data only. Run from the repository root:
    python tests/data/issue4_affine_root/generate_reference.py <dps> <out.json>
Consumes the shared research helper tests/data/issue4_research/mp_fixture_helper.py.
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


def main(dps, out):
    P.load_historical_fixture(REPO)
    _use_correctly_rounded_fixture_exp()
    P.set_dps(dps)
    st0 = P.base_state()[0]
    root, cond = P.affine_root(st0, 1e-5)
    import mpmath as mp
    data = dict(
        description="Extended-precision affine root for the synthetic _stage_state fixture (reference only)",
        dps=int(dps),
        root=[mp.nstr(x, 30) for x in root],
        cond2_note="see test file; infinity-norm cond recorded below",
        cond_inf=mp.nstr(cond, 10),
        provenance=dict(
            python=platform.python_version(), numpy=np.__version__, mpmath=mp.__version__,
            repo_commit=subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                                       cwd=str(REPO)).stdout.strip(),
            helper_sha256=P.HELPER_SHA256,
            step_for_exact_centred_jacobian=1e-5,
            map_chain="shared extended-precision port of _fixed_point_map_output -> physical K RHS -> "
                      "solve_cmc_discrete"))
    json.dump(data, open(out, "w"), indent=1)
    print("wrote", out, "mpmath", mp.__version__)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 30,
         sys.argv[2] if len(sys.argv) > 2 else str(HERE / "issue4_affine_root_reference.json"))
