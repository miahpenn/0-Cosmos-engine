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


def main(dps, out):
    P.load_historical_fixture(REPO)
    P.set_dps(dps)
    st0 = P.base_state()[0]
    root, cond = P.affine_root(st0, 1e-5)
    import mpmath as mp
    data = dict(
        description="Extended-precision affine root for the synthetic _stage_state fixture (reference only)",
        dps=int(dps),
        root=[mp.nstr(x, 30) for x in root],
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
