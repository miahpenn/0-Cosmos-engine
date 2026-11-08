"""Research-only guard against circular validation of the shared extended-precision port.

This file needs mpmath and therefore runs ONLY in the research diagnostic job (not in the blocking
contract job). It checks the shared port (tests/data/issue4_research/mp_fixture_helper.py) against the
production-tool float64 map, which is the code under test. A port transcription error would show here.

Predeclared threshold: max abs difference <= 1e-12 at working precision 30 digits (observed ~7e-15).
"""
import importlib.util
import pathlib
import sys

import numpy as np


import mpmath as mp  # hard import: a missing mpmath must fail, never silently skip

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "tests" / "data" / "issue4_research"))
import mp_fixture_helper as P  # noqa: E402

PORT_THRESHOLD = 1.0e-12


def test_shared_port_matches_float64_tool_map_at_base_state():
    T = P.load_historical_fixture(REPO)
    P.set_dps(30)
    base, _ = T._stage_state()
    st0 = P.base_state()[0]
    n = base["r"].size
    out64 = T._fixed_point_map_output(base, np.ones(n)).alpha
    outmp = P.map_output(st0, [mp.mpf(1)] * n)
    diff = max(abs(outmp[i] - mp.mpf(float(out64[i]))) for i in range(n))
    print(f"\nshared port vs float64 tool map (base): max abs = {mp.nstr(diff, 3)}")
    assert diff <= PORT_THRESHOLD


def test_shared_port_matches_float64_tool_map_at_nontrivial_lapse():
    T = P.load_historical_fixture(REPO)
    P.set_dps(30)
    base, _ = T._stage_state()
    st0 = P.base_state()[0]
    n = base["r"].size
    alpha = 1.0 - 0.03 * np.exp(-(base["r"] / 2.5) ** 2)
    out64 = T._fixed_point_map_output(base, alpha).alpha
    outmp = P.map_output(st0, [mp.mpf(float(x)) for x in alpha])
    diff = max(abs(outmp[i] - mp.mpf(float(out64[i]))) for i in range(n))
    print(f"\nshared port vs float64 tool map (nontrivial lapse): max abs = {mp.nstr(diff, 3)}")
    assert diff <= PORT_THRESHOLD
