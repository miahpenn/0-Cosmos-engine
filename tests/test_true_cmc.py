import numpy as np
import pytest

try:
    import scipy
except ImportError:
    pytest.skip("SciPy unavailable", allow_module_level=True)

from engine.true_cmc_reference import cmc_lapse, run
from engine import reference_pirk_unified as q


def test_true_cmc_lapse_is_finite_and_positive():
    g = q.Grid(24, 20.0)
    s, fields = q.make_initial(g, .003, 7.0)
    alpha, kd, amin, amax = cmc_lapse(g, s, fields)
    assert np.all(np.isfinite(alpha))
    assert amin > 0.0
    assert amax >= amin


@pytest.mark.slow
def test_true_cmc_smoke():
    out = run(24, 0.5, .06)
    assert out["ok"]
