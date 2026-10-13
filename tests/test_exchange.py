import numpy as np

from engine.v55_matter import exchange_pair


def test_dm_scalar_exchange_cancels():
    dphi_t = np.array([1.0, 2.0])
    dphi_r = np.array([2.0, -1.0])
    z = exchange_pair(
        -0.04,
        np.array([3.0, 4.0]),
        dphi_t,
        dphi_r,
    )
    assert np.allclose(z["net_t"], 0.0)
    assert np.allclose(z["net_r"], 0.0)
    assert np.allclose(z["phi_t"], -z["dm_t"])
    assert np.allclose(z["phi_r"], -z["dm_r"])
