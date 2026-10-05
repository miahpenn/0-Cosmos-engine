import numpy as np

from engine.v55_matter import exchange_pair


def test_dm_scalar_exchange_cancels():
    x = np.array([1.0, 2.0])
    z = exchange_pair(-0.04, np.array([3.0, 4.0]), x, z if False else x * 2.0)
    assert np.allclose(z["net_t"], 0.0)
    assert np.allclose(z["net_r"], 0.0)
    assert np.allclose(z["phi_t"], -z["dm_t"])
    assert np.allclose(z["phi_r"], -z["dm_r"])
