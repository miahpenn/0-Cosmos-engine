from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from analyze_centre_ricci_derivative_channels import (
    cell_derivative_fourth,
    cell_second_derivative_fourth,
    pinned_vendor_tf_ricci,
    provenance_metadata,
)


def test_pinned_fourth_derivative_is_exact_on_odd_linear_profile():
    dr = 0.2
    r = (np.arange(20, dtype=float) + 0.5) * dr
    derivative = cell_derivative_fourth(r, dr, parity=-1)
    np.testing.assert_allclose(derivative, 1.0, rtol=0.0, atol=2e-13)


def test_pinned_fourth_second_derivative_is_exact_on_even_quadratic_profile():
    dr = 0.125
    r = (np.arange(20, dtype=float) + 0.5) * dr
    derivative = cell_second_derivative_fourth(r**2, dr, parity=+1)
    np.testing.assert_allclose(derivative, 2.0, rtol=0.0, atol=3e-12)


def test_pinned_vendor_tf_ricci_derivative_and_algebraic_response_identity():
    n, rmax = 40, 80.0
    dr = rmax / n
    r = (np.arange(n, dtype=float) + 0.5) * dr
    q = 1.0 + 3e-5*r**2 + 2e-8*r**4
    a = q**(-4.0/3.0)
    b = q**(2.0/3.0)
    X = q**(1.0/3.0)
    d1a = cell_derivative_fourth(a, dr, +1)
    d1b = cell_derivative_fourth(b, dr, +1)
    metric_lambda = d1a/(2*a*a) - d1b/(a*b) + 2/r*(1/b-1/a)
    evolved_lambda = metric_lambda + 1e-5*r*np.exp(-(r/6.0)**2)

    tf_e, d4_e, alg_e = pinned_vendor_tf_ricci(a,b,X,evolved_lambda,r,dr)
    tf_m, d4_m, alg_m = pinned_vendor_tf_ricci(a,b,X,metric_lambda,r,dr)
    response = tf_m - tf_e
    predicted = (2.0/3.0)*X*X*(d4_m-d4_e) + (alg_m-alg_e)
    np.testing.assert_allclose(response, predicted, rtol=0.0, atol=1e-12)



def test_vendor_operator_provenance_labels_git_blob_sha1_unambiguously():
    provenance = provenance_metadata()
    digest = provenance["vendor_operator_git_blob_sha1"]
    assert len(digest) == 40
    assert all(character in "0123456789abcdef" for character in digest)
    assert "vendor_operator_blob_sha256" not in provenance
