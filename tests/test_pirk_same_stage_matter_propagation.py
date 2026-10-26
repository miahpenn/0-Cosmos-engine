from pathlib import Path
from types import SimpleNamespace
import hashlib
import json

import numpy as np

from tools.audit_pirk_same_stage_matter_propagation import (
    KAPPA,
    total_current_tangent,
    write_signed_artifacts,
)


class LinearGrid:
    def __init__(self, n=9, dr=0.2):
        self.centers = (np.arange(n, dtype=float) + 0.5) * dr
        self.dr = dr

    def cell_derivative_fourth(self, values, parity=1):
        """Simple fixed linear operator sufficient to test the chain rule."""
        u = np.asarray(values, dtype=float)
        out = np.empty_like(u)
        out[0] = (u[1] - parity * u[0]) / (2.0 * self.dr)
        out[1:-1] = (u[2:] - u[:-2]) / (2.0 * self.dr)
        out[-1] = (3.0 * u[-1] - 4.0 * u[-2] + u[-3]) / (2.0 * self.dr)
        return out


def _current(grid, geometry, scalars, matter):
    r = grid.centers
    root_gamma = np.sqrt(geometry.a) * geometry.b * r**2 / geometry.X**3
    d1 = lambda x: grid.cell_derivative_fourth(x, parity=1)
    total = (
        -scalars.PS * d1(scalars.S)
        -scalars.PD * d1(scalars.D)
        -scalars.Pi * d1(scalars.phi) / KAPPA
    )
    for name in ("dark_matter", "baryons", "radiation"):
        total += getattr(matter, name).momentum_r / (root_gamma * KAPPA)
    return total


def test_six_sector_current_tangent_matches_symmetric_directional_difference():
    grid = LinearGrid()
    r = grid.centers
    geometry = SimpleNamespace(
        a=1.0 + 0.03*r**2, b=1.0 + 0.02*r**2,
        X=1.0 + 0.01/(1.0+r**2),
    )
    scalars = SimpleNamespace(
        S=np.exp(-r), PS=0.02*np.exp(-0.5*r),
        D=0.3*np.exp(-0.7*r), PD=0.01*np.exp(-0.4*r),
        phi=0.1*np.exp(-0.3*r), Pi=0.005*np.exp(-0.6*r),
    )
    scalar_rates = SimpleNamespace(**{
        name: 1.0e-3*np.cos((index+1)*r)
        for index, name in enumerate(("S", "PS", "D", "PD", "phi", "Pi"))
    })
    matter = SimpleNamespace()
    matter_rates = {}
    for index, name in enumerate(("dark_matter", "baryons", "radiation")):
        momentum = (index+1)*1.0e-4*np.exp(-0.2*r)
        rate = 1.0e-5*np.sin((index+1)*r)
        setattr(matter, name, SimpleNamespace(momentum_r=momentum))
        matter_rates[name] = SimpleNamespace(momentum_r=rate)

    rates = {
        "a": 2.0e-4*np.sin(r),
        "b": -1.0e-4*np.cos(r),
        "X": 0.5e-4*np.sin(0.7*r),
    }
    state = SimpleNamespace(grid=grid, geometry=geometry, scalars=scalars, matter=matter)
    pieces, total_dot, volume_rate = total_current_tangent(
        state, rates, scalar_rates, matter_rates
    )

    eps = 2.0e-6
    plus_g = SimpleNamespace(**{
        name: getattr(geometry, name) + eps*rates[name] for name in ("a", "b", "X")
    })
    minus_g = SimpleNamespace(**{
        name: getattr(geometry, name) - eps*rates[name] for name in ("a", "b", "X")
    })
    plus_s = SimpleNamespace(**{
        name: getattr(scalars, name) + eps*getattr(scalar_rates, name)
        for name in ("S", "PS", "D", "PD", "phi", "Pi")
    })
    minus_s = SimpleNamespace(**{
        name: getattr(scalars, name) - eps*getattr(scalar_rates, name)
        for name in ("S", "PS", "D", "PD", "phi", "Pi")
    })
    plus_m, minus_m = SimpleNamespace(), SimpleNamespace()
    for name in ("dark_matter", "baryons", "radiation"):
        U = getattr(matter, name).momentum_r
        Udot = matter_rates[name].momentum_r
        setattr(plus_m, name, SimpleNamespace(momentum_r=U+eps*Udot))
        setattr(minus_m, name, SimpleNamespace(momentum_r=U-eps*Udot))

    finite_difference = (
        _current(grid, plus_g, plus_s, plus_m)
        - _current(grid, minus_g, minus_s, minus_m)
    ) / (2.0*eps)
    np.testing.assert_allclose(total_dot, finite_difference, rtol=2.0e-7, atol=2.0e-9)
    assert set(pieces) == {"S", "D", "COSMOS_phi", "dark_matter", "baryons", "radiation"}
    assert np.all(np.isfinite(volume_rate))


def test_signed_artifact_writer_preserves_profiles_and_manifest(tmp_path: Path):
    report = {
        "schema": "issue4-same-stage-matter-brown-propagation/v1",
        "status": "completed_rhs_only_diagnostic",
        "scientific_admission": False,
        "input_state_digest": "pinned-digest",
        "input_state_immutable": True,
    }
    arrays = {
        "r": np.array([0.25, 0.75, 1.25]),
        "accepted_minus_Brown_source_aware_residual": np.array([1.0e-4, -2.0e-5, 3.0e-7]),
        "connection_D1_Cr": np.array([-1.0e-3, 2.0e-4, 4.0e-5]),
    }

    sums = write_signed_artifacts(tmp_path, report, arrays)

    persisted = json.loads((tmp_path / "report.json").read_text(encoding="utf-8"))
    assert persisted["input_state_digest"] == "pinned-digest"
    assert persisted["scientific_admission"] is False
    with np.load(tmp_path / "signed_arrays.npz", allow_pickle=False) as saved:
        assert set(saved.files) == set(arrays)
        for name, expected in arrays.items():
            np.testing.assert_array_equal(saved[name], expected)
    manifest = (tmp_path / "SHA256SUMS.txt").read_text(encoding="utf-8")
    for name in ("report.json", "signed_arrays.npz"):
        expected = hashlib.sha256((tmp_path / name).read_bytes()).hexdigest()
        assert sums[name] == expected
        assert f"{expected}  {name}" in manifest

