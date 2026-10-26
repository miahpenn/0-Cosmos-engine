from pathlib import Path
from types import SimpleNamespace
import hashlib
import json

import numpy as np

from tools.audit_pirk_same_stage_matter_propagation import (
    KAPPA,
    total_current_tangent,
    geometric_constraint_tangent_terms,
    primary_l2_decomposition,
    solve_cmc_lapse_parity_center_pde_row,
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
        "stage_state_digest": "pinned-stage-digest",
        "stage_state_immutable": True,
    }
    arrays = {
        "r": np.array([0.25, 0.75, 1.25]),
        "accepted_minus_Brown_source_aware_residual": np.array([1.0e-4, -2.0e-5, 3.0e-7]),
        "connection_D1_Cr": np.array([-1.0e-3, 2.0e-4, 4.0e-5]),
    }

    sums = write_signed_artifacts(tmp_path, report, arrays)

    persisted = json.loads((tmp_path / "report.json").read_text(encoding="utf-8"))
    assert persisted["stage_state_digest"] == "pinned-stage-digest"
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

def test_parity_center_pde_row_reproduces_regular_even_quadratic_lapse():
    n = 9
    h = 0.2
    grid = SimpleNamespace(n=n, dr=h)
    r = (np.arange(n, dtype=float) + 0.5) * h
    curvature = -0.005
    alpha_exact = 1.0 + curvature * (r**2 - r[-1]**2)  # alpha[-1] = 1
    c = 2.0 / r
    m = np.zeros(n)
    rhs = np.zeros(n)
    rhs[0] = (
        (alpha_exact[1] - alpha_exact[0]) / h**2
        + c[0] * (alpha_exact[1] - alpha_exact[0]) / (2.0*h)
    )
    for i in range(1, n-1):
        rhs[i] = (
            (alpha_exact[i-1] - 2.0*alpha_exact[i] + alpha_exact[i+1]) / h**2
            + c[i] * (alpha_exact[i+1] - alpha_exact[i-1]) / (2.0*h)
        )
    solved, centre_residual = solve_cmc_lapse_parity_center_pde_row(grid, c, m, rhs)
    np.testing.assert_allclose(solved, alpha_exact, rtol=0.0, atol=2.0e-14)
    assert abs(centre_residual) < 1.0e-13
    assert not np.isclose(alpha_exact[0], alpha_exact[1], rtol=0.0, atol=1.0e-14)


def test_primary_l2_decomposition_reconstructs_pinned_operator_form():
    n = 7
    r = np.arange(n, dtype=float) + 0.5
    alpha = 0.9 + 0.01*r
    drdr_alpha = 0.002*np.cos(r)
    lap_alpha = -0.003*np.sin(r)
    rr = 0.004*np.cos(0.3*r)
    scalar_r = -0.001*np.sin(0.4*r)
    aa_lapse = -(drdr_alpha - lap_alpha/3.0)
    aa_ricci = alpha*(rr-scalar_r/3.0)
    aa_operator = aa_lapse + aa_ricci
    k_operator = -lap_alpha
    split = primary_l2_decomposition(
        alpha,
        {"DrDr_alpha": drdr_alpha, "lap_alpha": lap_alpha, "Rr_r": rr, "R": scalar_r},
        aa_operator,
        k_operator,
    )
    np.testing.assert_allclose(split["Aa_split_residual"], 0.0, rtol=0.0, atol=1.0e-18)
    np.testing.assert_allclose(split["K_split_residual"], 0.0, rtol=0.0, atol=1.0e-18)
    np.testing.assert_array_equal(split["Aa_reconstructed"], aa_operator)
    np.testing.assert_array_equal(split["K_reconstructed_minus_lap_alpha"], k_operator)


def test_geometric_constraint_tangent_terms_reconstruct_the_exact_chain_rule():
    grid = LinearGrid(n=9, dr=0.2)
    r = grid.centers
    state = SimpleNamespace(
        grid=grid,
        geometry=SimpleNamespace(
            Aa=1.2e-4*np.exp(-0.4*r),
            b=1.0+0.015*r**2,
            X=1.0+0.02/(1.0+r**2),
        ),
    )
    rates = {
        "Aa": 1.0e-5*np.sin(r),
        "K": -2.0e-5*np.cos(0.7*r),
        "b": 3.0e-5*np.sin(0.5*r),
        "X": -1.0e-5*np.cos(0.3*r),
    }
    terms, total = geometric_constraint_tangent_terms(state, rates)
    rebuilt = np.zeros_like(r)
    for values in terms.values():
        rebuilt += values
    np.testing.assert_array_equal(total, rebuilt)
    assert set(terms) == {
        "D1_Aa_dot",
        "minus_two_thirds_D1_K_dot",
        "six_Aa_dot_chi",
        "minus_three_Aa_D1_X_dot_over_X",
        "plus_three_Aa_D1_X_X_dot_over_X_squared",
        "three_Aa_dot_over_r",
        "one_point_five_Aa_dot_D1_b_over_b",
        "one_point_five_Aa_D1_b_dot_over_b",
        "minus_one_point_five_Aa_D1_b_b_dot_over_b_squared",
    }

