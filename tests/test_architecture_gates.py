"""Deterministic architecture gates; these tests do not evolve the model."""

import ast
from dataclasses import fields
import inspect
import sys

import pytest

import engine.production_kernel as production_kernel
import engine.run_production as run_production
from engine.cosmos import (
    CosmosParams,
    construct_present_day,
    cosmos_rhs,
    dV_dphi,
    potential,
)


def test_homogeneous_scalar_dm_exchange_cancels_in_continuity_identity():
    """The archive-locked beta source transfers energy internally, not to Q."""
    p = CosmosParams()
    s = construct_present_day(p)
    y = (
        s.a, s.H, s.phi, s.pi_phi,
        s.rho_dm, s.rho_b, s.rho_r, s.rho_shear,
    )
    dy = cosmos_rhs(y, p)
    _, H, phi, pi, rho_dm, *_ = y

    rho_phi = 0.5 * pi * pi + potential(phi, p)
    p_phi = 0.5 * pi * pi - potential(phi, p)
    dV = dV_dphi(phi, p)

    # Differentiate rho_phi = pi^2/2 + V(phi) using the actual RHS.
    scalar_balance = pi * dy[3] + dV * dy[2]
    scalar_balance += 3.0 * H * (rho_phi + p_phi)
    dm_balance = dy[4] + 3.0 * H * rho_dm
    expected_scalar_source = p.beta * rho_dm * pi

    # Cancellation subtracts terms much larger than the exchange itself.
    # The tolerance is a round-off bound scaled by those actual terms.
    roundoff_scale = max(
        abs(pi * dy[3]),
        abs(dV * dy[2]),
        abs(3.0 * H * (rho_phi + p_phi)),
        abs(dy[4]),
        abs(3.0 * H * rho_dm),
        abs(expected_scalar_source),
    )
    tolerance = 128.0 * sys.float_info.epsilon * roundoff_scale

    assert abs(scalar_balance - expected_scalar_source) <= tolerance
    assert abs(dm_balance + expected_scalar_source) <= tolerance
    assert abs(scalar_balance + dm_balance) <= 2.0 * tolerance


def test_unpromoted_homogeneous_exchange_q_fails_closed():
    p = CosmosParams()
    s = construct_present_day(p)
    y = (
        s.a, s.H, s.phi, s.pi_phi,
        s.rho_dm, s.rho_b, s.rho_r, s.rho_shear,
    )
    with pytest.raises(ValueError, match="exchange_Q is not a promoted physical source"):
        cosmos_rhs(y, p, exchange_Q=1.0e-12)


def _called_names(module):
    tree = ast.parse(inspect.getsource(module))
    names = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if isinstance(node.func, ast.Name):
            names.add(node.func.id)
        elif isinstance(node.func, ast.Attribute):
            names.add(node.func.attr)
    return names


def test_current_spatial_run_owns_one_state_and_does_not_silently_add_cosmos_lane():
    """Scope guard: changing this requires an approved ownership contract."""
    owned_fields = {field.name for field in fields(production_kernel.ProductionState)}
    assert {"grid", "geometry", "scalars", "matter", "t", "tau"} <= owned_fields
    assert "cosmos" not in owned_fields

    calls = _called_names(production_kernel) | _called_names(run_production)
    assert "geometry_driven_step" not in calls
    assert "cosmos_rhs_driven" not in calls
    assert "cosmos_rhs" not in calls
    assert "CosmosState" not in calls
