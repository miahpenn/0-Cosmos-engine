"""Evolved-state curvature consistency and local constraint-budget audit.

Compares the pinned vendor BSSN curvature with the repository's general-BSSN
ricci_terms() on identical states. Both use evolved Lambda in the Ricci
principal part. No polar-areal identity is used.

At t=0,1,2,3, also measures one-step Hamiltonian changes from the curvature,
extrinsic-curvature terms, S, D, COSMOS phi, and fluid energy. The regularity
projection's state-map offset is reported separately. Measurement only.
"""
import copy
import math
import pathlib
import sys
from types import SimpleNamespace

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from engine import reference_pirk_unified as ref
from engine import v55_pirk_adapter as adapter
from engine.discrete_consistent_initial import discrete_consistent_state
from engine.production_kernel import LAMBDA_M, V55ProductionKernel
from engine.scalar_system import ScalarFields, scalar_projection
from engine.stress_energy import assemble_total_stress_energy

_, vacuum, _ = adapter.vendor_modules()

N = 40
R_MAX = 40.0
AMPLITUDE = 0.01
D_AMPLITUDE = 1.0e-10
CFL = 0.03
TARGET_TIMES = (1.0, 2.0, 3.0)
PROBE_DT_FACTORS = (0.03, 0.015, 0.0075)
INITIAL_AGREEMENT_TOL = 1.0e-12


def reference_state(geometry):
    return SimpleNamespace(
        a=np.asarray(geometry.a, dtype=float),
        b=np.asarray(geometry.b, dtype=float),
        X=np.asarray(geometry.X, dtype=float),
        alpha=np.asarray(geometry.alpha, dtype=float),
        Lambda=np.asarray(geometry.Lambda, dtype=float),
    )


def compare_derivative_operators(state, reference_grid):
    """Compare every first/second derivative used by the curvature formula."""
    grid = state.grid
    geometry = state.geometry
    results = {}
    first_order = {
        "a": (geometry.a, +1),
        "b": (geometry.b, +1),
        "X": (geometry.X, +1),
        "Lambda": (geometry.Lambda, -1),
        "alpha": (geometry.alpha, +1),
    }
    second_order = {
        "a": (geometry.a, +1),
        "b": (geometry.b, +1),
        "X": (geometry.X, +1),
        "alpha": (geometry.alpha, +1),
    }

    for name, (values, parity) in first_order.items():
        vendor = np.asarray(grid.cell_derivative_fourth(values, parity=parity))
        reference = np.asarray(ref.D(reference_grid, values, parity))
        results[f"D1:{name}"] = vendor - reference

    for name, (values, parity) in second_order.items():
        vendor = np.asarray(grid.cell_second_derivative_fourth(values, parity=parity))
        reference = np.asarray(ref.D2(reference_grid, values, parity))
        results[f"D2:{name}"] = vendor - reference

    return results


def scalar_density_components(state):
    """Return isolated S, D, and phi energy densities for delta accounting.

    In the S-only and D-only isolated evaluations, phi=0 contributes the same
    constant potential offset at every time. That offset cancels in differences;
    these component values must not be summed as absolute total rho.
    """
    zero = np.zeros_like(state.scalars.S)
    components = {}
    names = ("S", "PS", "D", "PD", "phi", "Pi")
    for field_name in ("S", "D", "phi"):
        isolated = ScalarFields(*(
            getattr(state.scalars, name)
            if name == field_name
            or (field_name == "S" and name == "PS")
            or (field_name == "D" and name == "PD")
            or (field_name == "phi" and name == "Pi")
            else zero
            for name in names
        ))
        components[field_name] = np.asarray(
            scalar_projection(state.grid, state.geometry, isolated)[0],
            dtype=float,
        )
    return components


def hamiltonian_residual(state):
    raw = vacuum.constraints(state.grid, state.geometry)
    total = assemble_total_stress_energy(
        state.grid, state.geometry, state.scalars, state.matter
    )
    return np.asarray(raw["hamiltonian"]) - 16.0 * math.pi * np.asarray(total.rho)


def metric_connection_expression(grid, geometry):
    """Metric-derived conformal connection used by the vendor constraint."""
    radius = np.asarray(grid.centers)
    ap = np.asarray(grid.cell_derivative_fourth(geometry.a, parity=+1))
    bp = np.asarray(grid.cell_derivative_fourth(geometry.b, parity=+1))
    return (
        ap / (2.0 * np.asarray(geometry.a) ** 2)
        - bp / (np.asarray(geometry.a) * np.asarray(geometry.b))
        + 2.0 / radius * (
            1.0 / np.asarray(geometry.b) - 1.0 / np.asarray(geometry.a)
        )
    )


def momentum_components(state):
    """Signed terms in M_total = M_geometry - 8*pi*j, including source sectors."""
    grid, geometry = state.grid, state.geometry
    radius = np.asarray(grid.centers)
    aa = np.asarray(geometry.Aa)
    K = np.asarray(geometry.K)
    a = np.asarray(geometry.a)
    b = np.asarray(geometry.b)
    X = np.asarray(geometry.X)

    dAa = np.asarray(grid.cell_derivative_fourth(aa, parity=+1))
    dK = np.asarray(grid.cell_derivative_fourth(K, parity=+1))
    db = np.asarray(grid.cell_derivative_fourth(b, parity=+1))
    dX = np.asarray(grid.cell_derivative_fourth(X, parity=+1))
    chip = -0.5 * dX / X

    # Vendor convention: Ab=-Aa/2, so Aa-Ab=3 Aa/2.
    terms = {
        "M_DAa": dAa,
        "M_DK": -(2.0 / 3.0) * dK,
        "M_chi": 6.0 * aa * chip,
        "M_radial_Aa": 1.5 * aa * (2.0 / radius + db / b),
    }

    zero = np.zeros_like(state.scalars.S)
    names = ("S", "PS", "D", "PD", "phi", "Pi")
    scalar_js = {}
    for field_name in ("S", "D", "phi"):
        isolated = ScalarFields(*(
            getattr(state.scalars, name)
            if name == field_name
            or (field_name == "S" and name == "PS")
            or (field_name == "D" and name == "PD")
            or (field_name == "phi" and name == "Pi")
            else zero
            for name in names
        ))
        scalar_js[field_name] = np.asarray(
            scalar_projection(grid, geometry, isolated)[3], dtype=float
        )

    total = assemble_total_stress_energy(grid, geometry, state.scalars, state.matter)
    fluid_j = np.asarray(total.fluid_j, dtype=float)
    terms.update({
        "M_S_source": -8.0 * math.pi * scalar_js["S"],
        "M_D_source": -8.0 * math.pi * scalar_js["D"],
        "M_phi_source": -8.0 * math.pi * scalar_js["phi"],
        "M_fluid_source": -8.0 * math.pi * fluid_j,
    })

    raw = vacuum.constraints(grid, geometry)
    direct = np.asarray(raw["momentum"], dtype=float) - 8.0 * math.pi * np.asarray(total.j)
    reconstructed = np.sum(np.stack(list(terms.values())), axis=0)
    terms["M_total"] = direct
    terms["closure"] = reconstructed - direct
    return terms


def report_momentum_checkpoint(state, label):
    terms = momentum_components(state)
    print(f"[MOMENTUM_BUDGET {label}] t={state.t:.12g}")
    for name, values in terms.items():
        if name == "closure":
            print(
                f"  {name}: maxabs_cells0-4={np.max(np.abs(values[:5])):.3e}; "
                f"maxabs_all={np.max(np.abs(values)):.3e}"
            )
            continue
        print(
            f"  {name}: cell0={values[0]:+.6e}; "
            f"maxabs_cells0-4={np.max(np.abs(values[:5])):.3e}; "
            f"maxabs_all={np.max(np.abs(values)):.3e}"
        )


def constraint_components(state):
    """Term accounting for H; signs match the vendor constraint definition."""
    geometry = state.geometry
    raw = vacuum.constraints(state.grid, geometry)
    terms = vacuum.geometry_terms(state.grid, geometry)
    total = assemble_total_stress_energy(
        state.grid, geometry, state.scalars, state.matter
    )
    scalar_parts = scalar_density_components(state)
    R = np.asarray(terms["R"], dtype=float)
    Aa = np.asarray(geometry.Aa, dtype=float)
    K = np.asarray(geometry.K, dtype=float)
    A_term = -1.5 * Aa**2
    K_term = (2.0 / 3.0) * K**2
    H_vendor = np.asarray(raw["hamiltonian"]) - 16.0 * math.pi * np.asarray(total.rho)
    H_reconstructed = (
        R + A_term + K_term
        - 16.0 * math.pi * np.asarray(total.scalar_rho)
        - 16.0 * math.pi * np.asarray(total.fluid_rho)
    )
    return {
        "R": R,
        "A_term": A_term,
        "K_term": K_term,
        "rho_scalar": np.asarray(total.scalar_rho, dtype=float),
        "rho_fluid": np.asarray(total.fluid_rho, dtype=float),
        "S_rho": scalar_parts["S"],
        "D_rho": scalar_parts["D"],
        "phi_rho": scalar_parts["phi"],
        "H_vendor": H_vendor,
        "H_reconstructed": H_reconstructed,
        "closure": H_reconstructed - H_vendor,
    }


def report_curvature(state, reference_grid, label, initial_gate=False):
    geometry = state.geometry
    radius = np.asarray(state.grid.centers)
    vendor_R = np.asarray(vacuum.geometry_terms(state.grid, geometry)["R"])
    reference_R = np.asarray(ref.ricci_terms(reference_grid, reference_state(geometry))[0])
    curvature_diff = vendor_R - reference_R
    derivative_diffs = compare_derivative_operators(state, reference_grid)
    operator_maxima = {
        name: float(np.max(np.abs(values)))
        for name, values in derivative_diffs.items()
    }
    operator_max = max(operator_maxima.values(), default=0.0)
    center_curvature_max = float(np.max(np.abs(curvature_diff[:5])))
    all_curvature_max = float(np.max(np.abs(curvature_diff)))
    idx = int(np.argmax(np.abs(curvature_diff)))

    H = hamiltonian_residual(state)
    raw_constraints = vacuum.constraints(state.grid, geometry)
    total_stress = assemble_total_stress_energy(
        state.grid, geometry, state.scalars, state.matter
    )
    # The source normalization matches reference_pirk_unified.constraint():
    # M_total = M_geometry - 8*pi*j, with j assembled from all matter sectors.
    momentum_total = (
        np.asarray(raw_constraints["momentum"], dtype=float)
        - 8.0 * math.pi * np.asarray(total_stress.j, dtype=float)
    )
    constraint_profiles = {
        "connection": np.asarray(raw_constraints["connection"], dtype=float),
        "determinant": np.asarray(raw_constraints["determinant"], dtype=float),
        "momentum_geometric": np.asarray(raw_constraints["momentum"], dtype=float),
        "momentum_total_minus_8pi_j": momentum_total,
    }
    areal_departure = np.asarray(geometry.b) / np.asarray(geometry.X) ** 2 - 1.0

    print(f"[CURVATURE {label}] t={state.t:.12g}")
    print("  derivative operator max-absolute differences:")
    for name in sorted(operator_maxima):
        d = derivative_diffs[name]
        print(
            f"    {name}: all={operator_maxima[name]:.3e}; "
            f"cells0-4={np.max(np.abs(d[:5])):.3e}"
        )
    print(f"  curvature R vendor cells0-4 = {np.array2string(vendor_R[:5], precision=9)}")
    print(f"  curvature R ref    cells0-4 = {np.array2string(reference_R[:5], precision=9)}")
    print(f"  signed R difference cells0-4 = {np.array2string(curvature_diff[:5], precision=4)}")
    print(
        f"  |delta R| max: all={all_curvature_max:.3e} at cell {idx} "
        f"(r={radius[idx]:.6g}); cells0-4={center_curvature_max:.3e}; "
        f"cells>=5={np.max(np.abs(curvature_diff[5:])):.3e}"
    )
    print(
        f"  b/X^2 - 1: cells0-4={np.array2string(areal_departure[:5], precision=4)}; "
        f"max|.| all={np.max(np.abs(areal_departure)):.3e}"
    )
    print(
        f"  Hamiltonian residual: cells0-4={np.array2string(H[:5], precision=4)}; "
        f"max|H| all={np.max(np.abs(H)):.3e}"
    )
    print("  BSSN constraint profiles (no centre masking):")
    for name, values in constraint_profiles.items():
        print(
            f"    {name}: cells0-4={np.array2string(values[:5], precision=4)}; "
            f"cell0={values[0]:+.4e}; maxabs_cells0-4={np.max(np.abs(values[:5])):.3e}; "
            f"maxabs_all={np.max(np.abs(values)):.3e}"
        )

    if initial_gate:
        if operator_max > INITIAL_AGREEMENT_TOL:
            raise AssertionError(
                f"initial derivative-operator agreement failed: {operator_max:.3e} "
                f"> {INITIAL_AGREEMENT_TOL:.1e}"
            )
        if all_curvature_max > INITIAL_AGREEMENT_TOL:
            raise AssertionError(
                f"initial curvature agreement failed: {all_curvature_max:.3e} "
                f"> {INITIAL_AGREEMENT_TOL:.1e}"
            )
        print(
            f"  INITIAL_GATE=PASS (operator max {operator_max:.3e}, "
            f"curvature max {all_curvature_max:.3e}, tolerance "
            f"{INITIAL_AGREEMENT_TOL:.1e})"
        )


def metric_connection_time_derivative(grid, geometry, adot, bdot):
    """Directional derivative of the discrete metric-connection expression."""
    radius = np.asarray(grid.centers)
    a = np.asarray(geometry.a)
    b = np.asarray(geometry.b)
    da = np.asarray(grid.cell_derivative_fourth(a, parity=+1))
    db = np.asarray(grid.cell_derivative_fourth(b, parity=+1))
    dadot = np.asarray(grid.cell_derivative_fourth(adot, parity=+1))
    dbdot = np.asarray(grid.cell_derivative_fourth(bdot, parity=+1))
    adot = np.asarray(adot)
    bdot = np.asarray(bdot)
    return (
        dadot / (2.0 * a**2)
        - da * adot / a**3
        - dbdot / (a * b)
        + db * adot / (a**2 * b)
        + db * bdot / (a * b**2)
        + 2.0 / radius * (-bdot / b**2 + adot / a**2)
    )


def connection_rhs_budget(kernel, state):
    """Compare C_Lambda RHS with the momentum-driver source and metric rate."""
    rhs_state = copy.deepcopy(state)
    rhs_state.geometry = vacuum.enforce_algebraic_regularity(
        rhs_state.grid, rhs_state.geometry.copy()
    )
    rhs_state.geometry.alpha = kernel._solve_lapse(
        rhs_state.grid, rhs_state.geometry, rhs_state.scalars, rhs_state.matter
    )[0]
    rhs_state.geometry.beta.fill(0.0)
    rhs_state.geometry.B.fill(0.0)

    grid = rhs_state.grid
    geometry = rhs_state.geometry
    stage_terms = adapter.geometry_stage_terms(
        grid, geometry, rhs_state.scalars, rhs_state.matter, lambda_m=LAMBDA_M
    )
    adot = np.asarray(stage_terms["explicit"]["a"], dtype=float)
    bdot = np.asarray(stage_terms["explicit"]["b"], dtype=float)
    metric_connection_rate = metric_connection_time_derivative(
        grid, geometry, adot, bdot
    )
    lambda_rate = (
        np.asarray(stage_terms["lambda_l2"], dtype=float)
        + np.asarray(stage_terms["lambda_l3"], dtype=float)
    )
    c_rate = lambda_rate - metric_connection_rate

    # Remove only the configured lambda_m momentum-constraint channel, without
    # changing any state: lambda_l2 carries +lambda_m*alpha*M_geom/a and the
    # matter correction in lambda_l3 carries -8*pi*lambda_m*alpha*j/a.
    _, vacuum_local, _ = adapter.vendor_modules()
    lambda_rate_no_m = (
        np.asarray(vacuum_local.lambda_l2_rhs(grid, geometry, lambda_m=0.0))
        + np.asarray(vacuum_local.lambda_l3_rhs(grid, geometry))
    )
    c_rate_no_m = lambda_rate_no_m - metric_connection_rate
    raw = vacuum_local.constraints(grid, geometry)
    total = assemble_total_stress_energy(
        grid, geometry, rhs_state.scalars, rhs_state.matter
    )
    momentum_total = np.asarray(raw["momentum"]) - 8.0 * math.pi * np.asarray(total.j)
    momentum_channel = (
        LAMBDA_M * np.asarray(geometry.alpha) / np.asarray(geometry.a)
        * momentum_total
    )
    closure = c_rate - (c_rate_no_m + momentum_channel)

    print(f"[CONNECTION_RHS_BUDGET] t={state.t:.12g}")
    for name, values in (
        ("Cdot_full", c_rate),
        ("Cdot_no_lambda_m", c_rate_no_m),
        ("lambda_m_alpha_over_a_times_Mtotal", momentum_channel),
    ):
        print(
            f"  {name}: cell0={values[0]:+.6e}; "
            f"maxabs_cells0-4={np.max(np.abs(values[:5])):.6e}; "
            f"maxabs_all={np.max(np.abs(values)):.6e}"
        )
    print(
        f"  closure maxabs={np.max(np.abs(closure)):.3e}; "
        f"lambda-dot cell0={lambda_rate[0]:+.6e}; "
        f"metric-connection-dot cell0={metric_connection_rate[0]:+.6e}; "
        f"Mtotal cell0={momentum_total[0]:+.6e}"
    )


def one_step_budget(kernel, state):
    """Probe the state with short steps from a separately projected copy."""
    raw_H = hamiltonian_residual(state)
    raw_connection = np.asarray(
        vacuum.constraints(state.grid, state.geometry)["connection"], dtype=float
    )
    projected = copy.deepcopy(state)
    projected.geometry = vacuum.enforce_algebraic_regularity(
        projected.grid, projected.geometry.copy()
    )
    base_H = hamiltonian_residual(projected)
    projection_jump = base_H - raw_H
    base_connection = np.asarray(
        vacuum.constraints(projected.grid, projected.geometry)["connection"],
        dtype=float,
    )
    connection_projection_jump = base_connection - raw_connection
    base_metric_connection = metric_connection_expression(
        projected.grid, projected.geometry
    )
    print(f"[ONE_STEP_BUDGET] base_t={state.t:.12g}")
    print(
        f"  projection-only deltaH: cell0={projection_jump[0]:+.6e}; "
        f"maxabs_cells0-4={np.max(np.abs(projection_jump[:5])):.6e}; "
        f"maxabs_all={np.max(np.abs(projection_jump)):.6e}"
    )
    print(
        f"  projection-only deltaC_Lambda: cell0={connection_projection_jump[0]:+.6e}; "
        f"maxabs_cells0-4={np.max(np.abs(connection_projection_jump[:5])):.6e}; "
        f"maxabs_all={np.max(np.abs(connection_projection_jump)):.6e}"
    )

    base = constraint_components(projected)
    connection_rhs_budget(kernel, projected)
    momentum_base = momentum_components(projected)
    base_dr = projected.grid.dr
    for factor in PROBE_DT_FACTORS:
        dt = factor * base_dr
        advanced = kernel.step(copy.deepcopy(projected), dt)
        after = constraint_components(advanced)
        momentum_after = momentum_components(advanced)
        after_connection = np.asarray(
            vacuum.constraints(advanced.grid, advanced.geometry)["connection"],
            dtype=float,
        )
        delta_connection = after_connection - base_connection
        delta_lambda = (
            np.asarray(advanced.geometry.Lambda) - np.asarray(projected.geometry.Lambda)
        )
        after_metric_connection = metric_connection_expression(
            advanced.grid, advanced.geometry
        )
        delta_metric_connection = after_metric_connection - base_metric_connection
        connection_split_closure = (
            delta_connection - (delta_lambda - delta_metric_connection)
        )
        contributions = {
            "dR": after["R"] - base["R"],
            "dA_term": after["A_term"] - base["A_term"],
            "dK_term": after["K_term"] - base["K_term"],
            "S_term": -16.0 * math.pi * (after["S_rho"] - base["S_rho"]),
            "D_term": -16.0 * math.pi * (after["D_rho"] - base["D_rho"]),
            "phi_term": -16.0 * math.pi * (after["phi_rho"] - base["phi_rho"]),
            "fluid_term": -16.0 * math.pi * (after["rho_fluid"] - base["rho_fluid"]),
        }
        direct_delta = after["H_vendor"] - base["H_vendor"]
        reconstructed_delta = np.sum(np.stack(list(contributions.values())), axis=0)
        scalar_delta = -16.0 * math.pi * (after["rho_scalar"] - base["rho_scalar"])
        scalar_split_delta = (
            contributions["S_term"] + contributions["D_term"] + contributions["phi_term"]
        )
        print(f"  dt/dr={factor:.4g} dt={dt:.6g}: cell0 dH={direct_delta[0]:+.6e}")
        for name, values in contributions.items():
            print(
                f"    {name}: cell0={values[0]:+.6e}; "
                f"maxabs_cells0-4={np.max(np.abs(values[:5])):.6e}; "
                f"maxabs_all={np.max(np.abs(values)):.6e}"
            )
        print(
            f"    accounting closure: maxabs(dH_direct - sum_terms)="
            f"{np.max(np.abs(direct_delta - reconstructed_delta)):.3e}; "
            f"scalar split delta mismatch="
            f"{np.max(np.abs(scalar_delta - scalar_split_delta)):.3e}; "
            f"H reconstruction maxabs={np.max(np.abs(base['closure'])):.3e}"
        )
        print(
            f"    dC_Lambda: cell0={delta_connection[0]:+.6e}; "
            f"maxabs_cells0-4={np.max(np.abs(delta_connection[:5])):.6e}; "
            f"maxabs_all={np.max(np.abs(delta_connection)):.6e}"
        )
        print(
            f"    dLambda: cell0={delta_lambda[0]:+.6e}; "
            f"d(metric connection): cell0={delta_metric_connection[0]:+.6e}; "
            f"split closure={np.max(np.abs(connection_split_closure)):.3e}"
        )
        print("    momentum residual change by signed term:")
        for name, values in momentum_after.items():
            delta = values - momentum_base[name]
            if name == "closure":
                print(
                    f"      {name}: maxabs(delta)={np.max(np.abs(delta)):.3e}"
                )
            else:
                print(
                    f"      d{name}: cell0={delta[0]:+.6e}; "
                    f"maxabs_cells0-4={np.max(np.abs(delta[:5])):.3e}; "
                    f"maxabs_all={np.max(np.abs(delta)):.3e}"
                )



if __name__ == "__main__":
    kernel = V55ProductionKernel()
    state, _, _ = discrete_consistent_state(
        kernel,
        resolution=N,
        r_max=R_MAX,
        amplitude=AMPLITUDE,
        width=7.0,
        D_amplitude=D_AMPLITUDE,
        include_radiation=True,
    )
    reference_grid = ref.Grid(N, R_MAX)
    report_curvature(state, reference_grid, "initial slice", initial_gate=True)
    report_momentum_checkpoint(state, "initial")
    one_step_budget(kernel, state)

    dt = CFL * state.grid.dr
    for target in TARGET_TIMES:
        while state.t < target - 1.0e-12:
            state = kernel.step(state, min(dt, target - state.t))
        report_curvature(state, reference_grid, "evolved slice")
        report_momentum_checkpoint(state, f"t={target:g}")
        one_step_budget(kernel, state)
