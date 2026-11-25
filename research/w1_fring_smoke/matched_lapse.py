"""Test matched lapse stencils on frozen W1 corrector operands."""
import argparse
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import solve_banded

from research.w1_fring_smoke import stage_momentum as replay
from research.w1_fring_smoke import capture_lapse_stages as capture


WIDE = np.longdouble
LOWER, UPPER = 3, 2


def lapse_system(grid, geometry, potential, order):
    """Banded -laplacian + Q, with production parity and outer normalization."""
    dr = WIDE(grid.dr)
    a, b, X = (getattr(geometry, f).astype(WIDE) for f in ("a", "b", "X"))
    da, db, dX = (grid.cell_derivative_fourth(getattr(geometry, f), parity=1)
                  for f in ("a", "b", "X"))
    inverse = X**2 / a
    slope = -.5 * da / a + db / b - dX / X + 2 / grid.centers.astype(WIDE)
    if order == 4:
        offsets, first, second = range(-2, 3), (1, -8, 0, 8, -1), (-1, 16, -30, 16, -1)
        denominator = 12
    else:
        offsets, first, second = range(-1, 2), (-1, 0, 1), (1, -2, 1)
        denominator = 2
    band = np.zeros((LOWER + UPPER + 1, grid.n), dtype=WIDE)
    for i in range(grid.n - 1):
        if order == 4 and i == grid.n - 2:
            stencil = zip(range(-3, 2), (-1, 6, -18, 10, 3), (-1, 4, 6, -20, 11))
        else:
            stencil = zip(offsets, first, second)
        for offset, d1, d2 in stencil:
            j = i + offset
            if j < 0:
                j = -1 - j
            if j == i:
                continue
            second_denominator = 1 if order == 2 else 12
            weight = -inverse[i] * (
                WIDE(d2) / (second_denominator * dr**2)
                + slope[i] * WIDE(d1) / (denominator * dr)
            )
            band[UPPER + i - j, j] += weight
            band[UPPER, i] -= weight
        band[UPPER, i] += potential[i]
    band[UPPER, -1] = 1
    return band


def equation_error(band, departure, potential, rhs):
    # Subtract neighbouring departures before applying the derivative weights.
    # This avoids subtracting nearly equal O(1/dr**2) terms of the full lapse.
    result = potential.astype(WIDE) * departure - rhs
    for offset in (-3, -2, -1, 1, 2):
        start = max(0, -offset)
        end = min(len(rhs) - 1, len(rhs) - offset)
        i = slice(start, end)
        j = slice(start + offset, end + offset)
        result[i] += band[UPPER - offset, j] * (departure[j] - departure[i])
    result[-1] = departure[-1]
    return result


def fields_from(snapshot, name, vacuum):
    return vacuum.VacuumState(**{f: snapshot[f"{name}_{f}"].copy() for f in capture.FIELDS})


def matter_source(snapshot, name):
    rho, pr, pt = (snapshot[f"{name}_{field}"] for field in ("rho", "pr", "pt"))
    return rho + pr + 2.0 * pt, pr - pt


def l3_rhs(vacuum, grid, geometry, source):
    trace, anisotropy = source
    rhs = vacuum.primary_l3_rhs(grid, geometry)
    rhs["K"] += 4.0 * np.pi * geometry.alpha * trace
    rhs["Aa"] -= (16.0 * np.pi / 3.0) * geometry.alpha * anisotropy
    return rhs


def lapse_variant(grid, geometry, source, target, order):
    potential = 1.5 * geometry.Aa * geometry.Aa + geometry.K * geometry.K / 3.0
    potential += 4.0 * np.pi * source[0]
    target_rhs = WIDE(target) - potential.astype(WIDE)
    target_rhs[-1] = 0
    band = lapse_system(grid, geometry, potential, order)
    stored_departure = geometry.alpha.astype(WIDE) - 1
    baseline = equation_error(band, stored_departure, potential, target_rhs)
    band_float = band.astype(float)

    def correction(rhs):
        return solve_banded((LOWER, UPPER), band_float, rhs.astype(float)).astype(WIDE)

    departure = stored_departure - correction(baseline)
    for _ in range(3):
        residual = equation_error(band, departure, potential, target_rhs)
        departure -= correction(residual)
    residual = equation_error(band, departure, potential, target_rhs)
    lapse = (1 + departure).astype(float)
    if not np.all(np.isfinite(lapse)) or np.min(lapse) <= 0 or lapse[-1] != 1:
        raise ValueError("matched operator did not produce a normalized positive lapse")
    residual_float = equation_error(band, lapse.astype(WIDE) - 1, potential, target_rhs)
    return lapse, {
        "alpha_change_max": float(np.max(abs(lapse - geometry.alpha))),
        "alpha_change_band_rms": replay.local_rms(grid.centers, lapse - geometry.alpha),
        "own_equation_residual_before_max": float(np.max(abs(baseline))),
        "own_equation_residual_refined_max": float(np.max(abs(residual))),
        "own_equation_residual_stored_max": float(np.max(abs(residual_float))),
    }


def budget(r, snapshot, delta_parts):
    # Read the isolated geometry increments against the exact baseline
    # momentum change, leaving both matter contributions in the signed sum.
    start = snapshot["initial_lapse_M"]
    base = snapshot["corrector_before_boundary_M"]
    update = base - start + sum(delta_parts.values())
    return {
        "rms_delta_M": replay.local_rms(r, update),
        "rms_final_M": replay.local_rms(r, base + sum(delta_parts.values())),
        "signed_delta_M": update.tolist(),
        "signed_final_M": (base + sum(delta_parts.values())).tolist(),
        "changed_terms": {key: values.tolist() for key, values in delta_parts.items()},
        "changed_term_rms": {key: replay.local_rms(r, values)
                             for key, values in delta_parts.items()},
    }


def analyze(snapshot, grid, vacuum):
    dt = float(snapshot["dt"])
    n = replay.LOCAL_CELLS
    r = grid.centers
    start = fields_from(snapshot, "start", vacuum)
    primary = fields_from(snapshot, "primary", vacuum)
    final = fields_from(snapshot, "corrector_copy", vacuum)
    np.testing.assert_array_equal(final.alpha, primary.alpha)
    source0, source1 = (matter_source(snapshot, s) for s in ("start", "primary"))

    def combine(l2_start, l2_final, l3_start, l3_primary):
        return {
            key: getattr(start, key) + .5 * dt * (
                l2_start[key] + l2_final[key] + l3_start[key] + l3_primary[key]
            )
            for key in capture.PRIMARY_FIELDS
        }

    l2_start = vacuum.primary_l2_rhs(grid, start)
    l2_final = vacuum.primary_l2_rhs(grid, final)
    l3_start, l3_primary = (l3_rhs(vacuum, grid, state, source)
                            for state, source in ((start, source0), (primary, source1)))
    for name, block, rhs in (
        ("start", "primary_l2", l2_start), ("corrector_copy", "primary_l2", l2_final),
        ("start", "primary_l3", l3_start), ("primary", "primary_l3", l3_primary),
    ):
        for key in capture.PRIMARY_FIELDS:
            np.testing.assert_array_equal(rhs[key], snapshot[f"{name}_{block}_{key}"])
    base = combine(l2_start, l2_final, l3_start, l3_primary)
    for key in capture.PRIMARY_FIELDS:
        np.testing.assert_array_equal(base[key][:n], snapshot[f"corrector_before_boundary_{key}"])

    da = lambda values: grid.cell_derivative_fourth(values, parity=1)
    chi_slope = -.5 * da(final.X) / final.X
    b_slope = da(final.b) / final.b

    def linear_l2(geometry, change):
        # L2 is linear in lapse on fixed geometry. Differentiate the small
        # change directly to avoid subtracting twice-rounded full-lapse Hessians.
        terms = vacuum.geometry_terms(grid, geometry)
        ap = da(geometry.a) / geometry.a
        bp = da(geometry.b) / geometry.b
        xp = da(geometry.X) / geometry.X
        slope = da(change)
        curve = grid.cell_second_derivative_fourth(change, parity=1)
        inv = geometry.X**2 / geometry.a
        radial = inv * (curve - slope * (.5 * ap - xp))
        lap = inv * (curve - slope * (.5 * ap - bp + xp - 2 / r))
        return {"K": -lap, "Aa": -radial + lap / 3 + change * (terms["Rr_r"] - terms["R"] / 3)}

    def momentum_increments(dK, dAa):
        return {
            "K_gradient": (-2 / 3 * da(dK))[:n],
            "Aa_prime": da(dAa)[:n],
            "Aa_over_r_term": (3 * dAa / r)[:n],
            "metric_coupling": (6 * dAa * chi_slope + 1.5 * dAa * b_slope)[:n],
        }

    def variant(order):
        start_new, primary_new, final_new = (g.copy() for g in (start, primary, final))
        start_new.alpha, res0 = lapse_variant(
            grid, start, source0, snapshot["start_target"], order,
        )
        primary_new.alpha, res1 = lapse_variant(
            grid, primary, source1, snapshot["primary_target"], order,
        )
        final_new.alpha = primary_new.alpha.copy()
        modified = combine(
            vacuum.primary_l2_rhs(grid, start_new), vacuum.primary_l2_rhs(grid, final_new),
            l3_rhs(vacuum, grid, start_new, source0), l3_rhs(vacuum, grid, primary_new, source1),
        )
        dK = modified["K"] - base["K"]
        dAa = modified["Aa"] - base["Aa"]
        increments = momentum_increments(dK, dAa)

        change0 = start_new.alpha - start.alpha
        change1 = primary_new.alpha - primary.alpha
        first = linear_l2(start, change0)
        last = linear_l2(final, change1)
        # In zero shift with KO off, L3 is pointwise proportional to lapse.
        stable = {
            key: .5 * dt * (first[key] + last[key]
                            + l3_start[key] * (change0 / start.alpha)
                            + l3_primary[key] * (change1 / primary.alpha))
            for key in capture.PRIMARY_FIELDS
        }
        stable_parts = momentum_increments(stable["K"], stable["Aa"])
        k_stage = (-2 / 3 * da(modified["K"] - start.K))[:n]
        # Cross-check linear momentum response against the ordinary readback.
        base_geometry = final.copy()
        base_geometry.K, base_geometry.Aa = base["K"], base["Aa"]
        shifted_geometry = final.copy()
        shifted_geometry.K, shifted_geometry.Aa = modified["K"], modified["Aa"]
        difference = (vacuum.momentum_constraint(grid, shifted_geometry)
                      - vacuum.momentum_constraint(grid, base_geometry))[:n]
        closure = float(np.max(abs(difference - sum(increments.values()))))
        return {
            "start_solve": res0,
            "primary_solve": res1,
            "K_update_band_rms": replay.local_rms(r[:n], k_stage),
            "K_only": budget(r[:n], snapshot, {"K_gradient": increments["K_gradient"]}),
            "K_and_Aa": budget(r[:n], snapshot, increments),
            "linear_response": {
                "K_only": budget(r[:n], snapshot, {"K_gradient": stable_parts["K_gradient"]}),
                "K_and_Aa": budget(r[:n], snapshot, stable_parts),
            },
            "direct_momentum_response_closure": closure,
        }

    return {
        "production": budget(r[:n], snapshot, {}),
        "production_K_update_band_rms": replay.local_rms(
            r[:n], (-2 / 3 * da(base["K"] - start.K))[:n]),
        "three_point_control": variant(2),
        "fourth_order_match": variant(4),
        "r": r[:n].tolist(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    replay.provenance()
    go, vacuum, _ = replay.pins.ad.vendor_modules()
    if vacuum.KO_EPSILON != 0:
        raise ValueError("expected KO off")
    with np.load(args.snapshot, allow_pickle=False) as snapshot:
        if str(snapshot["engine"]) != replay.pins.SOURCE or str(snapshot["vendor"]) != replay.pins.VENDOR:
            raise ValueError("unexpected engine or vendor pin")
        n = len(snapshot["start_K"])
        grid = go.SphericalCellGrid(n, n * float(snapshot["dr"]))
        report = analyze(snapshot, grid, vacuum)
        report.update({key: snapshot[key].item() for key in (
            "dr", "dt", "t_start", "t_end", "input_sha256", "engine", "vendor", "capture_sha256",
        )})
    if report["capture_sha256"] != replay.sha(Path(capture.__file__)):
        raise ValueError("unexpected capture script pin")
    report.update(analysis_sha256=replay.sha(Path(__file__)),
                  numpy=np.__version__, scipy=scipy.__version__,
                  wide_mantissa_bits=np.finfo(WIDE).nmant + 1)
    args.output.write_text(json.dumps(report, allow_nan=False, separators=(",", ":")) + "\n")


if __name__ == "__main__":
    main()
