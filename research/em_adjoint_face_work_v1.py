"""Research-only adjoint map for the discrete radial-electric energy functional.

This module derives the face-charge work operator implied by
U = (dr/2) * sum_i w_i * ((Q_i + Q_{i+1})/2)**2.
It does not select a physical carrier, alter evolution, or include the
separate geometry-dependent metric tangent.
"""
from __future__ import annotations

import math

import numpy as np


FOUR_PI = 4.0 * math.pi
ROUND_OFF_ULPS = 128.0


def _finite_real_scalar(value, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value, (int, float, np.integer, np.floating)
    ):
        raise ValueError(f"{name} must be a finite real scalar")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be a finite real scalar")
    return result


def _finite_real_vector(values, name: str, *, min_size: int = 2) -> np.ndarray:
    try:
        raw = np.asarray(values)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real one-dimensional array") from exc
    if raw.ndim != 1 or raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be a finite real one-dimensional array")
    result = np.array(raw, dtype=float, copy=True)
    if result.size < min_size:
        raise ValueError(f"{name} must have at least {min_size} entries")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return result


def _validated_operands(dr, field_energy_weight, charge_faces):
    step = _finite_real_scalar(dr, "dr")
    if step <= 0.0:
        raise ValueError("dr must be positive")
    weight = _finite_real_vector(field_energy_weight, "field_energy_weight")
    if np.any(weight <= 0.0):
        raise ValueError("field_energy_weight must be strictly positive")
    charge = _finite_real_vector(
        charge_faces, "charge_faces", min_size=weight.size + 1
    )
    if charge.size != weight.size + 1:
        raise ValueError("charge_faces must have one more entry than cell weights")
    if charge[0] != 0.0:
        raise ValueError("regular-centre enclosed charge must be exactly zero")
    return step, weight, charge


def radial_electric_energy(dr, field_energy_weight, charge_faces) -> float:
    """Evaluate the existing cell-weighted energy functional from face Q.

    The cell charge is the adjacent-face average. The weight is passed in from
    the chosen geometry, so this function does not silently derive or modify
    the geometry convention.
    """
    step, weight, charge = _validated_operands(
        dr, field_energy_weight, charge_faces
    )
    charge_cell = 0.5 * (charge[:-1] + charge[1:])
    result = 0.5 * step * math.fsum(
        float(term) for term in weight * charge_cell**2
    )
    if not math.isfinite(result):
        raise ArithmeticError("field-energy functional is non-finite")
    return result


def charge_rate_from_rest_mass_flux(charge_to_rest_mass, rest_mass_flux_faces):
    """Map the declared constant-q/m rest-mass face flux to dQ/dt.

    RadialMaxwellState stores +4*pi*(q/m)*F as its charge-current flux and
    advances Q with a minus sign. This helper returns dQ/dt = -4*pi*(q/m)*F.
    It encodes that software convention only; it does not derive the carrier.
    """
    eta = _finite_real_scalar(charge_to_rest_mass, "charge_to_rest_mass")
    flux = _finite_real_vector(
        rest_mass_flux_faces, "rest_mass_flux_faces", min_size=3
    )
    if flux[0] != 0.0:
        raise ValueError("regular-centre rest-mass face flux must be exactly zero")
    rate = -FOUR_PI * eta * flux
    if not np.all(np.isfinite(rate)):
        raise ArithmeticError("face-charge rate is non-finite")
    return rate


def face_charge_energy_tangent(
    dr, field_energy_weight, charge_faces, charge_rate_faces
) -> dict[str, object]:
    """Return a direct cell tangent and its independently summed face adjoint.

    If q_i=(Q_i+Q_{i+1})/2 and U=(dr/2)*sum_i w_i*q_i**2, then
      dU/dQ_0 = dr*w_0*q_0/2,
      dU/dQ_k = dr*(w_{k-1}*q_{k-1}+w_k*q_k)/2 (interior),
      dU/dQ_n = dr*w_{n-1}*q_{n-1}/2.
    Thus dU/dt = sum_k (dU/dQ_k)*(dQ_k/dt). The endpoints are retained
    explicitly; for a regular centre dQ_0/dt=0, but the outer-face term need
    not vanish. The metric-weight tangent is deliberately not included.
    """
    step, weight, charge = _validated_operands(
        dr, field_energy_weight, charge_faces
    )
    rate = _finite_real_vector(
        charge_rate_faces, "charge_rate_faces", min_size=weight.size + 1
    )
    if rate.size != weight.size + 1:
        raise ValueError("charge_rate_faces must have one more entry than cells")
    if rate[0] != 0.0:
        raise ValueError("regular-centre face-charge rate must be exactly zero")

    charge_cell = 0.5 * (charge[:-1] + charge[1:])
    charge_cell_rate = 0.5 * (rate[:-1] + rate[1:])
    cell_terms = step * weight * charge_cell * charge_cell_rate

    face_gradient = np.empty(charge.size, dtype=float)
    face_gradient[0] = 0.5 * step * weight[0] * charge_cell[0]
    face_gradient[1:-1] = 0.5 * step * (
        weight[:-1] * charge_cell[:-1] + weight[1:] * charge_cell[1:]
    )
    face_gradient[-1] = 0.5 * step * weight[-1] * charge_cell[-1]
    face_terms = face_gradient * rate

    if not (
        np.all(np.isfinite(charge_cell))
        and np.all(np.isfinite(charge_cell_rate))
        and np.all(np.isfinite(cell_terms))
        and np.all(np.isfinite(face_gradient))
        and np.all(np.isfinite(face_terms))
    ):
        raise ArithmeticError("face-charge energy tangent operands overflowed")

    cell_tangent = math.fsum(float(term) for term in cell_terms)
    face_tangent = math.fsum(float(term) for term in face_terms)
    residual = cell_tangent - face_tangent
    scale = math.fsum(abs(float(term)) for term in cell_terms) + math.fsum(
        abs(float(term)) for term in face_terms
    )
    bound = ROUND_OFF_ULPS * np.finfo(float).eps * max(
        scale, np.finfo(float).tiny
    )
    ratio = 0.0 if residual == 0.0 else abs(residual) / bound
    if not math.isfinite(ratio) or ratio > 1.0:
        raise ArithmeticError(
            "cell/face adjoint identity exceeded operand-scaled roundoff: "
            f"residual={residual:.17e}, bound={bound:.17e}, ratio={ratio:.17e}"
        )

    return {
        "schema_version": 1,
        "field_energy": radial_electric_energy(step, weight, charge),
        "cell_charge": charge_cell.tolist(),
        "cell_charge_rate": charge_cell_rate.tolist(),
        "face_energy_gradient": face_gradient.tolist(),
        "charge_rate_faces": rate.tolist(),
        "cell_tangent_terms": cell_terms.tolist(),
        "face_tangent_terms": face_terms.tolist(),
        "cell_tangent": cell_tangent,
        "face_adjoint_tangent": face_tangent,
        "face_tangent_center": float(face_terms[0]),
        "face_tangent_interior": math.fsum(
            float(term) for term in face_terms[1:-1]
        ),
        "face_tangent_outer": float(face_terms[-1]),
        "closure_residual": residual,
        "operand_scale": scale,
        "roundoff_bound": bound,
        "bound_ratio": ratio,
        "closure_passed": True,
        "scope": (
            "charge/current tangent of the discrete radial-electric energy "
            "functional at fixed geometry; excludes the metric-weight tangent"
        ),
    }
