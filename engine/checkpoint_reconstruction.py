"""Dataset-scoped reconstruction of the archived N=320 Dbase checkpoints.

Diagnostic replay only: this module loads saved state and evaluates the
existing production diagnostics. It never evolves the state. Exact timestamp
matching is approved only for the archived artifact identified below.
"""
from __future__ import annotations

import math
import struct
from pathlib import Path

import numpy as np

from engine import v55_pirk_adapter as adapter
from engine.matter_system import ConservedSpecies
from engine.production_kernel import ProductionState, V55ProductionKernel
from engine.scalar_system import ScalarFields
from engine.v55_matter import V55MatterState

ARTIFACT_ID = 11638114258
ARTIFACT_SHA256 = "85b5fdfd379fe793cb0b4c70efd21fa43e598127ba45193f58b131380af61ba7"
SOURCE_COMMIT = "180324199196d7a5eacb7477187ce7e22ff6db64"
DATASET_SCOPE = "GitHub Actions artifact 11638114258 only; do not generalize exact timestamp matching"
TARGETS = (4.0, 8.0, 12.0, 16.0, 20.0, 22.5)
L2_FIELDS = (
    "hamiltonian_l2_inner",
    "hamiltonian_l2_outer",
    "momentum_l2_inner",
    "momentum_l2_outer",
)
# Existing full-grid outer-constraint bookkeeping criterion.
DECOMP_TOL = 1.0e-12
MATCHED, UNMATCHED, AMBIGUOUS = "MATCHED", "UNMATCHED", "AMBIGUOUS"


def load_checkpoint_npz(path: str | Path, r_max: float) -> ProductionState:
    """Reconstruct the saved ProductionState without changing any field."""
    grid_ops, vacuum, _ = adapter.vendor_modules()
    with np.load(Path(path), allow_pickle=False) as z:
        r = np.asarray(z["r"], dtype=float)
        grid = grid_ops.SphericalCellGrid(int(r.size), float(r_max))
        if not np.array_equal(np.asarray(grid.centers, dtype=float), r):
            raise ValueError("stored radial centers do not match rebuilt grid")

        def arr(key: str) -> np.ndarray:
            value = np.array(z[key], dtype=float, copy=True)
            if not np.all(np.isfinite(value)):
                raise FloatingPointError(f"non-finite values in checkpoint field {key}")
            return value

        geometry = vacuum.VacuumState(
            a=arr("a"), b=arr("b"), X=arr("X"), alpha=arr("alpha"),
            beta=arr("beta"), Aa=arr("Aa"), K=arr("K"),
            Lambda=arr("Lambda"), B=arr("B"),
        )
        scalars = ScalarFields(
            arr("S"), arr("PS"), arr("D"), arr("PD"), arr("phi"), arr("Pi"),
        )
        matter = V55MatterState(
            ConservedSpecies(arr("dm_rest"), arr("dm_energy_t"), arr("dm_momentum_r")),
            ConservedSpecies(arr("baryon_rest"), arr("baryon_energy_t"), arr("baryon_momentum_r")),
            ConservedSpecies(arr("radiation_rest"), arr("radiation_energy_t"), arr("radiation_momentum_r")),
        )
        t, tau, e_folds = (float(z[k]) for k in ("t", "tau", "e_folds"))
        if not all(math.isfinite(v) for v in (t, tau, e_folds)):
            raise FloatingPointError("non-finite checkpoint scalar metadata")

    return ProductionState(
        grid=grid, geometry=geometry, scalars=scalars, matter=matter,
        t=t, tau=tau, e_folds=e_folds,
    )


def reconstruct_l2(state: ProductionState) -> dict[str, float]:
    """Calculate the four frozen regional L2 values using production diagnostics."""
    diagnostics = V55ProductionKernel().diagnostics(state, profiles=False)
    result = {name: float(diagnostics[name]) for name in L2_FIELDS}
    result["hamiltonian_decomposition_error_max"] = float(
        diagnostics["hamiltonian_decomposition_error_max"]
    )
    if not all(math.isfinite(v) for v in result.values()):
        raise FloatingPointError("non-finite reconstructed diagnostic")
    return result


def match_ledger_row_by_checkpoint_time(checkpoint_t: float, rows: list[dict]) -> dict:
    """Match exact stored t equality for this archived dataset only.

    No tolerance, nearest-row, or nominal-target fallback is permitted.
    """
    hits = [
        (index, row) for index, row in enumerate(rows)
        if float(row["t"]) == float(checkpoint_t)
    ]
    if not hits:
        return {
            "status": UNMATCHED, "row_index": None, "row": None,
            "offset": None, "candidates": [],
        }
    if len(hits) != 1:
        return {
            "status": AMBIGUOUS, "row_index": None, "row": None,
            "offset": None, "candidates": [float(row["t"]) for _, row in hits],
        }
    index, row = hits[0]
    return {
        "status": MATCHED, "row_index": index, "row": row,
        "offset": float(row["t"]) - float(checkpoint_t),
        "candidates": [float(row["t"])],
    }


def compare_field(reconstructed: float, recorded: float | None) -> dict:
    """Test exact reconstruction; report absolute/relative errors without a fitted tolerance."""
    if recorded is None:
        return {
            "outcome": "UNDEFINED", "reason": "missing reference",
            "abs_err": None, "rel_err": None, "exact_equal": False,
            "float64_bitwise_equal": False,
        }
    reconstructed = float(reconstructed)
    recorded = float(recorded)
    if not math.isfinite(reconstructed) or not math.isfinite(recorded):
        return {
            "outcome": "FAIL", "reason": "non-finite value",
            "abs_err": None, "rel_err": None, "exact_equal": False,
            "float64_bitwise_equal": False,
        }
    if recorded == 0.0:
        return {
            "outcome": "UNDEFINED", "reason": "reference is exactly zero",
            "abs_err": abs(reconstructed - recorded), "rel_err": None,
            "exact_equal": reconstructed == recorded,
            "float64_bitwise_equal": struct.pack(">d", reconstructed) == struct.pack(">d", recorded),
        }
    abs_err = abs(reconstructed - recorded)
    rel_err = abs_err / abs(recorded)
    bitwise_equal = struct.pack(">d", reconstructed) == struct.pack(">d", recorded)
    exact_equal = reconstructed == recorded
    # This is reproduction of archived diagnostic values, not an accuracy
    # tolerance for physics. Any nonzero difference must be reported as FAIL.
    return {
        "outcome": "PASS" if bitwise_equal else "FAIL",
        "reason": None if bitwise_equal else "float64 values differ",
        "abs_err": abs_err,
        "rel_err": rel_err,
        "exact_equal": exact_equal,
        "float64_bitwise_equal": bitwise_equal,
    }
