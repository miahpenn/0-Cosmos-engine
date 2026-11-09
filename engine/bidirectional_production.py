"""GEAR V5.5 bidirectional production coupling: the interface ledger.

There is deliberately NO added interface source term.
All local/COSMOS feedback is carried by the single metric + total Tmunu.
This module records the derived exchange ledger only.

The ledger is computed on states evolved by ``V55ProductionKernel``: the same
initial data, CMC gauge, matter projection and worldtube extraction as
production. It formerly evolved its own copy of the V5.5 system
(``reference_pirk_unified`` + ``true_cmc_reference``), whose initializer put the
archive operating point into the wrong variables (phi <- a = 33.8983,
Pi <- phi = 0.179055). That stack is kept for provenance only; see
docs/PROVENANCE_OLDER_V55_PATHS.md.
"""
import argparse
from pathlib import Path

import numpy as np

from .production_kernel import V55ProductionKernel

# The production worldtube is extracted at the cell nearest coordinate r = 10.
WORLDTUBE_COORDINATE_R = 10.0

LEDGER_KEYS = (
    "H_eff", "tau_rate", "R_sigma", "M_MS", "chi_sigma", "flux_T", "work_pR",
    "flux_shift", "rho_outer", "p_outer", "j_outer",
)


def ledger_observables(kernel, state, Rcoord=WORLDTUBE_COORDINATE_R):
    """Ledger row of a production state, taken from the production diagnostics.

    rho_inner, phi_outer and Pi_outer are the means of the production total
    density and COSMOS scalar over the same inner (r <= worldtube cell) and
    outer (r >= 0.8 r_max) regions the production diagnostics use.
    """
    if float(Rcoord) != WORLDTUBE_COORDINATE_R:
        raise ValueError(
            "the production worldtube is fixed at coordinate r = 10; "
            f"Rcoord={Rcoord!r} is not supported"
        )
    d = kernel.diagnostics(state, profiles=True)
    row = {key: float(d[key]) for key in LEDGER_KEYS}
    r = np.asarray(state.grid.centers, dtype=float)
    k = int(np.argmin(np.abs(r - WORLDTUBE_COORDINATE_R)))
    inner = r <= r[k]
    outer = r >= 0.8 * state.grid.r_max
    row["rho_inner"] = float(np.mean(np.asarray(d["total_rho"])[inner]))
    row["phi_outer"] = float(np.mean(np.asarray(state.scalars.phi)[outer]))
    row["Pi_outer"] = float(np.mean(np.asarray(state.scalars.Pi)[outer]))
    return row


def run(N=40, T=24, cfl=.06, Rcoord=WORLDTUBE_COORDINATE_R, **initial_data):
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=N, r_max=40.0, **initial_data)
    dt = cfl * float(state.grid.dr)
    nsteps = max(1, round(T / dt))
    dt = T / nsteps
    hist = []
    for n in range(nsteps + 1):
        z = ledger_observables(kernel, state, Rcoord)
        z["t"] = float(state.t)
        hist.append(z)
        if n == nsteps:
            break
        state = kernel.step(state, dt)
    out = {k: np.asarray([z[k] for z in hist]) for k in hist[0]}
    dM = np.gradient(out["M_MS"], out["t"], edge_order=2)
    rhs = out["flux_T"] + out["work_pR"] + out["flux_shift"]
    mask = (out["t"] >= 2) & (out["t"] <= out["t"][-1] - 2)
    if not np.any(mask):
        mask = np.ones_like(out["t"], dtype=bool)
    meta = {
        "N": N, "T": float(out["t"][-1]), "dt": dt,
        "max_Misner_Sharp_residual": float(np.max(abs(dM[mask] - rhs[mask]))),
        "rms_Misner_Sharp_residual": float(np.sqrt(np.mean((dM[mask] - rhs[mask]) ** 2))),
        "H_final": float(out["H_eff"][-1]),
        "tau_final": float(state.tau),
        "no_added_source_term": True,
        "evolution": "V55ProductionKernel",
    }
    return out, meta


def save_ledger_archive(out, N, output_dir):
    """Save a ledger using its recorded final time in the artifact filename.

    A filename is part of the evidence provenance: never label an arbitrary
    time span as T24 just because that is the CLI's current default.
    """
    if "t" not in out:
        raise ValueError("ledger output must include its sampled 't' history")
    times = np.asarray(out["t"])
    if (
        times.ndim != 1
        or times.size == 0
        or times.dtype.kind not in "iuf"
        or not np.all(np.isfinite(times))
    ):
        raise ValueError("ledger output 't' must be a non-empty finite numeric vector")
    final_time = float(times[-1])
    if final_time < 0.0:
        raise ValueError("ledger final time must be nonnegative")
    time_tag = f"{final_time:g}".replace("-", "m").replace(".", "p")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    target = output_dir / f"GEAR_V5_5_BIDIRECTIONAL_PRODUCTION_LEDGER_N{N}_T{time_tag}.npz"
    np.savez_compressed(target, **out)
    return target


if __name__ == "__main__":
    import json

    parser = argparse.ArgumentParser(
        description="Run the bounded GEAR V5.5 bidirectional research ledger."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("artifacts/gear_v55"),
        help="directory for NPZ outputs (default: ./artifacts/gear_v55)",
    )
    args = parser.parse_args()

    for N in (40, 60):
        out, meta = run(N)
        output_path = save_ledger_archive(out, N, args.output_dir)
        print(json.dumps({**meta, "output_path": str(output_path)}, indent=2))
