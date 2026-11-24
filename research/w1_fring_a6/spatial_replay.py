"""Spatial-refinement replays to t 0.015 (Aster's README-34 setup, rebuilt; Aster's scripts are not in the repo).

Ring arm, open edge, eight-node areal initializer, KO off, edge tolerance 1.943e-7, the corrected W1 static slice
(committed native_static.npz for dr 0.05/0.025; prepare_spatial.py output for dr 0.0125). Engine and vendor must
match #584's source, as in origin_trace.

--recon (B): centre reconstruction inside the replay only, nothing written to engine or vendor files:
  * every vendor function that differentiates Aa as an even scalar (_d(grid, Aa, 1): the Aa advection term, the
    momentum constraint and the Lambda equation) uses the vendor's regularized_aa_prime(grid, Aa) instead;
  * after the initial data and after every accepted step, the vendor's enforce_algebraic_regularity is applied and
    only its first-cell values (a, b, Aa, beta, Lambda, B) are copied back.

Records per accepted step: max |M| over the grid and its cell, |M| at cell 0, max |M| over the band R 0.05-0.25,
max |H|. At t 0.015 it saves the full-grid M split (Aa', 3Aa/r, (3/2) Aa b'/b, 6 Aa chi', -(2/3) K', -8 pi j).

Usage (repo root): python -m research.w1_fring_a6.spatial_replay DR DT OUTDIR [--recon] [--slice NPZ]
"""
import argparse
import hashlib
import inspect
import json
import subprocess
from pathlib import Path

import numpy as np

from engine import scalar_system as ss, v55_initial as vi, v55_pirk_adapter as ad
from engine.adaptive_step import advance_with_radiation_admissibility_retries
from engine.production_kernel import V55ProductionKernel
from engine.run_production import assert_finite_state
from engine.stress_energy import assemble_total_stress_energy
from research.w1_fring_smoke import origin_trace as ot

T_END = .015
EDGE = 1.943e-7
H0_CHART = .0056689681106181795
STATIC = Path("results/research/w1_static_slice/native_static.npz")
BAND = (.05, .25)


def check_source():
    unchanged = subprocess.run(["git", "diff", "--quiet", ot.SOURCE, "--", "engine", "vendor"], check=False).returncode == 0
    vendor = subprocess.check_output(["git", "-C", "vendor/bb-palatini-unified-r0", "rev-parse", "HEAD"], text=True).strip()
    if not unchanged or vendor != ot.VENDOR:
        raise ValueError("engine or vendor differ from the smoke source")


def patch_aa_prime(vac):
    """Swap _d(grid, Aa, 1) for regularized_aa_prime(grid, Aa) in the vendor functions that use it."""
    patched = []
    for name, fn in list(vars(vac).items()):
        if inspect.isfunction(fn) and fn.__module__ == vac.__name__ and name != "regularized_aa_prime":
            src = inspect.getsource(fn)
            if "_d(grid, Aa, 1)" in src:
                exec(compile(src.replace("_d(grid, Aa, 1)", "regularized_aa_prime(grid, Aa)"), vac.__file__, "exec"),
                     vars(vac))
                patched.append(name)
    return patched


def project_first_cell(vac, grid, geometry):
    proj = vac.enforce_algebraic_regularity(grid, geometry)
    for k in ("a", "b", "Aa", "beta", "Lambda", "B"):
        getattr(geometry, k)[0] = getattr(proj, k)[0]


def split(grid, g, scal, mat, vac):
    d = lambda v, p: grid.cell_derivative_fourth(v, parity=p)
    r = np.asarray(grid.centers)
    Aap = vac.regularized_aa_prime(grid, g.Aa) if RECON else d(g.Aa, 1)
    parts = {"Aa_prime": Aap, "three_Aa_over_r": 3*g.Aa/r, "Aa_b_prime_over_b": 1.5*g.Aa*d(g.b, 1)/g.b,
             "six_Aa_chi_prime": 6*g.Aa*(-.5*d(g.X, 1)/g.X), "K_part": -(2/3)*d(g.K, 1),
             "j_part": -8*np.pi*assemble_total_stress_energy(grid, g, scal, mat).j}
    parts["M"] = vac.constraints(grid, g)["momentum"] + parts["j_part"]
    return parts


RECON = False


def main():
    global RECON
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("dr", type=float)
    p.add_argument("dt", type=float)
    p.add_argument("out")
    p.add_argument("--recon", action="store_true")
    p.add_argument("--slice", default=None)
    a = p.parse_args()
    RECON = a.recon
    check_source()
    go, vac, _ = ad.vendor_modules()
    if vac.KO_EPSILON != 0:
        raise ValueError("expected KO off")
    patched = patch_aa_prime(vac) if RECON else []
    n = round(60/a.dr)
    grid = go.SphericalCellGrid(n, 60.)
    key = f"dr{a.dr:g}"
    src = Path(a.slice) if a.slice else STATIC
    with np.load(src) as z:
        values = [np.array(z[key+"_"+k], dtype=float) for k in ("S", "PS", "D", "PD")]
    slice_sha = hashlib.sha256(src.read_bytes()).hexdigest()
    values.extend((np.full(n, .179055), np.full(n, .00341913)))
    label = (f"W1 corrected slice (A3/A4), static time, grid {key} from {src.name} sha256 {slice_sha}; "
             f"H0_phase={H0_CHART}; theta0=0; snapshot phase omega_u H0 r^2/2 per solver node; H0_chart={H0_CHART}")
    supplied = vi.SuppliedSlice(ss.ScalarFields(*values), np.asarray(grid.centers).copy(), H0_CHART, EDGE, label)
    kernel = V55ProductionKernel()
    state = kernel.initialize(resolution=n, r_max=60., outer_edge="open_edge", eight_node_areal=True,
                              initial_structure="supplied_slice", supplied_slice=supplied)
    r = np.asarray(grid.centers)
    band = (r >= BAND[0]) & (r <= BAND[1])
    rows = []

    def record():
        assert_finite_state(state)
        stress = assemble_total_stress_energy(grid, state.geometry, state.scalars, state.matter)
        raw = vac.constraints(grid, state.geometry)
        M = raw["momentum"] - 8*np.pi*stress.j
        H = raw["hamiltonian"] - 16*np.pi*stress.rho
        i = int(np.argmax(np.abs(M)))
        rows.append({"t": state.t, "max_abs_M": float(abs(M[i])), "max_cell": i, "M_cell0": float(M[0]),
                     "band_max_abs_M": float(np.max(np.abs(M[band]))), "max_abs_H": float(np.max(np.abs(H))),
                     "M_cells_0_5": M[:6].tolist()})

    if RECON:
        project_first_cell(vac, grid, state.geometry)
    record()
    steps = round(T_END/a.dt)
    if not np.isclose(steps*a.dt, T_END, atol=1e-17, rtol=0):
        raise ValueError("end time must be a multiple of dt")
    for _ in range(steps):
        adv = advance_with_radiation_admissibility_retries(kernel, state, a.dt)
        if adv.retry_count:
            raise RuntimeError("retry would break matched times")
        state = adv.state
        if RECON:
            project_first_cell(vac, grid, state.geometry)
        record()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tag = f"{'B' if RECON else 'A'}-{key}-dt{a.dt:g}"
    parts = split(grid, state.geometry, state.scalars, state.matter, vac)
    total = sum(parts[k] for k in ("Aa_prime", "three_Aa_over_r", "Aa_b_prime_over_b", "six_Aa_chi_prime",
                                    "K_part", "j_part"))
    np.savez(out/f"{tag}-split.npz", r=r, **parts)
    (out/f"{tag}.json").write_text(json.dumps({
        "dr": a.dr, "dt": a.dt, "recon": RECON, "patched": patched, "slice": str(src), "slice_sha256": slice_sha,
        "edge_tolerance": EDGE, "steps": steps, "t_end": state.t,
        "split_adds_back_max_abs": float(np.max(np.abs(total - parts["M"]))),
        "split_adds_back_rel_max_M": float(np.max(np.abs(total - parts["M"]))/np.max(np.abs(parts["M"]))),
        "rows": rows}, indent=1))
    print(tag, "done", rows[-1]["max_abs_M"], rows[-1]["band_max_abs_M"], abs(rows[-1]["M_cell0"]), flush=True)


if __name__ == "__main__":
    main()
