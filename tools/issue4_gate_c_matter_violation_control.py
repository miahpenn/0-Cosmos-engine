"""Gate C, matter side: locked full-Bianchi predictor on a matter-violating control.

On the accepted states every operand of the projected-Bianchi predictor is
truncation-sized, so terms proportional to constraint violations are tested only
weakly. This control makes the matter part of H_total and M_total O(eta) instead:

  1. Restore the exact #420 t = 0.60 endpoint profile (checksum-verified by the
     pinned stage-source helper), exactly as tools/audit_exact_profile_momentum_tangent.py.
  2. Perturb the matter only, keeping parity and primitive regularity:
       dark matter: rest, energy_t, momentum_r all scaled by (1 + eta_rho g(r))
                    (same velocities and specific energies, more density)
       baryons:     momentum_r += eta_j * energy_t * (r/w) * g(r)   (odd)
     with g(r) = exp(-(r/w)^2). The geometry is not changed, so H_total and M_total
     acquire matter-sourced violations of order eta.
  3. Re-solve the CMC lapse for the perturbed state with the production kernel's own
     solver (the production RHS refreshes the lapse at the start stage, so the stored
     lapse must match it). Skipped for the eta = 0 baseline, which must reproduce the
     archived accepted fixture identity.
  4. Run the unmodified exact-profile tangent pipeline (production RHS direction,
     two-sided tangent with CMC re-solve, independent six-sector total-T ledger).
  5. Feed the result to tools/issue4_full_bianchi_evaluator.py, unchanged, through the
     accepted-fixture loader with every check active except the archive identity pins,
     which are replaced by this control's own computed digests and reported as such.
     The N = 60-only accepted reference regression is not executed for the control.

Nothing is evolved, no production file is modified, and admission stays false.

Usage:
  python tools/issue4_gate_c_matter_violation_control.py --input-dir 60=<profile dir> \
      [--input-dir 120=... --input-dir 240=...] --diagnostic-source <4aec88d0> \
      --helper-baseline <7f47669a> --work-dir <dir> [--eta-rho 0.05 --eta-j 0.02 --width 3]
"""
import argparse
import hashlib
import io
import json
import math
import sys
import zipfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# The build phase (production pipeline, repo vendor copy) and the evaluate phase
# (locked evaluator, which requires vendor modules from the pinned diagnostic-source
# checkout) run in separate processes; imports are therefore made inside each phase.
FIXED_RADII = (1.0, 2.0, 3.0, 5.0, 8.0)
TERMS = ("normal_alpha_K_M", "normal_lapse_gradient_H", "normal_H_gradient", "connection_divergence",
         "H_total", "M_total", "H_geom", "M_geom", "rho_total", "j_total", "C_r")


def sha(b):
    return hashlib.sha256(b).hexdigest()


def regions(values, r, signal=None):
    x = np.asarray(values, float)
    inner = r <= 10.0
    annulus = (r >= 5.0) & (r <= 10.0)
    out = {
        "max_abs_r_le_10": float(np.max(np.abs(x[inner]))),
        "radius_of_max_r_le_10": float(r[inner][int(np.argmax(np.abs(x[inner])))]),
        "rms_annulus_5_10": float(np.sqrt(np.mean(x[annulus] ** 2))),
        "cells_0_1_2": [float(x[i]) for i in range(3)],
        "fixed_radius": {f"{t:g}": float(np.interp(t, r, x)) for t in FIXED_RADII},
    }
    if signal is not None:
        s = np.asarray(signal, float)
        out["relative_rms_annulus_5_10"] = out["rms_annulus_5_10"] / max(
            float(np.sqrt(np.mean(s[annulus] ** 2))), np.finfo(float).tiny)
    return out


def make_perturbation(S, eta_rho, eta_j, width, record, resolve_lapse):
    from engine.matter_system import ConservedSpecies
    from engine.production_kernel import V55ProductionKernel
    from engine.v55_matter import V55MatterState

    original = S._restore_profile_state

    def restore(profile_arrays, resolution):
        state = original(profile_arrays, resolution)
        r = np.asarray(state.grid.centers, float)
        g = np.exp(-(r / width) ** 2)
        dm = state.matter.dark_matter
        scale = 1.0 + eta_rho * g
        new_dm = ConservedSpecies(np.asarray(dm.rest) * scale, np.asarray(dm.energy_t) * scale,
                                  np.asarray(dm.momentum_r) * scale)
        ba = state.matter.baryons
        new_ba = ConservedSpecies(np.asarray(ba.rest).copy(), np.asarray(ba.energy_t).copy(),
                                  np.asarray(ba.momentum_r) + eta_j * np.asarray(ba.energy_t) * (r / width) * g)
        if eta_rho != 0.0 or eta_j != 0.0:
            state.matter = V55MatterState(dark_matter=new_dm, baryons=new_ba,
                                          radiation=state.matter.radiation.copy(), charged_dust=None)
        kernel = V55ProductionKernel()
        alpha_old = np.asarray(state.geometry.alpha, float).copy()
        alpha, _ = kernel._solve_lapse_for_state(state.grid, state.geometry, state.scalars, state.matter, None)
        record["stored_minus_cmc_resolved_lapse_max_abs"] = float(np.max(np.abs(alpha_old - np.asarray(alpha))))
        record["lapse_resolved"] = bool(resolve_lapse)
        if resolve_lapse:
            state.geometry.alpha = np.asarray(alpha, float).copy()
        record["lapse_change_max_abs"] = float(np.max(np.abs(np.asarray(state.geometry.alpha) - alpha_old)))
        record["dm_rest_change_max_abs"] = float(np.max(np.abs(new_dm.rest - np.asarray(dm.rest))))
        record["baryon_momentum_change_max_abs"] = float(np.max(np.abs(new_ba.momentum_r - np.asarray(ba.momentum_r))))
        return state

    return original, restore


def build(n, input_dir, work, eta_rho, eta_j, width, label, resolve_lapse):
    """Phase 1: restore, perturb, run the unmodified exact-profile pipeline, zip it."""
    import tools.audit_exact_profile_momentum_tangent as X
    import tools.issue4_accepted_state_fixture as F

    S = X.stage_source
    record = {"label": label, "N": n, "eta_rho": eta_rho, "eta_j": eta_j, "width": width}
    original, restore = make_perturbation(S, eta_rho, eta_j, width, record, resolve_lapse)
    S._restore_profile_state = restore
    try:
        out_dir = work / f"{label}-N{n}"
        report = X._run(Path(input_dir), n, out_dir)
    finally:
        S._restore_profile_state = original
    members = {name: (out_dir / name).read_bytes() for name in ("SHA256SUMS.txt", "report.json", "tangent_arrays.npz")}
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data in members.items():
            z.writestr(name, data)
    zip_bytes = buf.getvalue()
    zip_path = work / f"{label}-N{n}.zip"
    zip_path.write_bytes(zip_bytes)
    archived = dict(F.EXPECTED[n])
    record["identity"] = {
        "input_state_digest": report["input"]["input_state_digest"],
        "input_tangent_digest": report["input"]["input_tangent_digest"],
        "matches_archived_accepted_state": report["input"]["input_state_digest"] == archived["archive_state"],
        "matches_archived_tangent_state": report["input"]["input_tangent_digest"] == archived["helper_state"],
        "report_sha256": sha(members["report.json"]), "npz_sha256": sha(members["tangent_arrays.npz"]),
        "zip_sha256": sha(zip_bytes),
        "archive_pins": "replaced by this run's own computed digests for loading; all other loader checks active",
    }
    record["tangent_closures"] = [
        {"epsilon": t["epsilon"], "total_from_geom_plus_sector_current_ok": t["closure"]["total_from_geom_plus_sector_current_ok"]}
        for t in report["semi_discrete_Mdot_tangent_by_epsilon"]]
    record["total_T_projection_crosschecks"] = [
        t["projection_crosscheck_passed"] for t in report["independent_total_stress_energy_divergence_by_epsilon"]]
    record["zip_path"] = str(zip_path)
    record["register"] = {"zip": sha(zip_bytes), "report": sha(members["report.json"]),
                          "npz": sha(members["tangent_arrays.npz"]),
                          "archive_state": report["input"]["input_state_digest"],
                          "helper_state": report["input"]["input_tangent_digest"]}
    return record


_BASELINE = {}


def evaluate(record, diag, helper):
    """Phase 2 (fresh process): the unchanged locked evaluator on the control fixture."""
    import tools.issue4_accepted_state_fixture as F
    import tools.issue4_full_bianchi_evaluator as E

    n, label = record["N"], record["label"]
    zip_path = Path(record["zip_path"])
    if sha(zip_path.read_bytes()) != record["register"]["zip"]:
        raise RuntimeError("control fixture changed between build and evaluate")
    archived = dict(F.EXPECTED[n])
    F.EXPECTED[n] = dict(record["register"])
    original_regression = E.assert_reference_regression
    E.assert_reference_regression = lambda result: {
        "executed": False, "reason": f"{label}: control state, no archived reference regression"}
    try:
        res = E.evaluate_accepted_fixture(zip_path, n, diag, helper)
    finally:
        E.assert_reference_regression = original_regression
        F.EXPECTED[n] = archived
    r = np.asarray(res["grid_centers"], float)
    record["status"] = res["status"]
    record["production_steps_run"] = res["production_steps_run"]
    record["zero_step_tripwire"] = res["zero_step_tripwire"]
    record["ricci_offset_crosscheck"] = res["ricci_offset_crosscheck"]
    st = res["signed_terms"]
    record["terms"] = {k: regions(st[k], r) for k in TERMS}
    keep = {k: np.asarray(st[k], float).copy() for k in ("rho_total", "j_total", "H_total", "M_total")}
    keep["measured"] = np.asarray(res["channels"]["eps1e-05"]["measured_tangent"], float).copy()
    if label == "baseline":
        _BASELINE[n] = keep
    elif n in _BASELINE:
        record["minus_baseline"] = {k: regions(keep[k] - _BASELINE[n][k], r) for k in keep}
    record["channels"] = {}
    for eps, chans in res["channels"].items():
        meas = np.asarray(chans["measured_tangent"], float)
        record["channels"][eps] = {"measured_tangent": regions(meas, r)}
        for name in ("pointwise", "finite_volume"):
            ch = chans[name]
            record["channels"][eps][name] = {
                "stress_source": regions(ch["stress_source"], r),
                "predictor": regions(ch["predictor"], r),
                "residual": regions(ch["residual_tangent_minus_predictor"], r, signal=meas),
            }
    a, b = list(res["channels"])[:2]
    record["epsilon_sensitivity_max_abs"] = {
        name: float(np.max(np.abs(np.asarray(res["channels"][a][name]["residual_tangent_minus_predictor"])
                                  - np.asarray(res["channels"][b][name]["residual_tangent_minus_predictor"]))))
        for name in ("pointwise", "finite_volume")}
    return record


def order(values, ns):
    return [None if x == 0 or y == 0 else round(math.log(abs(x) / abs(y)) / math.log(m / k), 2)
            for x, y, k, m in zip(values, values[1:], ns, ns[1:])]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--phase", choices=("build", "evaluate"), required=True)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--input-dir", action="append", help="N=path (build)")
    p.add_argument("--diagnostic-source", type=Path, help="(evaluate)")
    p.add_argument("--helper-baseline", type=Path, help="(evaluate)")
    p.add_argument("--work-dir", type=Path, help="(build)")
    p.add_argument("--eta-rho", type=float, default=0.5)
    p.add_argument("--eta-j", type=float, default=0.2)
    p.add_argument("--width", type=float, default=3.0)
    p.add_argument("--output", type=Path)
    a = p.parse_args()
    if a.phase == "build":
        inputs = dict(item.split("=", 1) for item in a.input_dir)
        ns = sorted(int(k) for k in inputs)
        a.work_dir.mkdir(parents=True, exist_ok=True)
        built = []
        for n in ns:
            built.append(build(n, inputs[str(n)], a.work_dir, 0.0, 0.0, a.width, "baseline", False))
            built.append(build(n, inputs[str(n)], a.work_dir, 0.0, 0.0, a.width, "baseline_lapse_resolved", True))
            built.append(build(n, inputs[str(n)], a.work_dir, a.eta_rho, a.eta_j, a.width, "matter_stored_lapse", False))
            built.append(build(n, inputs[str(n)], a.work_dir, a.eta_rho, a.eta_j, a.width, "matter_lapse_resolved", True))
        a.manifest.write_text(json.dumps({"resolutions": ns, "perturbation":
                                          {"eta_rho": a.eta_rho, "eta_j": a.eta_j, "width": a.width},
                                          "records": built}, indent=1) + "\n")
        for rec in built:
            print(rec["label"], rec["N"], json.dumps(rec["identity"]), rec.get("lapse_change_max_abs"),
                  rec.get("stored_minus_cmc_resolved_lapse_max_abs"), rec["tangent_closures"],
                  rec["total_T_projection_crosschecks"])
        return
    man = json.loads(a.manifest.read_text())
    ns = man["resolutions"]
    a.eta_rho, a.eta_j, a.width = (man["perturbation"][k] for k in ("eta_rho", "eta_j", "width"))
    rows = [evaluate(rec, a.diagnostic_source, a.helper_baseline) for rec in man["records"]]
    eps = "eps1e-05"
    orders = {}
    for label in ("baseline", "baseline_lapse_resolved", "matter_stored_lapse", "matter_lapse_resolved"):
        sel = [row for row in rows if row["label"] == label]
        orders[label] = {}
        for name in ("pointwise", "finite_volume"):
            get = lambda f: [f(row["channels"][eps][name]["residual"]) for row in sel]  # noqa: E731
            orders[label][name] = {
                "rms_annulus_5_10": order(get(lambda s: s["rms_annulus_5_10"]), ns),
                **{f"r_{t:g}": order(get(lambda s, t=t: s["fixed_radius"][f"{t:g}"]), ns) for t in FIXED_RADII},
            }
    rep = {"schema": "issue4_gate_c_matter_violation_control_v1", "read_only": True,
           "scientific_admission": False, "acceptance_oracle_satisfied": False, "resolutions": ns,
           "perturbation": {"eta_rho": a.eta_rho, "eta_j": a.eta_j, "width": a.width},
           "rows": rows, "residual_orders_eps1e-05": orders}
    if a.output:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(json.dumps(rep, indent=1) + "\n")
    print("GATE_C_MATTER_REPORT_BEGIN")
    print(json.dumps(rep, separators=(",", ":")))
    print("GATE_C_MATTER_REPORT_END")
    for row in rows:
        t = row["terms"]
        s = row["channels"][eps]["pointwise"]["residual"]
        f = row["channels"][eps]["finite_volume"]["residual"]
        mb = row.get("minus_baseline")
        if mb:
            print(f"   {row['label']} N={row['N']} injected: d(rho)={mb['rho_total']['max_abs_r_le_10']:.2e} d(j)={mb['j_total']['max_abs_r_le_10']:.2e} "
                  f"d(H_tot)={mb['H_total']['max_abs_r_le_10']:.2e} d(M_tot)={mb['M_total']['max_abs_r_le_10']:.2e} d(meas)={mb['measured']['max_abs_r_le_10']:.2e}")
        print(f"{row['label']:22s} N={row['N']:>3} |H_tot|={t['H_total']['max_abs_r_le_10']:.2e} |M_tot|={t['M_total']['max_abs_r_le_10']:.2e} "
              f"|meas|={row['channels'][eps]['measured_tangent']['max_abs_r_le_10']:.2e} "
              f"res pw max={s['max_abs_r_le_10']:.2e} r2={s['fixed_radius']['2']:+.2e} r5={s['fixed_radius']['5']:+.2e} ann={s['rms_annulus_5_10']:.2e} "
              f"| fv r2={f['fixed_radius']['2']:+.2e} ann={f['rms_annulus_5_10']:.2e} "
              f"| state==archive {row['identity']['matches_archived_accepted_state']} lapse change {row['lapse_change_max_abs']:.2e} "
              f"(stored-cmc {row['stored_minus_cmc_resolved_lapse_max_abs']:.2e})")
    print("orders:", json.dumps(orders))


if __name__ == "__main__":
    main()
