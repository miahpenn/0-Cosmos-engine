#!/usr/bin/env python3
"""Recompute H/M/connection/determinant witnesses from a frozen replay ZIP.

This is a read-only audit of 27 saved accepted-state profiles. It imports no
engine code and advances no state. By default, it verifies the exact ZIP digest
and internal SHA256SUMS before recomputation.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import math
from pathlib import Path
import sys
import zipfile

import numpy as np

EXPECTED_ZIP_SHA256 = "5c3dee4cbe62fb410ed7bcdf97707385019e00bae4f925ab51d4242b42933d09"
SOURCE_BASELINE = "29d828b4d72b017c12a20b6dc90549f2511f025c"
VENDOR_PIN = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
CAPTURE_HEAD = "a854c3d6eafde1a1bae8066b7cc212bac7511ffd"
CAPTURE_SCHEMA = "radiation_explicit_predictor_stage_capture_v5_retry_identity_and_recomputable_independent_l3"
ROOT = "artifacts/current-main-radiation-frozen-state-replay"
KAPPA = 8.0 * math.pi
GEOMETRY_FIELDS = ("a", "b", "X", "Aa", "K", "Lambda", "alpha")
SCALARS = ("S", "PS", "D", "PD", "phi", "Pi")
SPECIES = ("dark_matter", "baryons", "radiation")
CONSERVED = ("rest", "energy_t", "momentum_r")


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def verify_zip(blob: bytes, *, allow_other_hash: bool = False):
    actual = sha256(blob)
    if not allow_other_hash and actual != EXPECTED_ZIP_SHA256:
        raise ValueError(f"wrong exact replay ZIP hash: {actual}")
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        names = z.namelist()
        def member(name):
            full = f"{ROOT}/{name}"
            if full not in names:
                raise ValueError(f"artifact member missing: {full}")
            return z.read(full)
        manifest = member("SHA256SUMS.txt").decode("utf-8")
        checked = {}
        for line in manifest.splitlines():
            line = line.strip()
            if not line:
                continue
            expected, name = line.split(maxsplit=1)
            name = name.lstrip("* ")
            hits = [n for n in names if n.rsplit("/", 1)[-1] == name]
            if len(hits) != 1:
                raise ValueError(f"manifest filename is ambiguous or missing: {name}")
            found = sha256(z.read(hits[0]))
            if found != expected:
                raise ValueError(f"internal manifest mismatch: {name}")
            checked[name] = found
        expected_names = {
            "report.json", "source_boundary_failure_accepted_state.npz",
            "source_boundary_failure_accepted_state.json", "stage_capture.json",
            "stage_capture_profiles.npz",
        }
        if set(checked) != expected_names:
            raise ValueError(f"unexpected manifest entries: {sorted(checked)}")
        report = json.loads(member("report.json"))
        capture = json.loads(member("stage_capture.json"))
        profile_bytes = member("stage_capture_profiles.npz")
    prov = report.get("provenance", {})
    if prov.get("candidate_head") != CAPTURE_HEAD:
        raise ValueError("observer head does not match the audited exact head")
    if prov.get("expected_base_main_commit") != SOURCE_BASELINE:
        raise ValueError("production source baseline changed")
    if prov.get("vendor_pin") != VENDOR_PIN or prov.get("expected_vendor_pin") != VENDOR_PIN:
        raise ValueError("vendor pin changed")
    if capture.get("schema") != CAPTURE_SCHEMA:
        raise ValueError("capture schema changed")
    if capture.get("accepted_state_count") != 27 or len(capture.get("accepted_state_history", [])) != 27:
        raise ValueError("expected exactly 27 saved accepted snapshots")
    if capture.get("geometry_source_term_record_count") != 703:
        raise ValueError("unexpected geometry-source capture count")
    return actual, checked, report, capture, np.load(io.BytesIO(profile_bytes), allow_pickle=False)


def d4(values, h: float, parity: int, order: int = 1):
    """Pinned vendor cell_derivative_fourth / cell_second_derivative_fourth."""
    u = np.asarray(values, dtype=float)
    n = u.size
    if u.ndim != 1 or n < 5 or parity not in (-1, 1):
        raise ValueError("invalid derivative input")
    ext = np.empty(n + 4, dtype=float)
    ext[2:-2] = u
    ext[1] = parity * u[0]
    ext[0] = parity * u[1]
    out = np.empty_like(u)
    for i in range(n - 2):
        s = ext[i:i + 5]
        if order == 1:
            out[i] = (s[0] - 8*s[1] + 8*s[3] - s[4]) / (12*h)
        elif order == 2:
            out[i] = (-s[0] + 16*s[1] - 30*s[2] + 16*s[3] - s[4]) / (12*h*h)
        else:
            raise ValueError("only first and second derivatives are supported")
    for i in (n - 2, n - 1):
        indices = np.arange(n - 5, n)
        offsets = indices.astype(float) - float(i)
        matrix = np.vstack([offsets**power for power in range(5)])
        target = np.zeros(5)
        target[order] = float(math.factorial(order))
        weights = np.linalg.solve(matrix, target)
        out[i] = float(weights @ u[indices]) / h**order
    return out


def geometry_terms(r, a, b, X, alpha, Lambda, h):
    """Same frozen reference-metric Ricci formulas as the pinned vendor."""
    ap = d4(a, h, 1)
    bp = d4(b, h, 1)
    app = d4(a, h, 1, 2)
    bpp = d4(b, h, 1, 2)
    Xp = d4(X, h, 1)
    Xpp = d4(X, h, 1, 2)
    chip = -0.5 * Xp / X
    chipp = -0.5 * (Xpp / X - (Xp / X)**2)
    Lp = d4(Lambda, h, -1)
    metric_lambda = (1.0 - a / b) / r**2
    inv = X**2 / a
    radial_bracket = (
        0.5 * app/a - a*Lp - 0.75*(ap/a)**2 + 0.5*(bp/b)**2
        - 0.5*Lambda*ap + ap/(r*b)
        + 2*metric_lambda*(1 + r*bp/b) + 4*chipp
        - 2*chip*(ap/a - bp/b - 2/r)
    )
    scalar_bracket = (
        0.5*app/a + bpp/b - a*Lp - (ap/a)**2 + 0.5*(bp/b)**2
        + 2*(3 - a/b)*bp/(r*b) + 4*metric_lambda
        + 8*(chipp + chip**2)
        - 8*chip*(0.5*ap/a - bp/b - 2/r)
    )
    return {"R": -inv*scalar_bracket, "Rr_r": -inv*radial_bracket}


def stress_projection(r, a, b, X, scalar, matter, h):
    """Reconstruct TotalStressEnergy using production primitive/projection maps."""
    n = len(r)
    inv = X**2 / a
    rho = np.zeros(n)
    pr = np.zeros(n)
    pt = np.zeros(n)
    current = np.zeros(n)

    # scalar_system.scalar_projection
    for field, momentum, potential, scale in (
        (scalar["S"], scalar["PS"], 0.5*scalar["S"]**2, 1.0),
        (scalar["D"], scalar["PD"], -0.5*scalar["D"]**2, 1.0),
        (scalar["phi"], scalar["Pi"],
         8.242522415500654e-5*(np.exp(-scalar["phi"]) - 0.10), KAPPA),
    ):
        fp = d4(field, h, 1)
        grad2 = inv * fp*fp
        rho += (0.5*momentum**2 + 0.5*grad2 + potential) / scale
        pr += (0.5*momentum**2 + 0.5*grad2 - potential) / scale
        pt += (0.5*momentum**2 - 0.5*grad2 - potential) / scale
        current += -momentum*fp / scale

    sqrt_gamma = np.sqrt(a) * b * r**2 / X**3
    gamma_rr = a / X**2
    gamma_rr_inv = X**2 / a

    # Production recover_dust -> species_projection, for DM and baryons.
    for name in ("dark_matter", "baryons"):
        conserved = matter[name]
        D = conserved["rest"] / sqrt_gamma
        S_r = conserved["momentum_r"] / sqrt_gamma
        s2 = gamma_rr_inv * S_r*S_r
        v2 = s2 / (D*D + s2)
        v_r = np.sign(S_r) * np.sqrt(np.maximum(v2 / gamma_rr, 0.0))
        W = 1.0 / np.sqrt(1.0 - v2)
        rest_density = D / W
        rho += rest_density*W*W / KAPPA
        pr += rest_density*W*W*v2 / KAPPA
        current += rest_density*W*W*gamma_rr*v_r / KAPPA

    # Production recover_radiation's pressure solve, including its exact stopping rule.
    conserved = matter["radiation"]
    E = conserved["energy_t"] / sqrt_gamma
    S_r = conserved["momentum_r"] / sqrt_gamma
    S2 = gamma_rr_inv * S_r*S_r
    pressure = np.empty(n)
    v_r = np.empty(n)
    for i in range(n):
        if S2[i] == 0.0:
            pressure[i] = E[i] / 3.0
            v_r[i] = 0.0
            continue
        lo, hi = 0.0, E[i] / 3.0
        def residual(p):
            v2 = S2[i] / (E[i] + p)**2
            correction = 4*p*v2/(1-v2)
            return (3*p - E[i]) + correction, v2
        f_lo, _ = residual(lo)
        f_hi, _ = residual(hi)
        if f_hi < 0.0:
            hi = np.nextafter(hi, math.inf)
            f_hi, _ = residual(hi)
        if not math.isfinite(f_hi) or f_lo > 0.0 or f_hi < 0.0:
            raise ValueError(f"radiation pressure bracket invalid at cell {i}")
        for _ in range(80):
            mid = 0.5*(lo + hi)
            f_mid, v2 = residual(mid)
            if abs(f_mid) <= 1.0e-12*max(E[i], 1.0):
                pressure[i] = mid
                v_r[i] = math.copysign(math.sqrt(max(v2, 0.0)/gamma_rr[i]), S_r[i])
                break
            if f_mid >= 0.0:
                hi = mid
            else:
                lo = mid
        else:
            mid = 0.5*(lo + hi)
            f_mid, v2 = residual(mid)
            if not math.isfinite(f_mid) or v2 >= 1.0:
                raise ValueError(f"radiation pressure inversion failed at cell {i}")
            pressure[i] = mid
            v_r[i] = math.copysign(math.sqrt(max(v2, 0.0)/gamma_rr[i]), S_r[i])
    rho_r = 3.0 * pressure
    enthalpy = rho_r + pressure
    W = 1.0 / np.sqrt(1.0 - gamma_rr*v_r*v_r)
    rho += (enthalpy*W*W - pressure) / KAPPA
    pr += (enthalpy*W*W*(gamma_rr*v_r*v_r) + pressure) / KAPPA
    pt += pressure / KAPPA
    current += enthalpy*W*W*gamma_rr*v_r / KAPPA
    return {"rho": rho, "pr": pr, "pt": pt, "j": current}


def peak(values, mask=None):
    x = np.asarray(values, dtype=float)
    ids = np.arange(x.size) if mask is None else np.flatnonzero(mask)
    i = int(ids[np.argmax(np.abs(x[ids]))])
    return {"max_abs": float(abs(x[i])), "signed": float(x[i]), "i": i, "r": float((i+0.5))}


def rms(values, mask):
    x = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(x[mask]*x[mask])))


def load_field(npz, prefix, name):
    key = f"{prefix}__{name}"
    if key not in npz:
        raise ValueError(f"required profile is missing: {key}")
    x = np.asarray(npz[key], dtype=float)
    if x.shape != (80,) or not np.all(np.isfinite(x)):
        raise ValueError(f"profile has wrong shape or non-finite values: {key}")
    return x


def audit(zip_path: Path, *, allow_other_hash: bool = False):
    blob = zip_path.read_bytes()
    zip_hash, manifest, report, capture, profiles = verify_zip(blob, allow_other_hash=allow_other_hash)
    records = []
    common_r = None
    for idx, state_meta in enumerate(capture["accepted_state_history"]):
        prefix = f"accepted_{idx:03d}"
        r = load_field(profiles, prefix, "r")
        a = load_field(profiles, prefix, "geometry__a")
        b = load_field(profiles, prefix, "geometry__b")
        X = load_field(profiles, prefix, "geometry__X")
        Aa = load_field(profiles, prefix, "geometry__Aa")
        K = load_field(profiles, prefix, "geometry__K")
        Lambda = load_field(profiles, prefix, "geometry__Lambda")
        alpha = load_field(profiles, prefix, "geometry__alpha")
        scalar = {n: load_field(profiles, prefix, f"scalars__{n}") for n in SCALARS}
        matter = {
            species: {f: load_field(profiles, prefix, f"matter__{species}__{f}") for f in CONSERVED}
            for species in SPECIES
        }
        if common_r is None:
            common_r = r.copy()
        elif not np.array_equal(r, common_r):
            raise ValueError("accepted-state grids do not match")
        h = float(r[1] - r[0])
        if not np.allclose(np.diff(r), h, rtol=0, atol=1e-14):
            raise ValueError("nonuniform grid encountered")
        if np.any(a <= 0) or np.any(b <= 0) or np.any(X <= 0):
            raise ValueError("non-positive metric profile encountered")

        projection = stress_projection(r, a, b, X, scalar, matter, h)
        geo = geometry_terms(r, a, b, X, alpha, Lambda, h)
        Ab = -0.5 * Aa
        H_R = geo["R"]
        H_A = -(Aa*Aa + 2.0*Ab*Ab)
        H_K = (2.0/3.0)*K*K
        H_rho = -16.0*math.pi*projection["rho"]
        H = H_R + H_A + H_K + H_rho

        chip = -0.5*d4(X, h, 1)/X
        M_A = d4(Aa, h, 1)
        M_K = -(2.0/3.0)*d4(K, h, 1)
        M_chi = 6.0*Aa*chip
        M_metric = (Aa - Ab)*(2.0/r + d4(b, h, 1)/b)
        M_j = -8.0*math.pi*projection["j"]
        M = M_A + M_K + M_chi + M_metric + M_j

        ap, bp = d4(a, h, 1), d4(b, h, 1)
        lambda_geom = ap/(2.0*a*a) - bp/(a*b) + 2.0/r*(1.0/b - 1.0/a)
        C_Lambda = Lambda - lambda_geom
        determinant = a*b*b - 1.0

        ann = r >= 5.0
        outer = r >= 75.0
        local_scale_H = np.abs(H_R) + np.abs(H_A) + np.abs(H_K) + np.abs(H_rho)
        local_scale_M = np.abs(M_A) + np.abs(M_K) + np.abs(M_chi) + np.abs(M_metric) + np.abs(M_j)
        local_scale_C = np.abs(Lambda) + np.abs(lambda_geom)
        scaled_H = np.divide(np.abs(H), local_scale_H, out=np.zeros_like(H), where=local_scale_H > 0)
        scaled_M = np.divide(np.abs(M), local_scale_M, out=np.zeros_like(M), where=local_scale_M > 0)
        scaled_C = np.divide(np.abs(C_Lambda), local_scale_C, out=np.zeros_like(C_Lambda), where=local_scale_C > 0)

        samples = {}
        for radius in (0.5, 1.5, 2.5, 5.5, 6.5, 10.5, 20.5, 40.5, 60.5, 63.5, 64.5, 70.5, 75.5, 79.5):
            i = int(np.argmin(np.abs(r-radius)))
            samples[f"{r[i]:.1f}"] = {
                "H": float(H[i]), "M": float(M[i]), "C_lambda": float(C_Lambda[i]),
                "determinant": float(determinant[i]), "rho_total": float(projection["rho"][i]),
                "j_total": float(projection["j"][i])
            }

        records.append({
            "index": idx, "t": float(state_meta["t"]), "tau": float(state_meta["tau"]),
            "radiation_outer_ratio": float(state_meta["input_cone"]["max_ratio"]),
            "H_full_peak": peak(H), "H_r_ge_5_peak": peak(H, ann), "H_outer_5_cell_peak": peak(H, outer),
            "H_rms_r_ge_5": rms(H, ann), "H_center_cells_0_to_2": [float(v) for v in H[:3]],
            "H_center_i0_terms": {"R": float(H_R[0]), "tracefree_A": float(H_A[0]), "K": float(H_K[0]), "matter": float(H_rho[0])},
            "H_max_abs_over_local_term_scale_r_ge_5": peak(scaled_H, ann),
            "M_full_peak": peak(M), "M_r_ge_5_peak": peak(M, ann), "M_outer_5_cell_peak": peak(M, outer),
            "M_rms_r_ge_5": rms(M, ann), "M_center_cells_0_to_2": [float(v) for v in M[:3]],
            "M_center_i0_terms": {"D_Aa": float(M_A[0]), "D_K": float(M_K[0]), "chi": float(M_chi[0]), "metric_ratio": float(M_metric[0]), "matter_current": float(M_j[0])},
            "M_max_abs_over_local_term_scale_r_ge_5": peak(scaled_M, ann),
            "C_lambda_full_peak": peak(C_Lambda), "C_lambda_r_ge_5_peak": peak(C_Lambda, ann), "C_lambda_outer_5_cell_peak": peak(C_Lambda, outer),
            "C_lambda_rms_r_ge_5": rms(C_Lambda, ann), "C_lambda_center_cells_0_to_2": [float(v) for v in C_Lambda[:3]],
            "C_lambda_center_i0_terms": {"Lambda": float(Lambda[0]), "minus_metric_connection": float(-lambda_geom[0])},
            "C_lambda_max_abs_over_local_term_scale_r_ge_5": peak(scaled_C, ann),
            "determinant_full_peak": peak(determinant), "determinant_rms": float(np.sqrt(np.mean(determinant*determinant))),
            "fixed_radius_samples": samples
        })

    # Cross-check our full stress projection against raw profiles recorded on the
    # exact same physical state by the independent stage observer.
    prefix = "accepted_026"
    term_prefix = "geometry_702"
    get = lambda n: load_field(profiles, term_prefix, f"metric__{n}")
    r = load_field(profiles, prefix, "r")
    a, b, X = (load_field(profiles, prefix, f"geometry__{n}") for n in ("a", "b", "X"))
    scalar = {n: load_field(profiles, prefix, f"scalars__{n}") for n in SCALARS}
    matter = {
        s: {f: load_field(profiles, prefix, f"matter__{s}__{f}") for f in CONSERVED}
        for s in SPECIES
    }
    h = float(r[1] - r[0])
    projection = stress_projection(r, a, b, X, scalar, matter, h)
    stress_diffs = {}
    for field in ("rho", "pr", "pt"):
        saved = load_field(profiles, term_prefix, f"source__l3_TotalStressEnergy_{field}")
        diff = projection[field] - saved
        stress_diffs[field] = {
            "max_abs": float(np.max(np.abs(diff))),
            "max_scaled_to_64eps": float(np.max(
                np.abs(diff) / (64.0*np.finfo(float).eps*(np.abs(projection[field])+np.abs(saved)) + np.finfo(float).tiny)
            ))
        }
        if stress_diffs[field]["max_scaled_to_64eps"] > 1.0:
            raise ValueError(f"stress projection failed roundoff comparison: {field}")

    # Cross-check copied Ricci operators through the captured actual primary-L2
    # trace-free Ricci contribution on the same terminal stage input.
    ga, gb, gX, galpha, gLambda = (get(n) for n in ("a", "b", "X", "alpha", "Lambda"))
    gR = geometry_terms(r, ga, gb, gX, galpha, gLambda, h)
    ricci_pred = galpha*(gR["Rr_r"] - gR["R"]/3.0)
    ricci_saved = load_field(profiles, term_prefix, "source__l2_Aa_ricci_tracefree")
    ricci_residual = ricci_pred - ricci_saved
    ricci_bound = 64.0*np.finfo(float).eps*(np.abs(ricci_pred) + np.abs(ricci_saved)) + np.finfo(float).tiny
    ricci_scaled = float(np.max(np.abs(ricci_residual)/ricci_bound))
    if ricci_scaled > 1.0:
        raise ValueError(f"pinned Ricci formula did not reproduce captured L2 term within 64 epsilon: ratio={ricci_scaled}")

    # Non-lapse geometry and matter are identical between accepted_026 and the
    # first RHS geometry stage after the CMC solve.
    for field in ("a", "b", "X", "Aa", "K", "Lambda"):
        if not np.array_equal(
            load_field(profiles, prefix, f"geometry__{field}"),
            load_field(profiles, term_prefix, f"metric__{field}")
        ):
            raise ValueError(f"terminal physical profile changed unexpectedly: {field}")

    summary = {
        "status": "completed_offline_constraint_audit",
        "diagnostic_only": True,
        "zip_sha256": zip_hash,
        "internal_manifest_files_verified": sorted(manifest),
        "source_baseline": SOURCE_BASELINE,
        "vendor_pin": VENDOR_PIN,
        "observer_head": CAPTURE_HEAD,
        "snapshots": len(records),
        "time_window": {"first_t": records[0]["t"], "last_t": records[-1]["t"], "first_tau": records[0]["tau"], "last_tau": records[-1]["tau"]},
        "terminal_stress_projection_roundoff_crosscheck": stress_diffs,
        "terminal_pinned_Ricci_L2_crosscheck": {
            "max_abs_residual": float(np.max(np.abs(ricci_residual))),
            "max_abs_residual_over_64eps_bound": ricci_scaled
        },
        "records": records,
        "interpretation_limits": [
            "This is only the saved terminal replay window; it does not show when the central constraints first became large earlier in the trajectory.",
            "All constraints are recomputed at one resolution (N=80, r_max=80) using the pinned production definitions/operators; this is not a continuum-convergence test.",
            "The inner max at r=0.5 is a first-cell witness, not a continuum value at r=0. It must be read alongside the r>=5 fixed-physical-region witnesses.",
            "The radiation predictor failure and the central constraint residuals are co-observed; this audit does not prove one causes the other.",
            "No state was advanced or modified. No production source, physics, vendor, boundary, CMC, or campaign setting was changed."
        ]
    }
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_zip", type=Path)
    parser.add_argument("--output-json", type=Path)
    parser.add_argument("--allow-other-zip-hash", action="store_true")
    args = parser.parse_args()
    report = audit(args.artifact_zip, allow_other_hash=args.allow_other_zip_hash)
    rendered = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output_json:
        args.output_json.write_text(rendered, encoding="utf-8")
    print(json.dumps({
        "status": report["status"], "snapshots": report["snapshots"],
        "first_t": report["time_window"]["first_t"], "last_t": report["time_window"]["last_t"],
        "stress_projection_roundoff_crosscheck": report["terminal_stress_projection_roundoff_crosscheck"],
        "ricci_l2_crosscheck": report["terminal_pinned_Ricci_L2_crosscheck"]
    }, indent=2, sort_keys=True))

if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, zipfile.BadZipFile, json.JSONDecodeError) as exc:
        print("AUDIT FAILED CLOSED: " + str(exc), file=sys.stderr)
        raise SystemExit(2)
