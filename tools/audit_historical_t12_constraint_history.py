#!/usr/bin/env python3
"""Recompute historical V5.5 H/M/connection/determinant constraints from an exact archived campaign ZIP.

This read-only verifier does not import the engine or advance any state. It
recomputes 18 full-field checkpoints, then checks the t=12 witnesses against
the original campaign summary and extracts threshold-crossing times from the
historical per-step ledgers.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile
import numpy as np

EXPECTED_ZIP_SHA256 = "dbbaa3be94242cdeecd437d0a41c00d4cf3f2d625e1f3aaaad9e1b45833f2047"
SOURCE_COMMIT = "6d8b3f2790644591326df57c55908863f994cb0c"
VENDOR_PIN = "d6052d605673ce9d82cdc99b0df79855fa2ea215"
RUN_ID = 38208253687
ARTIFACT_ID = 11728068723
KAPPA = 8.0 * math.pi
RESOLUTIONS = (40, 60, 80)
CHECKPOINT_TIMES = (2, 4, 6, 8, 10, 12)
ZONES = ((0.0, 2.0), (2.0, 5.0), (5.0, 10.0), (10.0, 20.0),
         (20.0, 40.0), (40.0, 60.0), (60.0, 75.0), (75.0, 80.000001)


def d4(values: np.ndarray, h: float, parity: int, order: int = 1) -> np.ndarray:
    """Pinned vendor cell_derivative_fourth / cell_second_derivative_fourth."""
    u = np.asarray(values, dtype=float)
    n = u.size
    if u.ndim != 1 or n < 5 or parity not in (-1, 1) or order not in (1, 2):
        raise ValueError("invalid fourth-order derivative input")
    ext = np.empty(n + 4)
    ext[2:-2] = u
    ext[1] = parity * u[0]
    ext[0] = parity * u[1]
    out = np.empty_like(u)
    for i in range(n - 2):
        s = ext[i:i+5]
        if order == 1:
            out[i] = (s[0] - 8*s[1] + 8*s[3] - s[4]) / (12*h)
        else:
            out[i] = (-s[0] + 16*s[1] - 30*s[2] + 16*s[3] - s[4]) / (12*h*h)
    for i in (n - 2, n - 1):
        idx = np.arange(n - 5, n)
        offsets = idx.astype(float) - float(i)
        matrix = np.vstack([offsets**power for power in range(5)])
        target = np.zeros(5)
        target[order] = float(math.factorial(order))
        out[i] = float(np.linalg.solve(matrix, target) @ u[idx]) / h**order
    return out


def geometry_terms(r, a, b, X, Lambda, h):
    ap, bp = d4(a, h, 1), d4(b, h, 1)
    app, bpp = d4(a, h, 1, 2), d4(b, h, 1, 2)
    Xp, Xpp = d4(X, h, 1), d4(X, h, 1, 2)
    chip = -0.5 * Xp / X
    chipp = -0.5 * (Xpp / X - (Xp / X)**2)
    Lambda_p = d4(Lambda, h, -1)
    metric_lambda = (1.0 - a/b) / r**2
    inv = X**2 / a
    radial = (0.5*app/a - a*Lambda_p - 0.75*(ap/a)**2 + 0.5*(bp/b)**2
              - 0.5*Lambda*ap + ap/(r*b) + 2*metric_lambda*(1+r*bp/b)
              + 4*chipp - 2*chip*(ap/a - bp/b - 2/r))
    scalar = (0.5*app/a + bpp/b - a*Lambda_p - (ap/a)**2 + 0.5*(bp/b)**2
              + 2*(3-a/b)*bp/(r*b) + 4*metric_lambda + 8*(chipp+chip**2)
              - 8*chip*(0.5*ap/a - bp/b - 2/r))
    return {"R": -inv*scalar, "Rr_r": -inv*radial}


def total_projection(r, a, b, X, scalar, matter, h):
    """Reconstruct archived scalar + DM + baryon + radiation projections."""
    n = len(r)
    inv = X**2 / a
    rho, pr, pt, current = (np.zeros(n) for _ in range(4))
    scalar_terms = (
        (scalar["S"], scalar["PS"], 0.5*scalar["S"]**2, 1.0),
        (scalar["D"], scalar["PD"], -0.5*scalar["D"]**2, 1.0),
        (scalar["phi"], scalar["Pi"],
         8.242522415500654e-5*(np.exp(-scalar["phi"])-0.10), KAPPA),
    )
    for field, momentum, potential, scale in scalar_terms:
        fp = d4(field, h, 1)
        grad_sq = inv * fp*fp
        rho += (0.5*momentum**2 + 0.5*grad_sq + potential) / scale
        pr += (0.5*momentum**2 + 0.5*grad_sq - potential) / scale
        pt += (0.5*momentum**2 - 0.5*grad_sq - potential) / scale
        current += -momentum*fp / scale

    sqrt_gamma = np.sqrt(a)*b*r**2/X**3
    gamma_rr, gamma_rr_inv = a/X**2, X**2/a
    for name in ("dark_matter", "baryons"):
        s = matter[name]
        D, Sr = s["rest"]/sqrt_gamma, s["momentum_r"]/sqrt_gamma
        s2 = gamma_rr_inv * Sr*Sr
        v2 = s2/(D*D+s2)
        vr = np.sign(Sr)*np.sqrt(np.maximum(v2/gamma_rr, 0.0))
        W = 1.0/np.sqrt(1.0-v2)
        rest_density = D/W
        rho += rest_density*W*W/KAPPA
        pr += rest_density*W*W*v2/KAPPA
        current += rest_density*W*W*gamma_rr*vr/KAPPA

    s = matter["radiation"]
    E, Sr = s["energy_t"]/sqrt_gamma, s["momentum_r"]/sqrt_gamma
    S2 = gamma_rr_inv*Sr*Sr
    pressure, vr = np.empty(n), np.empty(n)
    for i in range(n):
        if S2[i] == 0.0:
            pressure[i], vr[i] = E[i]/3.0, 0.0
            continue
        lo, hi = 0.0, E[i]/3.0
        def residual(p):
            v2 = S2[i]/(E[i]+p)**2
            return (3*p-E[i]) + 4*p*v2/(1-v2), v2
        f_lo, _ = residual(lo)
        f_hi, _ = residual(hi)
        if f_hi < 0.0:
            hi = np.nextafter(hi, math.inf)
            f_hi, _ = residual(hi)
        if not math.isfinite(f_hi) or f_lo > 0.0 or f_hi < 0.0:
            raise ValueError(f"radiation recovery bracket failed at cell {i}")
        for _ in range(80):
            mid = 0.5*(lo+hi)
            f_mid, v2 = residual(mid)
            if abs(f_mid) <= 1.0e-12*max(E[i], 1.0):
                pressure[i] = mid
                vr[i] = math.copysign(math.sqrt(max(v2, 0.0)/gamma_rr[i]), Sr[i])
                break
            if f_mid >= 0.0:
                hi = mid
            else:
                lo = mid
        else:
            raise ValueError(f"radiation recovery failed to converge at cell {i}")
    W = 1.0/np.sqrt(1.0-gamma_rr*vr*vr)
    enthalpy = 4.0*pressure
    rho += (enthalpy*W*W-pressure)/KAPPA
    pr += (enthalpy*W*W*(gamma_rr*vr*vr)+pressure)/KAPPA
    pt += pressure/KAPPA
    current += enthalpy*W*W*gamma_rr*vr/KAPPA
    return {"rho": rho, "pr": pr, "pt": pt, "j": current}


def peak(values, r, mask=None):
    values = np.asarray(values, dtype=float)
    ids = np.arange(values.size) if mask is None else np.flatnonzero(mask)
    if ids.size == 0:
        return {"max_abs": None, "signed": None, "index": None, "r": None}
    idx = int(ids[np.argmax(np.abs(values[ids]))])
    return {"max_abs": float(abs(values[idx])), "signed": float(values[idx]),
            "index": idx, "r": float(r[idx])}


def read_state(data):
    r = np.asarray(data["r"], dtype=float)
    h = float(r[1]-r[0])
    if not np.allclose(np.diff(r), h, rtol=0.0, atol=1e-12):
        raise ValueError("checkpoint grid is not uniform")
    a, b, X = (np.asarray(data[k], dtype=float) for k in ("a", "b", "X"))
    Aa, K, Lambda = (np.asarray(data[k], dtype=float) for k in ("Aa", "K", "Lambda"))
    if any(not np.all(np.isfinite(x)) for x in (r, a, b, X, Aa, K, Lambda)):
        raise ValueError("non-finite checkpoint geometry")
    if np.any(a <= 0) or np.any(b <= 0) or np.any(X <= 0):
        raise ValueError("non-positive checkpoint metric")
    scalar = {k: np.asarray(data[k], dtype=float) for k in ("S","PS","D","PD","phi","Pi")}
    matter = {
        name: {field: np.asarray(data[f"{prefix}_{field}"], dtype=float)
               for field in ("rest","energy_t","momentum_r")}
        for name, prefix in (("dark_matter","dm"),("baryons","baryon"),("radiation","radiation"))
    }
    proj = total_projection(r,a,b,X,scalar,matter,h)
    geo = geometry_terms(r,a,b,X,Lambda,h)
    Ab = -0.5*Aa
    H = geo["R"] - (Aa*Aa + 2*Ab*Ab) + (2.0/3.0)*K*K - 16*math.pi*proj["rho"]
    chi_r = -0.5*d4(X,h,1)/X
    M = (d4(Aa,h,1) - (2.0/3.0)*d4(K,h,1) + 6*Aa*chi_r
         + (Aa-Ab)*(2.0/r+d4(b,h,1)/b) - 8*math.pi*proj["j"])
    metric_connection = (d4(a,h,1)/(2*a*a) - d4(b,h,1)/(a*b)
                         + 2.0/r*(1.0/b-1.0/a))
    C = Lambda - metric_connection
    det = a*b*b-1.0
    outer_cut = float(r[-1] - 4.0*h)
    masks = {"full": np.ones_like(r,dtype=bool), "r_ge_5": r >= 5.0,
             "outer_five_cells": r >= outer_cut}
    output = {"N": int(len(r)), "t": float(data["t"]), "tau": float(data["tau"])}
    for cname, values in (("H",H),("M",M),("C_lambda",C),("determinant",det)):
        for mask_name, mask in masks.items():
            output[f"{cname}_{mask_name}_peak"] = peak(values,r,mask)
    output["radiation_outer_ratio"] = float(
        np.sqrt(max((X[-1]**2/a[-1]) *
                    (data["radiation_momentum_r"][-1] /
                     (np.sqrt(a[-1])*b[-1]*r[-1]**2/X[-1]**3))**2,0.0)) /
        (data["radiation_energy_t"][-1] /
         (np.sqrt(a[-1])*b[-1]*r[-1]**2/X[-1]**3))
    )
    # Fixed physical annuli use the actual cell centers, not resolution-dependent indices.
    faces = np.arange(len(r)+1,dtype=float)*h
    volumes = (faces[1:]**3-faces[:-1]**3)/3.0
    arrays = {"H":H,"M":M,"C_lambda":C,"determinant":det}
    output["fixed_physical_annuli"] = {}
    for lo, hi in ZONES:
        mask = (r >= lo) & (r < hi)
        if not np.any(mask):
            continue
        tag = f"{lo:g}-{hi:g}"
        output["fixed_physical_annuli"][tag] = {}
        for cname, values in arrays.items():
            p = peak(values,r,mask)
            p["volume_weighted_rms"] = float(
                np.sqrt(np.sum(volumes[mask]*values[mask]**2)/np.sum(volumes[mask]))
            )
            output["fixed_physical_annuli"][tag][cname] = p
    return output


def summarize_history(history):
    fields = ("hamiltonian_full_absmax","hamiltonian_inner_absmax","hamiltonian_max",
              "momentum_full_absmax","momentum_inner_absmax","momentum_max",
              "connection_full_absmax","connection_inner_absmax","connection_max",
              "determinant_constraint_max","radiation_outer_ratio")
    thresholds = (("hamiltonian_full_absmax",1e-2),("momentum_full_absmax",1e-3),
                  ("connection_full_absmax",1e-3),("hamiltonian_max",1e-4),
                  ("momentum_max",1e-4),("connection_max",1e-4))
    first = {}
    for key, threshold in thresholds:
        hit = next((x for x in history if float(x[key]) > threshold), None)
        first[f"{key}>{threshold:g}"] = None if hit is None else {
            "t": float(hit["t"]), "value": float(hit[key])
        }
    samples=[]
    for target in (0,.3,.6,1,2,4,6,8,10,12):
        row=min(history,key=lambda x:abs(float(x["t"])-target))
        samples.append({"target_t":target,"actual_t":float(row["t"]),
                        **{k:float(row[k]) for k in fields}})
    return {"records":len(history),"steps":len(history)-1,
            "first_step_t":float(history[1]["t"]),
            "initial":{k:float(history[0][k]) for k in fields},
            "first_threshold_crossings":first,"samples":samples}


def audit(artifact_zip: Path):
    blob=artifact_zip.read_bytes()
    digest=hashlib.sha256(blob).hexdigest()
    if digest != EXPECTED_ZIP_SHA256:
        raise ValueError(f"artifact ZIP SHA mismatch: {digest}")
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        names=set(z.namelist())
        expected={f"N{n}/checkpoint_t{t:05d}.000000.npz"
                  for n in RESOLUTIONS for t in CHECKPOINT_TIMES}
        expected |= {f"N{n}/ledgers.json" for n in RESOLUTIONS}
        expected.add("campaign_summary.json")
        missing=sorted(expected-names)
        if missing:
            raise ValueError(f"archive missing required members: {missing}")
        summary=json.loads(z.read("campaign_summary.json"))
        records=[]; checks=[]; histories=[]
        for n in RESOLUTIONS:
            for t in CHECKPOINT_TIMES:
                path=f"N{n}/checkpoint_t{t:05d}.000000.npz"
                data=np.load(io.BytesIO(z.read(path)),allow_pickle=False)
                q=read_state(data)
                records.append(q)
            ledger=json.loads(z.read(f"N{n}/ledgers.json"))
            history=ledger["history"]
            if len(history) != ledger["result"]["steps"]+1:
                raise ValueError(f"history/step count mismatch at N={n}")
            histories.append({"N":n,**summarize_history(history)})
            result=next(x for x in summary["results"] if int(x["resolution"])==n)
            diag=result["last_diagnostics"]
            q=next(x for x in records if x["N"]==n and abs(x["t"]-12.0)<1e-6)
            pairs=(
                ("hamiltonian_max","H_r_ge_5_peak"),
                ("hamiltonian_full_absmax","H_full_peak"),
                ("momentum_max","M_r_ge_5_peak"),
                ("momentum_full_absmax","M_full_peak"),
                ("connection_max","C_lambda_r_ge_5_peak"),
                ("connection_full_absmax","C_lambda_full_peak"),
                ("determinant_constraint_max","determinant_full_peak")
            )
            errors={}
            for source,target in pairs:
                reported=float(diag[source])
                recomputed=float(q[target]["max_abs"])
                errors[f"{source}_absolute_difference"]=abs(reported-recomputed)
            checks.append({"N":n,"errors":errors,
                           "max_absolute_difference":max(errors.values())})
    overall=max(x["max_absolute_difference"] for x in checks)
    report={
        "status":"completed_offline_historical_constraint_audit",
        "diagnostic_only":True,
        "artifact_zip_sha256":digest,
        "workflow_run_id":RUN_ID,"artifact_id":ARTIFACT_ID,
        "source_commit":SOURCE_COMMIT,"vendor_pin":VENDOR_PIN,
        "campaign":summary["campaign"],
        "checkpoint_count":len(records),
        "t12_validation":{"per_resolution":checks,"max_absolute_difference":overall,
                          "passed":overall < 1e-12},
        "per_step_ledger_constraint_onset":histories,
        "full_field_checkpoint_constraint_profiles":records,
        "interpretation_limits":[
            "The archived source is the historical workflow head, not current main; this is evidence about the archived run, not a current-main campaign.",
            "The t=12 full-field recomputation reproduces the source campaign diagnostic witnesses to floating-point precision, validating the constraint/source/operator port for the archived checkpoints.",
            "The campaign has one parameter setup and three resolutions; no convergence order is inferred from sign-changing center-cell maxima.",
            "First threshold-crossing times are the first recorded discrete history samples above illustrative levels; these are not physical event or admission criteria.",
            "The large center Hamiltonian residual exists at the initial slice; momentum/connection depart from roundoff during the first step. Correlation does not identify a unique faulty production equation.",
            "No state was advanced and no production source, physics, vendor, boundary, or campaign parameter was changed."
        ]
    }
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_zip",type=Path)
    parser.add_argument("--output-json",type=Path)
    args=parser.parse_args()
    report=audit(args.artifact_zip)
    rendered=json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n"
    if args.output_json:
        args.output_json.write_text(rendered,encoding="utf-8")
    print(json.dumps({"status":report["status"],"zip_sha256":report["artifact_zip_sha256"],
                      "checkpoint_count":report["checkpoint_count"],
                      "t12_validation":report["t12_validation"],
                      "ledger_counts":{str(x["N"]):x["records"] for x in report["per_step_ledger_constraint_onset"]}},
                     indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
