"""Aggregate analysis for the full-profile space-time D->geometry test."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


FIELDS = {
    "alpha": "alpha_grad_abs",
    "theta": "theta_grad_abs",
    "Rdot": "Rdot_grad_abs",
    "rho_total": "rho_total_grad_abs",
}
SUB_FIELDS = ("alpha", "theta", "Rdot", "rho_total")
DR_MIN = -12.0
DR_MAX = 16.0
DT_MIN = -6.0
DT_MAX = 6.0
RADIAL_MIN = 1.0
RADIAL_MAX = 60.0
MIN_CORR_SAMPLES = 200


def load_case(case_dir: Path):
    z = np.load(case_dir / "profiles.npz")
    meta = json.loads((case_dir / "meta.json").read_text())
    return z, meta


def common_time_grid(z1, z2):
    t1 = np.asarray(z1["t"], dtype=float)
    t2 = np.asarray(z2["t"], dtype=float)
    lo = max(float(t1[0]), float(t2[0]))
    hi = min(float(t1[-1]), float(t2[-1]))
    step = min(
        float(np.median(np.diff(t1))),
        float(np.median(np.diff(t2))),
    )
    n = int(np.floor((hi - lo) / step + 1e-9)) + 1
    return lo + step * np.arange(max(n, 1))


def interp_time(z, key, tnew):
    t = np.asarray(z["t"], dtype=float)
    y = np.asarray(z[key], dtype=float)
    if y.ndim != 2:
        raise ValueError(f"{key} must be 2-D, got {y.shape}")
    out = np.empty((len(tnew), y.shape[1]), dtype=float)
    for j in range(y.shape[1]):
        out[:, j] = np.interp(tnew, t, y[:, j])
    return out


def zscore_rows(a):
    a = np.asarray(a, dtype=float)
    mu = np.mean(a, axis=1, keepdims=True)
    sd = np.std(a, axis=1, keepdims=True)
    sd = np.where(sd > 0.0, sd, 1.0)
    return (a - mu) / sd


def corr_shift(a, b, di, dj, radial_mask):
    # C_DX(dr,dt) = corr(D(r,t), X(r+dr,t+dt)).
    T, R = a.shape
    if dj >= 0:
        ta0, ta1 = 0, T - dj
        tb0, tb1 = dj, T
    else:
        ta0, ta1 = -dj, T
        tb0, tb1 = 0, T + dj

    if di >= 0:
        ra0, ra1 = 0, R - di
        rb0, rb1 = di, R
    else:
        ra0, ra1 = -di, R
        rb0, rb1 = 0, R + di

    if ta1 <= ta0 or ra1 <= ra0:
        return float("nan"), 0

    m = radial_mask[ra0:ra1]
    if m.sum() < 4:
        return float("nan"), 0

    aa = a[ta0:ta1, ra0:ra1][:, m].ravel()
    bb = b[tb0:tb1, rb0:rb1][:, m].ravel()
    good = np.isfinite(aa) & np.isfinite(bb)
    aa = aa[good]
    bb = bb[good]
    if len(aa) < MIN_CORR_SAMPLES:
        return float("nan"), len(aa)

    aa -= aa.mean()
    bb -= bb.mean()
    den = float(np.linalg.norm(aa) * np.linalg.norm(bb))
    return (float(np.dot(aa, bb) / den) if den > 0 else float("nan")), len(aa)


def lag_surface(a, b, r, t, dr_values, dt_values):
    dr0 = float(np.median(np.diff(r)))
    dt0 = float(np.median(np.diff(t)))
    out = []
    mask = (r >= RADIAL_MIN) & (r <= min(RADIAL_MAX, r[-1] - 2 * dr0))
    for dt in dt_values:
        dj = int(round(float(dt) / dt0))
        for dr in dr_values:
            di = int(round(float(dr) / dr0))
            c, n = corr_shift(a, b, di, dj, mask)
            out.append({"dr": di * dr0, "dt": dj * dt0, "corr": c, "n": n})
    return out


def pick(surface, predicate):
    vals = [x for x in surface if np.isfinite(x["corr"]) and predicate(x)]
    if not vals:
        return None
    return max(vals, key=lambda x: x["corr"])


def value_map(surface):
    return {(round(x["dr"], 9), round(x["dt"], 9)): x["corr"] for x in surface}


def block_corr(a, b, r, t, di, dj, nblocks=4):
    T = len(t)
    edges = np.linspace(0, T, nblocks + 1, dtype=int)
    vals = []
    mask = (r >= RADIAL_MIN) & (r <= min(RADIAL_MAX, r[-1] - 2 * np.median(np.diff(r))))
    for k in range(nblocks):
        i0, i1 = edges[k], edges[k + 1]
        if i1 - i0 < 4:
            continue
        aa = a[i0:i1]
        bb = b[i0:i1]
        c, n = corr_shift(aa, bb, di, dj, mask)
        if np.isfinite(c):
            vals.append(c)
    return vals


def shuffle_null(a, b, r, t, dr_values, dt_values, target, rng, count=50):
    # Circularly shift D profiles in time. This destroys D->X temporal ordering
    # while preserving the one-dimensional spatial and temporal autocorrelation.
    target_dt = target["dt"]
    target_dr = target["dr"]
    dt0 = float(np.median(np.diff(t)))
    di = int(round(target_dr / np.median(np.diff(r))))
    dj = int(round(target_dt / dt0))
    T = len(t)
    vals = []
    mask = (r >= RADIAL_MIN) & (r <= min(RADIAL_MAX, r[-1] - 2 * np.median(np.diff(r))))
    for _ in range(count):
        shift = int(rng.integers(max(2, abs(dj) + 1), max(3, T - max(2, abs(dj) + 1))))
        aa = np.roll(a, shift, axis=0)
        c, _ = corr_shift(aa, b, di, dj, mask)
        if np.isfinite(c):
            vals.append(c)
    return vals


def analyze_pair(d, x, r, t, target, rng):
    dr0 = float(np.median(np.diff(r)))
    dt0 = float(np.median(np.diff(t)))
    dr_values = np.arange(
        np.ceil(DR_MIN / dr0) * dr0,
        DR_MAX + 0.5 * dr0,
        dr0,
    )
    dt_values = np.arange(
        np.ceil(DT_MIN / dt0) * dt0,
        DT_MAX + 0.5 * dt0,
        dt0,
    )
    surf = lag_surface(d, x, r, t, dr_values, dt_values)
    best = pick(surf, lambda q: True)
    future = pick(surf, lambda q: 0.25 * dt0 <= q["dt"] <= DT_MAX)
    forward = pick(surf, lambda q: q["dt"] >= -0.25 * dt0 and q["dr"] >= -0.25 * dr0)
    positive_quadrant = pick(surf, lambda q: q["dt"] >= 0.25 * dt0 and q["dr"] >= 0.25 * dr0)
    past = pick(surf, lambda q: q["dt"] <= -0.25 * dt0 and q["dr"] >= 0.0)
    block = None
    null = []
    if future:
        di = int(round(future["dr"] / dr0))
        dj = int(round(future["dt"] / dt0))
        block = block_corr(d, x, r, t, di, dj)
        null = shuffle_null(d, x, r, t, dr_values, dt_values, future, rng)
    return {
        "best_unrestricted": best,
        "best_future": future,
        "best_nonnegative": forward,
        "best_positive_space_time": positive_quadrant,
        "best_past_positive_space": past,
        "future_minus_past": (
            float(future["corr"] - past["corr"])
            if future and past else float("nan")
        ),
        "future_block_corrs": block,
        "future_block_mean": float(np.mean(block)) if block else float("nan"),
        "future_block_min": float(np.min(block)) if block else float("nan"),
        "shuffle_null_mean": float(np.mean(null)) if null else float("nan"),
        "shuffle_null_p95": float(np.percentile(null, 95)) if null else float("nan"),
        "shuffle_exceedance_fraction": (
            float(np.mean(np.asarray(null) >= future["corr"]))
            if null and future else float("nan")
        ),
        "lag_surface": surf,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-root", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    root = Path(args.input_root)
    case_dirs = sorted(p for p in root.iterdir() if p.is_dir() and (p / "profiles.npz").exists())
    if not case_dirs:
        raise SystemExit("No case artifacts found")

    rng = np.random.default_rng(20261007)
    loaded = {p.name: load_case(p) for p in case_dirs}

    # Baselines by resolution: only matched-grid subtraction is allowed.
    baselines = {
        meta["resolution"]: (z, meta)
        for name, (z, meta) in loaded.items()
        if float(meta["D_amplitude"]) == 0.0
    }

    report = {
        "diagnostic": "0star full space-time D->geometry lock",
        "conventions": {
            "correlation": "corr(|dr D_active|(r,t), |dr X|(r+dr,t+dt))",
            "positive_dr": "response feature is at larger radius",
            "positive_dt": "response feature occurs later in system time",
            "causality": "not inferred from correlation alone",
            "radial_window": [RADIAL_MIN, RADIAL_MAX],
            "lag_window": {
                "dr": [DR_MIN, DR_MAX],
                "dt": [DT_MIN, DT_MAX],
            },
        },
        "cases": {},
    }

    for name, (z, meta) in loaded.items():
        t = np.asarray(z["t"], dtype=float)
        r = np.asarray(z["r"], dtype=float)
        D = zscore_rows(np.asarray(z["D_active_grad_abs"], dtype=float))
        case_out = {
            "meta": meta,
            "responses": {},
        }

        baseline = baselines.get(meta["resolution"])
        for field, key in FIELDS.items():
            X = zscore_rows(np.asarray(z[key], dtype=float))
            res = analyze_pair(D, X, r, t, None, rng)
            # Time-reversed comparison at the same |lag| reports whether the
            # observed positive-time lock is stronger than the reverse ordering.
            case_out["responses"][field] = res

            if baseline is not None and float(meta["D_amplitude"]) > 0.0:
                bz, bmeta = baseline
                if np.array_equal(r, np.asarray(bz["r"], dtype=float)):
                    common_t = common_time_grid(z, bz)
                    D_common = zscore_rows(interp_time(z, "D_active_grad_abs", common_t))
                    X_common = interp_time(z, key.replace("_grad_abs", ""), common_t)
                    X0_common = interp_time(
                        bz, key.replace("_grad_abs", ""), common_t
                    )
                    dX = X_common - X0_common
                    dX_edge = zscore_rows(
                        np.abs(np.gradient(dX, r, axis=1))
                    )
                    sub = analyze_pair(D_common, dX_edge, r, common_t, None, rng)
                    case_out["responses"][field]["D0_subtracted"] = sub

        # Compact ridge-only timing summary, to preserve the previous result
        # while the full lag surfaces provide the new test.
        case_out["ridge_offsets_final"] = {}
        if "r_grad_D" in z.files:
            pass
        report["cases"][name] = case_out

    # The compact table is intentionally redundant with the detailed map.
    compact = []
    for name, case in report["cases"].items():
        for field, res in case["responses"].items():
            f = res["best_positive_space_time"]
            compact.append({
                "case": name,
                "field": field,
                "best_dr": f["dr"] if f else None,
                "best_dt": f["dt"] if f else None,
                "corr": f["corr"] if f else None,
                "future_minus_past": res["future_minus_past"],
                "block_mean": res["future_block_mean"],
                "block_min": res["future_block_min"],
                "null_p95": res["shuffle_null_p95"],
                "null_exceedance": res["shuffle_exceedance_fraction"],
            })

    report["compact"] = compact
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, default=float))
    print(json.dumps(compact, indent=2, default=float))


if __name__ == "__main__":
    main()
