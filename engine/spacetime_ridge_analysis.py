"""Propagating-ridge / comoving-lock diagnostic for the 0star campaign.

Diagnostic only: no evolution equations, source terms, or fitted coefficients.
The same radial window used by the existing space-time lock analysis is used.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

RMIN, RMAX = 1.0, 60.0
FIELDS = ("rho_total", "Rdot", "alpha", "theta")


def load_case(case_dir: Path):
    z = np.load(case_dir / "profiles.npz")
    meta = json.loads((case_dir / "meta.json").read_text())
    return z, meta


def ridge_track(z, key: str):
    r = np.asarray(z["r"], dtype=float)
    a = np.asarray(z[f"{key}_grad_abs"], dtype=float)
    mask = (r >= RMIN) & (r <= RMAX)
    idx = np.argmax(np.where(mask[None, :], a, -np.inf), axis=1)
    return r[idx], a[np.arange(a.shape[0]), idx]


def linear_speed(t, r):
    if len(t) < 2 or np.std(t) == 0:
        return float("nan")
    return float(np.polyfit(t, r, 1)[0])


def summarize_case(z, meta):
    t = np.asarray(z["t"], dtype=float)
    tracks = {}
    for key in ("D_active",) + FIELDS:
        tracks[key], _ = ridge_track(z, key)

    d = tracks["D_active"]
    out = {
        "meta": meta,
        "tracks": {
            "D_active": {
                "r_start": float(d[0]),
                "r_end": float(d[-1]),
                "net_delta_r": float(d[-1] - d[0]),
                "linear_speed": linear_speed(t, d),
            }
        },
        "relative": {},
    }

    vd = linear_speed(t, d)
    for key in FIELDS:
        x = tracks[key]
        delta = x - d
        corr = (
            float(np.corrcoef(d, x)[0, 1])
            if np.std(d) > 0 and np.std(x) > 0
            else float("nan")
        )
        out["tracks"][key] = {
            "r_start": float(x[0]),
            "r_end": float(x[-1]),
            "net_delta_r": float(x[-1] - x[0]),
            "linear_speed": linear_speed(t, x),
        }
        out["relative"][key] = {
            "offset_median": float(np.median(delta)),
            "offset_q25": float(np.percentile(delta, 25)),
            "offset_q75": float(np.percentile(delta, 75)),
            "offset_std": float(np.std(delta)),
            "co_ridge_position_corr": corr,
            "kinematic_time_equivalent": (
                float(np.median(delta) / vd) if np.isfinite(vd) and abs(vd) > 0 else float("nan")
            ),
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-root", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    root = Path(args.input_root)
    case_dirs = sorted(
        p for p in root.iterdir()
        if p.is_dir() and (p / "profiles.npz").exists()
    )
    if not case_dirs:
        raise SystemExit("No case artifacts found")

    report = {
        "diagnostic": "0star propagating-ridge / comoving-lock test",
        "conventions": {
            "ridge": "argmax of |dr(field)| in the existing r=1..60 analysis window",
            "co_motion": "geometry ridge position compared directly with D ridge position",
            "kinematic_time_equivalent": "spatial offset divided by fitted D ridge speed; not a causal delay",
            "evolution": "unchanged; diagnostic only",
        },
        "cases": {},
    }

    for case_dir in case_dirs:
        z, meta = load_case(case_dir)
        report["cases"][case_dir.name] = summarize_case(z, meta)

    # Matched-grid radiation ON/OFF check.
    on_dir = root / "D1e4-on-N160"
    off_dir = root / "D1e4-off-N160"
    if on_dir.exists() and off_dir.exists():
        zon, _ = load_case(on_dir)
        zof, _ = load_case(off_dir)
        if np.array_equal(zon["t"], zof["t"]) and np.array_equal(zon["r"], zof["r"]):
            cmp = {}
            for key in ("D_active",) + FIELDS:
                ron, _ = ridge_track(zon, key)
                roff, _ = ridge_track(zof, key)
                diff = ron - roff
                cmp[key] = {
                    "max_abs_ridge_difference": float(np.max(np.abs(diff))),
                    "median_ridge_difference": float(np.median(diff)),
                    "identical_all_samples": bool(np.array_equal(ron, roff)),
                }
            report["matched_on_off_N160"] = cmp

    # Same-physics resolution convergence.
    a_dir = root / "D1e4-on-N160"
    b_dir = root / "D1e4-on-N320"
    if a_dir.exists() and b_dir.exists():
        za, _ = load_case(a_dir)
        zb, _ = load_case(b_dir)
        cmp = {}
        ta = np.asarray(za["t"], dtype=float)
        tb = np.asarray(zb["t"], dtype=float)
        for key in ("D_active",) + FIELDS:
            ra, _ = ridge_track(za, key)
            rb, _ = ridge_track(zb, key)
            rb_i = np.interp(ta, tb, rb)
            diff = ra - rb_i
            cmp[key] = {
                "max_abs_ridge_difference": float(np.max(np.abs(diff))),
                "median_ridge_difference": float(np.median(diff)),
                "mean_abs_ridge_difference": float(np.mean(np.abs(diff))),
            }
        report["resolution_convergence"] = cmp

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, default=float))
    print(json.dumps(report, indent=2, default=float))


if __name__ == "__main__":
    main()
