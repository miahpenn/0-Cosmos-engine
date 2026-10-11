#!/usr/bin/env python3
"""Signed Hamiltonian-residual pole tracking on a saved trajectory; no evolution."""
import argparse, json, math
from pathlib import Path

N, RMAX, NSLICES = 160, 80.0, 49
INNER, RIDGE, TRACK_LIMIT = (1.25, 4.0), (4.0, 20.0), 1.25

def sign(x):
    return 1 if x > 0 else (-1 if x < 0 else 0)

def peak(r, h, lo, hi):
    i = max((i for i, x in enumerate(r) if lo <= x <= hi), key=lambda i: abs(h[i]))
    return {"cell": i, "r": r[i], "H": h[i], "abs_H": abs(h[i]), "sign": sign(h[i])}

def extrema(r, h):
    out = []
    for i in range(1, len(h)-1):
        if not 1.25 <= r[i] <= 20.0:
            continue
        if h[i] >= h[i-1] and h[i] > h[i+1]:
            kind = "max"
        elif h[i] <= h[i-1] and h[i] < h[i+1]:
            kind = "min"
        else:
            continue
        out.append({"kind": kind, "cell": i, "r": r[i], "H": h[i], "sign": sign(h[i])})
    return out

def sign_edges(r, h):
    return [{"r_left": r[i-1], "r_right": r[i], "H_left": h[i-1], "H_right": h[i]}
            for i in range(1, len(h))
            if r[i-1] >= 1.25 and r[i] <= 79.75 and sign(h[i-1])*sign(h[i]) < 0]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", type=Path)
    ap.add_argument("--output-dir", type=Path, default=Path("runs/strong-D-ridge-signed-pole-tracking"))
    args = ap.parse_args()
    rows = sorted((json.loads(s) for s in args.input.read_text().splitlines() if s.strip()), key=lambda x: float(x["t"]))
    if len(rows) != NSLICES:
        raise ValueError(f"expected {NSLICES} saved slices, got {len(rows)}")
    samples = []
    for j, row in enumerate(rows):
        r, h = list(map(float, row["radius"])), list(map(float, row["H_vendor"]))
        domain = row.get("domain", {})
        if len(r) != N or len(h) != N or int(domain.get("N", -1)) != N or not math.isclose(float(domain.get("r_max", -1)), RMAX):
            raise ValueError(f"unexpected profile/domain in slice {j}: {domain}")
        if not all(math.isfinite(x) for x in h):
            raise ValueError(f"nonfinite H in slice {j}")
        samples.append({"t": float(row["t"]), "r": r, "h": h, "extrema": extrema(r, h),
                        "inner": peak(r, h, *INNER), "ridge": peak(r, h, *RIDGE),
                        "offcentre": peak(r, h, 1.25, 79.75), "sign_edges": sign_edges(r, h)})

    def at(t):
        return min(samples, key=lambda x: abs(x["t"]-t))
    s6, s625, s85, s875 = [at(t) for t in (6.0, 6.25, 8.5, 8.75)]
    max6 = {x["r"]: x for x in s6["extrema"] if x["kind"] == "max"}
    max625 = {x["r"]: x for x in s625["extrema"] if x["kind"] == "max"}
    min6 = [x for x in s6["extrema"] if x["kind"] == "min" and 3.0 <= x["r"] <= 4.5]
    min625 = [x for x in s625["extrema"] if x["kind"] == "min" and 3.0 <= x["r"] <= 4.5]
    same_min = min(min6, key=lambda x: abs(x["r"]-3.75)) if min6 else None
    next_min = min(min625, key=lambda x: abs(x["r"]-same_min["r"])) if min625 and same_min else None
    checks = {
        "t6 ridge winner is r=8.25": abs(s6["ridge"]["r"]-8.25) < 1e-9,
        "separate r=5.75 local maximum coexists at t6": 5.75 in max6,
        "t6.25 ridge winner is r=5.75": abs(s625["ridge"]["r"]-5.75) < 1e-9,
        "r=8.25 local maximum is absent at t6.25": 8.25 not in max625,
        "same-radius local minimum changes negative to positive": bool(same_min and next_min and same_min["H"] < 0 < next_min["H"]),
        "radial sign-transition count changes 2 to 0": len(s6["sign_edges"]) == 2 and len(s625["sign_edges"]) == 0,
        "inner/ridge winners exchange dominance from t8.5 to t8.75":
            s85["inner"]["abs_H"] > s85["ridge"]["abs_H"] and s875["inner"]["abs_H"] < s875["ridge"]["abs_H"],
        "all inner/ridge winners at this crossover are positive":
            all(x[k]["H"] > 0 for x in (s85, s875) for k in ("inner", "ridge")),
    }
    # Match a local-minimum branch to the nearest minimum from each adjacent slice.
    min_track = []
    current = same_min
    if current:
        for s in reversed(samples[:samples.index(s6)]):
            candidates = [x for x in s["extrema"] if x["kind"] == "min" and abs(x["r"]-current["r"]) <= TRACK_LIMIT]
            if not candidates: break
            current = min(candidates, key=lambda x: abs(x["r"]-current["r"]))
            min_track.append({"t": s["t"], **current})
        min_track.reverse()
        min_track.append({"t": s6["t"], **same_min})
        current = same_min
        for s in samples[samples.index(s6)+1:]:
            candidates = [x for x in s["extrema"] if x["kind"] == "min" and abs(x["r"]-current["r"]) <= TRACK_LIMIT]
            if not candidates: break
            current = min(candidates, key=lambda x: abs(x["r"]-current["r"]))
            min_track.append({"t": s["t"], **current})
    crossings = []
    for a, b in zip(min_track, min_track[1:]):
        if a["H"]*b["H"] < 0:
            f = abs(a["H"])/(abs(a["H"])+abs(b["H"]))
            crossings.append({"t0": a["t"], "t1": b["t"], "r0": a["r"], "r1": b["r"],
                              "H0": a["H"], "H1": b["H"], "zero_t_linear": a["t"]+f*(b["t"]-a["t"])})
    crossovers = []
    for a, b in zip(samples, samples[1:]):
        d0 = a["ridge"]["abs_H"] - a["inner"]["abs_H"]
        d1 = b["ridge"]["abs_H"] - b["inner"]["abs_H"]
        if d0 < 0 <= d1:
            f = -d0/(d1-d0)
            crossovers.append({"t0": a["t"], "t1": b["t"], "equal_abs_H_t_linear": a["t"]+f*(b["t"]-a["t"])})

    report = {
        "schema": "strong_D_signed_pole_tracking_v1", "diagnostic_only": True,
        "source": {"artifact_id": 11706797173, "run_id": 38150925743,
                   "zip_sha256": "8b118d1f81e01f2a5d291d160b8bac42abc6a11051d8453efd9c19d0826d6d8b",
                   "N": N, "r_max": RMAX, "slice_count": len(samples),
                   "t_start": samples[0]["t"], "t_end": samples[-1]["t"]},
        "definitions": {"local_extrema": "three-point extrema of signed H_vendor in 1.25<=r<=20",
                        "identity": "nearest same-kind extremum between adjacent saved slices, |delta r|<=1.25",
                        "inner_competitor": "max |H| in 1.25<=r<4; r=1.25 is band edge, not necessarily a full-profile extremum",
                        "ridge_competitor": "max |H| in 4<=r<=20",
                        "sign_changes": "adjacent-cell sign transitions over 1.25<=r<=79.75; no deadband",
                        "guardrail": "tracked extrema are features of a constraint residual, not assumed physical objects"},
        "checks": checks, "all_checks_pass": all(checks.values()),
        "inner_outer_equal_amplitude_brackets": crossovers,
        "tracked_minimum_sign_crossings": crossings,
        "samples": [{"t": s["t"], "inner": s["inner"], "ridge": s["ridge"], "offcentre": s["offcentre"],
                     "sign_change_count": len(s["sign_edges"]), "sign_edges": s["sign_edges"],
                     "local_extrema": s["extrema"]} for s in samples],
        "tracked_minimum": min_track
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir/"signed_pole_tracking_report.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+"\n")
    lines = ["# Signed Hamiltonian-residual pole tracking", "",
        "Post-processing only: checksum-pinned R=80, N=160 saved production profiles. No evolution or production physics was changed.", "",
        f"Source artifact: 11706797173; run 38150925743; SHA-256: 8b118d1f81e01f2a5d291d160b8bac42abc6a11051d8453efd9c19d0826d6d8b.", "",
        "## Frozen checks", "", "| Check | Result |", "|---|---|"]
    lines += [f"| {k} | {'PASS' if v else 'FAIL'} |" for k, v in checks.items()]
    lines += ["", f"All checks passed: **{all(checks.values())}**.", "",
        "## Ridge-window winner: distinct local maxima", "",
        "| t | winning r | signed H | local maximum at r=5.75 |", "|---:|---:|---:|---:|"]
    for t in (5.75, 6., 6.25, 6.5, 7.5, 8.5, 8.75, 9., 12.):
        s = at(t)
        q = next((x for x in s["extrema"] if x["kind"] == "max" and abs(x["r"]-5.75) < 1e-8), None)
        lines.append(f"| {s['t']:.2f} | {s['ridge']['r']:.2f} | {s['ridge']['H']:+.6e} | {q['H']:+.6e} |" if q else f"| {s['t']:.2f} | {s['ridge']['r']:.2f} | {s['ridge']['H']:+.6e} | — |")
    lines += ["", "At t=6.00 the positive r=8.25 maximum and separate positive r=5.75 maximum coexist; at t=6.25 the r=8.25 maximum is gone and r=5.75 wins. The winning radius jumped because peak identity changed, not because one tracked maximum travelled between those radii.", "",
        "## Inner versus ridge-band dominance", "", "| t | inner winner (r, H) | ridge winner (r, H) |", "|---:|---:|---:|"]
    for t in (8., 8.25, 8.5, 8.75, 9., 12.):
        s = at(t)
        lines.append(f"| {s['t']:.2f} | ({s['inner']['r']:.2f}, {s['inner']['H']:+.6e}) | ({s['ridge']['r']:.2f}, {s['ridge']['H']:+.6e}) |")
    if crossovers:
        c = crossovers[0]
        lines += ["", f"Inner/ridge dominance crosses between t={c['t0']:.2f} and {c['t1']:.2f}; linear amplitude interpolation gives t≈{c['equal_abs_H_t_linear']:.4f}. Both competitors are positive here; this is not a sign reversal between the two winners.", ""]
    lines += ["## Tracked local minimum and profile sign changes", "",
        "| t | nearest tracked minimum r | signed H | profile sign-transition count |", "|---:|---:|---:|---:|"]
    for t in (0., 4., 5., 5.5, 5.75, 6., 6.25, 6.5, 7., 7.5, 9., 12.):
        s = at(t)
        entry = min(min_track, key=lambda x: abs(x["t"]-s["t"])) if min_track else None
        lines.append(f"| {s['t']:.2f} | {entry['r']:.2f} | {entry['H']:+.6e} | {len(s['sign_edges'])} |" if entry else f"| {s['t']:.2f} | — | — | {len(s['sign_edges'])} |")
    if crossings:
        c = crossings[0]
        lines += ["", f"The tracked minimum changes from {c['H0']:+.6e} at t={c['t0']:.2f}, r={c['r0']:.2f} to {c['H1']:+.6e} at t={c['t1']:.2f}, r={c['r1']:.2f}; linear interpolation estimates a zero at t≈{c['zero_t_linear']:.4f}. Profile sign-transition count falls from 2 to 0 over the same saved interval."]
    lines += ["", "## Conclusion", "",
        "Both observations occur in the saved Hamiltonian constraint residual: a local-minimum sign change and changes in peak identity/dominance. They are separately measured; timing alone does not show one causes the other. No material propagation, D-to-matter transfer, physical front, bounce or cycle is inferred.", ""]
    (args.output_dir/"signed_pole_tracking_report.md").write_text("\n".join(lines))
    print(json.dumps({"status": "PASS" if report["all_checks_pass"] else "FAIL", "checks": checks,
        "sign_crossings": crossings, "inner_outer_brackets": crossovers}, indent=2, sort_keys=True))
    if not report["all_checks_pass"]:
        raise RuntimeError("One or more frozen signed-pole checks failed")

if __name__ == "__main__":
    main()
