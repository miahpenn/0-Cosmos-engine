"""Compare stored-input PIRK rate arithmetic against a 60-digit replay."""
import argparse
from decimal import localcontext
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np

from research.w1_fring_smoke import capture_lapse_stages as capture
from research.w1_fring_smoke import matched_lapse, stage_momentum as replay
from research.w1_fring_smoke.live_lapse_audit import signs
from research.w1_fring_smoke.frozen_centre_derivative import shell_rms
from research.w1_fring_smoke.momentum_precision import wide

CELLS = replay.LOCAL_CELLS
SECOND_PARTS = ("ricci_a_second", "ricci_b_second", "ricci_X_second")
RATE_PARTS = ("hessian", *SECOND_PARTS, "ricci_ratio", "ricci_other", "L2_sum", "L3")
FIELDS = capture.PRIMARY_FIELDS


class FirstStepComplete(Exception):
    """Return from the pinned replay after its first accepted step."""


class CaptureFirst(capture.Capture):
    def __init__(self, vacuum, case):
        super().__init__(vacuum, case)
        self.end = .000375

    def _geometry_stage_terms(self, grid, geometry, scalars, matter, electromagnetic):
        call = self.geometry_calls
        terms = super()._geometry_stage_terms(
            grid, geometry, scalars, matter, electromagnetic
        )
        if self.capturing and call == 1:
            self.save_geometry("explicit", geometry)
            for key in FIELDS:
                self.data[f"explicit_primary_l2_{key}"] = terms["primary_l2"][key].copy()
        return terms


def capture_first(case, go, vacuum, meta):
    observer = CaptureFirst(vacuum, case)
    original = replay.advance
    accepted = {}

    def one_step(kernel, state, dt, watch):
        end, row = original(kernel, state, dt, watch)
        accepted.update(row)
        accepted["grid"] = end.grid
        raise FirstStepComplete

    with (patch.object(replay, "Observe", lambda _: observer),
          patch.object(replay, "advance", one_step)):
        try:
            replay.fine_replay(float(case), go, vacuum, meta)
        except FirstStepComplete:
            pass
    if accepted["t_end"] != observer.end:
        raise ValueError("first-step capture missed its target")
    return observer.data, accepted


def diff_at(values, i, h, parity=1, second=False):
    def at(j):
        return values[j] if j >= 0 else parity * values[-j-1]
    a, b, c, d, e = (at(j) for j in range(i-2, i+3))
    if second:
        return (-a + 16*b - 30*c + 16*d - e) / (12*h*h)
    return (a - 8*b + 8*d - e) / (12*h)


def centered_second(values, i, h):
    # Same even-parity stencil, with the constant cancelled before multiplying.
    a, b, c, d, e = (values[j if j >= 0 else -j-1] for j in range(i-2, i+3))
    return (16*((b-c) + (d-c)) - ((a-c) + (e-c))) / (12*h**2)


def l2_cell(r, state, first, second):
    # Use the same expression for the binary64 control and the widened readback.
    a, b, X, connection, alpha = state
    ap, bp, xp, connection_p, alphap = first
    app, bpp, xpp, alphapp = second
    half = type(r)("0.5")
    chip = -half * xp / X
    chipp = -half * (xpp/X - (xp/X)**2)
    metric_lambda = (1 - a/b) / r**2
    inverse = X**2 / a
    radial_bracket = (
        half*app/a - a*connection_p - (3*half/2)*(ap/a)**2 + half*(bp/b)**2
        - half*connection*ap + ap/(r*b)
        + 2*metric_lambda*(1 + r*bp/b) + 4*chipp
        - 2*chip*(ap/a - bp/b - 2/r)
    )
    scalar_bracket = (
        half*app/a + bpp/b - a*connection_p - (ap/a)**2 + half*(bp/b)**2
        + 2*(3 - a/b)*bp/(r*b) + 4*metric_lambda
        + 8*(chipp + chip**2) - 8*chip*(half*ap/a - bp/b - 2/r)
    )
    radial_hessian = inverse*(alphapp - alphap*(half*ap/a + 2*chip))
    lap = inverse*(alphapp - alphap*(half*ap/a - bp/b - 2*chip - 2/r))
    radial = -inverse*radial_bracket
    scalar = -inverse*scalar_bracket
    return -lap, -(radial_hessian - lap/3), alpha*(radial - scalar/3)


def same_input_l2(grid, vacuum, state, binary):
    fields = ("a", "b", "X", "Lambda", "alpha")
    parity = (1, 1, 1, -1, 1)
    stored = [getattr(state,f) for f in fields]
    first = [vacuum._d(grid,v,p) for v,p in zip(stored,parity)]
    second = [vacuum._d2(grid,stored[i],1) for i in (0,1,2,4)]
    y = [list(map(wide, v[:CELLS+2])) for v in stored]
    h = wide(grid.dr)
    out = {"K": {"hessian": []},
           "Aa": {key: [] for key in ("hessian", *SECOND_PARTS, "ricci_ratio", "ricci_other", "centered_second_control")}}
    for i in range(CELLS):
        r_float = grid.centers[i]
        r = wide(r_float)
        fp_state = [v[i] for v in stored]
        fp_first = [v[i] for v in first]
        fp_second = [v[i] for v in second]
        old = [binary["K"]["hessian"][i], binary["Aa"]["hessian"][i], binary["Aa"]["ricci"][i]]
        np.testing.assert_array_equal(l2_cell(r_float, fp_state, fp_first, fp_second), old)

        values = [v[i] for v in y]
        d1 = [diff_at(v,i,h,parity=p) for v,p in zip(y,parity)]
        d2 = [diff_at(y[j],i,h,second=True) for j in (0,1,2,4)]
        k, aa_hessian, aa_ricci = l2_cell(r,values,d1,d2)
        a, b, X, _, alpha = values
        app, bpp, xpp = (v-wide(f) for v,f in zip(d2[:3],fp_second[:3],strict=True))
        ratio = (1-a/b)/r**2 - wide((1-fp_state[0]/fp_state[1])/r_float**2)
        factor = -alpha*X**2/a
        second_delta = (factor*app/(3*a), -factor*bpp/(3*b), -2*factor*xpp/(3*X))
        stable_a, stable_b, stable_X = (
            wide(centered_second(v,i,grid.dr)) - wide(fp)
            for v,fp in zip(stored[:3],fp_second[:3],strict=True)
        )
        stable_delta = factor*(stable_a/(3*a) - stable_b/(3*b) - 2*stable_X/(3*X))
        ratio_delta = factor*ratio*(2*(1+r*d1[1]/b) - wide(4)/3)

        out["K"]["hessian"].append(k-wide(old[0]))
        out["Aa"]["hessian"].append(aa_hessian-wide(old[1]))
        for part, correction in zip(SECOND_PARTS, second_delta, strict=True):
            out["Aa"][part].append(correction)
        out["Aa"]["centered_second_control"].append(stable_delta)
        out["Aa"]["ricci_ratio"].append(ratio_delta)
        out["Aa"]["ricci_other"].append(aa_ricci-wide(old[2])-sum(second_delta)-ratio_delta)
    return out


def stage_l2(grid, vacuum, snapshot, name):
    state = matched_lapse.fields_from(snapshot, name, vacuum)
    binary_geometry = vacuum.geometry_terms(grid, state)
    binary = {
        "K": {"hessian": -binary_geometry["lap_alpha"]},
        "Aa": {
            "hessian": -(binary_geometry["DrDr_alpha"] - binary_geometry["lap_alpha"]/3),
            "ricci": state.alpha*(binary_geometry["Rr_r"] - binary_geometry["R"]/3),
        },
    }
    expanded = same_input_l2(grid, vacuum, state, binary)
    result = {}
    for key in FIELDS:
        original = snapshot[f"{name}_primary_l2_{key}"]
        np.testing.assert_array_equal(sum(binary[key].values()), original)
        result[key] = expanded[key]
        exact_binary_sum = [sum(map(wide, values)) for values in zip(*binary[key].values())]
        result[key]["L2_sum"] = [w - wide(b) for w,b in zip(exact_binary_sum,original)]
    return result


def stage_l3(grid, vacuum, snapshot, name):
    state = matched_lapse.fields_from(snapshot, name, vacuum)
    source = matched_lapse.matter_source(snapshot, name)
    baseline = matched_lapse.l3_rhs(vacuum, grid, state, source)
    result = {key: [] for key in FIELDS}
    pi = wide(np.pi)
    for key in FIELDS:
        np.testing.assert_array_equal(baseline[key], snapshot[f"{name}_primary_l3_{key}"])
    for i in range(CELLS):
        alpha, aa, k = map(wide, (state.alpha[i], state.Aa[i], state.K[i]))
        rho, pr, pt = (wide(snapshot[f"{name}_{f}"][i]) for f in ("rho", "pr", "pt"))
        new = {
            "Aa": alpha*k*aa - 16*pi/3*alpha*(pr-pt),
            "K": alpha*(3*aa**2/2 + k*k/3) + 4*pi*alpha*(rho + pr + 2*pt),
        }
        for key in FIELDS:
            result[key].append(new[key] - wide(baseline[key][i]))
    return result


def analyze(grid, vacuum, snapshot, stages):
    assert not np.any(snapshot["start_beta"]) and not np.any(snapshot["primary_beta"])
    h = wide(grid.dr)
    dt = wide(snapshot["dt"])
    n = CELLS-2
    r = grid.centers[:n]
    L2 = {s: stage_l2(grid, vacuum, snapshot, s)
          for s in ("start", "explicit", "corrector_copy")}
    L3 = {s: stage_l3(grid, vacuum, snapshot, s) for s in ("start", "primary")}
    outputs = {}
    for stage, final_l2, weights in (
        ("primary_before_boundary", "explicit", (dt, 0)),
        ("corrector_before_boundary", "corrector_copy", (dt/2, dt/2)),
    ):
        final = matched_lapse.fields_from(snapshot, final_l2, vacuum)
        dx, db = (vacuum._d(grid, f, 1)[:n] for f in (final.X, final.b))
        coefficient = 3/r - 3*dx/final.X[:n] + 1.5*db/final.b[:n]
        correction = {}
        for part in (*RATE_PARTS, "centered_second_control"):
            fields = {}
            for field in FIELDS:
                if part == "L3":
                    values = [weights[0]*a + weights[1]*b
                              for a,b in zip(L3["start"][field], L3["primary"][field])]
                else:
                    left = L2["start"][field].get(part, [0]*CELLS)
                    right = L2[final_l2][field].get(part, [0]*CELLS)
                    values = [dt*(a+b)/2 for a,b in zip(left,right)]
                diff = [diff_at(values,i,h) for i in range(n)]
                if field == "K":
                    fields[field] = np.array([float(-2*value/3) for value in diff])
                else:
                    fields[field] = np.array([float(d+wide(c)*v)
                                              for d,c,v in zip(diff,coefficient,values)])
            correction[part] = {**fields, "sum": sum(fields.values())}

        def stats(values):
            shell = shell_rms(r, values)
            return {"signed": values.tolist(), "band_rms": shell["rms"],
                    "rough_band_rms": shell["rough_rms"]}

        baseline = np.array(stages[stage]["rows"]["M"][:n])
        stable = correction.pop("centered_second_control")["sum"]
        second_only = sum(correction[p]["sum"] for p in SECOND_PARTS)
        all_parts = sum(v["sum"] for v in correction.values())
        shifted = baseline + all_parts
        hessian_only = baseline + correction["hessian"]["sum"]
        outputs[stage] = {
            "changes": {part: {key: stats(v) for key,v in pieces.items()}
                        for part,pieces in correction.items()},
            "all_parts": stats(all_parts), "production_M": stats(baseline),
            "shadow_M": stats(shifted), "production_shell_signs": signs(r,baseline),
            "shadow_shell_signs": signs(r,shifted),
            "hessian_only_M": stats(hessian_only),
            "hessian_only_signs": signs(r,hessian_only),
            "ricci_second_total": stats(second_only),
            "centered_second_control": stats(stable),
            "centered_minus_wide_second": stats(stable-second_only),
        }
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("0.025", "0.0125"))
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    meta = replay.provenance()
    go, vacuum, _ = replay.pins.ad.vendor_modules()
    if vacuum.KO_EPSILON != 0:
        raise ValueError("expected KO off")
    snapshot, accepted = capture_first(args.case, go, vacuum, meta)
    previous_path = Path("results/research/w1_fring_stage_momentum")/f"{args.case}.json"
    previous = json.loads(previous_path.read_text())
    if accepted["stages"] != previous["observed_steps"][0]["stages"]:
        raise ValueError("first-step stage replay differs from archived production")
    with localcontext() as ctx:
        ctx.prec = 60
        result = analyze(accepted["grid"],vacuum,snapshot,accepted["stages"])
    report = {"engine": replay.pins.SOURCE, "vendor": replay.pins.VENDOR,
              "input_sha256": previous["input_sha256"], "script_sha256": replay.sha(Path(__file__)),
              "reference_sha256": replay.sha(previous_path),
              "dt": snapshot["dt"], "dr": snapshot["dr"], "r": previous["r"][:CELLS-2],
              "stages": result}
    args.output.write_text(json.dumps(report,allow_nan=False,separators=(",",":"))+"\n")


if __name__ == "__main__":
    main()
