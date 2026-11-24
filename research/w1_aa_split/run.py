"""Split the W1 tracefree rate alpha*A_a into the coded A_a rate terms and attribute its curvature difference.

DESIGN-W1-A_a-rate-term-split-predeclared-Nov24-2026 (Drive 1tt1m0E8bkTvUGbsGsWG06DIiK2w3a5YEDEP-9XQRpTI).
Replays the conformal_ledger run (#587/#588) unchanged, captures every A_a rate the production step uses,
re-expresses each one as its coded terms (tags.py), checks the closure gate on cells 0-5 at every evaluation,
and carries each tag through the telescoped alpha*A_a into the ledger's tracefree tangent.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from engine import v55_matter as vm
from engine.v55_pirk_adapter import vendor_modules
from research.w1_aa_split import readings, tags as tagdef
from research.w1_fring_smoke import areal_trace, conformal_ledger, metric_balance, origin_trace, x_endpoint

WIDE = metric_balance.WIDE
LEDGER = Path("results/research/w1_fring_conformal_ledger/ledger.json")
CELLS = slice(0, 6)
OWN = [Path(__file__), Path(readings.__file__), Path(tagdef.__file__)]


def source_hashes():
    hashes = conformal_ledger.source_hashes()
    hashes.update({f"w1_aa_split/{p.name}": hashlib.sha256(p.read_bytes()).hexdigest() for p in OWN})
    hashes[LEDGER.name] = hashlib.sha256(LEDGER.read_bytes()).hexdigest()
    return hashes


def wide(x):
    return np.asarray(x, dtype=WIDE)


def run(divisor, directory):
    hashes = source_hashes()
    reference = json.loads(conformal_ledger.REFERENCE.read_text())["runs"][str(divisor)]["readback"]
    q = conformal_ledger
    archived = (q.quarter_stage_trace.RESULTS / "quarter.json" if divisor == 4 else
                q.quarter_trace.PREVIOUS / ("A-0.025-half.json" if divisor == 2 else "A-0.025.json"))
    prior = json.loads(archived.read_text())["quarter_rows" if divisor == 4 else "rows"]
    _, vacuum, _ = vendor_modules()
    n = round(60/x_endpoint.DR)
    dt = .03*x_endpoint.DR/divisor
    wdt = WIDE(dt)
    steps = 20*divisor
    names = tagdef.NAMES

    ledger_totals = np.zeros((2, n), dtype=WIDE)            # conformal_ledger "tangent" term, identical arithmetic
    tangent = {k: np.zeros(n, dtype=WIDE) for k in names + ("initial", "remainder")}
    F = None                                                 # alpha*A_a decomposition at the current start stage
    pending = None                                           # corrector pieces waiting for the next start alpha
    gate = {"evaluations": 0, "worst": {}, "failures": []}
    pieces = {"rate_closure": 0., "update_arithmetic": 0., "boundary": 0.}
    per_step = {}
    seen, captured, l2_calls, stages = [], [], [], {}
    state = {"initial": None, "endpoint": None, "f_init": None, "previous_Aa": None, "previous_ab": None, "step": 0}

    def check_gate(kind, parts, production):
        passed, ratio = readings.closure([v[CELLS] for v in parts.values()], production[CELLS])
        gate["evaluations"] += 1
        gate["worst"][kind] = max(gate["worst"].get(kind, 0.), ratio)
        if not passed:
            gate["failures"].append({"step": state["step"], "evaluation": kind, "gap_over_bound": ratio})

    def l2_wrapper(original):
        def wrapped(grid, s):
            out = original(grid, s)
            l2_calls.append((next(reversed(stages), None), grid, s.copy(), np.array(out["Aa"], copy=True)))
            return out
        return wrapped

    stage_terms = origin_trace.Observer._geometry_stage_terms

    def capture(self, grid, geometry, scalars, matter, electromagnetic):
        out = stage_terms(self, grid, geometry, scalars, matter, electromagnetic)
        if electromagnetic is not None:
            raise ValueError("the W1 ring has no electromagnetic sector")
        captured.append(dict(grid=grid, geometry=geometry.copy(), scalars=scalars, matter=matter,
                             l2=np.array(out["primary_l2"]["Aa"], copy=True),
                             l3=np.array(out["primary_l3"]["Aa"], copy=True)))
        return out

    def record(name, grid, geometry):
        nonlocal F, pending
        if name == "after_initial_lapse_refresh":
            stages.clear()
            captured.clear()
            l2_calls.clear()
        stages[name] = geometry
        if name != "corrected_after_boundary":
            return
        start, predictor = stages["after_initial_lapse_refresh"], stages["explicit_predictor"]
        primary_pre, primary = stages["primary_predictor_before_boundary"], stages["primary_predictor_after_lapse"]
        corrected_pre, corrected = stages["corrected_before_boundary"], geometry
        late = [call[1:] for call in l2_calls if call[0] == "primary_predictor_after_lapse"]
        if len(captured) != 3 or len(late) != 2:
            raise ValueError(f"unexpected number of rate evaluations in a step: {len(captured)} {len(late)}")
        for call, stage in zip(captured, (start, predictor, primary)):
            np.testing.assert_array_equal(call["geometry"].Aa, stage.Aa)
            np.testing.assert_array_equal(call["geometry"].alpha, stage.alpha)
        np.testing.assert_array_equal(late[0][2], captured[2]["l2"])
        final_grid, final_state, l2_final = late[1]
        for key in ("a", "b", "X"):
            np.testing.assert_array_equal(getattr(final_state, key), getattr(corrected, key))
        np.testing.assert_array_equal(final_state.alpha, primary.alpha)

        # Production A_a updates, reproduced bitwise from the captured rates.
        l2_0, l2_p, l3_0, l3_1 = captured[0]["l2"], captured[1]["l2"], captured[0]["l3"], captured[2]["l3"]
        np.testing.assert_array_equal(start.Aa + dt*(0.5*l2_0 + 0.5*l2_p + l3_0), primary_pre.Aa)
        np.testing.assert_array_equal(start.Aa + 0.5*dt*(l2_0 + l2_final + l3_0 + l3_1), corrected_pre.Aa)

        # Coded terms of every rate evaluation, each gated on cells 0-5.
        t_l2_0 = tagdef.l2_terms(vacuum, captured[0]["grid"], captured[0]["geometry"])
        t_l2_p = tagdef.l2_terms(vacuum, captured[1]["grid"], captured[1]["geometry"])
        t_l2_f = tagdef.l2_terms(vacuum, final_grid, final_state)
        t_l3_0 = tagdef.l3_terms(vacuum, vm, captured[0]["grid"], captured[0]["geometry"],
                                 captured[0]["scalars"], captured[0]["matter"])
        t_l3_1 = tagdef.l3_terms(vacuum, vm, captured[2]["grid"], captured[2]["geometry"],
                                 captured[2]["scalars"], captured[2]["matter"])
        for kind, parts, production in (("L2 start", t_l2_0, l2_0), ("L2 explicit predictor", t_l2_p, l2_p),
                                        ("L2 final", t_l2_f, l2_final), ("L3 start", t_l3_0, l3_0),
                                        ("L3 primary", t_l3_1, l3_1)):
            check_gate(kind, parts, production)

        def rates(*weighted):
            out = {k: np.zeros(n, dtype=WIDE) for k in names}
            for weight, parts in weighted:
                for k, v in parts.items():
                    out[k] += WIDE(weight)*v
            return out
        primary_rate = rates((.5, t_l2_0), (.5, t_l2_p), (1, t_l3_0))
        corrector_rate = rates((.5, t_l2_0), (.5, t_l2_f), (.5, t_l3_0), (.5, t_l3_1))
        production_primary = WIDE(.5)*wide(l2_0) + WIDE(.5)*wide(l2_p) + wide(l3_0)
        production_corrector = WIDE(.5)*(wide(l2_0) + wide(l2_final) + wide(l3_0) + wide(l3_1))

        a0, A0 = wide(start.alpha), wide(start.Aa)
        if F is None:
            state["initial"] = np.stack([getattr(start, k) for k in conformal_ledger.FIELDS]).astype(WIDE)
            state["f_init"] = a0*A0
            F = {k: np.zeros(n, dtype=WIDE) for k in names}
        else:
            np.testing.assert_array_equal(start.Aa, state["previous_Aa"])
            np.testing.assert_array_equal(np.stack([start.a, start.b]), state["previous_ab"])
            rate, dA, closure_part, arithmetic, boundary, A_prev, alpha_prev = pending
            for k in names[:-1]:
                F[k] = F[k] + a0*wdt*rate[k]
            F["lapse"] = F["lapse"] + A_prev*(a0 - alpha_prev)
            pieces["rate_closure"] = max(pieces["rate_closure"], float(np.max(abs(a0*closure_part)[CELLS])))
            pieces["update_arithmetic"] = max(pieces["update_arithmetic"], float(np.max(abs(a0*arithmetic)[CELLS])))
            pieces["boundary"] = max(pieces["boundary"], float(np.max(abs(a0*boundary)[CELLS])))
        f0 = a0*A0

        a1, A1 = wide(primary.alpha), wide(primary.Aa)
        F1 = {k: F[k] + a1*wdt*primary_rate[k] for k in names[:-1]}
        F1["lapse"] = F["lapse"] + A0*(a1 - a0)
        f1 = a1*A1
        for label, value in (("rate_closure", wdt*(production_primary - sum(primary_rate[k] for k in names[:-1]))),
                             ("update_arithmetic", wide(primary_pre.Aa) - A0 - wdt*production_primary),
                             ("boundary", A1 - wide(primary_pre.Aa))):
            pieces[label] = max(pieces[label], float(np.max(abs(a1*value)[CELLS])))

        # The ledger's tangent, with conformal_ledger's arithmetic, and its tag split.
        ledger_rates = np.array([[-2], [1]], dtype=WIDE)
        ledger_totals[...] += conformal_ledger.coords((WIDE(dt)/2)*ledger_rates*(f0 + f1))
        half = WIDE(dt)/2
        step_parts = {k: half*(F[k] + F1[k]) for k in names}
        step_parts["initial"] = half*(state["f_init"] + state["f_init"])
        reconstructed0 = state["f_init"] + sum(F.values())
        reconstructed1 = state["f_init"] + sum(F1.values())
        step_parts["remainder"] = half*((f0 - reconstructed0) + (f1 - reconstructed1))
        for k, v in step_parts.items():
            tangent[k] += v
        index = state["step"]
        if index in (0, steps - 1):
            per_step["first" if index == 0 else "last"] = {k: metric_balance.numbers(v[CELLS])
                                                           for k, v in step_parts.items()}

        # Corrector pieces; applied with the next step's refreshed lapse.
        dA = wide(corrected.Aa) - A0
        pending = (corrector_rate, dA,
                   wdt*(production_corrector - sum(corrector_rate[k] for k in names[:-1])),
                   wide(corrected_pre.Aa) - A0 - wdt*production_corrector,
                   wide(corrected.Aa) - wide(corrected_pre.Aa), A0, a0)
        state["previous_Aa"] = np.array(corrected.Aa, copy=True)
        state["previous_ab"] = np.stack([corrected.a, corrected.b])
        state["endpoint"] = metric_balance.metric_arrays(corrected, WIDE)
        state["step"] += 1

    def accepted(st, vac, row):
        index = len(seen)
        if index < len(prior) or index == steps:
            actual = areal_trace.observe(st, vac, row)
            expected = prior[index] if index < len(prior) else reference
            if actual != expected:
                raise ValueError(f"archived readback disagrees at step {index}")
        seen.append(st.t)

    original_l2 = vacuum.primary_l2_rhs
    vacuum.primary_l2_rhs = l2_wrapper(original_l2)
    origin_trace.Observer._geometry_stage_terms = capture
    try:
        result = origin_trace.run(x_endpoint.DR, "A", divisor, end_time=x_endpoint.END,
                                  accepted_observer=accepted, stage_observer=record, stage_step=None)
    finally:
        vacuum.primary_l2_rhs = original_l2
        origin_trace.Observer._geometry_stage_terms = stage_terms
    if state["step"] != steps or source_hashes() != hashes:
        raise ValueError("unexpected number of steps or changed inputs")
    split_sum = sum(tangent.values())
    np.testing.assert_allclose(split_sum, ledger_totals[1], atol=1e-24, rtol=0)
    metadata = {"attempts": result["attempts"], "time": seen[-1], "sha256": hashes, "gate": gate,
                "non_tag_pieces_max_cells_0_5": pieces, "per_step": per_step,
                "split_minus_ledger_max": float(np.max(abs(split_sum - ledger_totals[1])))}
    np.savez(directory/f"{divisor}.npz", endpoint=np.stack(state["endpoint"]), initial=state["initial"],
             ledger=ledger_totals, metadata=json.dumps(metadata, allow_nan=False),
             **{f"tag_{k}": v for k, v in tangent.items()})
    if gate["failures"]:
        raise RuntimeError(f"closure gate failed: {gate['failures'][:3]}")


def combine(directory):
    hashes = source_hashes()
    ledger = json.loads(LEDGER.read_text())
    grid = vendor_modules()[0].SphericalCellGrid(2400, 60.)
    keys = tagdef.NAMES + ("initial", "remainder")
    endpoints, tangents, ledgers, runs = {}, {}, {}, {}
    for divisor in (1, 2, 4):
        with np.load(directory/f"{divisor}.npz", allow_pickle=False) as data:
            runs[divisor] = json.loads(data["metadata"].item())
            if runs[divisor].pop("sha256") != hashes or runs[divisor]["attempts"] != 20*divisor:
                raise ValueError("inconsistent replay provenance")
            if runs[divisor]["gate"]["failures"]:
                raise ValueError("closure gate failed; no reading")
            endpoints[divisor] = data["endpoint"]
            ledgers[divisor] = data["ledger"]
            tangents[divisor] = {k: data[f"tag_{k}"] for k in keys}

    pairs = {}
    for left, right in ((1, 2), (2, 4)):
        lv, rv = endpoints[left], endpoints[right]
        x = (lv[2] + rv[2])/2
        pair = f"{left}-{right}"
        widths = {}
        for width in areal_trace.WIDTHS:
            base = metric_balance.curvature(grid, (*rv[:2], x), width)

            def response(tracefree):
                a = rv[0]*np.exp(0/3 - 2*tracefree)
                b = rv[1]*np.exp(0/3 + tracefree)
                return metric_balance.curvature(grid, (a, b, x), width) - base

            ledger_delta = ledgers[left][1] - ledgers[right][1]
            ledger_response = response(ledger_delta)
            archived = ledger["pairs"][pair]["curvature"][str(width)]["contributions"]["tracefree"]["tangent"]
            reproduced = metric_balance.numbers(ledger_response)
            if reproduced != archived:
                raise ValueError(f"D,T ledger tracefree tangent not reproduced for {pair} width {width}")
            delta = {k: tangents[left][k] - tangents[right][k] for k in keys}
            per_tag = {k: response(delta[k]) for k in tagdef.NAMES}
            others = {k: response(delta[k]) for k in ("initial", "remainder")}
            c_split = sum(per_tag.values()) - ledger_response
            coarse = {g: response(sum(delta[k] for k in members)) for g, members in tagdef.COARSE.items()}
            coarse_c = sum(coarse.values()) - ledger_response
            widths[width] = {
                "ledger_tracefree_tangent": reproduced,
                "c_split": metric_balance.numbers(c_split),
                "reproduces_target_within_10x_c_split": bool(
                    abs(sum(per_tag.values())[0] - archived[0]) <= 10*abs(c_split[0])),
                "tags": {k: metric_balance.numbers(v) for k, v in per_tag.items()},
                "non_tags": {k: metric_balance.numbers(v) for k, v in others.items()},
                "coarse": {k: metric_balance.numbers(v) for k, v in coarse.items()},
                "coarse_c_split": metric_balance.numbers(coarse_c),
            }
        pairs[pair] = widths

    out_readings = {}
    for width in areal_trace.WIDTHS:
        for label, key, ckey, origin in (("primary", "tags", "c_split", tagdef.ORIGIN),
                                         ("coarse", "coarse", "coarse_c_split", tagdef.COARSE_ORIGIN)):
            data = {p: {"tags": {k: v[0] for k, v in pairs[p][width][key].items()},
                        "total": pairs[p][width]["ledger_tracefree_tangent"][0],
                        "c_split": pairs[p][width][ckey][0]} for p in readings.PAIRS}
            out_readings[f"{label}/{width}"] = readings.classify(data, origin)
    print(json.dumps({"t": x_endpoint.END, "dr": x_endpoint.DR, "sha256": hashes, "source": origin_trace.SOURCE,
                      "vendor": origin_trace.VENDOR, "input_sha256": origin_trace.DIGEST,
                      "origin_tags": sorted(tagdef.ORIGIN), "coarse_groups": {k: list(v) for k, v in tagdef.COARSE.items()},
                      "runs": runs, "pairs": pairs, "readings": out_readings}, indent=2, allow_nan=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--divisor", type=int, choices=(1, 2, 4))
    args = parser.parse_args()
    if np.finfo(WIDE).nmant < 63:
        raise RuntimeError("at least 64 significand bits required")
    if args.divisor is None:
        combine(args.directory)
    else:
        run(args.divisor, args.directory)
