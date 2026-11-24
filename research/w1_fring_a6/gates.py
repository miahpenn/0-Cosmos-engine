"""Gates and readings for the f_ring A6 runs (A6 1Pc5Md4LCSX7_8ZHY3fAQDbsxl2soth5mPDMsZ69sOZ4 with A7
1dr_BKBvWuk5nh_wVMMFOPSiPtAAAfpstJg9NHypnR9U). Pure functions, fixed before any A6 step; synthetic.py exercises
each one on cases that must pass and cases that must fail. How each design line is read: PRERUN.md.
"""
import math

import numpy as np

FIRST_CHECK_TARGET = 8.5030e-11  # direct nine-point H, ring, dr 0.05, dt 1.5e-3, cell 0, t 0.015 (areal trace)
FIRST_CHECK_TOL = 1e-3           # relative
RELATIVE = 1e-3                  # A6 3(i)
FACTOR = 2.0                     # A6 3(ii): 2^(q-1), q = 2
FLOOR_MULT = 10                  # A7 2(c)


def first_check(value):
    """Direct nine-point H at t 0.015 must match 8.503e-11 within 1e-3 (relative) before anything is read."""
    gap = abs(value/FIRST_CHECK_TARGET - 1)
    return {"value": value, "target": FIRST_CHECK_TARGET, "relative_gap": gap, "pass": bool(gap <= FIRST_CHECK_TOL)}


def relative_size(value, terms):
    """A6 3(i): |C| <= 1e-3 x sum |terms| at each cell. ``terms`` is a list of arrays over the same cells."""
    value = np.abs(np.asarray(value, dtype=float))
    bound = RELATIVE*np.sum(np.abs(np.asarray(terms, dtype=float)), axis=0)
    ratio = np.divide(value, bound, out=np.full_like(value, np.inf), where=bound > 0)
    return {"ratio_to_bound": ratio.tolist(), "pass": bool(np.all(value <= bound))}


def resolution(value):
    return float(np.spacing(abs(float(value))))


def order_reading(dH_direct, write, increment, closure, H_fine):
    """A6 3(ii) read under A7 2(b,c) at cell 0 for one dt pair at one time.

    dH_direct: signed direct-curvature H difference, coarse minus fine dt. write, increment: the measured write and
    increment-arithmetic curvature responses for the pair. closure: closure magnitude of the split. H_fine: the
    finer run's direct H. Rate part = dH_direct - write - increment. Floor = 10 x max(closure, resolution).
    Halving dt must reduce H by >= 2x: (H_fine + rate)/H_fine >= 2, with H_fine + rate the coarse run's H
    without writes. Below the floor: "unresolved"."""
    rate = dH_direct - write - increment
    floor = FLOOR_MULT*max(abs(closure), resolution(rate))
    out = {"rate_part": rate, "floor": floor}
    if abs(rate) < floor:
        out.update(reading="unresolved", factor=None, q=None)
        return out
    factor = (H_fine + rate)/H_fine if H_fine != 0 else math.inf
    out["factor"] = factor
    out["q"] = 1 + math.log2(factor) if factor > 0 and math.isfinite(factor) else None
    out["reading"] = "pass" if factor >= FACTOR else "fail"
    return out


def growth(t, h):
    """A6 3(iii)/A5 1(iii) exponential test on one arm over the ring period.

    t, h: times (>= 0.015) and max |H| over cells 0-5. Fits ln h = a + g t (exponential) and ln h = a + p ln t
    (power law) by least squares. "exponential" when the exponential fit has the smaller residual AND the fitted
    growth over the span is at least a factor e (g x span >= 1)."""
    t, h = np.asarray(t, dtype=float), np.asarray(h, dtype=float)
    y = np.log(h)
    fits = {}
    for name, x in (("exponential", t), ("power", np.log(t))):
        A = np.vstack([np.ones_like(x), x]).T
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        fits[name] = {"slope": float(coef[1]), "rms": float(np.sqrt(np.mean((A @ coef - y)**2)))}
    span = float(t[-1] - t[0])
    exponential = (fits["exponential"]["rms"] < fits["power"]["rms"]
                   and fits["exponential"]["slope"]*span >= 1)
    return {"fits": fits, "span": span, "growth_over_span": float(h[-1]/h[0]),
            "growing": bool(h[-1] > h[0]), "exponential": bool(exponential)}


def instability(finite, growth_result, order_fail_at_end):
    """A6 3(iii): stop on non-finite values, exponential growth over the period, or no fall with dt refinement
    (a resolved order fail at the last save) while growing."""
    reasons = []
    if not finite:
        reasons.append("non-finite")
    if growth_result["exponential"]:
        reasons.append("exponential growth over the ring period")
    if order_fail_at_end and growth_result["growing"]:
        reasons.append("no fall with dt refinement while growing")
    return {"stop": bool(reasons), "reasons": reasons}
