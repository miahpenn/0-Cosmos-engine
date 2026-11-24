"""Show every A6 gate and reading can discriminate before any step (pass and fail cases). Exits 1 on a miss."""
import json

import numpy as np

from research.w1_fring_a6 import gates


def main():
    out, ok = {}, True

    def case(name, got, expected):
        nonlocal ok
        good = got == expected
        ok &= good
        out[name] = {"got": got, "expected": expected, "pass": good}

    T = gates.FIRST_CHECK_TARGET
    case("first check, 5e-4 off", gates.first_check(T*(1+5e-4))["pass"], True)
    case("first check, 2e-3 off", gates.first_check(T*(1-2e-3))["pass"], False)
    case("first check, sign flipped", gates.first_check(-T)["pass"], False)

    terms = [np.array([1e-8, 2e-8]), np.array([-1e-8, 3e-8]), np.array([5e-9, 1e-9])]
    total = np.sum(np.abs(terms), axis=0)
    case("3(i) at 5e-4 of terms", gates.relative_size(5e-4*total, terms)["pass"], True)
    case("3(i) at 2e-3 of terms at one cell", gates.relative_size(np.array([5e-4, 2e-3])*total, terms)["pass"], False)

    H = 2e-11
    case("3(ii) factor 4", gates.order_reading(3*H + 1e-13, 1e-13, 0., 1e-15, H)["reading"], "pass")
    case("3(ii) factor 1.5", gates.order_reading(.5*H, 0., 0., 1e-15, H)["reading"], "fail")
    case("3(ii) factor 4 hidden by writes", gates.order_reading(3*H, 2.6*H, 0., 1e-15, H)["reading"], "fail")
    case("3(ii) below floor", gates.order_reading(1e-14, 0., 0., 1e-14, H)["reading"], "unresolved")
    case("3(ii) no fall (sign reversed)", gates.order_reading(-.5*H, 0., 0., 1e-15, H)["reading"], "fail")

    t = np.linspace(.015, 3.291, 200)
    rng = np.random.default_rng(6)
    noise = np.exp(.02*rng.normal(size=t.size))
    case("3(iii) exponential e^t", gates.growth(t, 1e-11*np.exp(t)*noise)["exponential"], True)
    case("3(iii) power t^1.3", gates.growth(t, 1e-11*t**1.3*noise)["exponential"], False)
    case("3(iii) flat", gates.growth(t, 1e-11*noise)["exponential"], False)
    case("3(iii) slow exponential (factor 1.2 over the span)",
         gates.growth(t, 1e-11*np.exp(.0555*t)*noise)["exponential"], False)
    grow = gates.growth(t, 1e-11*t**1.3*noise)
    case("3(iii) stop: non-finite", gates.instability(False, grow, False)["stop"], True)
    case("3(iii) stop: order fail while growing", gates.instability(True, grow, True)["stop"], True)
    case("3(iii) no stop: power law, order passes", gates.instability(True, grow, False)["stop"], False)

    out["all_discriminate"] = ok
    print(json.dumps(out, indent=1))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
