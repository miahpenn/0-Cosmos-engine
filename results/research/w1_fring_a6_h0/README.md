# W1 starting-state H wiggle (t = 0, zero steps)

Jeremiah's go (07:24 UTC, via the coordinator): "Have Wright do it". No step is taken. `h0_build.py` rebuilds the t = 0 states with `kernel.initialize` only, and its H matches the archived replays' t = 0 H bitwise at dr .025 and .0125. `h0_mp.py` reads them in 60 digits.

## Result

The cell-to-cell H wiggle at dr .0125 is binary64 rounding in the stored geometry, which the stencils then amplify. It is not in the supplied S, PS, D, PD, and it is not truncation.

Band native cells r in [0.05, 0.25]:

| dr | H binary64 | H 60 digits on stored inputs | H 60 digits on rebuilt geometry |
|---|---|---|---|
| .05 | 0/3, RMS 2.18e-11 | 0/3, 2.19e-11 | 0/3, 2.18e-11 |
| .025 | 0/7, 1.66e-12 | 0/7, 1.47e-12 | 0/7, 1.38e-12 |
| .0125 | 12/15, 4.25e-12 | 11/15, 3.53e-12 | **0/15, 8.4e-14** |

- **Rebuilt geometry:** a, b, X and Lambda redone in 60 digits from the areal loop's last binary64 source. That means exact eight-node weights, an exact cumulative sum and exact powers. Aa, K and rho are kept as stored, since none of them is differentiated in H.
- **Stored against rebuilt:** a, b and X differ by at most 1 ulp (2.2e-16), and Lambda by 3e-14.
- **Scaling:** H on the rebuilt geometry is smooth and falls 16x per halving (2.2e-11, 1.4e-12, 8.4e-14). That is fourth-order truncation.
- **Rounding floor:** the rounding is about ε/dr² times the stencil weights. That gives about 7e-12 at .0125 and 2e-12 at .025.
  - At .0125 it is 40x above truncation, so it shows as alternation.
  - At .025 it is about the same size as truncation and doesn't flip signs.
  - At .05 it is hidden.
- **Evaluation rounding:** binary64 against 60 digits on the same inputs gives an RMS of 2.5e-12 at .0125 and 4.4e-13 at .025.

## Where it enters (code read)

`engine/v55_initial.py` `build_initial_data`, supplied-slice branch, `eight_node_areal=True`:
- **Lines 373–414:** the areal loop. `_eight_node_origin_integral` (lines 130–155) sums cell pieces with `np.cumsum`, then `B = integ / r`.
- **Lines 420–430:** `A = 1/sqrt(B)`, `a = A**(4/3)`, `b = A**(-2/3)`, `X = A**(-1/3)`, and `Lambda` from the 4th-order `_d1` of a and b.

Each of these is stored in binary64. H (vendor `constraints`) takes second derivatives of a, b and X and the first derivative of Lambda, so ulp-level rounding becomes about ε/dr². The scalar inputs and the HB solve are not the source: rebuilding from the same binary source removes the wiggle.

## Stage-1 controls do not use this path

`research/control_tests/ctl.py` `make()`, lines 31–38, calls `initialize(initial_structure="none")`. It then sets every vacuum field (`VAC.ALL_FIELDS`) to zero and a, b, X and alpha to exactly 1.0. 1b (`t0_checks.build_1b`, lines 17–19) sets K, X and alpha in closed form. No stage-1 geometry comes from the areal loop.

The generic effect, binary64 rounding amplified by ε/dr², is what the 1a floor arms and the B3/B5 bounds measure by design.
