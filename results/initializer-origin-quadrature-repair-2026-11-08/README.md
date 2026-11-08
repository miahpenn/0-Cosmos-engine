# Initializer origin-quadrature repair: evidence (2026-11-08)

Evidence for the production PR that makes the areal initial-data reconstruction's cumulative integral
fourth order at the origin (`engine/v55_initial.py`; diff in `initializer_fix.diff`, SHA-256
`2dda1891cba55c76f52610a4366959e3b4c21cfbc3a2786ee5b13d49ef423afc`).

## Mechanism

The initializer reconstructs B from (r B)' = src, integrating outward with a rectangle over [0, r0] and
plain trapezoids afterwards. Near the origin src ≈ 1 − c r² is even, and each of those first steps leaves an
O(h³) constant in the integral. A constant in r B is the homogeneous solution B ⊃ ε/r: smooth in the
continuum, but odd and singular at the origin, so the even-parity fourth-order stencils at cells 0–2 turn
ε ~ h³ into an O(1) Hamiltonian-constraint violation. Through the propagation identity's ⅙α D_rH − ⅓α′H
terms, that violation then drives a centre momentum constraint during evolution that does not converge.

Fixing only the first interval (exact integral of the even quadratic through cells 0–1) lowered cell 0 from
4.06e-3 to 7.0e-5 but left cell 1 at 3.8e-4, also non-converging, because the nearby trapezoid steps carry
the same kind of constant. The repair therefore makes the whole cumulative rule fourth order: the exact
integral of the even quartic through cells 0–2 over [0, r0], then trapezoids with the Euler–Maclaurin end
correction −(h²/12)[f′(r_i) − f′(r0)], f′ from the grid's even-parity derivative (f′(0) = 0 by evenness).

## Evolution comparison

`tools/initializer_origin_quadrature_evolution_comparison.py` (SHA-256
`2e35ec32510e27fc1db7c7c19064dc000378b442fc3d563a82439aa516a90497`), run from the fix branch (`main` at
`385e6dc0`, which already includes the CMC centre-row repair, plus this change; vendor `d6052d60`; Python
3.12.15, numpy 2.5.4, scipy 1.18.1). The legacy rule is swapped in by monkeypatching the helper; nothing
else differs. Default initial data, r_max = 40, dt = 0.06·dr, to t = 0.60. Output `compare.json`, SHA-256
`75d9c443a8c457a9eadaeb611d22b192e1de7af1cf6a3f3806ef92b91d3d16bb`.

Cells 0–2 max (5 ≤ r ≤ 10 max in brackets), legacy → repaired:

| | N = 60 | N = 120 | N = 240 |
|---|---|---|---|
| H_total at t = 0 | 4.03e-3 → 7.2e-6 (5.0e-6 → 5.3e-7) | 4.05e-3 → 5.1e-7 (1.2e-6 → 3.1e-8) | 4.06e-3 → 7.7e-8 (3.3e-7 → 2.1e-9) |
| M_total at t = 0.6 | 9.5e-4 → 8.8e-7 | 6.7e-4 → 5.0e-8 | 9.5e-4 → 5.4e-8 |
| H_total at t = 0.6 | 5.0e-4 → 8.2e-7 | 9.4e-3 → 1.2e-6 | 5.9e-3 → 6.3e-7 |
| C at t = 0.6 | 7.5e-4 → 1.3e-6 | 4.5e-4 → 6.1e-8 | 7.7e-5 → 2.4e-8 |

At N = 480 the repaired initial centre H is 2.0e-8 (legacy 4.06e-3). The initial data also becomes fourth
order away from the origin (annulus 5.3e-7 → 3.1e-8 → 2.1e-9 → 1.3e-10), where the legacy rule was second
order.

Not settled: the repaired centre momentum constraint at t = 0.6 is 5.0e-8 at N = 120 and 5.4e-8 at N = 240,
so it stops falling there; that floor was not investigated (time-stepping error is the first candidate, as
dt halves with dr here).

## Relation to earlier work

The October origin half-cell audits (`tools/audit_current_main_centre_proper_radial_ricci.py`,
`tools/audit_analytic_origin_half_cell_self_consistent_response.py` and related) already compared an
even-series half-cell integral against the production rectangle. They study the historical recurrence and
copy it verbatim, so in the repair PR their three tests run with the legacy rule through an explicit
`legacy_origin_recurrence` fixture rather than being rewritten.
