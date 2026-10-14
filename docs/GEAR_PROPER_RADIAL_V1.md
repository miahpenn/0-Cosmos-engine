# GEAR Proper-Radial Observer v1 — representation-aware diagnostic

**Status:** research-only post-processing utility; no evolution or physics equation changed  
**Source branch:** \`research/gear-proper-radial-observer-v1-2026-10-14\`  
**Archive anchor:** \`GEAR_proper_radial_space_experiment_v1_2026-10-03.md\`.

## Why representation matters

The archived radial engine defines proper distance by \(d\ell=a(r)\,dr\), because its \(a\) is the radial line-element factor. The current BSSN representation instead declares \(\gamma_{rr}=a/X^2\), so its proper-distance increment is \(d\ell=\sqrt{a}/X\,dr\). These are different variable conventions; reusing the same array by name would produce a wrong coordinate map.

This observer exposes both operations separately:

- \`archive_proper_radius(r, radial_scale_factor)\` for the historical radial-engine convention \(d\ell=a\,dr\).
- \`bssn_proper_radius(r, metric_a, X)\` for current BSSN data \(d\ell=\sqrt{a}/X\,dr\).
- \`profile_on_proper_grid(ell, profile, target_ell)\` for bounded diagnostic interpolation; it refuses extrapolation.

## Quadrature and contracts

For cell-centred radii \(r_i>0\), the first centre-to-origin segment uses the midpoint rule \(\ell_0=r_0 f_0\), where \(f\) is the relevant radial line factor. Subsequent centre intervals use trapezoidal integration:

\[
\ell_i=\ell_{i-1}+\frac{f_{i-1}+f_i}{2}(r_i-r_{i-1}).
\]

The implementation checks dimensions, finiteness, positive metric factors, strictly increasing areal radii, and a positive monotone resulting proper-radius map. Synthetic tests cover flat-space identity, constant-factor scaling, a nonconstant quadrature oracle, representation non-equivalence, profile remapping, and fail-closed invalid inputs.

Run:

\`\`\`bash
pytest -q tests/test_gear_proper_radial_v1.py
\`\`\`

## Interpretation boundary

This is a coordinate observer over an already-solved spatial slice. It does not replace areal radius \(r\), change the evolution grid, move flux faces, modify the metric equations, or add a new degree of freedom. A remapped profile is a diagnostic view of saved data, not another evolution state.

The tests validate the quadrature contract on synthetic inputs only. The next physical-data step—when separately authorized—is to apply both representation-appropriate maps to the *same pinned saved geometry*, compare with the archived proper-radius checkpoints, and report the radial quadrature/order and source hash. No new campaign or production change is part of this PR.
