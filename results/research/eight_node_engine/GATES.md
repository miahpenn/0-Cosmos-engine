# Eight-node areal integral into the engine: predeclared gates

Written 2026-11-24 before any run on this branch. Go: Jeremiah, "Eight-node first" (Lumen² channel, 01:19 UTC Nov 24), relayed by the coordinator.
Patch: Aster, B7-eight-node-areal.patch (Drive 1YHyR8_ckKVqo4gAlzPl73EXKu99Srfzk), applied to main 1e7bffc with conflicts resolved by keeping both sides (supplied_slice and eight_node_areal). Opt-in flag, default off.
Review: Lumen², REVIEW-B7-W1-origin-repair-eight-node (Drive 1cN0IOfLl_Z6tv3a1GBBAFNcQ3P5Lt4CjDPTI28AOHxY), item 1e asks for (f) and (g).
No gate below changes after numbers are seen. Corrections only as dated amendments; post hoc readings labelled post hoc.

| Gate | Pass condition | Why it can fail (discrimination) |
|---|---|---|
| (a) suites | `pytest tests` and root `pytest` (vendor submodule, sympy) show no failure that is absent on main 1e7bffc | The new tests assert exact polynomial integration on every stencil row; a mutated stencil fails them (checked by mutation before the run). |
| (b) default-off baselines | Flag off: every array of `build_initial_data` for none, rescaled_pocket (pocket_x 0.5), inherited_pocket, hand_floor_reference at n 1200 and n 2400, r_max 60, is bit-identical to main. Prediction: zero movement. Flag on: magnitudes of the change in B, a, X reported, not gated. | Flag on must give a nonzero change, else the comparison cannot see an integral change. |
| (c) G4c reproduction | Aster's full.py on native.npz, with her research rule replaced by the engine's `_eight_node_origin_integral`: official ring/PC2 ratio at dr 0.025 <= 10 and within 0.05 of 0.970. Also: engine integral vs Aster's `integrate(r,f,h,8)` on the W1 sources, max abs difference <= 1e-12. | The same script with the engine's default integral is run as a negative control; dr 0.05 gives 62.50x, outside the band. |
| (d) dr 0.05 failure | Recorded, not gated. Expected about 62.50x (FAIL on record). | n/a |
| (e) direct curvature | Nine- and seven-node direct areal ratios recorded only. Never a gate. | n/a |
| (f) order, non-polynomial | Integrand 1 + exp(-r^2), exact r + (sqrt(pi)/2) erf(r), r_max 4: observed order from dr 0.2 to 0.1 >= 6, and max error at dr 0.1 <= 1e-7. | The existing cumulative rule is measured in the same test and must show order < 6. |
| (g) outer-edge noise gain | Per-row noise gain sum_j abs(dI_i/df_j)/dr, from impulse responses at n 32: max over all rows <= 10. Derived prediction: about 4.67 at the last row (eight-point Adams-Moulton weights), about 1 in the interior. Recorded per row. | An extrapolating last-interval stencil (Adams-Bashforth eight-point, gain about 81.9) would fail. |

Pinned engine blobs for the allowlist are reported in the PR once the code is final. Our sessions do not edit the allowlist (workflow, audit tool, its test); Jeremiah makes that commit.
