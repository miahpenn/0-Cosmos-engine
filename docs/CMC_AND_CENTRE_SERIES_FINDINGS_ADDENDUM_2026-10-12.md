# CMC and centre-series findings addendum — 2026-10-12

This addendum reconciles the latest three diagnostic lines. All are frozen-state / read-only diagnostics; no production physics was changed and no long trajectory was run.

## 1. Signed production RHS split

Run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38294736014  
Artifact ID: 11768515284  
SHA-256: \`2128028ad4fcd80c4d603107d481a78a1add233d18b35eaa3d1d19a43cae2a91\`  
Full report: \`docs/PRODUCTION_AA_K_RHS_CENTRE_SOURCE_SPLIT_2026-10-12.md\`

Seventeen frozen states were audited. Signed RHS reconstruction closes to roundoff: maximum base Aa/K reconstruction error 2.17e-19; momentum-current reconstruction error 4.07e-20; directional-tangent reconstruction discrepancy 1.09e-14. This validates accounting of the implemented semi-discrete RHS, not continuum conservation.

At baseline N=60/80/120, first-cell Aa_dot is approximately +1.09e-3 and is dominated by trace-free Ricci (~+1.05e-3); the trace-free lapse-Hessian term is ~3.6–4.1e-5. Uniform fourth-order quadrature largely suppresses the Ricci contribution, but the lapse-Hessian remains ~3.87–4.15e-5 at N=80–160. Aa_dot/r0^2 rises from 1.53e-4 at N=80 to 6.21e-4 at N=160. The centre-series gate remains open.

On that uniform arm, K_dot at the first cell is about -1.38 to -1.52e-3 versus about +1.0 to +1.45e-4 at the next cell. The lapse-Laplacian contribution changes sign strongly between those cells; the matter source changes much less. D_r K_dot/r at the first cell grows from 5.11e-4 (N=40) to 8.68e-3 (N=160). This is a witness, not proof that the CMC row is faulty.

## 2. Read-only CMC centre-row/operator audit

Run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38296556473  
Artifact ID: 11768278320  
Artifact SHA-256: \`e02c2c0f384d20692fc02bb3ee771409d0743a2479dc39276ee3a14b1156101f\`  
Report: \`docs/CMC_LAPSE_CENTRE_ROW_AUDIT_2026-10-12.md\`  
Diagnostic commit: \`95606bb4db9260286191aec4a513a0430aed386a\`

Two archived t=0 states were checked: R80 shared-core N=160 and native R160 N=320. The CMC elliptic rows close to 7.32e-15 max; independently reconstructed production K-RHS rows close to 6.04e-15 max. Alpha's even-polynomial fit over the first eight cell centres has max residual 4.49e-8. The production lapse-geometric K term and independent CMC operator agree to about 2e-15 at the first cells and fixed radii r≈2,3,4,5.

This weakens the hypothesis that the archived-state CMC centre row itself is faulty. It does not test the uniform-quadrature initial arm and is not a formal refinement-order proof. Do not patch the centre row based on the uniform-arm jump alone.

## 3. CMC target-closure rank probe

Run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38294955148  
Artifact SHA-256: \`64cc1666a2010d73b462cac48c78706ac48be1a1d02fea7de343bd09172c4e24\`  
Report: \`docs/CMC_TARGET_CLOSURE_RANK_AUDIT_2026-10-12.md\`

Nine fixed Kdot targets yield output slope 0.999999999757 and closure residual slope -2.42686e-10; H residual is about -8.58e-14 and interior PDE rows close to 3.71e-15. The outer normalization cell is not an elliptic PDE row but contributes to the target average. This suggests near-redundancy of the target closure on this frozen slice, not authority to fit an adjustment or remove the normalization cell.

## 4. Independent blocker and next action

Matched quadrature one-step runs at N=240 and N=320 have first-cell radiation realizability failures (|S|/E about 1.0203 and 1.3603). This is a separate admission blocker. Do not clip primitives or loosen the cone.

Next: compare the exact CMC lapse operator with the production lapse term on the *same uniform-fourth-order quadrature initial states* N=40/60/80/120/160. Record alpha, central row residual, even-series fit, -Delta alpha, trace-free lapse Hessian, trace-free Ricci, Aa_dot, K_dot and matter source. If operators agree, inspect initializer geometry/matter construction and regular-centre identities for Aa/K/X/a/b; if not, isolate the exact mismatch before any repair.

No centre projection, damping, floors/clamps, source/gauge/boundary changes, or long t=160 campaign until a discrete defect is demonstrated and all admission gates pass. Production main remains unchanged. The current evidence is a reproducible numerical lead, not production readiness.
