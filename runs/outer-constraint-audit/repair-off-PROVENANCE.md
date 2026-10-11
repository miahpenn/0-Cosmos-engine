# Repair-OFF outer-constraint audit provenance

- Diagnostic-only workflow: [run 37855818725](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37855818725)
- Trigger/source commit: \`5049c7cd001b046bc37d727a763c3ef8f3f40b73\`; parent is ON-summary pin \`5c743f6b707ff085dad7d797767d6376dfc55f55\`.
- Artifact: \`outer-constraint-audit-off\`, ID \`11584329648\`.
- GitHub artifact SHA-256 and independently checked local ZIP SHA-256: \`f8a9ed2933c3254d0200f649ae82bc87858bcf43b58f5fc91bc32521c5cc737d\`.
- Extracted raw report: \`repair-off.json\`, SHA-256 \`d1617a1e2002d6c711b14eee2dd1ec357a49e5bfad85ad3614266e113cce8657\`.
- Admission passed: completed without numerical failure; frozen config matched; eight registered samples within 0.05 time tolerance; all 160 radial cells present at each sample; all profile values finite; full-grid Hamiltonian decomposition error exactly zero at all eight samples.
- Frozen configuration: D=1e-4, N=160, r_max=80, CFL=0.0075, radiation enabled, final time 22.5, targets 0, 0.01, 4, 8, 12, 16, 20, 22.5. Regions: r<=20, 20<r<64, r>=64.
- No physical Hamiltonian or momentum residual pass/fail threshold was applied. No physics code was changed.
- Physics source blobs match the Step 1 reference \`869a745\`: \`engine/production_kernel.py\` = \`d9a4740592f2521f7dc04d66a718a854c70a8cde\`; \`engine/step1_hamiltonian_terms.py\` = \`30254271cb24b621a1831711e985a6d6b4e5e7d9\`.
- ON counterpart: run 37852226397, artifact ID 11583388522, ZIP SHA-256 \`4e02f25548cf7ac6ea12fe14b031ef2413540da255e8665ccdab0e1dd53079d8\`, report SHA-256 \`5d2837ce52aa695df3655d6fc9ffbb50268597bb9a5ad5f72b191e0cbc25e945\`.
- See \`repair-on-off-comparison-summary.json\` for exact final-time center-cell H and M residuals, center-cell curvature, and full-grid/three-region H/M comparisons. At t=22.5 the center-cell momentum residual is opposite in sign and 38.56 times larger in magnitude with repair OFF. Middle and outer H/M regional maxima are numerically close, not identical.
- The comparison summary does not infer that all effects are confined to the center: it records the full-grid witnesses and complete radial profiles remain in the immutable Actions artifacts.
