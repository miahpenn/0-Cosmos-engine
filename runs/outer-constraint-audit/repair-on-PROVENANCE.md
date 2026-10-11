# Full-grid outer-constraint audit — repair ON provenance

Status: **admitted as a diagnostic artifact only**. This is not a physical-residual pass/fail result and does not change Step 1 admission.

- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37852226397
- Run ID: 37852226397; conclusion: success.
- Source commit: 92b383b8133ac747654c6d0cd6fb0656b35569bb.
- GitHub Actions artifact: [outer-constraint-audit-on, ID 11583388522](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37852226397/artifacts/11583388522)
- Artifact ZIP digest reported by GitHub and independently checked: `sha256:4e02f25548cf7ac6ea12fe14b031ef2413540da255e8665ccdab0e1dd53079d8`.
- Extracted report `repair-on.json`: 857,428 bytes; independently calculated SHA-256 `5d2837ce52aa695df3655d6fc9ffbb50268597bb9a5ad5f72b191e0cbc25e945`.
- Configuration: D=1e-4; N=160; r_max=80; CFL=0.0075; radiation enabled; final coordinate time=22.5; initial record at t=0 and seven registered targets.
- Admission gate passed: completed, config matched, eight records within sample-time tolerance, 160 cells in order per sample, all profile values finite, and full-grid Hamiltonian decomposition error <=1e-12 at every sample (exactly zero in all eight records).
- The gate applies no pass/fail threshold to physical H or M residuals.

## Cross-check against Step 1 repair ON

The production kernel and Step 1 recorder blobs are byte-identical to those in the Step 1 reference commit `869a7451523a2d6490acc40cc805f3c752909a62`. The sample times match. At t=22.5, the maximum of |H| over cells 2+ remains at cell 36, r=18.25 in both records: Step 1 `4.638237604395251e-4`; outer audit `4.638237605335193e-4`. The global maximum is cell 0, r=0.25 in both; values are respectively `1.1020277092482012e-3` and `1.1020277117824726e-3`. The small run-to-run numerical differences are preserved, not thresholded or silently rounded away.

At final time the global maximum and the regional summaries are also retained in `repair-on-summary.json`. The complete 160-cell profiles for all eight samples remain available from the linked Actions artifact, whose digest and extracted report hash are fixed above.

## Limits / sequencing

This is the repair-ON half only. Do not interpret its physical residuals in isolation, do not infer bounce/trapped-surface/singularity/shell behavior from this audit, and do not change physics. Repair OFF may be triggered only after this ON report is pinned and checked.
