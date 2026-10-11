# Strong-D constraint audit — pinned artifact

Status: diagnostic-only; not admitted into the Step 1 record.

- Workflow run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37845980166
- Run ID: 37845980166 (run 56), conclusion: success
- Run head / source commit: aae3ccd2d2d5145253e0f692ebfebd9d67030150
- Script: engine/constraint_resolution_audit.py (blob 4cd4611a18e028149dd92839f07d3b9008e0efdb)
- Workflow: .github/workflows/zero_star_constraint_resolution.yml
- Physics/kernel code was not changed by this pinning commit.
- Configuration: D_amplitude=1e-4, radiation enabled, final coordinate time t=24; nominal step is 0.03 * grid.dr.
- Archived cases set both resolution=N and r_max=N for N=160, 240, 320, 480.
- Original Actions artifact: 0star-constraint-resolution, artifact ID 11581237503.
- Artifact ZIP SHA-256: 208b6c8010e38e17588833f2ceed36e71b37d2ab6658e6680d44b90332e67dd0
- Extracted report SHA-256: 3d2f00b73559338dad34463d60d291aba7bffcac10d3dd9f1eaf1f82804edc57
- Extracted report Git blob SHA-1: 768ed837a14263b0a1d8a01b25115af881eb7f08

## Qualification before interpretation

The run varies resolution and radial domain together (resolution=N and r_max=N); it therefore does not isolate spatial-resolution convergence. Treat it as a mixed resolution/domain sweep until the grid-spacing/domain design is corrected or explicitly controlled in a follow-up. The trapping-related fields are preserved as raw diagnostics only and do not authorize bounce, trapped-surface, or singularity claims.

This report is a separate artifact. Do not merge its interpretation into Step 1 admission or change physics based on this run alone.
