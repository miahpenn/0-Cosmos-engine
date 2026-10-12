# Flat Geometry RHS Audit — Corrected Vacuum State

## Current result: PASS

- Workflow: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38202616225
- Test: `diagnostics/flat_vacuum_geometry_smoke.py`
- Commit: `8a40a51351e8906be5a78fcd04e65203179891ec`
- Artifact: `flat-vacuum-geometry-smoke`, ID `11725659439`

The corrected state sets (a=b=X=\alpha=1), (\beta=B=A_a=K=\Lambda=0), (S=P_S=D=P_D=\Pi=0), and (\phi=\ln 10). With the production potential
`V(phi)=V0*(exp(-phi)-0.1)`, this makes (V(\phi)=0). The matter API still requires tiny positive fluid seeds (`rest=energy_t=1e-30`), which are explicitly recorded in the test.

Observed:
- maximum absolute geometry RHS: (7.161375604511246\times10^{-15})
- threshold: (10^{-12})
- result: PASS

## Correction to the earlier failed test

The previous version set every scalar field to zero, including the COSMOS field `phi=0`. That was **not vacuum**: the potential at `phi=0` is (0.9V_0), so the geometry was correctly sourced and the K RHS magnitude was (7.41827\times10^{-5}). That failure was a test-state mistake, not evidence of a defect in the geometry equations. The test now uses the true zero-potential point and passes to roundoff.

## Scope

This validates the production geometry RHS blocks at a zero-potential flat state. It does not test the full production time stepper. A separate workflow is now testing one actual step of `V55ProductionKernel` from the repository's homogeneous archive cosmology operating point. That model has a rolling exponential scalar potential, so it is not the same as constant-Λ LCDM.
