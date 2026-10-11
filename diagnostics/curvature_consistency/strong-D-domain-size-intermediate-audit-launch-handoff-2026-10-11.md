# Strong-D intermediate-domain launch handoff

**Snapshot time:** 2026-10-11 05:19 UTC  
**Branch:** `diag/independent-metric-ricci-audit`  
**Purpose:** follow-up to `strong-D-domain-size-residual-audit-findings-2026-10-11.md`.

## Scientific reason for the follow-up

The completed ((N,r_{max})=(160,80)) and ((320,160)) profiles show a separate outward-moving maximum in the off-centre Hamiltonian-residual profile, but its physical meaning is unresolved. At fixed (\Delta r=0.5), the late outer-peak amplitude differs by about 8%, and the lapse profile on the shared radial interval differs substantially. Add ((N,r_{max})=(240,120)) to test whether the domain trend interpolates smoothly.

## Follow-up code and configuration

Diagnostic script: [strong_D_domain_size_intermediate_audit.py](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/independent-metric-ricci-audit/diagnostics/curvature_consistency/strong_D_domain_size_intermediate_audit.py)  
Script commit: `6c0eb68cdf342964ef5476cf01bd06c741d7ae0a`.

Two sequential trajectories are specified: untouched production initializer, then the regular-(F) cubic/Gauss shadow. Both use (N=240), (r_{max}=120), (\Delta r=0.5), CFL (=0.0075), (D=10^{-4}), amplitude (0.01), width (7), radiation enabled, final coordinate time (12), and full (H(r,t)) / lapse profiles every (0.25). Production physics, defaults, projection, and gauge are unchanged. This is a domain-size test, not spatial-resolution convergence.

## Live workflow state at snapshot

Both launches passed checkout, Python setup, compilation, fixed-spacing/sample/CFL validation, and dependency installation. Both had entered the numerical-evolution step; neither had published an artifact at this snapshot.

1. **First launch — duplicate to cancel if possible:** [workflow #38114304847](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38114304847), head commit `197900d98fe8ed3230f84b85deef88e03d9af76a`, job ID `114396025506`. This came from the newly added intermediate workflow.
2. **Second launch — stable registered route:** [workflow #38114367091](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38114367091), head commit `8b3bddf08c5a88e5abd02d65b0f67966e147d156`, job ID `114396208949`. This is the same intermediate script routed through the already-registered domain-audit workflow.

The second launch was triggered while repairing workflow routing because the first launch appeared delayed. Once both runs appeared, it was clear they were executing the same expensive two-trajectory campaign. No further campaign should be launched. The available GitHub connector exposes run reads, logs, artifacts and reruns, but does not expose a cancellation action; I could not stop the duplicate programmatically. To avoid duplicated computation, cancel **#38114304847** in GitHub Actions if you have access to its Cancel workflow control, and retain #38114367091.

## Next action after a terminal result

Inspect the artifact from the retained run, verify both trajectories reach (t=12), then compare (\alpha(r,t)), the inner and outer (|H|) peak radii/amplitudes, (\alpha_0(t)), and geometry health against both completed domain sizes at matched actual times. Do not infer a physical front, bounce, or turnaround from the residual ridge. Do not patch physics from this diagnostic alone.
