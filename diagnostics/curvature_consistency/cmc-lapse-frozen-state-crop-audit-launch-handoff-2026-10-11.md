# CMC Lapse Frozen-State Crop Audit — Launch Handoff

**Date:** 2026-10-11  
**Branch:** diag/independent-metric-ricci-audit  
**Purpose:** Separate domain dependence in the native elliptic CMC solve from differences accumulated by independently evolved domain runs.

## Completed prerequisite

Three-domain fixed-spacing diagnostic is complete and summarized in [strong-D-domain-size-three-domain-synthesis-2026-10-11.md](strong-D-domain-size-three-domain-synthesis-2026-10-11.md). All six trajectories at rmax=80,120,160, production and regular-F shadow initializers, reached t=12 without evolution failure. The outer residual feature varies smoothly with domain but cannot be promoted to a physical front. Lapse profiles differ materially between independently evolved domains, and central residuals remain initializer-dependent. No physics patch was supported.

## Active run

- Actions run: [#38117001920](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38117001920)
- Trigger commit: 06d23520b5dbe41443d914e2385e8a8752b4724a
- New diagnostic script: [cmc_lapse_frozen_state_crop_audit.py](cmc_lapse_frozen_state_crop_audit.py)
- Script introduction commit: e4210215a337ef161355c3d9326bd7d9904be760
- Latest observed state at handoff: job in progress; compile/configuration validation passed, dependency installation passed, numerical evolution active, artifact not yet uploaded.
- Artifact expected name: cmc-lapse-frozen-state-crop-audit.

No second campaign should be started. The run is a single trajectory, not another domain sweep.

## Fixed run configuration

One production evolution at N=320, rmax=160, dr=0.5; CFL=0.0075; D=1e-4; scalar amplitude=0.01, width=7; radiation enabled; target t=12. Full state is saved every 0.25 time units.

At each saved slice, the diagnostic:
1. Solves native CMC lapse on the full state.
2. Crops the same state to (N,rmax)=(160,80),(240,120),(320,160), preserving all shared interior geometry/scalar/conserved-matter arrays.
3. Solves the cropped native equation, using the cropped domain's native full-volume target and alpha(R)=1.
4. Solves counterfactuals using the full-domain target Kdot_CMC: first with cropped alpha(R)=1, then with the full-domain solved alpha at the cropped boundary.
5. Records coefficient arrays, per-profile difference statistics, state arrays, and the exact discrete residual of the assembled operator.

The full state is essential: previous domain sweeps recorded only H_vendor and alpha, which do not suffice to reconstruct the lapse operator/source or separate effects.

## How to analyze the artifact

First verify the run's conclusion is success, report.json status is completed, accepted_steps=3200, sample_count=49, final t=12, and cycle-event/handoff counts are zero. Verify the artifact SHA-256 against GitHub's published artifact digest before interpretation.

At each sample time and crop compare:
- Native crop vs full-resolved lapse: full domain dependence with identical interior input fields, including crop-induced change of target and boundary row/derivative closure.
- Crop solve with full-domain Kdot_CMC vs full-resolved lapse: removes the change in the global CMC target, retaining the crop operator and alpha(R)=1.
- Boundary-matched vs native-boundary-value counterfactual (both with fixed full-domain Kdot_CMC): isolates the Dirichlet-value effect for a fixed crop operator.
- Native crop target Kdot_CMC vs full-domain target: establishes how much the proper-volume projection shifts because the integration domain and sampled RHS differ.

Check every solver's stored discrete operator residual. The rmax=160 crop should reproduce the full-domain solution to numerical roundoff because it is the same grid/state/operator/boundary. A discrepancy there points first to the audit harness/assembly, not to new physics. Evaluate both shared-interior RMS and maximum differences; report the radius of each maximum.

Finally compare frozen-state crop results with the separately evolved rmax=80,120 lapse profiles. Those old-run differences include evolution-history changes and must not be mislabeled as pure elliptic boundary effects.

## Guardrails

- Do not alter production physics, defaults, gauge, projection, lapse equation, evolution boundary conditions, or evolution algorithm.
- Do not call the outer residual ridge a D-front, causal propagation, bounce, turnaround, or cycle.
- Do not patch the lapse until a specific defective implementation layer is isolated with a verified operator residual and matched-state comparisons.
- Do not start more runs while this single audit is active or before its artifact is checked.

## Related source-level finding

The native CMC lapse solver is in engine/cmc_gauge.py: its final cell enforces alpha(R)=1, while Kdot_CMC is computed as a proper-volume weighted mean of the current K RHS over the whole domain. This gives two immediate global paths for domain sensitivity: the outer Dirichlet normalization and the changed integration domain for the CMC projection. The evolution characteristic boundary can also generate different evolved states on different domains. This makes the frozen-state test a necessary diagnostic before causal attribution.

See [cmc-lapse-domain-sensitivity-source-audit-2026-10-11.md](cmc-lapse-domain-sensitivity-source-audit-2026-10-11.md).
