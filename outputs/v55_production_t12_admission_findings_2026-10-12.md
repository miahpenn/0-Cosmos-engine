# V5.5 bounded production campaign — admission findings and continuation

**Snapshot:** 12 October 2026  
**Repository:** `miahpenn/0-Cosmos-engine`  
**Diagnostic code commit:** `b3b50c0d087a9e01bbfe32808fce0742aa7dc148`  
**Branch:** `fix/strict-campaign-admission-2026-10-12`  
**Scope:** read-only interpretation of the completed N=40/60/80 production campaign; no physics modifications.

## 1. Exact run and artifact provenance

- Campaign workflow: [run 38207937002](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38207937002) — completed successfully.
- Strict admission smoke: [run 38207936911](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38207936911) — completed successfully.
- Campaign artifact: [artifact 11727843997](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38207937002/artifacts/11727843997), `v55-production-t12-n40-n60-n80`.
- Verified local ZIP SHA-256: `ffb65350bdc5af3f0826f46f8012110c3d5ec87e8ddb8412f64710b153e70550`.
- Configuration: resolutions `(40, 60, 80)`, `r_max=80`, CFL `0.03`, requested final coordinate time `t=12`, checkpoints every `2`.
- Repository tests passed in the production job. Every resolution reached `t=12`; the current artifact includes an explicit `t=0` row and 201, 301, and 401 history rows respectively (initial row plus 200, 300, and 400 accepted steps).

Do not confuse this artifact with the preceding campaign artifact (different commit, `a1e92d8...`). This note's values are from the exact `b3b50c0...` run and checksum above.

## 2. What passed — and what did not

The runner's process/admission record is `passed=true`, with no errors. This establishes that the requested resolutions are present, each trajectory reports `completed`, and each final coordinate time reaches the requested horizon within the configured floating-point tolerance. The new runner also checks state/diagnostic finiteness and positive geometry during evolution.

That is **execution admission, not physics/convergence admission**. `campaign_admission_errors()` currently gates completeness and final time; it does not reject a trajectory because Hamiltonian or momentum constraints grow or fail to converge. The numerical ledger therefore must be assessed separately.

### End-of-run summary

| Resolution | Accepted steps | Final proper time τ | Full-domain max |H| at t=0 | Full-domain max |H| at t=12 | Max |H| outside first two cells at t=0 | Same masked max at t=12 | Max full-domain |M| at t=12 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| N=40 | 200 | 11.9559993814 | 3.7692e-3 | 4.0417e-2 | 4.7327e-6 | 4.7781e-3 | 3.3249e-3 |
| N=60 | 300 | 11.9284259877 | 3.9316e-3 | 7.9698e-2 | 1.7195e-5 | 1.0721e-2 | 3.0020e-3 |
| N=80 | 400 | 11.9184114053 | 3.9886e-3 | 1.3494e-1 | 1.7741e-5 | 1.9625e-2 | 8.0814e-3 |

Definitions in `engine/production_kernel.py`:

- `H = vacuum.constraints(grid, geometry)["hamiltonian"] - 16*pi*total.rho`.
- `M = vacuum.constraints(grid, geometry)["momentum"] - 8*pi*total.j`.
- The full-domain maxima include the first two radial cells. The separate `hamiltonian_max` witness intentionally excludes those cells and reports its radius on that same masked domain.

The Hamiltonian residual is present on the initial slice near the center (about (3.8\times10^{-3}) to (4.0\times10^{-3})) and becomes substantially larger during evolution. More importantly, the off-centre/masked maximum grows from roughly (5\times10^{-6})–(2\times10^{-5}) at initialization to (4.8\times10^{-3})–(2.0\times10^{-2}) at (t=12). The full-domain maximum and the masked maximum both grow with resolution at the final time instead of displaying constraint convergence in this three-resolution campaign. This is not only a first-cell reporting issue.

Full-domain momentum begins at roundoff scale, then is (O(10^{-3})) by the final slice. By contrast, the determinant witness stays below (8\times10^{-10}) at (t=12), and the recorded last-interior CMC residual stays at approximately (10^{-8}) or below. Small determinant/CMC residuals do not close the Hamiltonian or momentum constraints.

Final τ values also differ across the resolution set (N=40 minus N=80 is about (0.03759)). They are useful observables, but a proper-time trend cannot substitute for constraint convergence.

## 3. Supported conclusion

**The bounded full-machine run executes cleanly to the target time, but its evolved constraints are not admitted as converged.** The initial slice already has a central Hamiltonian residual; the evolution then develops a much larger residual beyond the first two cells as well. Present evidence does not isolate whether the dominant growing contribution is geometric, stress-energy sourcing/projection, or a mismatch in the discrete operators/stage coupling.

Do not describe this as a completed production-validation pass. Equally, do not diagnose a particular physics equation from these scalar maxima alone. The production kernel's constraint diagnostic currently records the combined residual, not a matching full-profile budget of its geometric and matter pieces.

## 4. Next discriminating diagnostic — no physics edits

Use the *existing* checkpoint artifact as input; do not rerun the full evolution just to get another scalar maximum.

1. On the exact pinned kernel/vendor versions, reconstruct or independently evaluate the initial state and each saved `t=2,4,6,8,10,12` checkpoint for N=40/60/80.
2. At common physical radii, record full radial profiles of:
   - the raw geometric `vacuum.constraints(...)["hamiltonian"]` term;
   - the source term (16\pi\rho_{total}), plus sector contributions if readily available;
   - their signed difference (H);
   - the raw geometric momentum term, (8\pi j_{total}), and signed (M).
3. Keep full-domain and (r\geq 2\,dr) witnesses separate; compare at matched radii rather than matching grid indices. Preserve cell-centre arrays, pinned code/vendor identities, input artifact digest, and per-profile hashes.
4. Use the result to branch the investigation: if the raw geometric term is the growing piece, audit the metric/geometry-stage update and discrete curvature operators; if the source term is responsible, audit stress assembly and matter/scalar source projections; if both are individually large but their cancellation degrades, audit operator/source consistency. No floors, clamps, fitted terms, or changed evolution equations.
5. Add a separate scientific admission criterion only after the diagnostic defines the measured behavior and reference evidence. Do not invent or relax a tolerance simply to get a green workflow.

This should be an offline/checkpoint diagnostic with no new production trajectory. It is the highest-priority next step because it tells us which layer the full engine is actually failing in.

## 5. Parallel tracks — current cross-checks

### Strong-D ridge refinement

The fixed-domain N=160/N=320 comparison finished: [run 38165315602](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38165315602), artifact 11712817995, SHA-256 `d266c9aedf259222e720cd52067b97ab5cc44f7556f98177db5e00e656ee449f`. With (r_{max}=80), (D=10^{-4}), radiation ON and (t=12), late shell-window peak amplitudes agree closely (N=160 (2.829182\times10^{-4}) at (r=9.75), N=320 (2.827630\times10^{-4}) at (r=9.875); ratio 0.99945). The middle-time local peak identity switches, so this is promising late-time agreement, not a blanket convergence proof. The next conditional study in that track is a D/S/momentum/gradient source split once peak-branch identity has been checked.

### CMC global-target counterfactual

The full evolving R80 run with archived R160 (\dot K) target history closed 99.9941% of the final proper-time gap in that controlled setup. It is causal evidence of target-history sensitivity, not proof the imposed R160 target is uniquely physical or self-consistent on R80. The every-other-sample target-history sensitivity run [38207554072](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38207554072) was still reported in progress at the last check. Do not launch a duplicate; inspect its artifact when it completes.

### N=640 signed-pole path and canonical state

The N=640 signed-pole run [38183314088](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38183314088) is **cancelled**, with no admitted result/artifact. The existing `main` status note still says it is running and is stale on that point. Main was not changed in this audit.

## 6. Guardrails preserved

- No physics edits or promotion to `main` are justified by the current evidence.
- No artificial damping/floor/clamp, fitted feedback coefficient, manufactured bounce rule, or silent equation change.
- Distinguish workflow success from physical/campaign convergence.
- Keep the current CMC coarsening run single-instance; use the existing full-state artifacts where possible.
- After the source split, repair only the evidenced failing layer, then rerun the smallest admitted regression followed by the same bounded production campaign.

**Current disposition:** execution is successful; full-machine numerical constraint admission is open. The next action is a source-versus-geometry constraint budget on the saved trajectory, not a speculative physics change.
