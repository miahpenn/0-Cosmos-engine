# Issue #4 — shared-state momentum / total-stress calculation

Date: 2026-10-19  
Status: shared-state calculation executed successfully; full propagation closure remains open.

## Implementation

- PR: [#408](https://github.com/miahpenn/0-Cosmos-engine/pull/408), draft, diagnostic-only.
- Shared-state adapter: [audit_shared_state_momentum_total_stress_divergence.py](https://github.com/miahpenn/0-Cosmos-engine/blob/diagnostic/accepted-state-raw-momentum-witness/tools/audit_shared_state_momentum_total_stress_divergence.py).
- Workflow: [accepted-state-raw-momentum-witness.yml](https://github.com/miahpenn/0-Cosmos-engine/blob/diagnostic/accepted-state-raw-momentum-witness/.github/workflows/accepted-state-raw-momentum-witness.yml).
- Production source pin: `7f47669a1c0edfafef9fc5f27b2db80840d3f035`.
- Vendor pin: `d6052d605673ce9d82cdc99b0df79855fa2ea215`.
- Initial digest: `a256b57f274ee120b55c0f805ed488d91382a9799647c8f9811d7df82aa337c0`.
- Accepted digest: `c6acb3711d59a85c804a88208c9bc0356e8227a806c0684ae3a8afac72c9091a`.

The adapter advances the pinned production kernel exactly one step at resolution 60, then constructs one RHS-direction plus/minus pair at epsilon `1e-5`. It passes those same in-memory state objects to the signed momentum-constraint observer and the independent total-stress-energy divergence observer. It emits signed cell-centred and finite-volume divergence arrays plus the signed momentum tangent, and records state digests before/after observation.

A wrong-digest negative control must be rejected. Source/vendor provenance and both expected state digests fail closed.

## Scope and limits

This is the first shared-state calculation, not yet the full matter-coupled momentum-propagation identity. It reports residuals; it does not assert that they close. The existing six-sector observer still excludes charged dust and radial Maxwell and therefore cannot be promoted to an all-sector oracle. Brown's existing propagation audit remains vacuum-only.

No production equations, source/vendor pins, numerical bounds, or scientific-admission rules changed. `scientific_admission=false` and `acceptance_oracle_satisfied=false` remain mandatory. PR #408 stays draft until the exact-head workflow runs, signed results and controls are inspected, and the complete independent propagation RHS is assembled.

## Executed result

- Exact-head workflow: [run #39193664244](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/39193664244) — job completed successfully, including the shared-state calculation and artifact upload.
- Artifact: `accepted-state-raw-momentum-witness`, ID `12161012410`, archive SHA-256 `e5fab5a274d8e455a66ae5b037c5441444e4d0351bedef724857e875a9eeece2`.
- Shared report SHA-256: `0c95c52bbb2623a624bc07607a9b4ca96f94a76772f0a088d5bf110f6dd0e0e8`.
- Exact state pairing passed: initial and accepted digests matched the pins; plus/minus states had distinct digests; accepted state digest was unchanged after observation; the wrong-digest negative control was rejected.

### First signed measurements

At the two centre cells:

| Signed quantity | Cell 0 | Cell 1 |
|---|---:|---:|
| Measured momentum-constraint tangent | +3.933697986e-3 | -7.461987336e-4 |
| Total stress divergence, radial component (cell-centred) | +4.009516398e-6 | +8.605236521e-7 |
| Total stress divergence, radial component (finite volume) | +3.589113039e-8 | +2.995064137e-8 |
| Total stress divergence, time component (cell-centred) | +8.740353868e-8 | +2.995356864e-8 |
| Total stress divergence, time component (finite volume) | +3.045726886e-8 | +3.245581884e-8 |

Global maxima/RMS from the same report: cell-centred `div_T_r` max absolute `4.009516398e-6`, RMS `5.297798764e-7`; finite-volume `div_T_r` max absolute `4.199445192e-8`, RMS `1.103229817e-8`.

The six-sector raw divergence sum does **not** reconstruct the total time-component divergence directly (max gap `2.614321629e-16)); after the existing measured time-density-component accounting, the unexplained gap is `3.068892472e-20` maximum across channels. This is an accounting result, not a proof that the full momentum-propagation RHS closes.

### Decision

The same-state calculation and its negative control passed. The measured momentum tangent and independently observed total-`T` divergence are now available from one exact state pair, but they are distinct diagnostics with different units/roles and are not expected to be equated directly. The independent normalized Hamiltonian-gradient plus connection-constraint propagation assembly, complete enabled-sector/exchange coverage, and manufactured signed controls remain the next gate. Both admission flags remain false; PR #408 stays draft.


## Increment 2 — Brown structural cross-check on the same accepted state

- Commit: `11bbf8455a3cf6dec7d80d7adb217147877b4c4a`.
- Exact-head workflow: [run #39194181923](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/39194181923) — **success**, all 15 substantive steps passed.
- Artifact ID: `12161615073`; archive SHA-256: `20f257bbdb0dc432c0f504ba1be4a901538627650a87df7ffc35ae35241baba3`.
- Shared report SHA-256: `005764de4fea3f8aa4e4ab667df867f211e2c00923d047aa496e43e81c334fe3`.

The adapter now computes the measured geometric momentum tangent separately from the total momentum tangent, evaluates the independent Brown Eq. 11b vacuum structural predictor on the same accepted geometry, records the signed difference, and requires a wrong-sign Hamiltonian-gradient control to be detected. The exact initial/accepted digests still match; plus/minus states differ; the accepted state remains unchanged; the digest-mismatch control is rejected.

### Signed centre-cell comparison

| Quantity | Cell 0 | Cell 1 |
|---|---:|---:|
| Measured geometric momentum tangent | +3.8576351315e-3 | -9.4250659294e-4 |
| Brown vacuum structural prediction | -1.4702616878e-4 | -1.7039909435e-4 |
| Measured minus Brown prediction | +4.0046613003e-3 | -7.7210749859e-4 |
| Brown Hamiltonian-gradient term | -1.4565991903e-4 | -1.7112785312e-4 |

Global max/RMS for the structural residual: max absolute `4.0046613003e-3`, RMS `5.2717981054e-4`. The wrong-sign H-gradient mutation was detected. Total-T divergence results and six-sector reconstruction remain as recorded above.

### Interpretation and next action

The vacuum structural terms do **not** account for the observed matter-coupled geometric tangent at this state: the signed residual is much larger than the Brown prediction. That is the useful result, not a pass/fail verdict on production physics. The remaining residual mixes explicit matter/source contributions with finite-difference product-rule/operator effects. Next, separate the pinned production RHS source terms from the geometric stencil/discretization remainder, preserving signed cell-0/cell-1 and annular profiles, and add a manufactured control that distinguishes the two. No tuning or admission-gate changes.

The exact replay and strict JSON residual contract are green. Lightweight validation and strict campaign/radiation retry integration were still in progress at the time this note was updated.
