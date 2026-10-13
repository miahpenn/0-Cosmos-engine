# Radiation Outer-Metric Kinematic Identity Audit

**Date:** 2026-10-13 UTC  
**Scope:** Offline observer/artifact audit only. No simulation, production-source edit, boundary change, or physics change.  
**Coordination:** [Issue #5 cross-agent board](https://github.com/miahpenn/0-Cosmos-engine/issues/5#issuecomment-6137134957)  
**Production source baseline:** `29d828b4d72b017c12a20b6dc90549f2511f025c`  
**Pinned geometry/PIRK vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`  
**Current main at audit:** `23fc8403365d041bbc9323b36d1df1b5a7cb6209` (documentation-only divergence from latest tested code/test tree)

## Question

The checked frozen-state radiation replay localizes the last inadmissible explicit predictor to the outer cell, while the radiation matter RHS improves the admissibility margin and the explicit metric change worsens it. Is the stored `a/X` metric-rate ledger inconsistent with the metric kinematics encoded by `K`, `Aa`, and the CMC lapse?

For zero shift, `gamma_rr = a / X^2` and `K^r_r = K/3 + Aa`. The ADM kinematics therefore require

```
2*Xdot/X - adot/a = 2*alpha*(K/3 + Aa).
```

With the replay's zero-shift, zero-KO explicit block, the pinned vendor source separately gives `adot = -2*alpha*a*Aa` and `Xdot = alpha*X*K/3`. The cone margin used by this audit is

```
C = U_E - (X/sqrt(a))*abs(U_r)
```

so its explicit metric contribution must satisfy

```
dotC_metric
  = -(X/sqrt(a))*abs(U_r)*(Xdot/X - adot/(2*a))
  = -(X/sqrt(a))*abs(U_r)*alpha*(K/3 + Aa).
```

## Pinned evidence and method

- Canonical diagnostic lane: [PR #61](https://github.com/miahpenn/0-Cosmos-engine/pull/61), exact-head frozen replay [run 38337821862](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38337821862), [artifact 11785239489](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38337821862/artifacts/11785239489).
- Downloaded ZIP SHA-256: `5c3dee4cbe62fb410ed7bcdf97707385019e00bae4f925ab51d4242b42933d09`. The embedded five-entry `SHA256SUMS.txt` manifest was recomputed and every entry matched.
- Capture schema: `radiation_explicit_predictor_stage_capture_v5_retry_identity_and_recomputable_independent_l3`.
- Independent script: [`tools/audit_radiation_metric_kinematic_identity.py`](../blob/diag/radiation-metric-kinematic-identity-2026-10-13/tools/audit_radiation_metric_kinematic_identity.py). It reads the existing artifact ZIP, verifies provenance, and recomputes the stated identities; it does not import or run the engine.
- Source equations were cross-read from [production adapter](https://github.com/miahpenn/0-Cosmos-engine/blob/main/engine/v55_pirk_adapter.py), [production step](https://github.com/miahpenn/0-Cosmos-engine/blob/main/engine/production_kernel.py), and the pinned vendor's [`vacuum_bssn.py`](https://github.com/RonBibb/bb-palatini-unified-r0/blob/d6052d605673ce9d82cdc99b0df79855fa2ea215/vacuum_bssn.py).

The audit checks every saved geometry stage and all cells; for cone-term reconstruction it joins only exact matching **ordinal-0** RHS and geometry stage keys, then verifies those source profiles are identical. This avoids mixing predictor states that intentionally differ before/after the primary `Aa/K` update.

## Results

| Witness | Result |
|---|---:|
| Geometry stages | 703 |
| Geometry cells checked | 56,240 |
| Maximum absolute residual of `2*Xdot/X - adot/a = 2*alpha*(K/3 + Aa)` | `2.220446049250313e-16` |
| Maximum relative residual, normalized by sum of absolute operands | `2.431903378e-16` |
| Maximum absolute residual at outer cell | `5.551115123125783e-17` |
| Max residual in captured `adot = -2 alpha a Aa` | `0.0` |
| Max residual in captured `bdot = alpha b Aa` | `0.0` |
| Max residual in captured `Xdot = alpha X K/3` | `2.775557562891351e-17` |
| Exact ordinal-0 RHS↔geometry input pairs checked | 651 |
| Profile mismatches among those pairs | 0 |
| Max absolute difference between direct and reduced `dotC_metric`, all cells and ordinal-0 stages | `6.776263578034403e-21` |
| Max relative difference between direct and reduced `dotC_metric` | `2.289278791e-16` |

### Outer cell 79: first rejected predictor

At input `t=47.130000000001026`, attempt 1, `dt=0.03`:

- `C = +2.806443824920446e-7`.
- Matter RHS contribution to `dotC`: `+4.6733363346889415e-5`.
- Metric contribution, direct from `adot/Xdot`: `-5.8743458862801356e-5`.
- Metric contribution, reduced with `alpha*(K/3+Aa)`: the exact same recorded value.
- Net instantaneous tendency: `dotC = -1.201009551591194e-5`.

### Outer cell 79: terminal machine-resolution attempt

At input `t=47.15300258518951`, attempt 651, `dt=6.36411169718682e-15`:

- `C = +5.421010862427522e-20`.
- Matter RHS contribution to `dotC`: `+4.49658805832266e-5`.
- Metric contribution, both forms: `-6.01861977594991e-5`.
- Net instantaneous tendency: `dotC = -1.5220317176272503e-5`.

The strict probe remains short of target `t=47.160000000001034`; the archived campaign's `admission_passed` remains `false`.

## Interpretation

**The captured metric-rate identity passes at floating-point roundoff.** This rules out a mismatch in which the stored `adot/Xdot` values fail to satisfy the kinematic relation they are supposed to obey with the same-stage `alpha/K/Aa` profiles. It also confirms that the direct and reduced metric terms used in the radiation cone budget are equivalent.

It does **not** establish that the accumulated `K/Aa` trajectory, CMC lapse response, matter/geometry stage coupling, outer characteristic boundary history, or continuum limit is correct. The metric contribution is a real, internally consistent term in the captured semi-discrete evolution; a bad `a/X` rate formula is not supported as the immediate cause. At the two representative endpoints, the matter contribution is positive but smaller than the negative metric contribution.

No captured state was modified; no cone threshold, boundary operation, coefficient, source law, or physics equation was changed. No evolution was run.

## Next bounded discriminator

Do not patch the `a/X` explicit kinematic rates on this evidence. Continue offline from these exact stage profiles: cross-check the captured `Aa/K` L2/L3 rates and CMC-lapse inputs against the already-coordinated regular-centre, momentum/Hamiltonian constraint, and boundary witnesses. A production change should wait for a signed discrete identity that fails, a minimal profile/state reproducing it, and a bounded counterfactual that improves all implicated admission witnesses.
