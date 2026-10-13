# Historical Constraint-Onset Audit: Exact Archived t=12 Campaign

**Date:** 2026-10-13 UTC  
**Scope:** Reanalysis of an already-completed historical campaign; no new evolution, production edit, or parameter change.  
**Workflow:** [Run 38208253687](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208253687)  
**Artifact:** [11728068723 — v55-production-t12-n40-n60-n80](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38208253687/artifacts/11728068723)  
**Artifact ZIP SHA-256:** `dbbaa3be94242cdeecd437d0a41c00d4cf3f2d625e1f3aaaad9e1b45833f2047`  
**Historical source commit:** `6d8b3f2790644591326df57c55908863f994cb0c`  
**Pinned vendor submodule:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`

## Purpose

PR #66's late-window audit showed that large central H/M/connection constraint residuals coexist with the outer radiation-cone failure at (t\approx47.15), but its accepted-state capture began only at (t=47.13). This audit reuses the older complete production campaign checkpoints and per-step ledgers to answer an important chronology question: **were central constraint problems already present early, and do they appear to be the same spatial problem as the much later outer radiation failure?**

This is explicitly a historical-source analysis, not a claim about an exact current-main campaign. The archive's production source commit differs from current `main`; the vendor pin is the same. The recomputation is matched to the archived source through its original per-step ledger and the independently reconstructed (t=12) values below.

## Integrity and recomputation validation

The exact artifact ZIP SHA-256 matches the archived status record. The run contains:
- N=40, 200 accepted steps, 201 diagnostic-history records;
- N=60, 300 accepted steps, 301 diagnostic-history records;
- N=80, 400 accepted steps, 401 diagnostic-history records;
- 18 full-field checkpoints across (t\approx2,4,6,8,10,12).

I recomputed the production Hamiltonian constraint (H), momentum constraint (M), conformal-connection residual (C_\Lambda), and determinant witness (ab^2-1) from all 18 full-field checkpoints using the pinned vendor's cell-centred fourth-order derivative/Ricci operators and reconstructed scalar plus conservative-fluid stress projections.

**Cross-check:** the recomputed (t=12) max norms match the original campaign's `campaign_summary.json` at all three resolutions. The maximum difference across every compared H/M/(C_\Lambda)/determinant metric was (1.39\times10^{-17}). That matches the output diagnostics at floating-point precision and validates the copied operator/projection calculation on the archived state.

The portable offline verifier is [`tools/audit_historical_t12_constraint_history.py`](https://github.com/miahpenn/0-Cosmos-engine/blob/diag/historical-t12-constraint-onset-audit-2026-10-13/tools/audit_historical_t12_constraint_history.py). It fails closed if the ZIP digest is wrong, checks checkpoint shapes and required files, produces signed checkpoint/annulus summaries, and asserts that the recomputed terminal norms agree with the original campaign record. It does not import or execute the engine.

## Finding 1 — the Hamiltonian defect exists at initialization and worsens

At (t=0), the full-domain Hamiltonian witness is already about (4\times10^{-3}) in the first cell across all three resolutions. In contrast, the initial momentum and connection residuals are at roundoff, and the determinant error is also at roundoff.

| Resolution | Initial (|H|_\max), first cell | (t=12) (|H|_\max), first cell | (t=12) (|M|_\max), first cell | (t=12) (|C_\Lambda|_\max), first cell |
|---:|---:|---:|---:|---:|
| N=40, Δr=2 | (3.7692\times10^{-3}) | (4.0417\times10^{-2}) | (3.3249\times10^{-3}) | (2.0684\times10^{-2}) |
| N=60, Δr=4/3 | (3.9316\times10^{-3}) | (7.9698\times10^{-2}) | (3.0020\times10^{-3}) | (2.5164\times10^{-2}) |
| N=80, Δr=1 | (3.9886\times10^{-3}) | (1.3494\times10^{-1}) | (8.0814\times10^{-3}) | (3.0858\times10^{-2}) |

The first-cell (H) is not monotonic: it initially decreases before growing later. But its (t=12) magnitude increases substantially under refinement, rather than approaching a common value. No convergence order is fitted to these centre-cell samples.

## Finding 2 — the momentum and connection constraints depart from roundoff in the first ordinary step

The per-step ledger records the first accepted step at (t=0.06,0.04,0.03) for N=40/60/80 respectively. By that first step:
- full-domain (|M|_\max\approx1.56\times10^{-4}) across all three resolutions;
- full-domain (|C_\Lambda|_\max\approx9.9\times10^{-6}), (6.4\times10^{-6}), and (4.8\times10^{-6}) for N=40/60/80;
- the determinant witness remains small, (O(10^{-11})) or below.

This is consistent with a nonzero initial semi-discrete constraint tendency, not an error first appearing at the (t\approx47) radiation endpoint. The ledger alone does not identify which geometric/source term fails to compensate the matter tangent.

The first discrete step where full-domain witnesses pass illustrative thresholds is:

| Resolution | (|M|_\max>10^{-3}) | (|C_\Lambda|_\max>10^{-3}) | (|H|_\max>10^{-2}) |
|---:|---:|---:|---:|
| N=40 | (t=0.42) | (t=0.66) | (t=4.92) |
| N=60 | (t=0.32) | (t=0.60) | (t=3.48) |
| N=80 | (t=0.24) | (t=0.54) | (t=1.59) |

These thresholds are descriptive markers for the saved time series, not scientific-admission tolerances or physical event definitions.

## Finding 3 — the early constraint pattern is central, not at the outer numerical boundary

By (t=12), the full-domain peaks sit at the first cell, but the (r\ge5) production witness remains nonzero and peaks in the (r=5)–6 region:

| Resolution | H exterior max (radius) | M exterior max (radius) | (C_\Lambda) exterior max (radius) | Max determinant residual |
|---:|---:|---:|---:|---:|
| N=40 | (4.7781\times10^{-3}) (r=5) | (6.4924\times10^{-4}) (r=5) | (7.1661\times10^{-3}) (r=5) | (2.7571\times10^{-9}) |
| N=60 | (5.0775\times10^{-5}) (r=6) | (6.6747\times10^{-5}) (r=6) | (1.3309\times10^{-3}) (r=6) | (4.3519\times10^{-10}) |
| N=80 | (6.9634\times10^{-5}) (r=5.5) | (9.8804\times10^{-5}) (r=5.5) | (1.8500\times10^{-4}) (r=5.5) | (7.7068\times10^{-10}) |

In the last physical band (75\le r<80), all three resolutions instead have (t=12) maxima of roughly (7.5\times10^{-8}) for (H), (7.3\text{–}7.7\times10^{-9}) for (M), and (9.5\times10^{-8}\text{–}1.05\times10^{-7}) for (C_\Lambda). This separates the early central/near-(r=5) constraint behavior from the far outer boundary.

At (t=12), the outer-cell radiation cone ratio is only 0.0288–0.0337 across N=40/60/80. The archived early campaign has no outer cone crossing; the later cone failure is a separate late-time issue, not the same instantaneous boundary state.

## Interpretation

1. There is an initial Hamiltonian residual concentrated in the first cell, already present at (t=0).
2. Momentum/connection constraints, initially roundoff-sized, acquire a visible tendency immediately during evolution.
3. Central residuals grow strongly by (t=12), with their magnitude depending substantially on resolution. The exterior witness near (r=5) is smaller and has mixed resolution behavior, while the far outermost band remains very small in this early campaign.
4. This does not prove a unique faulty equation, a continuum cause, or that the central constraints cause the later radiation-cone failure. It makes the central evolution/regularity/constraint-tangent gate an early and separate production blocker, rather than a late by-product of the (t\approx47) radiation event.

## Next discriminator

Use the current owners' signed centre-stage and fixed-slice CMC audits to choose one exact discrete identity to test on these already-saved early full-field checkpoints. Keep the historical and current-main source pins separate. Do not repeat a long campaign until the source/stage identity at fault is identified and an isolated counterfactual improves all implicated admission witnesses.

No production code, vendor code, physics coefficient, boundary condition, or campaign setting was changed; no evolution was performed by this audit.
