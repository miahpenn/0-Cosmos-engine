# Brown matter-source bridge and absolute propagation residual — 2026-10-25

**Status:** research diagnostic only. This addendum refines the interpretation of the v7 accepted-state crosswalk. It does not change a production equation or set scientific admission.

**Source run:** [40159027760](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40159027760), all N=60/120/240 jobs passed.  
**PR:** [#420](https://github.com/miahpenn/0-Cosmos-engine/pull/420).  
**Observer commit:** ddc560964067c0e2d1569e401ff6cc92fe969adf.  
**Input accepted state:** t=0.60; source 4aec88d0e3ff028e0a65a8ac36b30948625df792; vendor d6052d605673ce9d82cdc99b0df79855fa2ea215.

## 1. The S-sector source-map residual is balanced by Brown's matter-dependent terms

Define:

- \(B_{\rm tot}\): the published Brown \(M_r\) source evaluated on matter-inclusive H/M constraints;
- \(B_{\rm geom}\): the same Brown construction using the geometric-only constraints;
- \(S_{\rm act}\): the existing geometric matter-injection plus current-tangent composite;
- \(S_Q=-8\pi\alpha Q_r\): the provisional stress-divergence source;
- \(\Delta S=S_{\rm act}-S_Q\): the source-map residual;
- \(D_{\rm tot}=\dot M_{\rm measured}-(B_{\rm tot}+S_Q)\);
- \(D_{\rm geom}=\dot M_{\rm measured}-(B_{\rm geom}+S_{\rm act})\).

The same accepted-state arrays satisfy the signed bridge
\[
B_{\rm tot}-B_{\rm geom}\simeq \Delta S,
\qquad
(B_{\rm tot}-B_{\rm geom})-\Delta S
= D_{\rm geom}-D_{\rm tot}.
\]

The first relation explains why the source-map residual can be large while the full matter-inclusive Brown-plus-stress-divergence candidate is small: the Brown matter terms are part of that source accounting. They cannot be dropped and then interpreted as a missing physical source.

For the full-continuous-profile arm, the maximum absolute gap between the Brown source difference and the source-map residual decreases across resolution:

| N | max abs of \((B_{\rm tot}-B_{\rm geom})-\Delta S\) |
|---:|---:|
| 60 | 3.625e-10 |
| 120 | 3.009e-11 |
| 240 | 2.015e-12 |

The equivalent signed bookkeeping relation \((B_{\rm tot}-B_{\rm geom})-\Delta S-(D_{\rm geom}-D_{\rm tot})=0\) closes to maximum absolute 5.08e-21, 7.46e-21, and 1.06e-20 for N=60/120/240 respectively. This is an algebraic reconciliation of the current Brown decompositions on the saved tangent arrays, not an independent physical-source proof.

At fixed r=2, the source-map residual and Brown source difference are:

| N | \(\Delta S\) | \(B_{\rm tot}-B_{\rm geom}\) |
|---:|---:|---:|
| 60 | +5.8090e-5 | +5.8090e-5 |
| 120 | +5.7798e-5 | +5.7798e-5 |
| 240 | +5.8351e-5 | +5.8351e-5 |

At fixed r=5, they are likewise +6.1904e-5/+6.1904e-5, +6.2535e-5/+6.2535e-5, and +6.2464e-5/+6.2464e-5. Thus the persistent source-map residual is almost exactly represented by the difference between Brown evaluated on total versus geometric-only constraints. Reading the source-map residual in isolation—as though it must vanish before the total Brown formula can close—was incomplete.

## 2. The full matter-inclusive Brown candidate has a convergent absolute residual

The direct signed residual is \(D_{\rm tot}=\dot M_{\rm measured}-(B_{\rm tot}+S_Q)\), evaluated on the same directional tangent and pointwise divergence source. Its absolute values for the full-profile arm are:

| Witness | N=60 | N=120 | N=240 |
|---|---:|---:|---:|
| max absolute, all cells | 8.975e-6 | 9.317e-7 | 4.133e-8 |
| RMS, all cells | 1.522e-6 | 1.407e-7 | 7.375e-9 |
| signed residual at r=1 | +2.310e-6 | +4.368e-7 | +2.947e-8 |
| signed residual at r=2 | +8.975e-6 | +6.122e-7 | +4.030e-8 |
| signed residual at r=5 | +9.004e-7 | +1.778e-8 | +1.045e-9 |
| signed residual at r=8 | -1.194e-6 | -7.328e-8 | -4.591e-9 |

At r=2 the absolute residual decreases by factors 14.66 and 15.19 under successive halvings of dr, corresponding to observed orders 3.87 and 3.93. At r=8 it decreases by factors 16.30 and 15.96, orders 4.03 and 4.00. The global maximum also falls from 8.98e-6 to 9.32e-7 to 4.13e-8, although the radius of the maximum shifts as the grid is refined.

This is a substantive correction to the earlier reading based on the candidate's normalized relative RMS: a high relative RMS is not, by itself, evidence that the signed absolute residual is failing to converge when the reference derivative becomes very small. The correct next metric is the signed absolute profile and its fixed-physical-radius order, while retaining relative RMS as a separate conditioning diagnostic.

For comparison, the partial-joint arm has global max absolute residual 1.127e-4, 8.287e-4, 3.927e-3; the production baseline has 5.088e-4, 7.550e-3, 3.249e-2. The full-profile arm is the only one of these three with a clearly decreasing global absolute residual across all three resolutions.

## 3. What this does and does not close

This crosswalk strengthens the numerical evidence for the **full matter-inclusive Brown \(M_r\) propagation candidate** on the accepted \(t=0.60\) full-profile shadow. It also shows why the sector source-map residual cannot be treated alone as a failed conservation identity: the Brown matter-dependent terms balance it to the same order as the direct candidate residual.

Still, this is not yet a scientific-admission result. The exact convention/prefactors and regular-origin limits must stay tied to the archive; the total stress-energy sector conservation ledger is a distinct gate; the weak COSMOS-phi exchange comparator still requires its archived normalization/tolerance review; and the near-centre finite-radius evidence does not alone prove the complete analytical identity under every configured sector and gauge condition. Keep the partial first-half-cell arm as a negative control and the full continuous profile as a research shadow.

## Provenance and guardrails

Input artifacts from workflow 40159027760 were SHA-verified by the observer:
- N=60 artifact 12584544612, SHA-256 045cd07a28a1e5d55e88650fb06df5771796999f3887d8c3cc49b050715b7645
- N=120 artifact 12585144988, SHA-256 69c42bb0e072449d085c63c45add2f9d75d8c886a9a09f2b3707537b1500ced4
- N=240 artifact 12585018461, SHA-256 d817c8e9fb63e7d753fa3d688a8195a79b160b0ce1df628f95ec26df7351ebd3

Protected engine/vendor/dependency paths are unchanged. No production initializer/equation was modified. Full scientific admission remains false until the complete archive-backed propagation and sector-conservation gates are explicitly closed.
