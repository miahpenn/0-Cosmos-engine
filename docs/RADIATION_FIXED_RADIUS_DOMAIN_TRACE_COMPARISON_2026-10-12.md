# Fixed-Physical-Radius Radiation Comparison Across Existing Domain Traces

**Date:** 12 October 2026 (UTC)  
**Classification:** Offline comparison of existing artifacts only. No trajectory was launched and no production equation or boundary law was changed.  
**Diagnostic branch:** `diag/radiation-fixed-radius-domain-compare-2026-10-12`  
**Production tree reference:** `2035acc1210db18da0cededea18750b7ed813260` (production files unchanged).

## Finding

At a matched physical radius \(r=79.5\), the existing \(R_{\max}=120\) trace is substantially more radiation-admissible around central proper time \(\tau\approx21.913\) than three \(R_{\max}=80\) traces. In the expanded-domain trace, the near-null state is instead concentrated at the new outermost radius \(r=119.5\). This is evidence of a strong **domain-sensitive / outer-region effect**, and argues against claiming that the recorded event is already demonstrated to occur at a unique, fixed physical radius independent of the computational domain.

It does not prove that the outer boundary formula is the cause. The traces differ in source snapshot and/or recovery policy, as explicitly recorded below. They are a cross-run forensic comparison, not a controlled same-code, same-policy convergence experiment.

## Compared artifacts and provenance

1. **Expanded domain \(R_{\max}=120\)** — [run 38028306184](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38028306184), [artifact 11661332148](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38028306184/artifacts/11661332148). Configuration \(N=120\), \(\Delta r=1\), CFL \(=0.03\), radiation on, accepted-metric recovery off, adaptive admissibility retries on. Frozen production source commit `a7170c12ed08b4fdb3c52ca98f02e106a3f39b82`; `trace.json` SHA-256 `f418a023c44dcd51802696e6468e3e09691d01ea44ecbb5a9463881be6b14e76`; artifact ZIP SHA-256 `f4f06a54bf66f0e5eed597efe0b6b0e605e31375edd7303bcb4edd4cdb97fac5`.

2. **\(R_{\max}=80\), accepted-metric recovery control** — [run 38023932307](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38023932307), [artifact 11660061575](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38023932307/artifacts/11660061575). Configuration \(N=80\), \(\Delta r=1\), CFL \(=0.03\), recovery policy on, adaptive admissibility retries on. Source commit `1e50e24acd9d75a0c855acbcaa086b73fb60c63e`; `trace.json` SHA-256 `70740b6449613029258ce3cb0d9b4451f42b75f2c96a18e552c0ce7b211e74fe`; artifact ZIP SHA-256 `2cedb6279bf476405d8d0465a31c689af93f5d39f3415415081608f0aa5c6fba`.

3. **\(R_{\max}=80\), rate-corrected trace** — [run 38310237842](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38310237842), [artifact 11775230036](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38310237842/artifacts/11775230036). Configuration \(N=80\), \(\Delta r=1\), CFL \(=0.03\), accepted-metric recovery off, no retries. Frozen source commit `1c63090d4f42396f899421a2379f70320689d5ad`; `trace.json` SHA-256 `53a1f9f188e0eea4575da4458d09891f50144b9f89f945c995cd9980221fb619`; artifact ZIP SHA-256 `f64298e2819e0dd9c4c2df712a4d1250e4855a5b0e467430fe7aa419a533e49f`.

4. **\(R_{\max}=80\), stacked RHS/rest-current source** — [run 38321456272](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38321456272), [artifact 11779635426](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38321456272/artifacts/11779635426). Configuration \(N=80\), \(\Delta r=1\), CFL \(=0.03\), accepted-metric recovery off, retries off. Frozen source commit `cdcec8282bff0443f2fe4674b15929822ff38736`; `trace.json` SHA-256 `06138c2b31db6d6bd2f7807a4cdb655d535eb48ec392c22aaae897ff855391ce`; artifact ZIP SHA-256 `448f8184571a7d0d0e4603e53531ef3148f07c0b565acfa250e31a230a8377d0`. This follow-up records both predictor and completed-state cone-margin budgets with zero stored reconstruction-error counts and zero instrumentation errors.

The first two traces (R80 recovery control and R120 domain control) use the **same SHA-256 values** for production `engine/production_kernel.py`, `engine/matter_system.py`, `engine/matter_rhs.py`, `engine/valencia.py`, `engine/v55_matter.py`, `engine/v55_pirk_adapter.py`, `engine/cmc_gauge.py`, `engine/stress_energy.py`, and the pinned vendor boundary file. Their explicit control difference is the recovery policy, while their domain/resolution differ. The two newer R80 traces include later source changes and are used as corroborating observations, not exact matched controls.

## Same-radius, same-central-clock comparison

For each artifact, I selected the stored `accepted_start_CMC_lapse` snapshot nearest \(\tau=21.913\) that contains cell \(i=79\), \(r=79.5\). No interpolation or evolution was performed.

| Existing trace | \(t\) | \(\tau\) | \(C\) at \(r=79.5\) | \(|S|/E\) at \(r=79.5\) |
|---|---:|---:|---:|---:|
| \(R_{\max}=80\), accepted-metric recovery control | 46.620000000001006 | 21.910498546691237 | \(+3.4913097\times10^{-6}\) | 0.983917532 |
| \(R_{\max}=80\), rate-corrected trace | 46.620000000001006 | 21.910482973554274 | \(+4.1883593\times10^{-6}\) | 0.981123974 |
| \(R_{\max}=80\), stacked RHS/rest-current trace | 46.68000000000101 | 21.91233698555736 | \(+4.3783676\times10^{-6}\) | 0.977667218 |
| \(R_{\max}=120\), expanded-domain trace | 61.04625000000155 | 21.912982335453805 | \(+1.0553626\times10^{-5}\) | 0.838596494 |

\(C=U_E-\sqrt{\gamma^{rr}}|U_r|\) is evaluated from each trace's stored conservative variables and its stored metric at that stage. All four selected states remain inside the radiation cone.

At the same R120 accepted snapshot, the outermost cell \(i=119\), \(r=119.5\), is already nearly null: \(C=+1.2077498\times10^{-8}\), \(|S|/E=0.999993703\). The R120 run then reaches a genuine near-null rejection at that outer cell at \(t=61.046503364247606\), \(\tau=21.9130694389674\), during the explicit predictor. The final tiny retry has ratio exactly \(1.0\) to stored precision. By contrast, at comparable proper time the \(R_{\max}=80\) traces show their \(r=79.5\) state still inside the cone with \(|S|/E\approx0.978\)–0.984; that physical-radius cell only approaches \(|S|/E\approx1\) later, near \(\tau\approx22.02\).

The effect is not a simple monotonic ratio profile across the full R120 outer region. At \(\tau=21.91298\), the R120 trace's stored ratios include approximately 0.7725 at \(r=89.5\), 0.0885 at \(r=99.5\), 0.8521 at \(r=109.5\), and 0.999994 at \(r=119.5\). So this is a coupled, spatially structured solution, not simply a scalar front translated to the new edge.

## Interpretation

**Supported:** in the existing artifacts, the near-null radiation event moves to the new outermost region when \(R_{\max}\) is expanded. At a common physical radius \(r=79.5\), the expanded-domain trace is much more admissible at comparable central proper time. This is consistent with a domain/global-CMC/outer-region influence in the current discretization and argues against treating the R80 event alone as a fixed-radius, domain-independent physical event.

**Not established:** this is not a strict domain-convergence experiment. The accepted-metric recovery policy differs between the closest shared-code pair, and the other R80 traces use newer matter/RHS code than the R120 run. Therefore the ratio difference cannot be attributed uniquely to the domain boundary, gauge, recovery method, or physical model. A fixed-radius comparison of current source would require a matched pair of domain runs, which this report deliberately does not launch.

This result complements the completed-boundary and primary-predictor \(A_a\) counterfactuals: neither single-call boundary-only suppression rescues the captured failure, and the domain comparison shows the location/severity is domain-sensitive. The next cause-level question is the joint spatial-domain/CMC/predictor-stage response, not a boundary-law rollback.

## Disposition

No production change is justified. Keep \(E\geq |S|\) exact. No boundary disable, floor, clamp, limiter, damping, fitted correction, or forced retry is proposed. Reuse these existing fixed-radius and stage-budget witnesses in the next equation-level review; do not launch another long run until a preregistered same-source domain question distinguishes a specific discrete identity or stage incompatibility.

No production main, pinned vendor, canonical GEAR archive, or physical equation was changed by this report.
