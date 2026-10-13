# Uniform-quadrature CMC centre-row one-step comparison

**Date:** 13 October 2026 UTC  
**Diagnostic run:** [38327602102](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38327602102) — SUCCESS  
**Diagnostic head:** `b7ad7172a6078655b656d1bb7ce8dadd4d1a06ce`  
**Production source:** `6d8b3f2790644591326df57c55908863f994cb0c`  
**Pinned vendor:** `d6052d605673ce9d82cdc99b0df79855fa2ea215`  
**Artifact:** [11781305984](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38327602102/artifacts/11781305984)  
**Artifact ZIP SHA-256:** `95484ab33fc5bffe27cb5c251dc0201132e8fa030df586979c9f950a921a1225`  
**Report JSON SHA-256:** `a53e4ab8bfec5b138763b1de96cd1040f6a67afa6aaf7a76152010dc43e57eee`

## Scope

This diagnostic reruns the frozen uniform-fourth-order-quadrature initializer at N=40/60/80/120/160. It injects one of three centre-row treatments into the diagnostic kernel instance's lapse solver and then executes exactly one ordinary `V55ProductionKernel.step()` at `dt = 0.03*dr`. The chosen row is used for the initializer and all four internal CMC solves (g0, explicit predictor, primary predictor and completed final slice).

The three alternatives are:

1. Current row: `alpha[0] = alpha[1]`.
2. Second-order even-reflection ghost PDE row.
3. Direct fourth-order even-parity PDE row at cell 0; interior CMC PDE rows and outer `alpha[-1]=1` remain unchanged.

Production source, vendor, physical equations/sources, matter definitions, outer boundary reconstruction, step ordering and stage validators are byte-identical to the pinned reference. No damping, projection, clipping, retries, physical thresholds or parameter tuning were introduced.

## Exact one-step results

Every arm completed the one ordinary step and passed its finite/positive geometry check: **15/15 arms passed**. Each arm invoked the selected lapse solver five times (initializer plus four actual production step solves). No initial/final diagnostic field was non-finite.

| N | Centre row | Initial `Aa_rhs/r₀²` | Max `|K_rhs−target|` over rows 0..N−2 | Final `H` full abs max | Final `M` full abs max | Final `C_lambda` abs max | Final radiation `|S|/E` outer |
|---:|---|---:|---:|---:|---:|---:|---:|
| 40 | Current | +1.4570e−4 | 1.4404e−3 | 3.8469e−4 | 1.5298e−5 | 9.32e−7 | 2.5500e−4 |
| 40 | Even ghost | +1.0512e−4 | 8.4817e−5 | 3.8366e−4 | 2.0684e−5 | 1.27e−6 | 2.5507e−4 |
| 40 | Direct D4 | +1.0738e−4 | 5.1100e−5 | 3.8372e−4 | 2.0384e−5 | 1.25e−6 | 2.5507e−4 |
| 60 | Current | +1.1664e−4 | 1.5197e−3 | 1.0089e−4 | 5.8501e−6 | 2.39e−7 | 1.6860e−4 |
| 60 | Even ghost | +2.2780e−5 | 4.1978e−5 | 9.9847e−5 | 3.3671e−6 | 1.38e−7 | 1.6862e−4 |
| 60 | Direct D4 | +2.5303e−5 | 3.3961e−5 | 9.9875e−5 | 3.2195e−6 | 1.32e−7 | 1.6862e−4 |
| 80 | Current | +1.5347e−4 | 1.5501e−3 | 3.7756e−5 | 6.1229e−6 | 1.87e−7 | 1.2593e−4 |
| 80 | Even ghost | −1.4692e−5 | 2.4533e−5 | 3.6701e−5 | 4.1910e−7 | 1.30e−8 | 1.2593e−4 |
| 80 | Direct D4 | −1.2072e−5 | 2.1991e−5 | 3.6717e−5 | 3.3275e−7 | 1.03e−8 | 1.2593e−4 |
| 120 | Current | +3.3563e−4 | 1.5728e−3 | 1.0191e−5 | 6.1851e−6 | 1.26e−7 | 8.3603e−5 |
| 120 | Even ghost | −4.4750e−5 | 1.1208e−5 | 1.0561e−5 | 3.3403e−7 | 6.75e−9 | 8.3604e−5 |
| 120 | Direct D4 | −4.2059e−5 | 1.0808e−5 | 1.0555e−5 | 3.7352e−7 | 7.56e−9 | 8.3604e−5 |
| 160 | Current | +6.2149e−4 | 1.5810e−3 | 5.9256e−6 | 6.1663e−6 | 9.44e−8 | 6.2572e−5 |
| 160 | Even ghost | −5.5983e−5 | 6.3658e−6 | 6.8759e−6 | 2.9649e−7 | 4.51e−9 | 6.2572e−5 |
| 160 | Direct D4 | −5.3266e−5 | 6.2973e−6 | 6.8720e−6 | 3.1893e−7 | 4.85e−9 | 6.2572e−5 |

At N=80–160, both alternatives remove the large positive centre `Aa_rhs/r₀²) witness and reduce the CMC K-target mismatch by roughly 1–2 orders of magnitude. The full-domain momentum and connection residuals also fall substantially, while full-domain Hamiltonian residuals stay in the same order of magnitude (and rise modestly in the alternatives at N=120/160). Radiation remains admissible for this single step in all three variants. The Misner–Sharp mass and `flux_T + work_pR` values are essentially identical across row choices at each N over this one-step interval.

The alternatives therefore repair important *discrete centre/gauge witnesses* on the frozen profile and in one real step, but do not show long-horizon constraint preservation, conservation or convergence. At N=40/60 some central quotient remains positive; resolution dependence is not fully closed.

## Interpretation and next gate

This evidence strengthens the diagnosis that imposing equal lapse samples at two different positive radii is a major source of the initial centre-series and K-target defect. It does **not** yet authorise either centre row for production. The direct-D4 centre residual is roundoff by construction; the even-ghost alternative avoids that particular tautology but still needs refinement and full operator checks.

Next: preserve these results and compare the two candidate centre rows against the same-resolution constraint/centre-series trajectory through several accepted steps (only if the bounded stage checks remain clean), including H/M/C_lambda, determinant, direct-D2 lapse-Hessian/Ricci split, radiation admissibility at all stages, proper time, and Misner–Sharp current balance. Keep the current row as the counterfactual control. Do not start the expensive long-horizon campaign or promote the row before these separate gates support it.
