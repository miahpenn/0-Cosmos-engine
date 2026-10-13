# Native-cell Misner–Sharp mass-gradient identity — 13 October 2026

**Disposition:** bounded diagnostic only; no scientific admission.  
**Branch:** \`diag/current-main-native-ms-gradient-identity-2026-10-13\`  
**Exact-head workflow:** [38345293215](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38345293215)  
**Artifact:** [11788637675](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38345293215/artifacts/11788637675)

## Scope and provenance

The observer initializes the exact current-main production machine at each resolution and does not take any evolution steps. It evaluates the native-cell signed radial identity

\[
\mathcal E_M(r)=D_r M_{\rm MS}
-4\pi R^2\left(\rho_{\rm total}D_rR+j_{\rm total}\frac{\dot R}{\alpha}\right)
\]

at native cell centres bracketing areal radius \(R=10\), alongside the full sampled \(8\le R\le12\) annulus and the first five centre cells. It also records the same-state Hamiltonian/momentum/connection/determinant constraints and the central \(A_a\) RHS divided by \(r_0^2\). The production Misner–Sharp mass derivative uses odd radial parity. A separately coded centered fourth-order stencil is used only as a derivative-representation cross-check at interior cells.

- Current-main source commit: \`6df02e1f91dd4aa414bd760b7ee0d7d86f1605e3\` (the output-preflight driver fix was already merged).
- Diagnostic head: \`b86c08461da920f35f6eecd42c1b4988820c3af2\`.
- Production kernel blob: \`2954e7a79656849b5f088053ea0f4ac97d4bb443\`.
- Initializer blob: \`4cfbfe4be6ac3ce6ed492f4e435186997d94116a\`.
- Worldtube blob: \`aea19090cb0fd2ebe5481a8b7af0d26f119a9ad6\`.
- Stress-energy blob: \`c31c647f74fab6893254b3dea1366eada15b3110\`.
- V5.5 adapter blob: \`f671b838d9aefd8926702d9c2235209f52cd6815\`.
- Vendor submodule: \`d6052d605673ce9d82cdc99b0df79855fa2ea215\`.
- Workflow provenance gate confirms all these pins, no production engine/vendor/dependency diffs, and no state mutation by the observer.
- Full repository suite: **239 passed**. All five resolutions initialized at \(t=0\), with **zero time steps**.

Artifact ZIP SHA-256: \`d0fea1f808708f22f4a16b97df4aa926424a5c00286db6fbbb656837c009c6de\`. The embedded \`report.json\` SHA-256 is \`ada0b9ba62a1c3b1d6ccca5db3f5718bd618441d55e9edd827c4ed2b4bced4b1\`; both hashes were independently recomputed from the downloaded artifact and its manifest.

## Results

### Native radial mass-gradient identity

| Resolution \(N\) | Max \(|\mathcal E_M|\), \(8\le R\le12\) | RMS \(\mathcal E_M\), \(8\le R\le12\) | Signed \(\mathcal E_M\), lower R=10 bracket | Signed \(\mathcal E_M\), upper R=10 bracket |
|---:|---:|---:|---:|---:|
| 80 | \(2.95857\times10^{-5}\) | \(2.42320\times10^{-5}\) | \(+2.95857\times10^{-5}\) at \(r=9.5\) | \(+2.52128\times10^{-5}\) at \(r=10.5\) |
| 120 | \(1.21691\times10^{-5}\) | \(9.99461\times10^{-6}\) | \(+1.21691\times10^{-5}\) at \(r=9.66667\) | \(+1.11673\times10^{-5}\) at \(r=10.33333\) |
| 160 | \(6.64987\times10^{-6}\) | \(5.47061\times10^{-6}\) | \(+6.64987\times10^{-6}\) at \(r=9.75\) | \(+6.27637\times10^{-6}\) at \(r=10.25\) |
| 240 | \(2.89233\times10^{-6}\) | \(2.38374\times10^{-6}\) | \(+2.89074\times10^{-6}\) at \(r=9.83333\) | \(+2.79101\times10^{-6}\) at \(r=10.16667\) |
| 320 | \(1.62237\times10^{-6}\) | \(1.33153\times10^{-6}\) | \(+1.61201\times10^{-6}\) at \(r=9.875\) | \(+1.57159\times10^{-6}\) at \(r=10.125\) |

The maximum annulus identity error decreases smoothly. The pairwise convergence exponents from \(N=80\to120\to160\to240\to320\) are approximately **2.19, 2.10, 2.05, 2.01**, giving end-to-end order about **2.09** for this initialized-state family. Because these are native-cell measurements, this confirms that the gap is not solely caused by interpolating to an off-grid \(R=10\) tube. It remains compatible with a second-order spatial discretization/identity error and does not by itself demonstrate a continuum or physical conservation failure.

### Interior stencil cross-check

| \(N\) | Max production-vs-independent centered-4 \(D_rM_{\rm MS}\) difference | Operand-scaled roundoff bound | Result |
|---:|---:|---:|---|
| 80 | 0 | \(4.27\times10^{-13}\) | Pass |
| 120 | 0 | \(6.62\times10^{-13}\) | Pass |
| 160 | 0 | \(8.96\times10^{-13}\) | Pass |
| 240 | 0 | \(1.37\times10^{-12}\) | Pass |
| 320 | 0 | \(1.84\times10^{-12}\) | Pass |

The independently coded interior five-point derivative matches the production derivative exactly on these stored float64 arrays. This is useful confirmation of stencil representation, **not** an independent derivation of the stress-derived radial identity and not a conservation pass.

### Central constraints remain a separate problem

| \(N\) | Initial central \(H\) | Full-grid max \(|M|\) | Full-grid max \(|C_\Lambda|\) | \(A_{a,\rm RHS}/r_0^2\) |
|---:|---:|---:|---:|---:|
| 80 | \(+3.98859\times10^{-3}\) | \(3.60\times10^{-18}\) | 0 | \(+4.35652\times10^{-3}\) |
| 120 | \(+4.02811\times10^{-3}\) | \(9.40\times10^{-18}\) | 0 | \(+9.79154\times10^{-3}\) |
| 160 | \(+4.04156\times10^{-3}\) | \(4.66\times10^{-18}\) | 0 | \(+1.74021\times10^{-2}\) |
| 240 | \(+4.05101\times10^{-3}\) | \(2.00\times10^{-17}\) | 0 | \(+3.91474\times10^{-2}\) |
| 320 | \(+4.05429\times10^{-3}\) | \(7.58\times10^{-18}\) | 0 | \(+6.95911\times10^{-2}\) |

The initialized Hamiltonian peak remains near \(4.0\times10^{-3}\); momentum and connection are at roundoff in this state. At the same time, the central \(A_a\) RHS quotient grows by roughly 16× between \(N=80\) and \(N=320\). This is an independent warning that the centre-series/constraint issue is not resolved by the outer radial identity's smooth decay.

## Interpretation and next gate

The evidence supports three narrow conclusions:

1. The current-main native-cell radial identity gap exists at \(t=0\) and decays approximately as \(N^{-2}\) over the tested resolutions.
2. The production and independently coded interior derivative operators agree; that removes a straightforward interior-stencil representation mismatch from the current hypotheses but does not establish the physical source identity.
3. The centre-sector Hamiltonian and \(A_a\) RHS witnesses do not converge in this initial-state check and must remain separate from the radial mass-gradient gate.

No production/source/vendor code, physics, coefficients, boundary conditions, or campaign settings were changed. No evolution or long campaign was run. Do **not** substitute the numerical \(D_rM_{\rm MS}\) for the stress-derived term in the physical flux law, fit weights to the residual, or promote an initializer based on this spatial identity alone. The next useful gate is a source-pinned multi-time observer comparing signed native-cell radial identity, fixed-areal-radius temporal current/flux/work balance, and the central constraint witnesses at matched resolutions. Scientific admission remains false.
