# Exact-stage Lambda current and connection-constraint propagation crosswalk — 2026-10-25

**Status:** research diagnostic only. This is a stage-provenance/geometry-constraint witness, not the full matter-coupled momentum-propagation identity or a scientific-admission decision.

**PR:** [#420](https://github.com/miahpenn/0-Cosmos-engine/pull/420)  
**Run:** [40158222918](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40158222918), all N=60/120/240 jobs passed.  
**Source anchor:** \`4aec88d0e3ff028e0a65a8ac36b30948625df792\`; **vendor:** \`d6052d605673ce9d82cdc99b0df79855fa2ea215\`.  
**Inputs:** SHA-verified accepted \(t=0.60\) baseline/partial-joint and full-profile artifacts. Each arm is probed with one copied kernel step only; the saved accepted state is unchanged and the probe is not continued.

## 1. The exact stage currents are now captured

The diagnostic retains the total current \(j\) at three concrete stages: the step-local CMC-refreshed start, explicit predictor, and primary predictor. On each of these exact stages it compares the canonical \`total_matter_projection.j\` with the independent total projection and compares the matter part of the actual \(\Lambda_{L3}\) adapter with the direct stage expression
\[
L3_{\mathrm{matter},\Lambda}=-8\pi\lambda_m\,\alpha j/a,\qquad \lambda_m=2.
\]

The independent current projections agree to roundoff: the RMS difference over \(2\le r\le10\) stays between \(0\) and \(4.5\times10^{-22}\) across all arms/resolutions. The actual adapter matter source and the direct expression agree at every captured stage; the largest recorded absolute difference over \(r\le10\) is \(1.36\times10^{-20}\). These checks confirm that the source and \(j\) arrays are paired to the correct stage in this implementation.

The start and primary-predictor values must not be conflated. At cell 0 for the production baseline:

| N | start \(j_0\) | primary-predictor \(j_0\) | start \(L3_{\rm matter,\Lambda,0}\) | predictor source |
|---:|---:|---:|---:|---:|
| 60 | \(-1.32287\times10^{-6}\) | \(-1.34981\times10^{-6}\) | \(+6.64087\times10^{-5}\) | \(+6.76516\times10^{-5}\) |
| 120 | \(-6.72374\times10^{-7}\) | \(-6.79508\times10^{-7}\) | \(+3.37569\times10^{-5}\) | \(+3.40853\times10^{-5}\) |
| 240 | \(-3.37466\times10^{-7}\) | \(-3.39298\times10^{-7}\) | \(+1.69414\times10^{-5}\) | \(+1.70260\times10^{-5}\) |

The cell-0 current changes by about 2.04%, 1.06%, and 0.54% from start to primary predictor as N doubles. This matches the fact that the step size halves with refinement. Replacing the start/predictor pair with a lone endpoint current would erase a real stage distinction.

## 2. The implemented Lambda update itself pairs to the captured stage budget

For all three arms at all three resolutions, the accepted \(\Lambda\) rate and its exact PIRK stage budget agree to at most \(2.04\times10^{-17}\) absolute on \(r\le10\). The observer also verifies the identity-rate arithmetic on the exact same one-step probe. This is an implementation pairing check, not a physical conservation result: the stage source comparison is expected to reproduce the source expression already installed in the adapter.

## 3. The connection-constraint propagation witness separates the metric arms

For each saved accepted state the probe computes
\[
\dot C_\Lambda \approx
\frac{(\Lambda-\Lambda_{\rm metric})_{\rm probe}
      -(\Lambda-\Lambda_{\rm metric})_{\rm accepted}}{\Delta t},
\]
where \(\Lambda_{\rm metric}\) is independently reconstructed from the metric profile at each endpoint. This is the change in the connection-constraint residual over one copied production-kernel step, not a new equation or a retuned acceptance oracle.

The maximum \(|\dot C_\Lambda|\) in the central region \(r\le3\) is:

| N | Production baseline | Partial half-cell + regular CMC | Full continuous-profile + regular CMC |
|---:|---:|---:|---:|
| 60 | \(2.4839\times10^{-3}\) | \(1.7581\times10^{-4}\) | \(1.3108\times10^{-5}\) |
| 120 | \(2.2517\times10^{-3}\) | \(3.6518\times10^{-4}\) | \(9.2070\times10^{-7}\) |
| 240 | \(1.1310\times10^{-3}\) | \(5.7666\times10^{-4}\) | \(6.8357\times10^{-8}\) |

For the full-profile arm, the central maximum drops by factors of 14.24 and 13.47 over successive halvings of \(\Delta r\), corresponding to observed orders about 3.83 and 3.75. The fixed-annulus RMS (\(2\le r\le10\)) falls \(6.2860\times10^{-6}\to4.0278\times10^{-7}\to3.1517\times10^{-8}\), also approximately fourth-order over this bounded sequence.

By contrast, the partial arm's central maximum grows from \(1.76\times10^{-4}\) to \(5.77\times10^{-4}\). At N=240, the full-profile central maximum is about **8,400 times smaller** than the partial arm and about **16,500 times smaller** than the baseline. The cellwise signs are retained in the machine-readable profile artifact; the comparison is not based on absolute peaks alone.

## 4. Interpretation and limits

This is a stronger mechanism discriminator than the previous endpoint-only current check:

- It directly archives start and primary-predictor \(j\) and \(\Lambda_{L3}\) source arrays from the stages used by the real kernel step.
- It confirms that the implemented matter-source insertion is paired consistently with the stage current, and that canonical/independent stage current projections agree to roundoff.
- On the same one-step probe, the full continuous metric-profile shadow makes the central \(\dot C_\Lambda\) witness converge rapidly; the partial first-half-cell correction does not.

It does **not** prove the complete matter-coupled momentum-constraint propagation identity. In particular, it does not resolve the source-aware normalized Hamiltonian-gradient term, reference-metric connection-response coupling convention, the scalar spatial-operator balance, or total stress-energy sector conservation/exchange. An adapter reproducing its own source formula is implementation provenance, not an independent physical source-law derivation.

## Provenance / guardrails

- Workflow: [run 40158222918](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40158222918).
- Artifact ZIP digests: N60 \`c5a067131937d3607311218f607579778b83e08b7e026d9caea794392a28ea44\`; N120 \`e3fe876888788419830f41d56324431b662249d32d0d039ce63205d39c507079\`; N240 \`510351f235adcab607e559a2ad59da822bebe7ad9b38026eccba588eb47fce67\`.
- Source/vendor pins pass; the protected engine/vendor/dependency diff is empty; accepted input-state digests are unchanged during each probe; all captured profiles are finite.
- No production initializer/equation, engine, vendor or dependency file changed. \`full_matter_coupled_propagation_identity_closed=false\`, \`total_stress_energy_sector_conservation_closed=false\`, \`scientific_admission=false\`.

## Disposition

Keep the full continuous radial metric profile as a useful geometric shadow and the partial correction as a negative control. The next physics gate remains the independent, convention-pinned full matter-coupled propagation and sector conservation identity—not a production edit to the already stage-consistent \(\Lambda_{L3}\) source.


## 5. Component attribution of the central connection-constraint rate

The v2 rerun explicitly records the start and final-primary Lambda-L2 arrays, plus the start/primary Lambda-L3 vacuum and matter components. This removes the need to infer the final-primary Lambda-L2 term by rearranging the total rate. The decomposition is
\[
\dot C_\Lambda =
\underbrace{\tfrac12(L2_{\rm start}+L2_{\rm final})}_{\text{averaged }L2}
+\underbrace{\tfrac12(L3_{{\rm vac,start}}+L3_{{\rm vac,primary}})}_{\text{vacuum }L3}
+\underbrace{\tfrac12(L3_{{\rm matter,start}}+L3_{{\rm matter,primary}})}_{\text{matter }L3}
-\dot\Lambda_{\rm metric}.
\]

The four components close to the directly measured \(\dot C_\Lambda\) within \(2.16\times10^{-17}\) absolute on \(r\le10\). At cell 0:

| N | Arm | averaged \(L2\) | vacuum \(L3\) | matter \(L3\) | \(-\dot\Lambda_{\rm metric}\) | resulting \(\dot C_\Lambda\) |
|---:|---|---:|---:|---:|---:|---:|
| 60 | baseline | \(-3.701e-4\) | \(+1.322e-6\) | \(+6.703e-5\) | \(+2.786e-3\) | \(+2.484e-3\) |
| 60 | partial joint | \(+1.542e-6\) | \(-6.210e-8\) | \(+6.676e-5\) | \(-2.440e-4\) | \(-1.758e-4\) |
| 60 | full profile | \(+1.501e-6\) | \(-6.915e-9\) | \(+6.677e-5\) | \(-5.787e-5\) | \(+1.039e-5\) |
| 120 | baseline | \(-7.578e-4\) | \(-1.282e-7\) | \(+3.392e-5\) | \(-1.528e-3\) | \(-2.252e-3\) |
| 120 | partial joint | \(+1.633e-7\) | \(+1.989e-8\) | \(+3.378e-5\) | \(+3.312e-4\) | \(+3.652e-4\) |
| 120 | full profile | \(+1.977e-7\) | \(-1.055e-9\) | \(+3.378e-5\) | \(-3.339e-5\) | \(+5.862e-7\) |
| 240 | baseline | \(-1.526e-3\) | \(-2.203e-9\) | \(+1.698e-5\) | \(+3.776e-4\) | \(-1.131e-3\) |
| 240 | partial joint | \(+2.146e-8\) | \(-8.405e-9\) | \(+1.691e-5\) | \(-5.936e-4\) | \(-5.767e-4\) |
| 240 | full profile | \(+2.064e-8\) | \(-1.346e-10\) | \(+1.691e-5\) | \(-1.692e-5\) | \(+1.492e-8\) |

**Attribution:** the matter-\(L3\) source is nearly unchanged between partial and full-profile arms at each resolution. The decisive changed term is \(-\dot\Lambda_{\rm metric}\): in the full-profile arm it nearly cancels the matter-\(L3\) contribution, with the small averaged-\(L2\) and vacuum-\(L3\) terms left over as the residual. In the partial arm the metric-identity rate has the wrong magnitude (and at N=120 the opposite sign), so it overwhelms the same source and leaves a large \(\dot C_\Lambda\).

This is strong evidence that the full-profile metric shadow repairs the **kinematic Lambda/metric-identity pairing** at the origin in the one-step probe; it is not evidence that the physical source/exchange law is fully derived or the momentum constraint propagates. The baseline also shows a large independent averaged-\(L2\) contribution, so all four terms must remain in the ledger.

The exact-head v2 run is [40158404689](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40158404689); all three resolution jobs passed. Updated artifact hashes and component arrays are recorded in the machine-readable summary. Guardrails remain unchanged.
