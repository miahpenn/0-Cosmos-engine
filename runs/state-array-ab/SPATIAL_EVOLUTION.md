# Repair ON/OFF state-array spatial evolution

Diagnostic-only post-processing of the immutable state-array A/B artifact. No simulation was run; no state-difference or physical-significance threshold was applied.

## Provenance
- Source run: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37859756564
- Artifact: https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37859756564/artifacts/11588005336
- Code commit: `7b9bac434f402df597641cfe7b24cf1a126a02e1`
- Artifact ZIP SHA-256: `349c239c8f0fd394831e8dee2bb77e922d771e92a5ccc80599d480920dc64abf`
- Raw JSON SHA-256: `4e53d8c62a369ebd5517a72307966730a49fcfab243bf39b9849d410fcc33b07`
- Reproducible script: [engine/state_array_spatial_evolution.py](../../blob/0star-central-clock/engine/state_array_spatial_evolution.py)

Method: absolute OFF-minus-ON difference; regional maximum and unweighted RMS for each selected numeric state array. Four bands partition the grid: r<20, 20<=r<40, 40<=r<64, r>=64. The ON/OFF states already differ at t=0, so this is not an identical-initial-state experiment.

## Regional maximum absolute difference versus time

Values are maxima in each radial band, not physical residuals or a significance measure.

### geometry.Lambda
| t | r<20 | 20-40 | 40-64 | r>=64 |
|---:|---:|---:|---:|---:|
| 0 | 3.5085e-4 | 0 | 0 | 0 |
| 0.01125 | 3.5085e-4 | 2.9816e-19 | 4.9806e-19 | 3.0832e-19 |
| 4.00125 | 3.5086e-4 | 2.2659e-15 | 3.7432e-15 | 5.0497e-15 |
| 7.99875 | 3.5082e-4 | 3.9750e-13 | 1.1216e-14 | 6.6380e-15 |
| 12 | 3.5152e-4 | 9.7238e-13 | 3.1301e-14 | 2.4093e-14 |
| 16.00125 | 3.5153e-4 | 8.1439e-12 | 4.6282e-14 | 5.6086e-14 |
| 19.99875 | 3.5153e-4 | 9.6606e-12 | 3.3426e-13 | 1.1199e-13 |
| 22.5 | 3.5153e-4 | 1.0118e-11 | 8.0903e-13 | 2.3091e-13 |

### geometry.Aa
| t | r<20 | 20-40 | 40-64 | r>=64 |
|---:|---:|---:|---:|---:|
| 0 | 5.1653e-15 | 0 | 0 | 0 |
| 0.01125 | 1.1605e-5 | 1.1543e-16 | 2.3990e-17 | 8.8471e-17 |
| 4.00125 | 9.2971e-5 | 1.0540e-12 | 1.3166e-13 | 3.6659e-14 |
| 7.99875 | 1.3346e-4 | 9.5885e-12 | 1.1932e-12 | 2.9373e-13 |
| 12 | 3.3537e-5 | 4.6272e-12 | 7.5709e-13 | 1.8700e-13 |
| 16.00125 | 7.8236e-6 | 3.1529e-12 | 1.0796e-13 | 4.9018e-14 |
| 19.99875 | 2.2013e-5 | 6.7550e-12 | 3.1070e-13 | 4.8348e-14 |
| 22.5 | 3.6485e-5 | 4.4855e-12 | 9.9012e-13 | 1.3930e-13 |

### geometry.alpha
| t | r<20 | 20-40 | 40-64 | r>=64 |
|---:|---:|---:|---:|---:|
| 0 | 1.1490e-9 | 1.3984e-12 | 5.1159e-13 | 1.4033e-13 |
| 0.01125 | 1.1207e-9 | 1.7468e-12 | 5.9508e-13 | 1.5388e-13 |
| 4.00125 | 2.1264e-8 | 6.2942e-12 | 6.2134e-12 | 3.3941e-12 |
| 7.99875 | 1.5845e-7 | 5.5360e-10 | 8.6952e-11 | 1.6590e-11 |
| 12 | 2.5822e-8 | 2.9436e-10 | 3.0814e-11 | 2.8748e-11 |
| 16.00125 | 4.0595e-9 | 1.2678e-10 | 2.8055e-11 | 2.3784e-11 |
| 19.99875 | 7.8095e-10 | 1.0215e-10 | 2.7224e-11 | 4.2959e-12 |
| 22.5 | 3.6749e-10 | 1.1062e-10 | 4.8908e-11 | 6.1651e-12 |

### scalars.D
| t | r<20 | 20-40 | 40-64 | r>=64 |
|---:|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 | 0 |
| 0.01125 | 2.2014e-14 | 4.3344e-22 | 3.2741e-33 | 5.6077e-56 |
| 4.00125 | 2.9114e-9 | 3.7896e-16 | 6.7670e-27 | 2.0825e-47 |
| 7.99875 | 1.7736e-8 | 3.1144e-13 | 2.0680e-22 | 1.5132e-42 |
| 12 | 5.8917e-9 | 4.6730e-12 | 7.0950e-20 | 1.9429e-38 |
| 16.00125 | 9.3459e-9 | 1.1969e-11 | 2.6826e-17 | 8.1017e-34 |
| 19.99875 | 9.9988e-9 | 9.2317e-12 | 6.1595e-15 | 1.2019e-29 |
| 22.5 | 1.0141e-8 | 8.2347e-12 | 1.3336e-13 | 5.5361e-27 |

The reproducible script also calculates these regional max/RMS profiles for geometry.a, geometry.b, and matter.dark_matter.energy_t.

## Interpretation boundary
- The largest differences in Lambda, Aa, a, b, and D remain strongly concentrated in r<20 across all recorded checkpoints.
- Some noncentral differences increase in absolute magnitude over time. For Lambda, the 20-40 band grows from zero at t=0 to 1.0118e-11 at t=22.5, while the inner maximum remains about 3.5153e-4.
- Alpha differs outside the inner region already at t=0 and changes non-monotonically in time.
- Therefore, differences are not literally confined to the center, and some noncentral amplitudes evolve. This does not establish physically meaningful outward spreading, causal propagation, or a propagation speed. No threshold was defined.
- The audit does not identify a mechanism, isolate an R-only effect, certify constraint satisfaction or convergence, or support bounce/trapped-surface claims.

At t=22.5, tau_on=6.983326299808656, tau_off=6.983325414935567, and tau_off minus tau_on=-8.848730894683854e-7.

Reproduction: extract state-array-ab.json from the source artifact, then run `python engine/state_array_spatial_evolution.py state-array-ab.json spatial-evolution.json`. This is post-processing only; it does not run the engine or modify physics.
