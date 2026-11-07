# Full-evolution regular-origin initializer shadow — 25 October 2026

**Status:** successful diagnostic-only full evolutions to accepted \(t=0.60\); not production authorization and not a scientific-admission pass.  
**Draft PR:** [#416](https://github.com/miahpenn/0-Cosmos-engine/pull/416)  
**Workflow:** [run 40153734946](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40153734946)

## Question and design

Does the one-step improvement from the combined regular-origin shadow survive a fully coupled trajectory? Four arms were evolved for each \(N=60,120,240\), each to \(t=0.60\), with the same pinned production step/RHS/PIRK/boundary ordering: production baseline; copied half-cell integral correction only; copied regular-origin CMC central row only; joint half-cell integral correction plus CMC row.

The initializer candidate changes only the half-cell integral \(I(0\rightarrow r_0)\) to the continuous-reference integral over \([0,r_0]\). Later radial intervals retain the pinned trapezoid recurrence, the six source/B iterations are retained, and matter is reinitialized against the candidate geometry. The alternate CMC centre row is a copied diagnostic solve only. No production engine or vendor file was changed.

## Provenance and reproduction

- Production/source anchor: 4aec88d0e3ff028e0a65a8ac36b30948625df792
- Vendor gitlink and checkout: d6052d605673ce9d82cdc99b0df79855fa2ea215
- Runtime: Python 3.12.15, NumPy 2.5.4, SciPy 1.18.1.
- All arms: \(r_{\max}=80\), amplitude \(0.01\), width \(7\), \(D_{\mathrm{amplitude}}=10^{-10}\), radiation enabled, \(\Delta t=0.03\,\Delta r\); 15/30/60 accepted steps for N=60/120/240.
- All three matrix jobs passed source/vendor/dependency fail-closed checks, compilation, full evolution, and artifact upload.
- Baseline initial geometry matches saved initial geometry to a maximum per-field discrepancy of \(3.3\times10^{-14}\) (N=60), \(5.4\times10^{-14}\) (N=120), and zero (N=240).
- Baseline full saved final state matches the archived state within \(1.45\times10^{-12}\), \(7.08\times10^{-12}\), and \(7.58\times10^{-13}\). The \(A_a\) field itself agrees within \(1.80\times10^{-16}\), \(2.63\times10^{-15}\), and \(1.34\times10^{-15}\). Largest differences are in lapse/proper-time fields, not \(A_a\).

Artifact ZIP digest and internal checksums:

| N | Artifact ID/link | ZIP SHA-256 | report.json SHA-256 | profiles.npz SHA-256 |
|---:|---|---|---|---|
| 60 | [12584100185](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40153734946/artifacts/12584100185) | 826707e47b0824bae1a83bd58a929fb1bb6f4acf055a2dc792592ebcfcb42f93 | 6a95b16bc4b17f20a312ee780bac76186a37e573c3cbd73fab874211c29fe0c6 | af8bc011e18a55607ef9cb52552e5c84fd2dcc777c6b41898656d40aa4ff6781 |
| 120 | [12583184844](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40153734946/artifacts/12583184844) | 7213829442684a67518688e5cbd2a32f3fb70f28b0b95acaae7a46b038149a9c | 3483d665c1450427196a5b02db5ff13e5c41519d8018110656c16a4dce0b0bad | cb150833a077280488f92344cf98964df7d6d1f1bd2d3745d8fa3f7922d61bc2 |
| 240 | [12583826100](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/40153734946/artifacts/12583826100) | dba72873e5618a50493956f290608bbe7deda4d6167dd30636af853a6c440000 | 0081e962917ab075d340d759a65d1d942719f303d3ead45a5c8f7fc555b9e474 | cef1c10d68ee253a594ca902860764af38872eba29e5c65848307ac99bb0d1e7 |

Each downloaded artifact's internal SHA256SUMS file verified.

## Accepted \(t=0.60\): centre-region residuals

Values are signed cell-zero \(A_a/r_0^2\), followed by maximum absolute \(H\), matter-coupled momentum residual \(M\), and connection residual \(C_\Lambda=\Lambda-\Lambda_{\rm metric}\), each maximum taken over \(r\leq3\).

| N | Arm | \(A_a/r_0^2\) | \(\max|H|\) | \(\max|M|\) | \(\max|C_\Lambda|\) |
|---:|---|---:|---:|---:|---:|
| 60 | Baseline | +8.18e-4 | 7.67e-4 | 1.27e-3 | 1.07e-3 |
| 60 | Half-cell only | -1.54e-5 | 4.36e-4 | 1.10e-4 | 1.17e-4 |
| 60 | CMC row only | +7.61e-4 | 6.16e-4 | 1.29e-3 | 1.10e-3 |
| 60 | Joint | -7.24e-5 | 2.89e-4 | 9.20e-5 | 8.91e-5 |
| 120 | Baseline | -7.05e-4 | 2.60e-4 | 1.03e-3 | 7.12e-4 |
| 120 | Half-cell only | +4.27e-4 | 9.12e-4 | 1.13e-4 | 9.94e-5 |
| 120 | CMC row only | -9.71e-4 | 4.99e-4 | 9.74e-4 | 7.49e-4 |
| 120 | Joint | +1.62e-4 | 4.28e-4 | 1.72e-4 | 6.21e-5 |
| 240 | Baseline | +1.66e-4 | 1.53e-2 | 6.51e-4 | 8.59e-4 |
| 240 | Half-cell only | -1.98e-3 | 5.42e-3 | 1.64e-3 | 3.61e-4 |
| 240 | CMC row only | +1.60e-3 | 9.35e-3 | 7.18e-4 | 4.47e-4 |
| 240 | Joint | -5.49e-4 | 5.50e-4 | 2.71e-4 | 5.16e-5 |

### Interpretation

- **Connection residual:** the joint arm lowers the centre maximum at all three resolutions, by approximately \(12\times\), \(11.5\times\), and \(16.6\times\) relative to baseline.
- **Momentum residual:** the joint arm lowers the centre maximum at all three resolutions, by approximately \(13.8\times\), \(6.0\times\), and \(2.4\times\). The improvement weakens at N=240.
- **Hamiltonian residual:** not monotonically improved. The joint arm is lower than baseline at N=60 and N=240 but higher at N=120. At N=240, the baseline peak is \(1.53\times10^{-2}\) at cell 0; the joint value is \(5.50\times10^{-4}\), but this does not produce a monotonic centre-limit sequence.
- **Centre regularity:** \(|A_a|/r_0^2\) for the joint arm is \(7.24\times10^{-5}\), \(1.62\times10^{-4}\), and \(5.49\times10^{-4}\) at N=60/120/240. It grows with refinement and does not settle to a common limit. Therefore the one-step regularity signal does **not** survive to accepted \(t=0.60\) as a convergence claim.
- **Single-factor controls:** half-cell-only often reduces \(M\) and \(C_\Lambda\), but its N=240 \(M\) and \(H\) residuals are worse than baseline. CMC-row-only does not repair the centre connection residual. The interaction matters; neither control alone suffices.

## Matched-radius profiles

Linear interpolation of accepted \(A_a\) and \(H\) profiles at shared radii \(r=1,\ldots,5\) is shown for baseline and joint arms. N=60 has first centre \(r_0=2/3\), so its \(r=1\) interpolation spans the very wide \([2/3,2]\) interval. Treat the inner region \(r\lesssim2\) as unresolved by N=60; do not interpret that interpolation as a centre estimate.

### \(A_a(r)\): baseline / joint

| r | N=60 baseline / joint | N=120 baseline / joint | N=240 baseline / joint |
|---:|---:|---:|---:|
| 1 | +2.5651e-4 / -2.5698e-5 | -3.4335e-5 / +4.9196e-6 | +2.2447e-5 / -1.4695e-5 |
| 2 | -6.4842e-5 / -6.2920e-6 | -2.2823e-5 / -4.7716e-5 | -3.7201e-5 / -3.8056e-5 |
| 3 | -6.3256e-5 / -7.7834e-5 | -7.0691e-5 / -7.2788e-5 | -7.0668e-5 / -7.0892e-5 |
| 4 | -8.9938e-5 / -1.0933e-4 | -1.0193e-4 / -1.0272e-4 | -1.0200e-4 / -1.0210e-4 |
| 5 | -1.2054e-4 / -1.2090e-4 | -1.2390e-4 / -1.2428e-4 | -1.2347e-4 / -1.2352e-4 |

### \(H(r)\): baseline / joint

| r | N=60 baseline / joint | N=120 baseline / joint | N=240 baseline / joint |
|---:|---:|---:|---:|
| 1 | +5.2895e-4 / +2.6225e-4 | -4.4529e-5 / +2.4979e-4 | +3.7936e-4 / +5.8717e-5 |
| 2 | -1.8492e-4 / +1.8287e-4 | +9.9805e-5 / +1.3142e-5 | +4.3741e-6 / +3.6754e-6 |
| 3 | -1.6908e-5 / -2.8554e-6 | -3.6149e-6 / -5.9351e-6 | -1.7306e-6 / -1.7301e-6 |
| 4 | +2.0685e-6 / -4.0677e-5 | -7.5369e-6 / -7.6529e-6 | -2.0506e-6 / -2.0507e-6 |
| 5 | -2.7391e-5 / -1.3915e-5 | -5.0372e-6 / -5.0378e-6 | -1.3615e-6 / -1.3615e-6 |

The \(A_a\) profiles agree closely at \(r=3\)–5 for N=120 and N=240, while the largest \(H\) and \(A_a\) differences occur in the inner region. This localizes sensitivity but does not establish central convergence.

## Fixed physical support fit \(r\leq5\)

The same local even-polynomial fits were applied to each accepted \(A_a(r)\) profile at all resolutions. The support contains 4, 8, and 15 cell centres at N=60, 120, and 240. Model A is \(A_a=c_0+c_2r^2\); model B is \(A_a=c_0+c_2r^2+c_4r^4\). RMS and maximum residuals are measured at the cell centres used in each fit.

| N | Arm | Model | \(c_0\) | \(c_2\) | \(c_4\) | RMS residual | Max residual |
|---:|---|---|---:|---:|---:|---:|---:|
| 60 | Baseline | quadratic | +1.864e-4 | -1.678e-5 | — | 1.376e-4 | 1.847e-4 |
| 60 | Baseline | quartic | +3.081e-4 | -6.395e-5 | +2.075e-6 | 9.597e-5 | 1.504e-4 |
| 60 | Joint | quadratic | -1.798e-5 | -4.961e-6 | — | 2.256e-5 | 3.153e-5 |
| 60 | Joint | quartic | -9.808e-6 | -8.126e-6 | +1.393e-7 | 2.157e-5 | 3.379e-5 |
| 120 | Baseline | quadratic | -3.395e-5 | -3.774e-6 | — | 2.539e-5 | 5.544e-5 |
| 120 | Baseline | quartic | -3.584e-5 | -3.095e-6 | -2.832e-8 | 2.533e-5 | 5.566e-5 |
| 120 | Joint | quadratic | -8.670e-6 | -5.369e-6 | — | 1.757e-5 | 2.727e-5 |
| 120 | Joint | quartic | +7.942e-6 | -1.134e-5 | +2.492e-7 | 9.965e-6 | 2.213e-5 |
| 240 | Baseline | quadratic | +2.767e-6 | -6.320e-6 | — | 2.496e-5 | 7.714e-5 |
| 240 | Baseline | quartic | +1.919e-5 | -1.293e-5 | +3.101e-7 | 2.030e-5 | 6.516e-5 |
| 240 | Joint | quadratic | -1.714e-5 | -5.030e-6 | — | 6.617e-6 | 1.363e-5 |
| 240 | Joint | quartic | -1.021e-5 | -7.820e-6 | +1.308e-7 | 2.505e-6 | 5.408e-6 |

The intercept depends on model and resolution, including a sign flip for the joint N=120 fit between models. N=60 has only four points and is particularly sparse; residuals remain large relative to the profile amplitude. These are fixed-support diagnostics, **not** valid centre extrapolations or a convergence-order estimate.

## Conclusion and open gates

The half-cell integral is a strong initial-geometry lever and robustly improves the centre connection residual. The CMC row has a separate interaction with it, and the joint arm improves several central constraint signals. However, \(A_a/r_0^2\) remains nonconvergent at accepted \(t=0.60\), the Hamiltonian residual is resolution-dependent, and the inner-region fits cannot establish a regular centre limit. The research result narrows the cause without licensing a patch.

No production source, initializer, vendor, dependency pin, strict bound, or admission flag was changed. The full matter-coupled momentum-propagation identity and total stress-energy conservation/exchange gates remain open. Scientific admission remains false.
