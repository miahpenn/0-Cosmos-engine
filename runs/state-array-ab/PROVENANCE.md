# State-array repair ON/OFF audit — provenance

## Record identity

- Diagnostic run: [GitHub Actions run 37859756564](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37859756564)
- Workflow: `.github/workflows/zero_star_state_array_ab.yml`
- Code commit used by the run: `7b9bac434f402df597641cfe7b24cf1a126a02e1`
- Branch: `0star-central-clock`
- Artifact: [state-array-ab-comparison](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37859756564/artifacts/11588005336), artifact ID `11588005336`
- Raw file inside artifact: `state-array-ab.json`
- Artifact ZIP size: 531,182 bytes
- Artifact ZIP SHA-256: `349c239c8f0fd394831e8dee2bb77e922d771e92a5ccc80599d480920dc64abf`
- Raw JSON size: 2,196,349 bytes
- Raw JSON SHA-256: `4e53d8c62a369ebd5517a72307966730a49fcfab243bf39b9849d410fcc33b07`
- Summary: [repair-on-off-comparison-summary.json](./repair-on-off-comparison-summary.json)

The artifact checksum was independently calculated after downloading the Actions artifact and matches the SHA-256 reported by GitHub for that artifact.

## Audit implementation

The diagnostic implementation at the run's code commit is:

- `engine/state_array_ab_audit.py` — captures read-only snapshots of numeric state arrays and calculates ON/OFF differences.
- `tests/test_state_array_ab_audit.py` — checks the audit implementation.
- `.github/workflows/zero_star_state_array_ab.yml` — runs repository tests, then the sequential state-array comparison, and uploads the resulting JSON.

The run executed the repair-ON case first. The repair-OFF case used the existing `NoRepairKernel` with only the center repair overridden. No production evolution equation or physics code was changed by this audit. Snapshots include geometry, scalar, matter, and grid arrays, as well as recorded time, `tau`, and e-fold values.

## Frozen configuration and admission

- Initial D amplitude: `1e-4`
- Resolution: `N=160`
- Domain: `r_max=80`
- CFL: `0.0075`
- Radiation: enabled
- Requested final time: `22.5`
- Requested samples: `t=0, 0.01, 4, 8, 12, 16, 20, 22.5`

Both ON and OFF cases completed without numerical failure and independently passed the audit's admission checks: frozen configuration matched, eight registered samples were present within the configured 0.05 time tolerance, recorded arrays covered all 160 cells, and all recorded array values were finite. Repository tests and the `state-array-ab` check passed on the exact run commit.

Admission intentionally applies **no state-difference threshold and no physical-residual threshold**. It confirms that the comparison is complete and well-formed; it does not certify physical correctness.

## Results

### 1. Difference exists at initialization

At (t=0), before the positive-time evolution samples, the largest ON/OFF absolute differences are in the innermost recorded cell, (r=0.25):

| Array | Maximum absolute difference | Cells not bitwise equal |
|---|---:|---:|
| `geometry.Lambda` | (3.50848\times10^{-4}) | 1 |
| `geometry.a` | (2.92375\times10^{-5}) | 4 |
| `geometry.b` | (1.46179\times10^{-5}) | 21 |
| `geometry.alpha` | (1.14900\times10^{-9}) | 159 |
| `geometry.Aa` | (5.16535\times10^{-15}) | 1 |

So this is not an initially identical-state experiment. The initialization difference must be accounted for when interpreting later differences.

### 2. Differences evolve in multiple arrays

By the first positive-time sample ((t=0.01125)), 21 of the 26 recorded arrays contain at least one non-bitwise-equal cell. By (t\approx4), many key geometry and scalar arrays differ at all 160 cells. Exact inequality across the grid is not by itself evidence of a resolved physical signal: the audit has no difference threshold, and the amplitude and spatial profile must also be considered.

The maximum absolute difference in `geometry.Aa` changes with time: approximately (5.2\times10^{-15}) initially, (1.16\times10^{-5}) at (t\approx0.011), (9.30\times10^{-5}) at (t\approx4), and (1.33\times10^{-4}) at (t\approx8); it then decreases to (3.35\times10^{-5}) at (t\approx12), (7.82\times10^{-6}) at (t\approx16), and (3.65\times10^{-5}) at (t=22.5). The maximum `geometry.Lambda` difference stays close to (3.51\times10^{-4}) across the recorded checkpoints. Maximum differences in `geometry.a` and `geometry.b` remain about (6.0\times10^{-5}) and (3.0\times10^{-5}), respectively, in the later samples.

At (t=22.5), recorded physical times match exactly. The recorded `tau` values are 6.983326299808656 (ON) and 6.983325414935567 (OFF), giving OFF minus ON of (-8.84873\times10^{-7}).

### 3. Spatial distribution at (t=22.5)

The grid regions below use the audit's existing partition: inner (r\le20) (40 cells), middle (20<r<64) (88 cells), and outer (r\ge64) (32 cells).

| Array | Inner max | Middle max | Outer max |
|---|---:|---:|---:|
| `geometry.Lambda` | (3.51534\times10^{-4}) | (1.01180\times10^{-11}) | (2.30907\times10^{-13}) |
| `geometry.Aa` | (3.64852\times10^{-5}) | (4.48552\times10^{-12}) | (1.39302\times10^{-13}) |
| `geometry.a` | (6.02154\times10^{-5}) | (1.56969\times10^{-10}) | (4.95060\times10^{-12}) |
| `geometry.b` | (3.01064\times10^{-5}) | (7.22054\times10^{-11}) | (1.79479\times10^{-12}) |
| `scalars.D` | (1.01413\times10^{-8}) | (8.23473\times10^{-12}) | (5.53608\times10^{-27}) |
| `geometry.alpha` | (3.67486\times10^{-10}) | (1.10618\times10^{-10}) | (6.16507\times10^{-12}) |
| `matter.dark_matter.energy_t` | (7.62055\times10^{-11}) | (5.08475\times10^{-12}) | (3.11126\times10^{-13}) |

The largest persistent differences in these geometry/scalar arrays are strongly concentrated in the inner region. However, the comparison is **not strictly center-only**: some arrays are bitwise different throughout the grid, and their small outer-region differences remain nonzero. Also, `geometry.alpha` has inner and middle maxima of the same order, so it should not be described as having the same degree of central concentration as `geometry.Lambda`.

## Interpretation boundary

This audit establishes what the recorded ON/OFF state arrays do under the frozen configuration; it does not identify a physical mechanism. In particular, it does **not** establish causal propagation, a propagation speed, that the repair acts only through (R), absence of structure elsewhere, constraint satisfaction or convergence, a bounce, or a trapped surface. No threshold was preregistered to translate array differences into a physical pass/fail judgment.

The full snapshots and per-field comparisons remain in the immutable Actions artifact listed above. This provenance file and the JSON summary are the durable, human-readable record; the large raw array payload is not duplicated into Git.
