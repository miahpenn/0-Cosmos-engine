# W1 corrected slice: 'f under Amendments 5 and 6 (S3 separated)

**Result:** the corrected slice passes every input check on both grids. G0.5'f now reads within 4e-5 of 1 on the
admission grid (dr 0.025), and within 2.2e-4 on dr 0.05.

The setup was declared in `PRERUN.md` (committed 34e7677, before the build). The go was Jeremiah's "Go", Nov 24
03:47 UTC. Zero-step only: no evolution, no f_ring run.

## Reading on the corrected slice (`native_static.npz`, #580)

| grid | G0.1-4, entry, 'a-'e | sum-back S / PS (tolerance) | 'f D | 'f S | band |
|---|---|---|---|---|---|
| dr 0.025 (admission) | PASS | 2.2e-19 / 4.4e-20 (1.94e-13) | 0.999961 | 0.999973 | [0.999, 1.001] PASS |
| dr 0.05 (record only) | PASS | 2.2e-19 / 2.0e-20 (1.95e-13) | 1.000050 | 0.999782 | [0.9985, 1.0015] PASS |

**Admission reading (A3 5): PASS.**

**How 'f was computed:**
- **Mean:** kappa-weighted, over whole local wavelengths, at the producer's discrete k_d. The windows the script
  computed match PRERUN exactly: D 17 / 16 wavelengths, S 15 / 15.
- **Derivatives:** the engine's 4th-order operator.
- **S field:** taken from the S3 harmonic, built from the stored vectors through the producer's unchanged `mapped()`
  with the S1 entries zeroed. S1 + S3 sums back to the slice at round-off.

The A3 version of 'f is kept in `report.json` as a record, not a gate.

The ~0.6% excess in #580's dr 0.05 S reading does not appear here. That fits A6's account of it: continuum k against
the producer's discrete k_d, plus the unweighted end term.

## Fault cases (`checks.json`)

Each case's S split follows its own construction.

| case | 'f D / S, dr 0.025 | 'f D / S, dr 0.05 | fails |
|---|---|---|---|
| uncorrected (admitted) | 1.2362 / 1.2844 | 1.2258 / 1.2811 | 'f |
| v^2 piece missing | **0.9867** / 0.9783 | 0.9918 / 0.9785 | 'a, 'c, 'f |
| wrong-sign phase | 1.4767 / 1.5654 | 1.4650 / 1.5668 | 'f |
| D tail reversed | -0.53 / 1.0000 | -0.54 / 0.9998 | 'd, 'e, 'f |
| S tail reversed | not read (sum-back stop) | not read | 'e, plus the sum-back STOP |

The v^2-missing D case, the narrowest margin, reads 0.9867 on the admission grid. That is 13x the band's half-width
away from 1.

## Reproduce

```
python3 research/w1_c9_g05p/run.py native_static.npz members.json OUT --static build.json --harmonics V10_DIR native.py
python3 research/w1_static_slice/checks.py BUILD_DIR native.npz members.json OUT --harmonics V10_DIR native.py
```

The inputs are the same as for `results/research/w1_static_slice`.

## Next, on Jeremiah's go

The zero-step C9 set (G1-G4) on the new bytes, then the f_ring run.
