# Pointwise CMC target/lapse closure audit — 2026-10-12

**Status:** completed, successful. **Diagnostic-only:** no production source or kernel changed.

- Workflow: [38210163221](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38210163221)
- Artifact: [11729296911](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38210163221/artifacts/11729296911)
- Artifact ZIP SHA-256: `cfb1b7a71d4267c15f2b07dd123ecac2f83195b9be012a4e6a84cf1f86966c33`
- Report SHA-256: `cb2c576a093975e263279c33d7bc72a060273c4553501c27695c54162286c5fb`
- Per-cell JSONL SHA-256: `83129b94087efc03806d95cc162bcd17e580c6f465d99892aa1e9d3eddde71db`
- Code: `diagnostics/cmc_target_lapse_pointwise_residual.py`
- Trigger commit: [ff7a14d2e30c9f99a0dadbafe0914bef43511bfd](https://github.com/miahpenn/0-Cosmos-engine/commit/ff7a14d2e30c9f99a0dadbafe0914bef43511bfd)

## Purpose

The previous scalar closure scan found no sign change in (H(k)=\mathrm{target\_kdot}(\alpha_{solve}(k))-k) over the tested bracket. This audit asks whether the roughly (-8.58\times10^{-14}) offset comes from the interior lapse PDE solve or from the way the pointwise K-RHS is averaged.

## Results

Three fixed target values were tested on the archived t=0 R80 state: native R80 target, midpoint between R80 and R160 targets, and archived R160 target. The lapse solver returned finite, positive lapse profiles in all cases.

| Trial | Kdot trial | Recomputed target | (H(k)) | Weighted mean pointwise residual |
|---|---:|---:|---:|---:|
| Native R80 | (6.198841550010287\times10^{-5}) | (6.198841541429242\times10^{-5}) | (-8.581045487891792\times10^{-14}) | (-8.581045884335190\times10^{-14}) |
| Midpoint | (6.216183207871493\times10^{-5}) | (6.216183199288647\times10^{-5}) | (-8.582846618750833\times10^{-14}) | (-8.582846974608617\times10^{-14}) |
| Archived R160 | (6.233524865732698\times10^{-5}) | (6.233524857152139\times10^{-5}) | (-8.580558952166889\times10^{-14}) | (-8.580558790364570\times10^{-14}) |

- The 159 elliptic PDE rows excluding the outermost normalization cell close at approximately (1.3\!-!1.5\times10^{-15}) RMS and (3.6\!-!5.4\times10^{-15}) maximum absolute residual.
- For each trial, (H(k)) minus the weighted mean pointwise residual is only (1.6\!-!4.0\times10^{-21}), consistent to roundoff.
- The largest all-cell residual is at the outermost normalization cell: about (-4.5882\times10^{-12}).
- That outermost cell prescribes the boundary lapse; it is not an elliptic PDE row. Its pointwise K-RHS nevertheless enters the global target average.

## Interpretation

The interior elliptic solve is not the source of this closure offset. The numerical offset is accounted for by the boundary-row/global-average seam: a cell whose lapse is imposed by normalization participates in the target average despite not being an interior PDE row.

This does **not** prove the production target definition is wrong. Inclusion/exclusion of this cell must be checked against the authoritative CMC derivation before any production patch. The alternative target averages proposed below are diagnostic probes only, not substitute physics.

## Next diagnostic

On the same frozen R80 slice, calculate target averages with the same weights for:
1. all cells (current production definition);
2. PDE rows only;
3. all cells except the outermost lapse-normalization cell.

Verify the target shift equals the exact outer-cell weighted contribution. Then compare (H(k)) under these alternate diagnostic definitions. Do not edit `engine.cmc_gauge.target_kdot`, add damping, fit parameters, or modify production equations. Resolve the intended summation domain from the authoritative CMC derivation before proposing a source change.
