# Lock-analysis audit bundle registry — 2026-10-09

Audit summary: [LOCK_ANALYSIS_AUDIT_ADDENDUM_2026-10-09.md](https://github.com/miahpenn/0-Cosmos-engine/blob/0star-central-clock/GEAR_group/LOCK_ANALYSIS_AUDIT_ADDENDUM_2026-10-09.md)

**Decision:** INCONCLUSIVE for all D-bearing cases; Gate 4 frame validity is 35/90 (0.389) for three D cases and 33/90 (0.367) for Dhalf, against the frozen 0.85 requirement. No machine/evolution run was launched.

## Exact companion-output hashes

The full ZIP below contains the audit addendum, aggregate summary JSON, the checkpoint-level SHA-256 manifest (540 checkpoint members) and six per-case JSON outputs. The archive is attached in the accompanying ChatGPT response; it is not a GitHub Actions artifact.

- Companion ZIP: `LOCK_ANALYSIS_AUDIT_BUNDLE_2026-10-09.zip`
  - SHA-256: `f6a67eb6dc1e84d44c260549bebb6cd1774b105f3b562f399b245623249e1e3e`
- Full audit addendum: `LOCK_ANALYSIS_AUDIT_ADDENDUM_2026-10-09.md`
  - SHA-256: `3e6b1903d10dccf7d052d3fa16c35f8bce6b2c7e0214c239ace80a9ec866e2aa`
- Aggregate JSON summary: `lock_analysis_audit_summary_2026-10-09.json`
  - SHA-256: `b6aa3aaf738ce511ac3af61a9a2c3e172d3ab904b49ea58f2b0377abd0a992b9`
- Checkpoint manifest: `lock_analysis_checkpoint_manifest_2026-10-09.json`
  - SHA-256: `d5e767276a8b54749e20ce9442a99d8a3b19066a67ba753a7112287cf86334af`

### Per-case JSONs

| File | SHA-256 |
|---|---|
| `lock_case_r160_D0_2026-10-09.json` | `44db402cbea0fa11eea71d37a92cd79ff25c9ec6864a1806334050bded2a709a` |
| `lock_case_r160_Dbase_2026-10-09.json` | `fbb6c4ef04efebbf6ed5c3a38e725143667fdf682adb11c9dd79f153a5ff57a5` |
| `lock_case_r160_Dhalf_2026-10-09.json` | `8a0114b922a5765f8c989632985a7ee3bc218e8bb82c5dc716714994525e4642` |
| `lock_case_r80_N160_D0_2026-10-09.json` | `626d9eb9fbcb878d45b09630e77ab4d6b85540cc2316b7e93a59a0366648c5f2` |
| `lock_case_r80_N160_Dbase_2026-10-09.json` | `f80e374dd2e4215ea0b73dbc18b635e7f59d2865ea6f5966783c7a1a686c70ec` |
| `lock_case_r80_N320_Dbase_2026-10-09.json` | `6f8d5c474f8640902f4899f85a08760f19afb6621d74eb7ab2065f7c1438ff13` |

The checkpoint hashes refer to exact checkpoint member bytes inside the four original GitHub Actions artifact ZIPs. Those original ZIP digests were verified against GitHub artifact API metadata. See the linked addendum for the scientific decision and limitations.

**Important limitations:** The scoping memo fixed a PRNG seed but not the PRNG implementation; the independent reconstruction used NumPy default_rng/PCG64. The source report's `lock_full_results.json` was not present in the retrieved artifacts, so the six case JSONs are reconstructed. The memo did not operationalize a control-uncertainty estimator/threshold; descriptive pair-to-pair spread is not a formal uncertainty decision. No physical claim is entered and no new run is authorized.
