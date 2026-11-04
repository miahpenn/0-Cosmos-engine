# Issue #4 zero-step observer — implementation audit log

Date: 2026-11-04  
Status: diagnostic tooling only; no numerical state evaluation.

## Accepted-source blob audit

The accepted comparison source anchor is
`4aec88d0e3ff028e0a65a8ac36b30948625df792`. Its Git tree was inspected
directly, and the observer's complete source-blob table was compared against
the actual Git blob identities at that exact commit.

The initial observer table contained five stale values that belonged to the
historical helper-baseline source identity rather than the accepted diagnostic
anchor. They were corrected in the observer module. The five corrected paths
are:

| Path | Accepted-anchor blob |
|---|---|
| `engine/production_kernel.py` | `1720615628a98a8d7cfc9fcadcc984fe28c88be1` |
| `engine/matter_system.py` | `0969f4d996e4f8e0fc5b70a7dff147f024907302` |
| `engine/matter_rhs.py` | `0153c232ee15b5a9d3dbac37ea858e4d2d730011` |
| `engine/scalar_system.py` | `e2c2dae9a016b4d6454cf74a9eae8b0167d5733d` |
| `engine/v55_pirk_adapter.py` | `be28ceea63fa1b59efc5cc5c460b9df19ae1cbc3` |

All 11 required engine source blobs now match the accepted anchor. Both helper
blobs also match the accepted anchor exactly:

- tangent helper: `3f00a4e0a2063b86d2e007c5a825f9edfff61a4c`;
- total-stress helper: `c10f9e8b9e654d5195fd55c493c3fe06ba5c7082`.

The vendor identity was independently confirmed as a Git submodule gitlink at
`vendor/bb-palatini-unified-r0`, pointing to
`d6052d605673ce9d82cdc99b0df79855fa2ea215`.

The historical helper baseline remains separately pinned to
`7f47669a1c0edfafef9fc5f27b2db80840d3f035`. It is not used as the accepted
runtime source anchor. The correction removes the accidental mixing of five
historical blob identities into the accepted-source blob table.

## Tests added

The foundation test suite now includes mocked provenance tests that require:

1. the accepted diagnostic checkout and historical helper baseline to resolve
   to their distinct exact commits;
2. all 11 accepted-source blobs and both helper blobs to match their pinned
   identities;
3. the vendor submodule to match its exact commit;
4. the historical baseline commit to be rejected if supplied as the accepted
   diagnostic checkout.

The tests are committed but have **not yet been executed** in a repository
checkout. No workflow was dispatched and no numerical observer was run.

## Remaining implementation boundary

The observer is still only a provenance/digest foundation. The canonical
full-Bianchi evaluator, accepted-state archive regression, production-step
tripwire, and manufactured / initial-state adapters are not implemented.
This log does not authorize or report any state-class evaluation.
