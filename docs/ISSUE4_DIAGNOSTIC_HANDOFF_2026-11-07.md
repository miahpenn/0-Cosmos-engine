# Issue #4 diagnostic handoff — 2026-11-07

**Status:** diagnostic-only handoff; no production or scientific-admission change.  
**Related issue:** [#4 — discrete momentum-constraint propagation at the central stencil](https://github.com/miahpenn/0-Cosmos-engine/issues/4)  
**Published diagnostic branch:** [`issue4-integrated-test-diagnostics`](https://github.com/miahpenn/0-Cosmos-engine/tree/issue4-integrated-test-diagnostics)  
**Exact branch tip:** `8965a1985c99d644137321016dc5f172d5695213`

This record consolidates the recovered November 5–6 split-derivative work, its November 7 publication, and the exact-head CI evidence. It is deliberately a handoff, not a claim that Issue #4 is scientifically closed.

## 1. Recovered branch and immutable provenance

Claude's original bundle was verified before publication:

- Bundle: `issue4-integrated-test-diagnostics.bundle`
- SHA-256: `8544465262bfb4806f14dd18b477387fa6391c930928ef165c0eeab3683d6e1e`
- Required base: `f8d5a65ea0e9bb1064a55f6b0691f860cfaa748a`
- Branch: `issue4-integrated-test-diagnostics`
- Exact tip: `8965a1985c99d644137321016dc5f172d5695213`

The bundle passed `git bundle verify`; its sole prerequisite was the stated base commit. The published branch contains exactly these two commits:

1. `2937ef96dd1a7cc90e154b6dc151173758477de5` — `Test-only: integrated Issue #4 contract and research diagnostics`
2. `8965a1985c99d644137321016dc5f172d5695213` — `CI: separate blocking contract tests from visible non-blocking witnesses (prepared, not pushed)`

The second commit's historical subject is retained verbatim. The branch was pushed and the remote ref was verified at the exact tip above. No merge or production-branch update is implied by publication. Do not rewrite either commit or rebase this recovery branch.

The earlier Issue #4 recovery-audit comment stating that these commits were not yet on GitHub predates the successful November 7 push; the verified remote ref above supersedes that earlier publication status.

## 2. Exact-head CI evidence

Workflow run: [42407237209](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/42407237209)  
Checked-out SHA in the research job: `8965a1985c99d644137321016dc5f172d5695213`

The research job used CPython 3.12.15, NumPy 2.5.4, pytest 9.1.1, and mpmath 1.3.0.

Job outcomes:

| Job | Result | What this establishes |
|---|---|---|
| `independent-cmc` | success | The blocking affine-root and implicit-derivative contracts ran successfully; the visible historical witnesses executed; the zero-production-call contract step completed. |
| `accepted-state-crosswalk` | success | The immutable accepted-state artifact was downloaded, extracted, and the six-sector/pinned-K RHS crosswalk completed. |
| `research-diagnostics` | failure | Port validation passed (2 tests), but exact regeneration comparison of the two committed numerical references failed. |

The failing research step regenerated both references to temporary runner files. Its exact comparisons were:

- Affine-root: committed `root` and `cond_inf` must equal regenerated values.
- Implicit derivative: every committed `derivative[epsilon]` value must equal the regenerated value for that epsilon.

The job printed:

```
affine-root reference identical: False
derivative reference identical: False
```

Thus this is a **numerical reference reproducibility mismatch**, not merely a metadata mismatch. The log does not print the numerical deltas, so this run alone does not quantify their magnitude or identify their cause. Both generators reported mpmath 1.3.0 in this run. The committed files were only read; regenerated outputs were written under the runner temporary directory. No committed reference was overwritten.

Do not describe the entire workflow as green: the `research-diagnostics` job failed even though the blocking-contract and accepted-state-crosswalk jobs succeeded.

## 3. Interpretation and limits

- Passing blocking contracts supports the specific test contracts exercised; it does not establish the correctness of a full matter-coupled continuum propagation identity.
- The accepted-state crosswalk passing is not, by itself, proof of total stress-energy conservation or scientific admission.
- The research reference mismatch remains open. The run establishes non-identity under the pinned CI environment, but does not establish whether the source is environment, generator provenance, arithmetic evaluation, or another reproducibility dependency.
- The mismatch is not resolved by pointing only to a version difference: both reference generators reported mpmath 1.3.0 in run 42407237209.
- The available investigation did not recover a verified original source file for `mp_fixture_port.py`. Do not claim its implementation or provenance has been recovered until the actual file and hash are located.
- JSON metadata differences (including fields such as `cond2_note` or `helper_sha256`) must be reported separately from the numerical fields actually compared by CI. They are not the comparison predicate in the workflow step above.

## 4. Related November 5–7 research notes

The following Library handoffs informed this record and should be preserved as supporting context, not silently treated as repository evidence until their source artifacts are copied or independently verified:

- `PROPOSAL_issue4_implicit_derivative_diagnostic.md` and its duplicate upload: proposes a pinned mpmath 1.3.0 regeneration procedure while explicitly preserving the committed reference. Its synthetic-fixture report gives maximum discrepancy `7.913e-13` against high-precision exact-root finite differences and implicit residual `1.6e-17`. These are synthetic-fixture results, not production-physics evidence.
- The November 7 Hamiltonian-slot correction/options notes withdraw the earlier assertion that metric-only and BSSN Hamiltonian slots differ only by the Lambda term. Neither slot is selected as preferred; retain the two-slot reporting approach and investigate the unexplained non-Lambda difference. This is a separate Issue #4 derivation lane and must not be conflated with the CI regeneration mismatch.
- The Brown local reproduction note is a separate reproduction record; its claims should remain tied to its own source, command, inputs, and outputs.

## 5. Required next actions — in order

1. **Preserve current evidence.** Keep the published recovery branch, both committed reference JSON files, the generator scripts, and the failing CI result unchanged.
2. **Capture actual regenerated outputs.** Add an opt-in diagnostic mode or a separate script that writes regenerated root and derivative JSON to named artifacts, then reports per-field numerical differences against the committed references. Do not replace the committed references.
3. **Trace exact provenance.** Record hashes for both committed JSON files, both generator scripts, imported helper modules, and any external fixture/port source. Locate the original `mp_fixture_port.py` or explicitly record that it is unavailable. Record Python, NumPy, mpmath, platform, and relevant BLAS/runtime details.
4. **Reproduce in the pinned CI environment.** Use CPython 3.12.15, NumPy 2.5.4, pytest 9.1.1, and mpmath 1.3.0. Compare numerical fields separately from metadata. Preserve every run and its deltas; do not average run-to-run variation away.
5. **Only after the source of the mismatch is identified**, propose a narrowly scoped diagnostic/test correction. Do not relax tolerances, regenerate/replace committed references, or alter the blocking contract merely to make CI green.
6. Keep Issue #4 and the shared coordination/task board current with the exact branch/commit, run URL, artifacts, result, and remaining limitation.

## 6. Scientific and implementation guardrails

- No production equations, evolution code, initializer, vendor source, solver tolerances, coefficients, gauge, boundary conditions, floors/clamps, or damping are changed by this handoff.
- No coefficient fitting, hidden hand-tuning, residual-based operator selection, or choosing between Hamiltonian slots by whichever residual is smaller.
- Research adapters and diagnostic fixtures remain separate from production.
- A successful software/contract check does not close the full matter-coupled momentum-constraint propagation gate.
- Keep `scientific_admission=false` and all conservation/closure admission flags false until their explicit evidence gates are satisfied.
- Preserve exact source SHA, runner environment, and artifact provenance for every result.

## 7. Explicitly unresolved

- Numerical deltas between the committed and regenerated affine-root and implicit-derivative reference values.
- The complete dependency/provenance chain for the reference generators and the missing `mp_fixture_port.py` source.
- Whether the CI mismatch is deterministic across repeated clean runs in the same pinned environment.
- The separate Hamiltonian-slot non-Lambda discrepancy.
- The full frozen, convention-pinned matter-coupled momentum-constraint propagation identity required by Issue #4.

**Disposition:** retain the branch and evidence as diagnostic-only. No merge or scientific-admission claim follows from this record.
