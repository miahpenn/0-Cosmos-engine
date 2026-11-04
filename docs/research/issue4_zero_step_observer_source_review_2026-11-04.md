# Issue #4 — zero-step observer source review

Date: 2026-11-04  
Reviewed branch: `research/issue4-zero-step-bianchi-observer-2026-11-04`  
Reviewed head: `eaad31c21a122a4cb739057f41115fb4783a7a9c`  
Scope: source-level review only. No evaluator execution, workflow trigger, production step, production/vendor/initializer edit, or scientific-admission change.

## Review disposition

The accepted-state evaluator has a coherent diagnostic structure and preserves the principal scientific distinctions:

- archive ZIP, report, NPZ, and per-array checksums are validated before state construction;
- upstream archive identity and tangent-helper canonical state digest remain distinct;
- the accepted source, protected helper baseline, source blobs, and vendor commit are separately pinned;
- metric-only Ricci is computed on a disposable geometry copy, while the evolved-connection response is reconstructed separately;
- (\chi=X^2) is kept distinct from (\chi_r=-D_1X/(2X)), with odd parity for (C^r) and even parity for metric fields;
- pointwise and finite-volume total-stress divergence channels remain separate;
- the regression gate checks signed cell witnesses, annular and interior maxima, and the directional epsilon sweep;
- the evaluator does not construct a production kernel instance or call an RHS/step.

This is a source review, not a runtime pass. The committed evaluator has not yet been run against the archives, and the latest Ricci-offset gate has not yet been executed.

## Findings and remediation status

### 1. Imported-module provenance is not fail-closed

**Location:** `tools/issue4_full_bianchi_evaluator.py`, `_load_vendor`.

The function verifies that the vendor checkout exists and is pinned, then prepends its path to `sys.path` and performs ordinary imports. Python may nevertheless reuse same-named modules already present in `sys.modules`, or resolve a module from a different path. The source checkout is verified, but the origin of the modules actually used by the evaluator is not.

**Remediated on the research branch:** `_load_vendor` now checks each loaded module's resolved `__file__` against the exact pinned vendor root. A preloaded same-named module from outside that root is refused rather than silently reused. A regression test covers that refusal. This verifies the three explicitly imported vendor modules; transitive dependencies are not yet independently enumerated.

### 2. Evolution tripwires are not restored

**Location:** `evaluate_accepted_fixture`, where `V55ProductionKernel.step` and the vendor PIRK function are replaced.

The replacements correctly refuse accidental calls, but the original attributes are not restored. In a long-lived interpreter or test process this leaves globally imported classes/modules modified after evaluation, even when evaluation succeeds or fails.

**Remediated on the research branch:** the `_zero_step_guard` context manager saves both original methods, arms the tripwires, and restores them in `finally`. The full protected evaluation is inside the context. Regression tests exercise restoration after a normal exit, a tripwire exception, and an injected exception.

### 3. File output appends a literal escape sequence

**Location:** `main`, both `args.output.write_text` branches.

The source currently appends `"\\n"`, which writes the two literal characters backslash and n rather than a line terminator. Standardize both branches on an actual newline (`"\n"`) and test the resulting file ending.

This is a serialization defect, not a numerical one.


## Remediation commits and verification boundary

The three tooling fixes and their regression tests were committed on the research branch:

- `3cf79e2`: validate actual vendor module origins;
- `b1c191e`: add the restoring tripwire context manager and correct JSON output newlines;
- `d4d2e51`: add module-origin, tripwire-restoration, and newline regression tests;
- `3e73fd2`: finalize the two newline writes after source inspection.

These changes have received source-level inspection only. No test suite, evaluator execution, workflow, or production step has been run after these commits. The next gate is a local, zero-step execution of the focused tests and then the accepted-state evaluator under exact pins; do not use the existing lightweight-validation workflow because it contains a production-step smoke.

## Scope boundary still outstanding

The evaluator currently exposes only `evaluate_accepted_fixture` and accepts the three archived accepted-state resolutions (N=60/120/240) at (t=0.60). The foundation CLI's `--state-class` option is only reported as not evaluated; it does not implement manufactured or frozen-production-(t=0) adapters.

Therefore the original three-class observer objective remains incomplete:

1. manufactured regular profiles;
2. frozen production (t=0);
3. accepted (t=0.60) state.

The first two classes must be added under the same pinned source/vendor contract, same predictor conventions, signed term accounting, and zero-step guarantees. They must not be made to appear equivalent by changing their initializer, production equations, or source state. If a class cannot satisfy the same contract, the observer should report the specific missing operand/assumption and refuse comparison.

## Test-coverage observations

The current tests establish useful helper-level contracts and static source properties. They do not yet execute the committed evaluator against an accepted archive, exercise module-origin refusal, assert tripwire restoration, or validate the output-file newline. The existing reference-regression test uses synthetic arrays to test the gate; it is not a substitute for evaluating archived data.

The existing lightweight-validation workflow includes a production-kernel one-step smoke. Do not use it to run these checks. Before any CI invocation, inspect the complete workflow and ensure the selected route cannot advance a production state. A local zero-step execution path is preferable once archive artifacts and exact pinned checkouts are available.

## Scientific status

No physical term is selected or rejected by this review. No residual is interpreted as proof of a physical cause. The accepted-state annular near-fourth-order behavior and localized centre-cell remainder remain as previously documented. No manufactured or (t=0) comparison has been performed.

All scientific-admission and acceptance flags remain false. No production step was run during this review.
