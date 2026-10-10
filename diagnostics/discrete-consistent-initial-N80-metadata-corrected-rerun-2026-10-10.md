# N=80 matched baseline/candidate decomposition — metadata-corrected rerun

**Purpose:** Re-run the completed N=80 pair after correcting a report metadata defect. The previous numerical trajectory used N=80 (200 accepted steps at dt=0.015 over t=3), but its JSON `settings.resolution` field incorrectly said 40 because the runner serialized the default settings rather than the selected settings. This rerun is required for clean, auditable evidence.

- Trigger marker: `[run-n80-decomposition-pair]`
- N=80, r_max=40, dr=0.5, CFL=0.03, dt=0.015, final coordinate time t=3
- Baseline and candidate run sequentially from the same source commit.
- Tests run before the trajectories.
- Both JSON reports must record resolution=80.
- Diagnostic only; no production physics or defaults change.

The earlier run completed successfully and showed candidate central H(t=3)=-1.1238565881e-5 and baseline central H(t=3)=+1.3997218214e-3. Treat these as provisional until the metadata-corrected rerun is admitted.
