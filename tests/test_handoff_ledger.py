from engine.handoff import CycleLedger, handoff_from_ledger


def test_cycle_ledger_records_only_actual_crossings():
    ledger = CycleLedger()
    assert ledger.observe(1, 0.1, 0.1, 1.0, 0.9) is None
    event = ledger.observe(2, 0.2, 0.2, 0.1, 0.0)
    assert event is not None
    assert event.kind == "turnaround"
    assert len(ledger.events) == 1


def test_production_handoff_policy_is_event_level():
    # This regression is intentionally structural: HandoffPoint construction
    # is driven by CycleLedger events, while per-step diagnostics belong in
    # state.history.
    source = open("engine/production_kernel.py", encoding="utf-8").read()
    assert "if event is not None:" in source
    assert "candidate.handoffs.append(hp)" in source
