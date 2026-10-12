# Checkpoint constraint profile audit protocol

This diagnostic reads the saved N=40, N=60 and N=80 production checkpoints from the completed t=12 run 38208253687. It restores each checkpoint into the exact production source and pinned vendor numerical kernel and recomputes the signed Hamiltonian, momentum and connection profiles. The state is not evolved.

It reports, at each saved time t=2,4,6,8,10,12:
- peak absolute residual and its radius inside r<5;
- peak absolute residual and its radius on the fixed exterior region r>=5;
- full-domain peak and radius;
- signed profile samples at fixed target radii.

The audit is intended to distinguish a central regularity defect from a residual that propagates into the exterior. It is post-processing only; the profiles share the production stress-energy assembly and pinned constraint operator, so this is not an independent derivation of the Einstein constraints and does not authorize a physical interpretation.
