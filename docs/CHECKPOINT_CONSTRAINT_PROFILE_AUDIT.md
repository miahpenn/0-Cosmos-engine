# Checkpoint constraint profile audit protocol

This diagnostic reads the saved N=40, N=60 and N=80 production checkpoints from completed run 38208253687. It restores each checkpoint into the exact production source and pinned vendor numerical kernel and recomputes signed Hamiltonian, momentum and connection profiles. It does not evolve the saved state.

The audit also recreates the production initial slice at each resolution without advancing time. Comparing t=0 with saved checkpoints helps distinguish an initial-data/operator mismatch from a residual that develops during evolution.

For each initial state and scheduled checkpoint (t=2,4,6,8,10,12), it records:
- peak absolute residual and its radius inside r<5;
- peak absolute residual and its radius on the fixed exterior region r>=5;
- full-domain peak and radius;
- signed profile samples at fixed target radii;
- geometric Hamiltonian and momentum terms separately from their total stress-energy source terms, including scalar/fluid Hamiltonian-source split;
- the scheduled checkpoint time separately from the actual saved integrator time.

The audit is post-processing only. Its reconstructed residuals share the production stress-energy assembly and pinned constraint operator; it is not an independent derivation of the continuum Einstein constraints and does not itself authorize a physical interpretation.
