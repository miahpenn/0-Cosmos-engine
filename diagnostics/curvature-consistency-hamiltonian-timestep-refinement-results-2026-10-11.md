# Fixed-grid Hamiltonian timestep refinement — results

Date: 2026-10-11  
Branch: `diag/independent-metric-ricci-audit`  
Workflow: [run 38097445638](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38097445638)  
Commit: `2834643ce53ac2a2f1084f911edde9e455029745`  
Artifact: [independent-metric-ricci-hamiltonian-time-refine](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/38097445638/artifacts/11687031435)  
Artifact ZIP SHA-256: `c8d1a2d622e5838e35a0270b62d7a3392240b782caf90e900d7b31295e65133f`

## Design

- Fixed N=80, (r_{\max}=40), (\Delta r=0.5)
- One discrete-consistent initial state, deep-copied for each trajectory
- CFL values 0.03, 0.015, 0.0075; nominal timesteps 0.015, 0.0075, 0.00375
- Samples (t=0,1,2,3), radiation on, amplitude 0.01, width 7, D amplitude (10^{-10})
- On each identical evolved state, report vendor H and counterfactual H with only the curvature scalar replaced by the independent metric-derived value.

No production physics, gauge, projection, or defaults changed.

## Central residual at t=3

| CFL | Signed vendor H | Signed metric-counterfactual H | (R_{vendor}-R_{metric}) |
|---:|---:|---:|---:|
| 0.0300 | (-1.12360\times10^{-5}) | (-1.46766\times10^{-5}) | (+3.44066\times10^{-6}) |
| 0.0150 | (-1.06604\times10^{-5}) | (-1.39482\times10^{-5}) | (+3.28772\times10^{-6}) |
| 0.0075 | (-1.05181\times10^{-5}) | (-1.37678\times10^{-5}) | (+3.24966\times10^{-6}) |

At (t=3), the central cell is also the maximum-(|H|) cell in all three timestep runs. The curvature gap is of order (3.3\times10^{-6}), while the vendor H residual is of order (1.1\times10^{-5}).

## Timestep sensitivity

For the maximum all-grid difference between trajectories as CFL is halved:

| Sample time | Field | max difference, CFL .03 vs .015 | max difference, CFL .015 vs .0075 | Ratio |
|---:|---|---:|---:|---:|
| 1 | (H_{vendor}) | (1.1680\times10^{-7}) | (2.9606\times10^{-8}) | 3.945 |
| 2 | (H_{vendor}) | (2.9980\times10^{-7}) | (7.4928\times10^{-8}) | 4.001 |
| 3 | (H_{vendor}) | (5.7553\times10^{-7}) | (1.4232\times10^{-7}) | 4.044 |
| 3 | (H_{metric}) | (7.2846\times10^{-7}) | (1.8038\times10^{-7}) | 4.038 |
| 3 | curvature gap | (1.5293\times10^{-7}) | (3.8062\times10^{-8}) | 4.018 |

The ratios near four are consistent with second-order timestep sensitivity in these measured differences. They are empirical indicators, not a stand-alone formal temporal convergence proof.

## Interpretation

1. Timestep size measurably affects the constraint trajectory and the curvature gap.
2. At (t=3), reducing CFL from 0.03 to 0.0075 changes the signed central vendor residual from about (-1.124\times10^{-5}) to (-1.052\times10^{-5}), approximately a 6.4% change. The residual does not vanish under timestep refinement.
3. The curvature gap itself also remains near (3.25\times10^{-6}), so timestep sensitivity does not erase the spatial/formulation discrepancy.
4. Earlier fixed-(\Delta t) resolution runs saw the central vendor residual fall from (6.41\times10^{-5}) (N=40) to (1.05\times10^{-5}) (N=80) to (1.87\times10^{-6}) (N=160). That is strong resolution dependence, but the effective CFL changed with N.
5. The next comparison should keep CFL fixed across N. This cross-check separates the resolution trend under a common dimensionless timestep rule. Do not infer causal ownership of any equation or change the vendor formulation based on this diagnostic alone.

## Gates

Workflow successful; diagnostic gate failures: none. Vendor Hamiltonian decomposition and curvature-replacement identity closed at floating-point tolerance. The N=80 initializer's vendor-defined H is at solver tolerance at t=0; the metric-derived counterfactual has a nonzero t=0 residual by construction, which should not be interpreted as evidence that either curvature route is physically preferred.
