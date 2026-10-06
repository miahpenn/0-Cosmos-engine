# Campaign Log

This page records the engineering campaigns that materially changed the diagnosis. It is intentionally concise; raw checkpoints and ledgers remain in GitHub Actions artifacts.

## Worldtube / boundary campaign

- [Run 28 — characteristic outer boundary](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37404042983) — reduced the outer geometry catastrophe and exposed the next matter-shell layer.
- [Run 29 — Rmax 80 worldtube scaling](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37405247985) — moved the worldtube outward; the first remaining failure became radiation conservative-state admissibility at the outer cell.
- [Run 35 — CMC boundary residual audit](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37409512861) — showed the CMC elliptic equation is well solved until the outer shell itself becomes contaminated.
- [Run 41 — deep causal instrumentation](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37410453267) — showed late radiation remains outgoing and the late matter/CMC layer is not an incoming radiation reflection.

## Worldtube-radius scaling

- [Rmax=160 campaign](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37413250578) — moving the worldtube outward moved the numerical turnaround substantially later while preserving healthy evolution for a longer interval.

This established that the earlier post-turnaround lapse-gradient layer is not simply a fixed local bulk feature. The remaining question is whether the observed timing dependence is a coordinate/time-interface effect, a genuine global CMC coupling, or both.

## 0* / central-clock investigation

- [Central proper-time experimental PR](https://github.com/miahpenn/0-Cosmos-engine/pull/1) — isolated archive-derived central-clock CMC coordinate.
- [First central-clock campaign](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37493375153) — lapse-only rescaling failed before turnaround at an inner radiation cell. This exposed the need for a consistent transformation of every component that interprets the evolution time parameter.
- [Current central-clock campaign](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37497656695) — tests the corrected CMC normalization/reparameterization path.
- [Low-CFL central-clock campaign](https://github.com/miahpenn/0-Cosmos-engine/actions/runs/37497656802) — independent lower-CFL consistency test.

## Interpretation rule

Failures are not promoted to physical claims until the invariant witnesses, resolution behavior, and worldtube dependence agree. Borrowed numerical tools remain pinned and are treated as diagnostic/numerical infrastructure; adoption into the physical model requires the project 0* and archive-consistency checks.

## Standing workflow

**Run the complete machine → expose the actual defect → repair only that layer → rerun.**