# Cosmos Engine — Project Intent

## Purpose

This repository is an attempt to build a high-fidelity numerical laboratory for the coupled dynamical evolution of spacetime and multi-sector matter under general relativity.

The central goal is not merely to evolve fields, but to keep geometry, constraints, and stress-energy consistent with one another as the system evolves. In spherical symmetry, the engine is being developed around a reference-metric BSSN formulation coupled to multi-sector matter sources. The major continuum identities governing long-term accuracy—especially the matter-coupled momentum-constraint propagation identity, Hamiltonian-constraint propagation and control, reference-metric connection response, and total stress-energy / sector-exchange balance—must be derived from the accepted equations, frozen under explicit conventions, and subjected to rigorous diagnostic verification.

## Working philosophy

- Continuum identities are derived before residuals are interpreted.
- Conventions—including index placement, density weight, divergence measure, source normalization, connection coupling, and time split—are explicitly frozen and recorded before claims are made.
- Diagnostics are fail-closed, provenance-controlled, and operate on SHA-identified states with provenance-verified inputs. Hash identity is necessary for traceability, but does not alone establish scientific correctness or reproducibility.
- Scientific admission is withheld until the corresponding continuum and discrete behavior has been closed under those conventions.
- Production changes are authorized only after the relevant identities have been verified.
- Archive-backed research and accepted source pins govern the work. External formulations are comparators until their assumptions and mappings have been reconciled.

## Scope

Spherical symmetry is adopted as the controlled setting in which the consistency problem can be examined with transparency and precision. It is a deliberate scope choice: the exact spherical branch does not by itself represent the archive's full Bianchi-I shear. The diagnostic and provenance framework is intended to remain auditable and reproducible as the system extends to more ambitious dynamical regimes.

## Validation boundary

Implementation progress, software checks, diagnostic closure, continuum derivation, and scientific admission are separate claims.

An implemented evolution pathway, a passing test suite, a small or convergent residual, an algebraic cancellation, a frozen-state counterfactual, or a short smoke test does not independently establish the full matter-coupled propagation identity or total stress-energy conservation. In particular, the full matter-coupled momentum-constraint propagation relation and total six-sector stress-energy / sector-exchange closure remain gated by their explicit derivations and discrete verification. The current truth of those gates is recorded in [docs/STATUS.md](STATUS.md), the dedicated [momentum-constraint propagation issue](https://github.com/miahpenn/0-Cosmos-engine/issues/4), and the corresponding research pull requests.

Scientific claims should be limited to what has been explicitly verified under the stated conventions and the exact source, vendor, input-state and workflow provenance of that evidence.

## Long-term ambition

The longer-term ambition is a numerical system in which spacetime and matter evolve together according to the equations rather than according to accumulated constraint or conservation error. The goal is to build a numerical-relativity machine whose internal consistency is sufficiently well understood and controlled that it can serve as a credible foundation for simulating the dynamical evolution of a universe.

This is an ongoing research effort, not a claim that the complete coupled system has already met that standard.
