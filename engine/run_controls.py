"""Executable cheap control runner for the canonical engine.

This deliberately runs the current reference RHS modules without inventing
an exchange/source closure. It is a diagnostic harness, not the production
coupled integrator.
"""
from dataclasses import asdict
from engine.cosmos import CosmosParams, cosmos_rhs, cosmos_rho_p
from engine.local import LocalParams, local_rhs, local_energy
from engine.interface import exchange_from_flux
from engine.diagnostics import friedmann_constraint


def euler_step(y, rhs, dt):
    dy = rhs(y)
    return tuple(v + dt * dv for v, dv in zip(y, dy))


def run_controls():
    cp = CosmosParams()
    lp = LocalParams()

    # COSMOS control state. It is deliberately not forced to satisfy the
    # Friedmann constraint; the residual is reported rather than hidden.
    cosmos = (1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    rho, _ = cosmos_rho_p(cosmos, cp)
    friedmann0 = friedmann_constraint(cosmos[1], rho)

    # Local GEAR-03 control.
    local = (0.02, 0.0, 1.0e-8, 0.0)
    local_rhs0 = local_rhs(local, lp, H=0.0)
    local_e0 = local_energy(local)

    # Interface algebraic conservation control.
    exchange = exchange_from_flux(1.0)

    return {
        "physics_locks": {
            "kappa": lp.kappa,
            "g": lp.g,
            "G": lp.G,
            "alpha_D": lp.alpha_D,
            "V0": cp.V0,
        },
        "cosmos": {
            "rho0": rho,
            "friedmann_constraint0": friedmann0,
            "rhs0": cosmos_rhs(cosmos, cp),
        },
        "local": {
            "state0": local,
            "rhs0": local_rhs0,
            "energy0": local_e0,
        },
        "interface": exchange,
        "status": "CONTROL_ONLY_NO_RECIPROCAL_SOURCE",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(run_controls(), indent=2))
