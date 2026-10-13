from engine.cosmos import CosmosParams
from engine.local import LocalParams
from engine.interface import exchange_from_flux
from engine.diagnostics import friedmann_constraint

def test_physics_locks():
    p=LocalParams()
    assert p.kappa==0.0 and p.g==1.0 and p.G==1.0
    assert p.alpha_D==1.0

def test_corrected_v0():
    assert CosmosParams().V0>0 and CosmosParams().V0!=1.0

def test_exchange_zero_sum():
    x=exchange_from_flux(3.0)
    assert x["net"]==0.0

def test_friedmann_control():
    assert friedmann_constraint(1.0,3.0)==0.0
