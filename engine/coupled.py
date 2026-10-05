"""Top-level coupled state and explicit exchange boundary."""
from dataclasses import dataclass
from .cosmos import CosmosParams,CosmosState
from .local import LocalParams,LocalState
from .interface import InterfaceState

@dataclass
class CoupledState:
    cosmos: CosmosState
    local: LocalState
    interface: InterfaceState

@dataclass(frozen=True)
class CoupledParams:
    cosmos: CosmosParams = CosmosParams()
    local: LocalParams = LocalParams()

def pack(state):
    return {"cosmos":state.cosmos.__dict__.copy(),"local":state.local.__dict__.copy(),"interface":state.interface.__dict__.copy()}
