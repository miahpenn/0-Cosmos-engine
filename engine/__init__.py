"""Canonical 0-Cosmos engine package.

The package root exposes the production-facing state graph and leaves the
older CMC/reference driver available only through its explicit module.
"""
from .cosmos import (
    CosmosParams, CosmosState, cosmos_rhs, cosmos_rhs_driven,
    construct_present_day, geometry_driven_step,
)
from .local import LocalParams, LocalState, local_rhs
from .interface import InterfaceState, exchange_from_flux
from .handoff import CycleLedger, HandoffPoint
from .production_kernel import ProductionState, V55ProductionKernel

__all__ = [
    "CosmosParams",
    "CosmosState",
    "construct_present_day",
    "cosmos_rhs",
    "cosmos_rhs_driven",
    "geometry_driven_step",
    "LocalParams",
    "LocalState",
    "local_rhs",
    "InterfaceState",
    "exchange_from_flux",
    "CycleLedger",
    "HandoffPoint",
    "ProductionState",
    "V55ProductionKernel",
]
