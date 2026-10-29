"""Product electron swarm solver comparison package.

Public product configs use schema v2 and canonical solver ids:
``two_term``, ``multi_term``, ``monte_carlo``, and ``propagator``.
"""

from .core.config import (
    CrossSectionFileConfig,
    CrossSectionsConfig,
    GasComponent,
    RequestedSolverConfig,
    SolverId,
    SwarmConfig,
)
from .core.config_parser import load_config
from .core.numerics import eedf_from_f0, eepf_from_eedf, widths_from_centers
from .core.results import (
    EnergyAngleDistributionResult,
    RateResult,
    SwarmCaseResult,
)

__all__ = [
    "GasComponent",
    "CrossSectionFileConfig",
    "CrossSectionsConfig",
    "EnergyAngleDistributionResult",
    "RateResult",
    "RequestedSolverConfig",
    "SolverId",
    "SwarmCaseResult",
    "SwarmConfig",
    "eedf_from_f0",
    "eepf_from_eedf",
    "load_config",
    "run",
    "run_from_config",
    "widths_from_centers",
]

__version__ = "0.12.0"


def __getattr__(name: str):
    if name in {"run", "run_from_config"}:
        from importlib import import_module

        value = getattr(import_module(".runner", __name__), name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
