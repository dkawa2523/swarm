"""Product electron swarm solver comparison package.

Public product configs use schema v2 and canonical solver ids:
``two_term``, ``multi_term``, and ``monte_carlo``.
"""

from .core.config import load_config, SwarmConfig
from .runner import run_from_config, run

__all__ = ["load_config", "SwarmConfig", "run_from_config", "run"]

__version__ = "0.2.0"
