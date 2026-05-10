"""Unified electron swarm calculation package.

This package adds a common configuration/results layer and multiple solver slots:

* particle Monte Carlo: delegated to an existing repository implementation through
  a command/API adapter;
* Boltzmann two-term approximation: production native BOLSIG-like
  energy-space solver with an optional BOLOS reference backend.
* multi-term Boltzmann entry point: moment-closure by default, lmax=1 as a
  two-term reference adapter, and lmax>1 as a reference-anchored closure for
  integral-cross-section inputs.

The code is intentionally additive so it can be copied into an existing swarm
repository without forcing a rewrite of the current Monte Carlo code.
"""

from .core.config import load_config, SwarmConfig
from .runner import run_from_config, run

__all__ = ["load_config", "SwarmConfig", "run_from_config", "run"]

__version__ = "0.2.0"
