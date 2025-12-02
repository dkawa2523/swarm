# Copyright (c) 2020-2021 ETH Zurich
"""
Lightweight caches used to speed up repeated simulations.
"""

from __future__ import annotations

from typing import Dict, Tuple

from swarm_mc.config import Config
from swarm_mc.gas_mixture import GasMixture


class GasMixtureCache:
    """
    Cache GasMixture instances keyed by the immutable parts of the configuration.

    Re-using the same GasMixture across sweeps avoids repeatedly parsing cross
    section files and constructing interpolations, which dominates wall-clock time
    for many swarm studies.
    """

    def __init__(self):
        self._cache: Dict[Tuple, GasMixture] = {}

    @staticmethod
    def _key_from_config(cfg: Config) -> Tuple:
        def _round_iter(values):
            return tuple(round(float(v), 12) for v in values)
        return (
            tuple(cfg.gases),
            tuple(cfg.paths_to_cross_section_files),
            _round_iter(cfg.fractions),
            round(float(cfg.max_cross_section_energy), 12),
        )

    def get(self, cfg: Config) -> GasMixture:
        """
        Retrieve a GasMixture for the given configuration, building and caching
        it on-demand.
        """
        key = self._key_from_config(cfg)
        if key not in self._cache:
            self._cache[key] = GasMixture(
                cfg.gases,
                cfg.paths_to_cross_section_files,
                cfg.fractions,
                cfg.max_cross_section_energy,
            )
        return self._cache[key]

    def clear(self) -> None:
        """Clear cached mixtures (useful for long-running interactive sessions)."""
        self._cache.clear()
