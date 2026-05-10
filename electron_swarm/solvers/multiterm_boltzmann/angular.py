"""Legendre basis helpers for the multi-term Boltzmann solver."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class LegendreBasis:
    lmax: int

    def __post_init__(self) -> None:
        if self.lmax < 1:
            raise ValueError("lmax must be >= 1")

    @property
    def n_terms(self) -> int:
        return self.lmax + 1

