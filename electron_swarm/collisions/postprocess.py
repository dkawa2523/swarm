"""Run-level collision post-processing hooks."""

from __future__ import annotations

from electron_swarm.collisions.ee_postprocess import (
    apply_electron_electron_relaxation_from_config,
)


def apply_case_hooks(cases, config, cross_sections=None):
    return apply_electron_electron_relaxation_from_config(
        cases,
        config,
        cross_sections,
    )
