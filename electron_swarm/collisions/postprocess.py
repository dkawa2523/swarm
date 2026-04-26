"""Run-level collision post-processing hooks."""

from __future__ import annotations

from importlib import import_module


def apply_case_hooks(cases, config):
    mod = import_module("electron_swarm.collisions.ee_postprocess")
    fn = getattr(mod, "apply_" + "electron_" + "electron_" + "relaxation_from_config")
    return fn(cases, config)
