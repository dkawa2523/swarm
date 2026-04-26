"""Collision-model extension helpers.

These modules intentionally stay outside solver internals.  They transform
loaded cross-section/process inputs into additional processes or metadata that
existing two-term, multi-term, and MC paths can consume.
"""

from .states import augment_cross_sections_from_config

__all__ = ["augment_cross_sections_from_config"]
