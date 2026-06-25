"""External reference result readers for product benchmarks."""

from .bolsig import load_bolsig_reference
from .common import (
    ExternalReferenceConfig,
    ReferenceCaseResult,
    load_external_reference_configs,
    load_reference_cases,
    parse_external_reference_configs,
)
from .mcig import load_mcig_reference

__all__ = [
    "ExternalReferenceConfig",
    "ReferenceCaseResult",
    "load_external_reference_configs",
    "load_bolsig_reference",
    "load_mcig_reference",
    "load_reference_cases",
    "parse_external_reference_configs",
]
