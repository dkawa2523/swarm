"""External reference result readers for product benchmarks."""

from .bolsig import load_bolsig_reference
from .common import ReferenceCaseResult, load_reference_cases
from .mcig import load_mcig_reference

__all__ = [
    "ReferenceCaseResult",
    "load_bolsig_reference",
    "load_mcig_reference",
    "load_reference_cases",
]
