"""Build provenance-bound canonical COMSOL input bundles."""

from .bundle import export_comsol_bundle
from .contracts import ComsolExportError, ComsolExportSummary

__all__ = (
    "ComsolExportError",
    "ComsolExportSummary",
    "export_comsol_bundle",
)
