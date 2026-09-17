"""Model-independent closure selection and anchor-fallback planning."""

from .anchor_fallback import (
    ANCHOR_FALLBACK_SCHEMA,
    ANCHOR_FALLBACK_SOURCE,
    ANCHOR_FALLBACK_STAGE,
    COMPOSITE_QUALITY_COLUMNS,
    COMPOSITE_QUALITY_SCHEMA,
    AnchorFallbackError,
    AnchorFallbackSummary,
    build_mc_anchor_fallback_plan,
)
from .closure import (
    MC_QUALIFICATION_TABLE,
    PHYSICAL_CONTEXT_KEY,
    ClosureSelectionError,
    compatible_physics,
    contained_file,
    context_from_metadata,
    file_sha256,
    physical_context,
    read_json_object,
    read_selection,
    table_evidence,
    validate_selected_tables,
)

__all__ = (
    "ANCHOR_FALLBACK_SCHEMA",
    "ANCHOR_FALLBACK_SOURCE",
    "ANCHOR_FALLBACK_STAGE",
    "COMPOSITE_QUALITY_COLUMNS",
    "COMPOSITE_QUALITY_SCHEMA",
    "MC_QUALIFICATION_TABLE",
    "PHYSICAL_CONTEXT_KEY",
    "AnchorFallbackError",
    "AnchorFallbackSummary",
    "ClosureSelectionError",
    "build_mc_anchor_fallback_plan",
    "compatible_physics",
    "contained_file",
    "context_from_metadata",
    "file_sha256",
    "physical_context",
    "read_json_object",
    "read_selection",
    "table_evidence",
    "validate_selected_tables",
)
