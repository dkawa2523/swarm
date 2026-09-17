"""Shared contracts for reproducible GEC CCP plots."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from swarm_workflow.comsol.models.gec_ccp.contracts import GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE


__all__ = [
    "CANONICAL_RESTRICTED_TRANSPORT_CLOSURE",
    "COMSOL_DIFFERENCE_COLORS",
    "COMSOL_RAINBOW_LIGHT_COLORS",
    "CURRENT_SOLVER_COMPARISON_CLOSURES",
    "FUNCTION_EEDF_REACTION_MODELS",
    "GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE",
    "GROUNDED_BOUNDARY_COLOR",
    "GecCcpPlotError",
    "GecCcpPlotSummary",
    "HYBRID_FUNCTION_EEDF_REACTION_MODEL",
    "POWERED_BOUNDARY_COLOR",
]

POWERED_BOUNDARY_COLOR = "#9b3c7d"
GROUNDED_BOUNDARY_COLOR = "#65752e"
COMSOL_RAINBOW_LIGHT_COLORS = (
    "#3b4cc0",
    "#4f8bd6",
    "#55b9c5",
    "#72c878",
    "#d6d65c",
    "#eda24c",
    "#c94741",
)
COMSOL_DIFFERENCE_COLORS = ("#3b67b0", "#f7f7f7", "#c94741")
FUNCTION_EEDF_REACTION_MODELS = frozenset(
    {"function_eedf", "function_eedf_preintegrated_inelastic"}
)
HYBRID_FUNCTION_EEDF_REACTION_MODEL = (
    "function_eedf_preintegrated_inelastic"
)
CANONICAL_RESTRICTED_TRANSPORT_CLOSURE = "comsol_specify_all_restricted"
CURRENT_SOLVER_COMPARISON_CLOSURES = {
    "two_term": {
        "electron_transport": "swarm_mobility_einstein",
        "reaction_model": "function_eedf",
    },
    "monte_carlo": {
        "electron_transport": "swarm_mobility_einstein",
        "reaction_model": "function_eedf_preintegrated_inelastic",
    },
}


# Chart contract:
# - Coefficient plots answer which DC-swarm closure COMSOL receives.
# - Spatial plots compare the built-in reference when it was explicitly run.
# - The spatial-cut overview shows four period-averaged fields along both the
#   axial centerline and radial midplane on linear axes; the phase map shows
#   three fields at 20 axial nodes x 51 RF phases for the accepted Swarm-table
#   solution.
# - Full-domain maps show the GEC-CCP axisymmetric gas/plasma region and its
#   powered-electrode, grounded-electrode, chamber-wall, and symmetry-axis
#   boundaries.  Built-in and Swarm panels share row-wise color limits.
# - Full-domain surface plots approximate COMSOL 6.x RainbowLight styling:
#   white graphics background, equal axis scale, thin neutral boundaries,
#   right-side legends, and less-saturated blue-cyan-green-yellow-red colors.
#   Compared maps remain explicitly labeled, use linear color legends, and
#   share limits.
# - PNG and SVG are the final QA surfaces.
# - Waveform plots compare electrode voltage/current over an RF period.
# - COMSOL plots are skipped, never fabricated, when result CSVs are absent.


class GecCcpPlotError(RuntimeError):
    """Raised when GEC CCP plots cannot be produced from supplied data."""


@dataclass(frozen=True, slots=True)
class GecCcpPlotSummary:
    output_directory: Path
    figures: tuple[Path, ...]
    manifest: Path
    skipped: tuple[str, ...]
