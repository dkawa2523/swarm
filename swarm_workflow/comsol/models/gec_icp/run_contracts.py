"""Execution contracts for the repository GEC-ICP COMSOL adapter."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

from swarm_workflow.comsol.input.eedf_audit import ComsolEedfAuditPlan
from swarm_workflow.comsol.runtime import ComsolExecutionSummary

from .contracts import GecIcpMapping, GecIcpMphContract


class GecIcpWorkflowError(RuntimeError):
    """Raised when an ICP run cannot be prepared or completed honestly."""


class GecIcpQualityError(GecIcpWorkflowError):
    """Raised when a completed ICP solve fails its acceptance gates."""


@dataclass(frozen=True, slots=True)
class GecIcpFunctionEedfSpec:
    table: str
    function_tag: str
    interpolation: Literal["structured_spreadsheet_linear_projection"]
    extrapolation: Literal["constant"]


@dataclass(frozen=True, slots=True)
class GecIcpBundleSpec:
    path: Path
    expected_source: Literal["two_term", "monte_carlo", "propagator", "composite"]
    expected_field_type: Literal["dc"]
    expected_transport_definition: str
    expected_mc_qualification_profile: str | None


@dataclass(frozen=True, slots=True)
class GecIcpClosureSpec:
    source_field: Literal["steady_dc"]
    electron_transport: Literal["swarm_mobility_comsol_einstein"]
    reaction_model: Literal["function_eedf"]
    elastic_energy_loss_model: Literal["comsol_cross_section_integral"]
    chemistry_owner: Literal["comsol_embedded_cross_sections"]
    rf_owner: Literal["comsol_frequency_transient"]
    function_eedf: GecIcpFunctionEedfSpec


@dataclass(frozen=True, slots=True)
class GecIcpConvergenceSpec:
    electron_inventory_relative_change: float
    ion_inventory_relative_change: float
    metastable_inventory_relative_change: float
    mean_energy_relative_change: float
    absorbed_power_relative_change: float
    coil_power_relative_error: float
    mean_energy_support_relative_guard: float


@dataclass(frozen=True, slots=True)
class GecIcpRunSpec:
    output_mph: Path
    power_W: float
    frequency_Hz: float
    gas_temperature_K: float
    pressure_Pa: float
    final_time_s: float
    output_points_per_decade: int
    clear_saved_solution: bool
    convergence: GecIcpConvergenceSpec


@dataclass(frozen=True, slots=True)
class GecIcpRunMapping:
    path: Path
    root: Path
    model_mapping_path: Path
    model_mapping: GecIcpMapping
    bundle: GecIcpBundleSpec
    closure: GecIcpClosureSpec
    run: GecIcpRunSpec
    output_directory: Path
    log_path: Path


@dataclass(frozen=True, slots=True)
class GecIcpClosureReadbackPlan:
    """Inputs and artifacts for the licensed post-apply closure audit."""

    mapping: GecIcpRunMapping
    java_path: Path
    report_path: Path
    transport_path: Path
    transport_sha256: str
    anchor_count: int


@dataclass(frozen=True, slots=True)
class GecIcpPlan:
    mapping: GecIcpRunMapping
    mph_contract: GecIcpMphContract
    bundle_evidence: dict[str, Any]
    output_directory: Path
    plan_json: Path
    closure_java: Path
    closure_support_java: tuple[Path, ...]
    closure_readback: GecIcpClosureReadbackPlan
    native_eedf_audit: ComsolEedfAuditPlan
    solve_java: Path
    expected_result_files: tuple[Path, ...]


@dataclass(frozen=True, slots=True)
class GecIcpRunSummary:
    plan: GecIcpPlan
    executions: tuple[ComsolExecutionSummary, ...]
    status_json: Path


__all__ = [
    "GecIcpBundleSpec",
    "GecIcpClosureReadbackPlan",
    "GecIcpClosureSpec",
    "GecIcpConvergenceSpec",
    "GecIcpFunctionEedfSpec",
    "GecIcpPlan",
    "GecIcpQualityError",
    "GecIcpRunMapping",
    "GecIcpRunSpec",
    "GecIcpRunSummary",
    "GecIcpWorkflowError",
]
