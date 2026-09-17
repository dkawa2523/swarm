"""Typed contracts shared by the GEC-CCP mapping and execution layers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from swarm_workflow.comsol.input.eedf_audit import ComsolEedfAuditPlan
from swarm_workflow.comsol.runtime import ComsolExecutionSummary


AVOGADRO_PER_MOL = "6.02214076e23[1/mol]"
FUNCTION_EEDF_TABLE = "eedf_f0_comsol_2d.csv"
FUNCTION_EEDF_INTERPOLATION = "structured_spreadsheet_linear_projection"
FUNCTION_EEDF_PREINTEGRATED_INELASTIC = "function_eedf_preintegrated_inelastic"
GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE = "swarm_hybrid_einstein_de"
GEC_RESTRICTED_TRANSPORT_CLOSURE = "comsol_specify_all_restricted"
GEC_EXTERNAL_ELASTIC_ENERGY_LOSS = "external_solver_native"
GEC_COMSOL_ELASTIC_ENERGY_LOSS = "comsol_mratio"
GEC_RESULT_ROLES = {
    "physical_target",
    "diagnostic_control",
    "ablation",
    "historical_validation",
}
GEC_LOOKUP_INTERPOLATION = "differentiable_log_piecewise_cubic"


class GecCcpWorkflowError(RuntimeError):
    """Raised when the GEC CCP workflow cannot be prepared or completed."""


class GecCcpQualityError(GecCcpWorkflowError):
    """Raised when a completed solve fails physics-quality acceptance."""


@dataclass(frozen=True, slots=True)
class GecModelSpec:
    input_mph: Path
    baseline_output_mph: Path | None
    output_mph: Path
    component: str
    physics: str
    plasma_feature: str
    study: str
    time_periodic_feature: str
    external_study: str
    external_time_periodic_feature: str
    external_solution: str
    conversion_study: str
    expected_original_eedf: str


@dataclass(frozen=True, slots=True)
class GecReactionSpec:
    name: str
    feature: str
    process_type: str


@dataclass(frozen=True, slots=True)
class GecFunctionEedfSpec:
    """COMSOL two-argument Function-EEDF import contract.

    Every solver source stores the same rectangular physical ``f0(E, meanE)``
    tensor in one COMSOL 6.4-compatible spreadsheet representation.
    """

    table: str
    function_tag: str
    interpolation: Literal["structured_spreadsheet_linear_projection"]
    extrapolation: Literal["constant"]


@dataclass(frozen=True, slots=True)
class GecClosureSpec:
    source_field: Literal["steady_dc", "rf_time_periodic"]
    electron_transport: Literal[
        "comsol",
        "swarm_mobility_einstein",
        "swarm_hybrid_einstein_de",
        "comsol_specify_all_restricted",
    ]
    zero_field_isotropization_Td: float | None
    thermal_diffusion_model: Literal[
        "off_restricted_diagonal",
        "comsol_grad_diffusivity",
    ]
    gradient_response_policy: Literal[
        "standard_local_energy",
        "require_full",
    ]
    reaction_model: Literal[
        "comsol_eedf",
        "external_rates",
        "function_eedf",
        "function_eedf_preintegrated_inelastic",
    ]
    external_rate_processes: tuple[str, ...]
    elastic_energy_loss_model: Literal[
        "comsol_mratio",
        "external_solver_native",
    ]
    lookup_jacobian: Literal["exact"]
    interpolation: Literal["differentiable_log_piecewise_cubic"]
    function_eedf: GecFunctionEedfSpec | None


@dataclass(frozen=True, slots=True)
class GecRunSpec:
    power_parameter: str
    power_W: float
    include_builtin_reference: bool
    nonlinear_globalization: Literal[
        "model_defined",
        "automatic_newton_no_recovery",
        "double_dogleg",
    ]
    source_stabilization: bool
    reaction_source_stabilization: bool
    axis_dataset: str
    radial_dataset: str
    period_dataset: str | None
    external_period_dataset: str
    phase_dataset: str
    baseline_waveform_dataset: str | None
    external_waveform_dataset: str
    support_policy: Literal["strict", "qualified_guard_sensitivity"]
    low_energy_guard_bundle: Path | None
    validation_mean_energy_floor_eV: float | None


@dataclass(frozen=True, slots=True)
class GecPathSpec:
    path: Path


@dataclass(frozen=True, slots=True)
class GecBundleSpec:
    path: Path
    expected_source: str
    expected_field_type: Literal["dc"]
    expected_transport_definition: str
    expected_mc_transport_estimator_schema: str | None
    expected_mc_tail_estimator_schema: str | None
    expected_mc_qualification_profile: str | None
    expected_two_term_transport_kernel_schema: str | None
    require_solver_selection: bool = False


@dataclass(frozen=True, slots=True)
class GecResultSpec:
    output_directory: Path
    role: Literal[
        "physical_target",
        "diagnostic_control",
        "ablation",
        "historical_validation",
    ]


@dataclass(frozen=True, slots=True)
class GecCcpMapping:
    path: Path
    root: Path
    model: GecModelSpec
    bundle: GecBundleSpec
    reactions: tuple[GecReactionSpec, ...]
    closure: GecClosureSpec
    run: GecRunSpec
    results: GecResultSpec
    logs: GecPathSpec


@dataclass(frozen=True, slots=True)
class MphHeavySpeciesContract:
    tag: str
    species_type: str
    enabled: bool
    molar_mass_kg_mol: float
    from_mass_constraint: bool
    thermal_diffusion: float
    charge_number: float


@dataclass(frozen=True, slots=True)
class MphSurfaceReactionContract:
    tag: str
    enabled: bool
    formula: str


@dataclass(frozen=True, slots=True)
class MphContract:
    comsol_version: str | None
    physics_operation: str
    physics_tag: str
    original_eedf: str
    plasma_feature: str
    reaction_features: tuple[str, ...]
    studies: tuple[str, ...]
    datasets: tuple[str, ...]
    heavy_species_molar_masses_kg_mol: tuple[tuple[str, float], ...]
    common_heavy_species_molar_mass_kg_mol: float | None
    heavy_species: tuple[MphHeavySpeciesContract, ...]
    heavy_species_selection: str
    axisymmetric: bool
    heavy_species_formulation: str
    heavy_species_diffusion_model: str
    heavy_species_migration: bool
    heavy_species_convection: bool
    mixture_diffusion_correction: bool
    ion_tensor_properties: bool
    ion_electric_field_time_model: str
    surface_reactions: tuple[MphSurfaceReactionContract, ...]


@dataclass(frozen=True, slots=True)
class GecCcpPlan:
    mapping: GecCcpMapping
    contract: MphContract
    output_directory: Path
    plan_json: Path
    baseline_run_java: Path | None
    baseline_export_java: Path | None
    apply_java: Path
    external_run_java: Path
    external_export_java: Path
    native_eedf_audit: ComsolEedfAuditPlan | None
    expected_result_files: tuple[Path, ...]


@dataclass(frozen=True, slots=True)
class GecCcpRunSummary:
    plan: GecCcpPlan
    executions: tuple[ComsolExecutionSummary, ...]
    status_json: Path
