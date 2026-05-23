"""Product YAML configuration model and validation.

Product mode uses schema v2 only. Legacy v1 runner-mode YAML is rejected at
load time with an explicit migration message. Existing numerical solvers are
fed through internal implementation configs derived from the product schema.
"""

from __future__ import annotations

from dataclasses import dataclass, field as dc_field
from pathlib import Path
from typing import Any, Literal, cast

import yaml

SolverId = Literal["two_term", "multi_term", "monte_carlo"]
UnsupportedPolicy = Literal["fail", "skip_solver", "approximate"]
DegradedPolicy = Literal["fail", "warn", "record_only"]

CANONICAL_SOLVER_IDS: tuple[str, ...] = ("two_term", "multi_term", "monte_carlo")
OBSOLETE_PUBLIC_NAMES = {"both", "all", "boltzmann_two_term", "multiterm_boltzmann"}
OBSOLETE_TOP_LEVEL_SOLVER_KEYS = {
    "boltzmann_two_term",
    "multiterm_boltzmann",
    "monte_carlo",
}
OUTPUT_FIELDS = {"directory", "base_name", "float_format"}
TWO_TERM_SOLVER_FIELDS = {
    "backend",
    "nonconservative_model",
    "min_momentum_cross_section_m2",
}
MULTI_TERM_SOLVER_FIELDS = {
    "formulation",
    "method",
    "lmax",
    "lmax_convergence_tolerance",
    "field_coupling_scale",
    "dense_threshold",
}
MONTE_CARLO_SOLVER_FIELDS = {
    "particles",
    "max_collisions",
    "timeout_s",
    "command",
    "python_api",
    "working_directory",
    "environment",
    "output_summary_csv",
    "output_eedf_csv",
    "passthrough",
}
OBSOLETE_OUTPUT_FIELDS = {
    "compatibility",
    "write_summary",
    "write_eedf",
    "write_rates",
    "write_comparison",
    "write_plots",
}
MIGRATION_ERROR = (
    "schema v2 is required; use schema_version: 2 and run.solvers with "
    "two_term, multi_term, or monte_carlo"
)

BoltzmannBackend = Literal["auto", "native_sg", "native_bolsig", "internal", "bolos"]
InternalBoltzmannBackend = Literal["auto", "native_bolsig", "internal", "bolos"]


@dataclass(slots=True)
class GasComponent:
    species: str
    fraction: float
    mass_amu: float


@dataclass(slots=True)
class ConditionsConfig:
    gas_temperature_K: float = 300.0
    pressure_Pa: float | None = None
    gas_number_density_m3: float | None = None
    length_scale_m: float | None = None
    gas_mixture: list[GasComponent] = dc_field(default_factory=list)


@dataclass(slots=True)
class CrossSectionFileConfig:
    path: Path
    species: str | None = None
    format: str | None = None


@dataclass(slots=True)
class CrossSectionsConfig:
    format: str = "csv"
    files: list[CrossSectionFileConfig] = dc_field(default_factory=list)
    high_energy_extrapolation: Literal["zero", "hold", "error"] = "zero"


@dataclass(slots=True)
class EnergyGridRefinementConfig:
    enabled: bool = False
    threshold_padding_eV: float = 0.15
    points_per_threshold: int = 8
    max_extra_points: int = 120


@dataclass(slots=True)
class EnergyGridConfig:
    min_eV: float = 1.0e-4
    max_eV: float = 100.0
    n: int = 600
    spacing: Literal["linear", "quadratic", "log"] = "quadratic"
    refine: EnergyGridRefinementConfig = dc_field(
        default_factory=EnergyGridRefinementConfig
    )


@dataclass(slots=True)
class MultiTermEnergyGridConfig:
    min_eV: float = 1.0e-3
    max_eV: float = 200.0
    n: int = 220
    spacing: Literal["linear", "log", "log_linear"] = "log_linear"
    linear_until_eV: float = 2.0
    refine: EnergyGridRefinementConfig = dc_field(
        default_factory=EnergyGridRefinementConfig
    )


@dataclass(slots=True)
class AdaptiveGridConfig:
    enabled: bool = True
    max_cycles: int = 4
    mean_energy_multiplier: float = 15.0
    min_max_eV: float = 20.0
    max_max_eV: float = 2000.0
    tail_probability: float = 1.0e-8
    tail_cells_fraction: float = 0.05
    edge_to_peak: float = 1.0e-10


@dataclass(slots=True)
class ConvergenceConfig:
    max_iterations: int = 120
    tolerance: float = 1.0e-8
    eigenvalue_tolerance: float = 1.0e-8
    residual_tolerance: float = 1.0e-7
    relaxation: float = 0.7
    clip_negative: bool = True


@dataclass(slots=True)
class BoltzmannTwoTermConfig:
    enabled: bool = True
    backend: InternalBoltzmannBackend = "native_bolsig"
    energy_grid: EnergyGridConfig = dc_field(default_factory=EnergyGridConfig)
    adaptive_grid: AdaptiveGridConfig = dc_field(default_factory=AdaptiveGridConfig)
    convergence: ConvergenceConfig = dc_field(default_factory=ConvergenceConfig)
    initial_electron_temperature_eV: float = 2.0
    nonconservative_model: Literal["growth", "ignore"] = "growth"
    ionization_energy_sharing: Literal[
        "equal", "primary_secondary", "loss_only"
    ] = "equal"
    secondary_electron_energy_eV: float = 0.0
    min_momentum_cross_section_m2: float = 1.0e-24


@dataclass(slots=True)
class MultiTermBoltzmannConfig:
    enabled: bool = True
    lmax: int = 4
    method: Literal["moment_closure"] = "moment_closure"
    hydrodynamic: bool = False
    dense_threshold: int = 280
    eedf_shape: Literal["maxwellian", "druyvesteyn"] = "maxwellian"
    lmax_convergence_tolerance: float = 0.03
    field_coupling_scale: float = 1.0
    product_method: Literal[
        "pn_closure_surrogate", "pn_closure_direct", "pn_dcs"
    ] = "pn_closure_surrogate"
    formulation: Literal["PN"] = "PN"
    energy_grid: MultiTermEnergyGridConfig = dc_field(
        default_factory=MultiTermEnergyGridConfig
    )


@dataclass(slots=True)
class MonteCarloAdapterConfig:
    enabled: bool = True
    command: str | None = None
    python_api: str | None = None
    working_directory: Path | None = None
    timeout_s: float | None = None
    environment: dict[str, str] = dc_field(default_factory=dict)
    output_summary_csv: Path | None = None
    output_eedf_csv: Path | None = None
    passthrough: dict[str, Any] = dc_field(default_factory=dict)
    particles: int | None = None
    max_collisions: int | None = None


@dataclass(slots=True)
class InternalSolverConfigs:
    two_term: BoltzmannTwoTermConfig = dc_field(default_factory=BoltzmannTwoTermConfig)
    multi_term: MultiTermBoltzmannConfig = dc_field(
        default_factory=MultiTermBoltzmannConfig
    )
    monte_carlo: MonteCarloAdapterConfig = dc_field(
        default_factory=MonteCarloAdapterConfig
    )


@dataclass(slots=True)
class TwoTermProductConfig:
    backend: BoltzmannBackend = "native_sg"
    nonconservative_model: Literal["growth", "ignore"] = "growth"
    min_momentum_cross_section_m2: float = 1.0e-24


@dataclass(slots=True)
class MultiTermProductConfig:
    formulation: Literal["PN"] = "PN"
    method: Literal[
        "pn_closure_surrogate", "pn_closure_direct", "pn_dcs"
    ] = "pn_closure_surrogate"
    lmax: int = 4
    lmax_convergence_tolerance: float = 0.03
    field_coupling_scale: float = 1.0
    dense_threshold: int = 280


@dataclass(slots=True)
class MonteCarloProductConfig:
    particles: int | None = None
    max_collisions: int | None = None
    timeout_s: float | None = None
    command: str | None = None
    python_api: str | None = None
    working_directory: Path | None = None
    environment: dict[str, str] = dc_field(default_factory=dict)
    output_summary_csv: Path | None = None
    output_eedf_csv: Path | None = None
    passthrough: dict[str, Any] = dc_field(default_factory=dict)


@dataclass(slots=True)
class SolversConfig:
    two_term: TwoTermProductConfig = dc_field(default_factory=TwoTermProductConfig)
    multi_term: MultiTermProductConfig = dc_field(default_factory=MultiTermProductConfig)
    monte_carlo: MonteCarloProductConfig = dc_field(
        default_factory=MonteCarloProductConfig
    )


@dataclass(slots=True)
class RequestedSolverConfig:
    id: SolverId
    enabled: bool = True
    label: str | None = None


@dataclass(slots=True)
class RunConfig:
    solvers: list[RequestedSolverConfig] = dc_field(default_factory=list)
    e_over_n_Td: list[float] = dc_field(default_factory=lambda: [100.0])
    case_prefix: str = "case"


@dataclass(slots=True)
class MagneticFieldConfig:
    enabled: bool = False
    B_T: float = 0.0
    angle_EB_deg: float = 0.0


@dataclass(slots=True)
class FieldConfig:
    type: Literal["dc", "rf", "time_dependent"] = "dc"
    magnetic_field: MagneticFieldConfig = dc_field(default_factory=MagneticFieldConfig)


@dataclass(slots=True)
class AngularScatteringConfig:
    model: Literal["isotropic", "momentum_power", "maxent_p1"] = "isotropic"
    higher_moment_closure: Literal["zero", "power", "maxent"] = "zero"


@dataclass(slots=True)
class ElectronElectronConfig:
    enabled: bool = False
    model: str = "none"
    relaxation_fraction: float = 0.05
    conserve_mean_energy: bool = True
    fallback_temperature_eV: float = 2.0


@dataclass(slots=True)
class IonizationConfig:
    energy_sharing: Literal["equal", "primary_secondary", "loss_only"] = "equal"
    secondary_electron_energy_eV: float = 0.0


@dataclass(slots=True)
class EnergyGridPolicyConfig:
    adaptive: bool = True
    threshold_refinement: bool = True
    tail_probability_target: float = 1.0e-8
    tail_rate_fraction_target: float = 1.0e-3
    tail_metrics: bool = True
    tail_threshold_eV: float | None = None
    tail_rate_warning_fraction: float = 0.05
    max_eV_limit: float = 2000.0


@dataclass(slots=True)
class FiniteKConfig:
    enabled: bool = False
    k_m_inv: float | None = None


@dataclass(slots=True)
class PhysicsConfig:
    field: FieldConfig = dc_field(default_factory=FieldConfig)
    angular_scattering: AngularScatteringConfig = dc_field(
        default_factory=AngularScatteringConfig
    )
    electron_electron: ElectronElectronConfig = dc_field(
        default_factory=ElectronElectronConfig
    )
    ionization: IonizationConfig = dc_field(default_factory=IonizationConfig)
    energy_grid_policy: EnergyGridPolicyConfig = dc_field(
        default_factory=EnergyGridPolicyConfig
    )
    finite_k: FiniteKConfig = dc_field(default_factory=FiniteKConfig)


@dataclass(slots=True)
class ComparisonConfig:
    enabled: bool = False
    reference_solver: SolverId | None = None
    candidate_solvers: list[SolverId] = dc_field(default_factory=list)
    compare_eedf: bool = True
    required: bool = False


@dataclass(slots=True)
class FeaturePolicyConfig:
    unsupported: UnsupportedPolicy = "fail"
    degraded: DegradedPolicy = "warn"
    allow_unsupported_fallback: bool = False


@dataclass(slots=True)
class OutputConfig:
    directory: Path = Path("outputs")
    base_name: str = "swarm"
    float_format: str = "%.10e"


@dataclass(slots=True)
class SwarmConfig:
    schema_version: int = 2
    run: RunConfig = dc_field(default_factory=RunConfig)
    conditions: ConditionsConfig = dc_field(default_factory=ConditionsConfig)
    cross_sections: CrossSectionsConfig = dc_field(default_factory=CrossSectionsConfig)
    physics: PhysicsConfig = dc_field(default_factory=PhysicsConfig)
    solvers: SolversConfig = dc_field(default_factory=SolversConfig)
    comparison: ComparisonConfig = dc_field(default_factory=ComparisonConfig)
    feature_policy: FeaturePolicyConfig = dc_field(default_factory=FeaturePolicyConfig)
    output: OutputConfig = dc_field(default_factory=OutputConfig)
    source_path: Path | None = None
    internal: InternalSolverConfigs = dc_field(default_factory=InternalSolverConfigs)


def _as_path(value: str | Path | None, base: Path | None) -> Path | None:
    if value is None:
        return None
    path = Path(value)
    if not path.is_absolute() and base is not None:
        path = (base / path).resolve()
    return path


def _list_float(value: Any) -> list[float]:
    if isinstance(value, (int, float)):
        return [float(value)]
    if isinstance(value, list):
        return [float(v) for v in value]
    raise TypeError("run.e_over_n_Td must be a number or list of numbers")


def _read_mapping(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as fp:
        data = yaml.safe_load(fp) or {}
    if not isinstance(data, dict):
        raise ValueError("YAML root must be a mapping")
    return data


def _validate_literal(value: str, allowed: set[str], field_name: str) -> str:
    if value not in allowed:
        raise ValueError(f"{field_name} must be one of {sorted(allowed)}")
    return value


def _reject_unknown_fields(raw: dict[str, Any], allowed: set[str], field_name: str) -> None:
    unknown = set(raw) - allowed
    if unknown:
        raise ValueError(f"Unsupported {field_name} fields: {sorted(unknown)}")


def _as_mapping_section(raw: dict[str, Any], key: str, field_name: str) -> dict[str, Any]:
    value = raw.get(key, {}) or {}
    if not isinstance(value, dict):
        raise ValueError(f"{field_name} must be a mapping")
    return value


def _bool_field(raw: dict[str, Any], key: str, default: bool, field_name: str) -> bool:
    value = raw[key] if key in raw else default
    if not isinstance(value, bool):
        raise ValueError(f"{field_name} must be a boolean")
    return value


def _validate_schema(raw: dict[str, Any]) -> None:
    if raw.get("schema_version") != 2:
        raise ValueError(MIGRATION_ERROR)
    if "run" in raw and isinstance(raw["run"], dict) and "mode" in raw["run"]:
        raise ValueError(MIGRATION_ERROR)
    present = OBSOLETE_TOP_LEVEL_SOLVER_KEYS.intersection(raw)
    if present:
        raise ValueError(
            f"{MIGRATION_ERROR}; move {sorted(present)} under solvers.*"
        )
    out_raw = raw.get("output", {}) or {}
    if isinstance(out_raw, dict):
        unknown = set(out_raw) - OUTPUT_FIELDS
        if unknown:
            obsolete = sorted(unknown & OBSOLETE_OUTPUT_FIELDS)
            if obsolete:
                raise ValueError(
                    "schema v2 writes fixed canonical outputs; remove output fields "
                    f"{obsolete}"
                )
            raise ValueError(f"Unsupported output fields: {sorted(unknown)}")


def _canonical_solver(value: Any, field_name: str) -> SolverId:
    solver = str(value)
    if solver in OBSOLETE_PUBLIC_NAMES:
        raise ValueError(
            f"{field_name} uses obsolete solver id {solver!r}; use two_term, "
            "multi_term, or monte_carlo"
        )
    return cast(
        SolverId,
        _validate_literal(solver, set(CANONICAL_SOLVER_IDS), field_name),
    )


def _grid_refinement_config(
    raw: dict[str, Any], field_name: str
) -> EnergyGridRefinementConfig:
    value = raw.get("refine", {}) or {}
    if isinstance(value, bool):
        value = {"enabled": value}
    if not isinstance(value, dict):
        raise TypeError(f"{field_name}.refine must be a mapping or boolean")
    cfg = EnergyGridRefinementConfig(
        enabled=_bool_field(value, "enabled", False, f"{field_name}.refine.enabled"),
        threshold_padding_eV=float(value.get("threshold_padding_eV", 0.15)),
        points_per_threshold=int(value.get("points_per_threshold", 8)),
        max_extra_points=int(value.get("max_extra_points", 120)),
    )
    if cfg.threshold_padding_eV < 0.0:
        raise ValueError(f"{field_name}.refine.threshold_padding_eV must be >= 0")
    if cfg.points_per_threshold < 2:
        raise ValueError(f"{field_name}.refine.points_per_threshold must be >= 2")
    if cfg.max_extra_points < 0:
        raise ValueError(f"{field_name}.refine.max_extra_points must be >= 0")
    return cfg


def _parse_run(raw: dict[str, Any]) -> RunConfig:
    run_raw = raw.get("run", {}) or {}
    solvers_raw = run_raw.get("solvers")
    if not isinstance(solvers_raw, list) or not solvers_raw:
        raise ValueError("schema v2 requires non-empty run.solvers")
    solvers: list[RequestedSolverConfig] = []
    for index, item in enumerate(solvers_raw):
        if isinstance(item, str):
            item = {"id": item}
        if not isinstance(item, dict) or "id" not in item:
            raise ValueError("run.solvers entries must be mappings with id")
        solvers.append(
            RequestedSolverConfig(
                id=_canonical_solver(item["id"], f"run.solvers[{index}].id"),
                enabled=_bool_field(item, "enabled", True, f"run.solvers[{index}].enabled"),
                label=str(item["label"]) if item.get("label") is not None else None,
            )
        )
    run = RunConfig(
        solvers=solvers,
        e_over_n_Td=_list_float(
            run_raw.get("e_over_n_Td", run_raw.get("E_over_N_Td", [100.0]))
        ),
        case_prefix=str(run_raw.get("case_prefix", "case")),
    )
    if any(v <= 0 for v in run.e_over_n_Td):
        raise ValueError("All run.e_over_n_Td values must be positive")
    return run


def _parse_conditions(raw: dict[str, Any]) -> ConditionsConfig:
    cond_raw = raw.get("conditions", {}) or {}
    gas_mixture = [
        GasComponent(
            species=str(item["species"]),
            fraction=float(item.get("fraction", 1.0)),
            mass_amu=float(item.get("mass_amu", item.get("mass", 39.948))),
        )
        for item in (cond_raw.get("gas_mixture", []) or [])
    ]
    if not gas_mixture:
        species = str(cond_raw.get("species", "Ar"))
        gas_mixture = [
            GasComponent(
                species=species,
                fraction=1.0,
                mass_amu=float(cond_raw.get("mass_amu", 39.948)),
            )
        ]
    fraction_sum = sum(g.fraction for g in gas_mixture)
    if fraction_sum <= 0.0:
        raise ValueError(
            "conditions.gas_mixture fractions must sum to a positive value"
        )
    gas_mixture = [
        GasComponent(g.species, g.fraction / fraction_sum, g.mass_amu)
        for g in gas_mixture
    ]
    conditions = ConditionsConfig(
        gas_temperature_K=float(cond_raw.get("gas_temperature_K", 300.0)),
        pressure_Pa=(
            float(cond_raw["pressure_Pa"])
            if cond_raw.get("pressure_Pa") is not None
            else None
        ),
        gas_number_density_m3=(
            float(cond_raw["gas_number_density_m3"])
            if cond_raw.get("gas_number_density_m3") is not None
            else None
        ),
        length_scale_m=(
            float(cond_raw["length_scale_m"])
            if cond_raw.get("length_scale_m") is not None
            else None
        ),
        gas_mixture=gas_mixture,
    )
    if conditions.pressure_Pa is None and conditions.gas_number_density_m3 is None:
        raise ValueError(
            "Either conditions.pressure_Pa or conditions.gas_number_density_m3 must be set"
        )
    if conditions.length_scale_m is not None and conditions.length_scale_m <= 0.0:
        raise ValueError("conditions.length_scale_m must be positive when set")
    return conditions


def _parse_cross_sections(
    raw: dict[str, Any], base: Path
) -> CrossSectionsConfig:
    xs_raw = raw.get("cross_sections", {}) or {}
    xs_format = str(xs_raw.get("format", "csv"))
    xs_files: list[CrossSectionFileConfig] = []
    for item in xs_raw.get("files", []):
        xs_files.append(
            CrossSectionFileConfig(
                path=_as_path(item["path"], base) or Path(item["path"]),
                species=item.get("species"),
                format=item.get("format", xs_format),
            )
        )
    if not xs_files:
        raise ValueError(
            "cross_sections.files must contain at least one cross-section file"
        )
    high_energy_extrapolation = cast(
        Literal["zero", "hold", "error"],
        _validate_literal(
            str(xs_raw.get("high_energy_extrapolation", "zero")),
            {"zero", "hold", "error"},
            "cross_sections.high_energy_extrapolation",
        ),
    )
    return CrossSectionsConfig(
        format=xs_format,
        files=xs_files,
        high_energy_extrapolation=high_energy_extrapolation,
    )


def _parse_physics(raw: dict[str, Any], base: Path) -> PhysicsConfig:
    phys_raw = _as_mapping_section(raw, "physics", "physics")
    physics_unknown = set(phys_raw) - {
        "field",
        "angular_scattering",
        "electron_electron",
        "ionization",
        "energy_grid_policy",
        "finite_k",
    }
    if physics_unknown:
        raise ValueError(f"Unsupported physics fields: {sorted(physics_unknown)}")
    field_raw = _as_mapping_section(phys_raw, "field", "physics.field")
    magnetic_raw = _as_mapping_section(
        field_raw, "magnetic_field", "physics.field.magnetic_field"
    )
    ee_raw = _as_mapping_section(
        phys_raw, "electron_electron", "physics.electron_electron"
    )
    angular_raw = _as_mapping_section(
        phys_raw, "angular_scattering", "physics.angular_scattering"
    )
    ion_raw = _as_mapping_section(phys_raw, "ionization", "physics.ionization")
    grid_raw = _as_mapping_section(
        phys_raw, "energy_grid_policy", "physics.energy_grid_policy"
    )
    finite_k_raw = _as_mapping_section(phys_raw, "finite_k", "physics.finite_k")
    _reject_unknown_fields(
        field_raw, {"type", "magnetic_field"}, "physics.field"
    )
    _reject_unknown_fields(
        magnetic_raw,
        {"enabled", "B_T", "angle_EB_deg"},
        "physics.field.magnetic_field",
    )
    _reject_unknown_fields(
        ion_raw,
        {"energy_sharing", "secondary_electron_energy_eV"},
        "physics.ionization",
    )
    _reject_unknown_fields(finite_k_raw, {"enabled", "k_m_inv"}, "physics.finite_k")
    grid_unknown = set(grid_raw) - {
        "adaptive",
        "threshold_refinement",
        "tail_probability_target",
        "tail_rate_fraction_target",
        "tail_metrics",
        "tail_threshold_eV",
        "tail_rate_warning_fraction",
        "max_eV_limit",
    }
    if grid_unknown:
        raise ValueError(
            "Unsupported physics.energy_grid_policy fields: "
            f"{sorted(grid_unknown)}"
        )
    ee_unknown = set(ee_raw) - {
        "enabled",
        "model",
        "relaxation_fraction",
        "conserve_mean_energy",
        "fallback_temperature_eV",
    }
    if ee_unknown:
        raise ValueError(
            "Unsupported physics.electron_electron fields: "
            f"{sorted(ee_unknown)}"
        )
    angular_unknown = set(angular_raw) - {"model", "higher_moment_closure"}
    if angular_unknown:
        raise ValueError(
            "Unsupported physics.angular_scattering fields: "
            f"{sorted(angular_unknown)}"
        )
    angular_name = str(angular_raw.get("model", "isotropic")).lower()
    angular_closure_defaults = {
        "isotropic": "zero",
        "momentum_power": "power",
        "maxent_p1": "maxent",
    }
    if angular_name not in angular_closure_defaults:
        raise ValueError(
            "physics.angular_scattering.model must be isotropic, momentum_power, or maxent_p1"
        )
    closure_raw = angular_raw.get(
        "higher_moment_closure", angular_closure_defaults[angular_name]
    )
    higher_moment_closure = cast(
        Literal["zero", "power", "maxent"],
        _validate_literal(
            str(closure_raw).lower(),
            {"zero", "power", "maxent"},
            "physics.angular_scattering.higher_moment_closure",
        ),
    )
    if higher_moment_closure != angular_closure_defaults[angular_name]:
        raise ValueError(
            "physics.angular_scattering model and higher_moment_closure mismatch"
        )
    ee_enabled = _bool_field(ee_raw, "enabled", False, "physics.electron_electron.enabled")
    ee_model = str(ee_raw.get("model", "none")).lower()
    if ee_enabled and ee_model != "relaxation_postprocess":
        raise ValueError(
            "physics.electron_electron.model must be relaxation_postprocess when enabled"
        )
    if not ee_enabled and ee_model != "none":
        raise ValueError(
            "physics.electron_electron.model must be none when electron_electron is disabled"
        )
    ee_relaxation_fraction = float(ee_raw.get("relaxation_fraction", 0.05))
    if not (0.0 <= ee_relaxation_fraction <= 1.0):
        raise ValueError("physics.electron_electron.relaxation_fraction must be in [0, 1]")
    ee_fallback_temperature = float(ee_raw.get("fallback_temperature_eV", 2.0))
    if ee_fallback_temperature <= 0.0:
        raise ValueError(
            "physics.electron_electron.fallback_temperature_eV must be positive"
        )
    tail_threshold = (
        float(grid_raw["tail_threshold_eV"])
        if grid_raw.get("tail_threshold_eV") is not None
        else None
    )
    if tail_threshold is not None and tail_threshold < 0.0:
        raise ValueError("physics.energy_grid_policy.tail_threshold_eV must be >= 0")
    tail_rate_warning_fraction = float(
        grid_raw.get("tail_rate_warning_fraction", 0.05)
    )
    if not (0.0 <= tail_rate_warning_fraction <= 1.0):
        raise ValueError(
            "physics.energy_grid_policy.tail_rate_warning_fraction must be in [0, 1]"
        )
    tail_metrics_raw = grid_raw.get("tail_metrics", True)
    if not isinstance(tail_metrics_raw, bool):
        raise ValueError("physics.energy_grid_policy.tail_metrics must be a boolean")
    return PhysicsConfig(
        field=FieldConfig(
            type=cast(
                Literal["dc", "rf", "time_dependent"],
                _validate_literal(
                    str(field_raw.get("type", "dc")),
                    {"dc", "rf", "time_dependent"},
                    "physics.field.type",
                ),
            ),
            magnetic_field=MagneticFieldConfig(
                enabled=_bool_field(
                    magnetic_raw,
                    "enabled",
                    False,
                    "physics.field.magnetic_field.enabled",
                ),
                B_T=float(magnetic_raw.get("B_T", 0.0)),
                angle_EB_deg=float(magnetic_raw.get("angle_EB_deg", 0.0)),
            ),
        ),
        angular_scattering=AngularScatteringConfig(
            model=cast(
                Literal["isotropic", "momentum_power", "maxent_p1"], angular_name
            ),
            higher_moment_closure=higher_moment_closure,
        ),
        electron_electron=ElectronElectronConfig(
            enabled=ee_enabled,
            model=ee_model,
            relaxation_fraction=ee_relaxation_fraction,
            conserve_mean_energy=_bool_field(
                ee_raw,
                "conserve_mean_energy",
                True,
                "physics.electron_electron.conserve_mean_energy",
            ),
            fallback_temperature_eV=ee_fallback_temperature,
        ),
        ionization=IonizationConfig(
            energy_sharing=cast(
                Literal["equal", "primary_secondary", "loss_only"],
                _validate_literal(
                    str(ion_raw.get("energy_sharing", "equal")),
                    {"equal", "primary_secondary", "loss_only"},
                    "physics.ionization.energy_sharing",
                ),
            ),
            secondary_electron_energy_eV=float(
                ion_raw.get("secondary_electron_energy_eV", 0.0)
            ),
        ),
        energy_grid_policy=EnergyGridPolicyConfig(
            adaptive=_bool_field(
                grid_raw, "adaptive", True, "physics.energy_grid_policy.adaptive"
            ),
            threshold_refinement=_bool_field(
                grid_raw,
                "threshold_refinement",
                True,
                "physics.energy_grid_policy.threshold_refinement",
            ),
            tail_probability_target=float(
                grid_raw.get("tail_probability_target", 1.0e-8)
            ),
            tail_rate_fraction_target=float(
                grid_raw.get("tail_rate_fraction_target", 1.0e-3)
            ),
            tail_metrics=tail_metrics_raw,
            tail_threshold_eV=tail_threshold,
            tail_rate_warning_fraction=tail_rate_warning_fraction,
            max_eV_limit=float(grid_raw.get("max_eV_limit", 2000.0)),
        ),
        finite_k=FiniteKConfig(
            enabled=_bool_field(finite_k_raw, "enabled", False, "physics.finite_k.enabled"),
            k_m_inv=(
                float(finite_k_raw["k_m_inv"])
                if finite_k_raw.get("k_m_inv") is not None
                else None
            ),
        ),
    )


def _parse_solvers(
    raw: dict[str, Any], base: Path, physics: PhysicsConfig
) -> tuple[SolversConfig, InternalSolverConfigs]:
    solvers_raw = raw.get("solvers", {}) or {}
    if not isinstance(solvers_raw, dict):
        raise ValueError("solvers must be a mapping in schema v2")
    unknown = set(solvers_raw) - set(CANONICAL_SOLVER_IDS)
    if unknown:
        for key in unknown:
            if key in OBSOLETE_PUBLIC_NAMES:
                _canonical_solver(key, f"solvers.{key}")
        raise ValueError(f"Unsupported solver config sections: {sorted(unknown)}")

    tt_raw = solvers_raw.get("two_term", {}) or {}
    mt_raw = solvers_raw.get("multi_term", {}) or {}
    mc_raw = solvers_raw.get("monte_carlo", {}) or {}
    _reject_unknown_fields(tt_raw, TWO_TERM_SOLVER_FIELDS, "solvers.two_term")
    _reject_unknown_fields(mt_raw, MULTI_TERM_SOLVER_FIELDS, "solvers.multi_term")
    _reject_unknown_fields(mc_raw, MONTE_CARLO_SOLVER_FIELDS, "solvers.monte_carlo")

    public_backend = cast(
        BoltzmannBackend,
        _validate_literal(
            str(tt_raw.get("backend", "native_sg")),
            {"auto", "native_sg", "native_bolsig", "internal", "bolos"},
            "solvers.two_term.backend",
        ),
    )
    internal_backend: InternalBoltzmannBackend = (
        "native_bolsig" if public_backend == "native_sg" else cast(InternalBoltzmannBackend, public_backend)
    )
    two_term = TwoTermProductConfig(
        backend=public_backend,
        nonconservative_model=cast(
            Literal["growth", "ignore"],
            _validate_literal(
                str(tt_raw.get("nonconservative_model", "growth")),
                {"growth", "ignore"},
                "solvers.two_term.nonconservative_model",
            ),
        ),
        min_momentum_cross_section_m2=float(
            tt_raw.get("min_momentum_cross_section_m2", 1.0e-24)
        ),
    )

    method = cast(
        Literal["pn_closure_surrogate", "pn_closure_direct", "pn_dcs"],
        _validate_literal(
            str(mt_raw.get("method", "pn_closure_surrogate")),
            {"pn_closure_surrogate", "pn_closure_direct", "pn_dcs"},
            "solvers.multi_term.method",
        ),
    )
    lmax = int(mt_raw.get("lmax", 4))
    multi_term = MultiTermProductConfig(
        formulation=cast(
            Literal["PN"],
            _validate_literal(
                str(mt_raw.get("formulation", "PN")),
                {"PN"},
                "solvers.multi_term.formulation",
            ),
        ),
        method=method,
        lmax=lmax,
        lmax_convergence_tolerance=float(
            mt_raw.get("lmax_convergence_tolerance", 0.03)
        ),
        field_coupling_scale=float(mt_raw.get("field_coupling_scale", 1.0)),
        dense_threshold=int(mt_raw.get("dense_threshold", 280)),
    )
    if multi_term.lmax < 1:
        raise ValueError("solvers.multi_term.lmax must be >= 1")

    monte_carlo = MonteCarloProductConfig(
        particles=(
            int(mc_raw["particles"]) if mc_raw.get("particles") is not None else None
        ),
        max_collisions=(
            int(mc_raw["max_collisions"])
            if mc_raw.get("max_collisions") is not None
            else None
        ),
        timeout_s=(
            float(mc_raw["timeout_s"]) if mc_raw.get("timeout_s") is not None else None
        ),
        command=mc_raw.get("command"),
        python_api=mc_raw.get("python_api"),
        working_directory=_as_path(mc_raw.get("working_directory"), base),
        environment={
            str(k): str(v) for k, v in (mc_raw.get("environment", {}) or {}).items()
        },
        output_summary_csv=_as_path(mc_raw.get("output_summary_csv"), base),
        output_eedf_csv=_as_path(mc_raw.get("output_eedf_csv"), base),
        passthrough=mc_raw.get("passthrough", {}) or {},
    )

    refine = EnergyGridRefinementConfig(
        enabled=physics.energy_grid_policy.threshold_refinement
    )
    adaptive = AdaptiveGridConfig(
        enabled=physics.energy_grid_policy.adaptive,
        max_max_eV=physics.energy_grid_policy.max_eV_limit,
        tail_probability=physics.energy_grid_policy.tail_probability_target,
    )
    boltzmann = BoltzmannTwoTermConfig(
        backend=internal_backend,
        energy_grid=EnergyGridConfig(refine=refine),
        adaptive_grid=adaptive,
        nonconservative_model=two_term.nonconservative_model,
        ionization_energy_sharing=physics.ionization.energy_sharing,
        secondary_electron_energy_eV=physics.ionization.secondary_electron_energy_eV,
        min_momentum_cross_section_m2=two_term.min_momentum_cross_section_m2,
    )
    mt_internal = MultiTermBoltzmannConfig(
        lmax=multi_term.lmax,
        method="moment_closure",
        hydrodynamic=False,
        dense_threshold=multi_term.dense_threshold,
        lmax_convergence_tolerance=multi_term.lmax_convergence_tolerance,
        field_coupling_scale=multi_term.field_coupling_scale,
        product_method=multi_term.method,
        formulation=multi_term.formulation,
        energy_grid=MultiTermEnergyGridConfig(refine=refine),
    )
    mc_internal = MonteCarloAdapterConfig(
        command=monte_carlo.command,
        python_api=monte_carlo.python_api,
        working_directory=monte_carlo.working_directory,
        timeout_s=monte_carlo.timeout_s,
        environment=monte_carlo.environment,
        output_summary_csv=monte_carlo.output_summary_csv,
        output_eedf_csv=monte_carlo.output_eedf_csv,
        passthrough=monte_carlo.passthrough,
        particles=monte_carlo.particles,
        max_collisions=monte_carlo.max_collisions,
    )
    return (
        SolversConfig(two_term, multi_term, monte_carlo),
        InternalSolverConfigs(
            two_term=boltzmann,
            multi_term=mt_internal,
            monte_carlo=mc_internal,
        ),
    )


def _parse_comparison(raw: dict[str, Any]) -> ComparisonConfig:
    cmp_raw = _as_mapping_section(raw, "comparison", "comparison")
    _reject_unknown_fields(
        cmp_raw,
        {"enabled", "reference_solver", "candidate_solvers", "compare_eedf", "required"},
        "comparison",
    )
    reference = cmp_raw.get("reference_solver")
    candidates = cmp_raw.get("candidate_solvers", []) or []
    return ComparisonConfig(
        enabled=_bool_field(cmp_raw, "enabled", False, "comparison.enabled"),
        reference_solver=(
            _canonical_solver(reference, "comparison.reference_solver")
            if reference is not None
            else None
        ),
        candidate_solvers=[
            _canonical_solver(item, "comparison.candidate_solvers")
            for item in candidates
        ],
        compare_eedf=_bool_field(cmp_raw, "compare_eedf", True, "comparison.compare_eedf"),
        required=_bool_field(cmp_raw, "required", False, "comparison.required"),
    )


def _parse_feature_policy(raw: dict[str, Any]) -> FeaturePolicyConfig:
    policy_raw = raw.get("feature_policy", {}) or {}
    _reject_unknown_fields(
        policy_raw,
        {"unsupported", "degraded", "allow_unsupported_fallback"},
        "feature_policy",
    )
    return FeaturePolicyConfig(
        unsupported=cast(
            UnsupportedPolicy,
            _validate_literal(
                str(policy_raw.get("unsupported", "fail")),
                {"fail", "skip_solver", "approximate"},
                "feature_policy.unsupported",
            ),
        ),
        degraded=cast(
            DegradedPolicy,
            _validate_literal(
                str(policy_raw.get("degraded", "warn")),
                {"fail", "warn", "record_only"},
                "feature_policy.degraded",
            ),
        ),
        allow_unsupported_fallback=_bool_field(
            policy_raw,
            "allow_unsupported_fallback",
            False,
            "feature_policy.allow_unsupported_fallback",
        ),
    )


def _parse_output(raw: dict[str, Any], base: Path) -> OutputConfig:
    out_raw = raw.get("output", {}) or {}
    return OutputConfig(
        directory=_as_path(out_raw.get("directory", "outputs"), base)
        or (base / "outputs"),
        base_name=str(out_raw.get("base_name", "swarm")),
        float_format=str(out_raw.get("float_format", "%.10e")),
    )


def load_config(path: str | Path) -> SwarmConfig:
    """Load and validate a schema v2 product YAML configuration file."""

    cfg_path = Path(path).resolve()
    base = cfg_path.parent
    raw = _read_mapping(cfg_path)
    _validate_schema(raw)
    run = _parse_run(raw)
    conditions = _parse_conditions(raw)
    cross_sections = _parse_cross_sections(raw, base)
    physics = _parse_physics(raw, base)
    solvers, internal = _parse_solvers(raw, base, physics)
    comparison = _parse_comparison(raw)
    feature_policy = _parse_feature_policy(raw)
    output = _parse_output(raw, base)
    return SwarmConfig(
        schema_version=2,
        run=run,
        conditions=conditions,
        cross_sections=cross_sections,
        physics=physics,
        solvers=solvers,
        comparison=comparison,
        feature_policy=feature_policy,
        output=output,
        source_path=cfg_path,
        internal=internal,
    )
