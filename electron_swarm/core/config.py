"""Product YAML configuration model and validation.

Product mode uses schema v2 only. Legacy v1 runner-mode YAML is rejected at
load time with an explicit migration message.
"""

from __future__ import annotations

from dataclasses import dataclass, field as dc_field
from math import isfinite as np_isfinite
from pathlib import Path
from typing import Any, Literal, cast

import yaml

from electron_swarm.core.legacy_schema import (
    MIGRATION_ERROR,
    OBSOLETE_PUBLIC_NAMES,
    reject_removed_schema,
)
from electron_swarm.core.solver_registry import CANONICAL_SOLVER_IDS

SolverId = Literal["two_term", "multi_term", "monte_carlo"]
UnsupportedPolicy = Literal["fail", "skip_solver"]
DegradedPolicy = Literal["fail", "record"]
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
    "population_model",
    "seed",
    "particles",
    "warmup_collisions",
    "max_collisions",
}
REMOVED_MONTE_CARLO_SOLVER_FIELDS = {
    "backend",
    "angular_scattering",
    "timeout_s",
    "command",
    "python_api",
    "working_directory",
    "environment",
    "output_summary_csv",
    "output_eedf_csv",
    "passthrough",
}
TOP_LEVEL_FIELDS = {
    "schema_version",
    "run",
    "conditions",
    "cross_sections",
    "physics",
    "solvers",
    "comparison",
    "feature_policy",
    "output",
}

BoltzmannBackend = Literal["auto", "native_sg", "native_bolsig", "internal", "bolos"]


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
class TwoTermProductConfig:
    backend: BoltzmannBackend = "native_sg"
    nonconservative_model: Literal["growth", "ignore"] = "growth"
    min_momentum_cross_section_m2: float = 1.0e-24


@dataclass(slots=True)
class MultiTermProductConfig:
    formulation: Literal["PN"] = "PN"
    method: Literal["pn_closure_direct", "pn_dcs"] = "pn_closure_direct"
    lmax: int = 4
    lmax_convergence_tolerance: float = 0.03
    field_coupling_scale: float = 1.0
    dense_threshold: int = 280


@dataclass(slots=True)
class MonteCarloProductConfig:
    population_model: Literal[
        "fixed_particle_single_daughter", "weighted_branching"
    ] = "fixed_particle_single_daughter"
    seed: int | None = None
    particles: int | None = None
    warmup_collisions: int | None = None
    max_collisions: int | None = None


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
class MomentTableConfig:
    path: Path
    format: Literal["normalized_legendre_moments"] = "normalized_legendre_moments"
    provenance: Literal["dcs_derived", "model_derived", "unknown"] = "unknown"
    extrapolation: Literal["error"] = "error"


@dataclass(slots=True)
class AngularScatteringConfig:
    model: Literal["isotropic", "momentum_power", "maxent_p1", "moment_table"] = "isotropic"
    higher_moment_closure: Literal["zero", "power", "maxent", "table"] = "zero"
    moment_table: MomentTableConfig | None = None


@dataclass(slots=True)
class ElectronElectronConfig:
    enabled: bool = False
    model: Literal["none", "relaxation_postprocess", "fp_energy"] = "none"
    strength_model: Literal["simple_relaxation", "density_based"] = "simple_relaxation"
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
    degraded: DegradedPolicy = "record"


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


def _float_field(raw: dict[str, Any], key: str, default: float, field_name: str) -> float:
    try:
        return float(raw.get(key, default))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must be a finite number") from exc


def _parse_finite_k(raw: dict[str, Any]) -> FiniteKConfig:
    enabled = _bool_field(raw, "enabled", False, "physics.finite_k.enabled")
    if raw.get("k_m_inv") is None:
        if enabled:
            raise ValueError("physics.finite_k.k_m_inv is required when finite_k is enabled")
        return FiniteKConfig(enabled=enabled, k_m_inv=None)
    try:
        k_m_inv = float(raw["k_m_inv"])
    except (TypeError, ValueError) as exc:
        raise ValueError("physics.finite_k.k_m_inv must be finite and > 0") from exc
    if not np_isfinite(k_m_inv) or k_m_inv <= 0.0:
        raise ValueError("physics.finite_k.k_m_inv must be finite and > 0")
    return FiniteKConfig(enabled=enabled, k_m_inv=k_m_inv)


def _validate_schema(raw: dict[str, Any]) -> None:
    reject_removed_schema(raw, output_fields=OUTPUT_FIELDS)
    _reject_unknown_fields(raw, TOP_LEVEL_FIELDS, "top-level")
    out_raw = raw.get("output", {}) or {}
    if not isinstance(out_raw, dict):
        raise ValueError("output must be a mapping")
    _reject_unknown_fields(out_raw, OUTPUT_FIELDS, "output")


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


def _parse_run(raw: dict[str, Any]) -> RunConfig:
    run_raw = raw.get("run", {}) or {}
    if not isinstance(run_raw, dict):
        raise ValueError("run must be a mapping")
    _reject_unknown_fields(run_raw, {"solvers", "e_over_n_Td", "case_prefix"}, "run")
    solvers_raw = run_raw.get("solvers")
    if not isinstance(solvers_raw, list) or not solvers_raw:
        raise ValueError("schema v2 requires non-empty run.solvers")
    solvers: list[RequestedSolverConfig] = []
    for index, item in enumerate(solvers_raw):
        if not isinstance(item, dict) or "id" not in item:
            raise ValueError("run.solvers entries must be mappings with id")
        _reject_unknown_fields(item, {"id", "enabled"}, f"run.solvers[{index}]")
        solvers.append(
            RequestedSolverConfig(
                id=_canonical_solver(item["id"], f"run.solvers[{index}].id"),
                enabled=_bool_field(item, "enabled", True, f"run.solvers[{index}].enabled"),
            )
        )
    run = RunConfig(
        solvers=solvers,
        e_over_n_Td=_list_float(run_raw.get("e_over_n_Td", [100.0])),
        case_prefix=str(run_raw.get("case_prefix", "case")),
    )
    if any(v <= 0 for v in run.e_over_n_Td):
        raise ValueError("All run.e_over_n_Td values must be positive")
    return run


def _parse_conditions(raw: dict[str, Any]) -> ConditionsConfig:
    cond_raw = raw.get("conditions", {}) or {}
    if not isinstance(cond_raw, dict):
        raise ValueError("conditions must be a mapping")
    _reject_unknown_fields(
        cond_raw,
        {
            "gas_temperature_K",
            "pressure_Pa",
            "gas_number_density_m3",
            "length_scale_m",
            "gas_mixture",
            "species",
            "mass_amu",
        },
        "conditions",
    )
    gas_mixture_raw = cond_raw.get("gas_mixture", []) or []
    if not isinstance(gas_mixture_raw, list):
        raise ValueError("conditions.gas_mixture must be a list")
    for index, item in enumerate(gas_mixture_raw):
        if not isinstance(item, dict):
            raise ValueError("conditions.gas_mixture entries must be mappings")
        _reject_unknown_fields(
            item,
            {"species", "fraction", "mass_amu"},
            f"conditions.gas_mixture[{index}]",
        )
        if "species" not in item:
            raise ValueError(f"conditions.gas_mixture[{index}].species is required")
    gas_mixture = [
        GasComponent(
            species=str(item["species"]),
            fraction=float(item.get("fraction", 1.0)),
            mass_amu=float(item.get("mass_amu", 39.948)),
        )
        for item in gas_mixture_raw
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
    if not isinstance(xs_raw, dict):
        raise ValueError("cross_sections must be a mapping")
    _reject_unknown_fields(
        xs_raw,
        {"format", "files", "high_energy_extrapolation"},
        "cross_sections",
    )
    xs_format = str(xs_raw.get("format", "csv"))
    xs_files: list[CrossSectionFileConfig] = []
    files_raw = xs_raw.get("files", [])
    if not isinstance(files_raw, list):
        raise ValueError("cross_sections.files must be a list")
    for index, item in enumerate(files_raw):
        if not isinstance(item, dict):
            raise ValueError("cross_sections.files entries must be mappings")
        _reject_unknown_fields(
            item,
            {"path", "species", "format"},
            f"cross_sections.files[{index}]",
        )
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
        "strength_model",
        "relaxation_fraction",
        "conserve_mean_energy",
        "fallback_temperature_eV",
    }
    if ee_unknown:
        raise ValueError(
            "Unsupported physics.electron_electron fields: "
            f"{sorted(ee_unknown)}"
        )
    angular_name = str(angular_raw.get("model", "isotropic")).lower()
    if angular_name == "dcs_table" or "dcs_table" in angular_raw:
        raise ValueError(
            "physics.angular_scattering.dcs_table input is not implemented; "
            "use model=moment_table"
        )
    angular_unknown = set(angular_raw) - {
        "model",
        "higher_moment_closure",
        "moment_table",
    }
    if angular_unknown:
        raise ValueError(
            "Unsupported physics.angular_scattering fields: "
            f"{sorted(angular_unknown)}"
        )
    moment_table_raw = _as_mapping_section(
        angular_raw, "moment_table", "physics.angular_scattering.moment_table"
    )
    angular_closure_defaults = {
        "isotropic": "zero",
        "momentum_power": "power",
        "maxent_p1": "maxent",
        "moment_table": "table",
    }
    if angular_name not in angular_closure_defaults:
        raise ValueError(
            "physics.angular_scattering.model must be isotropic, momentum_power, maxent_p1, or moment_table"
        )
    if angular_name != "moment_table" and moment_table_raw:
        raise ValueError(
            "physics.angular_scattering.moment_table is only valid with model=moment_table"
        )
    closure_raw = angular_raw.get(
        "higher_moment_closure", angular_closure_defaults[angular_name]
    )
    higher_moment_closure = cast(
        Literal["zero", "power", "maxent", "table"],
        _validate_literal(
            str(closure_raw).lower(),
            {"zero", "power", "maxent", "table"},
            "physics.angular_scattering.higher_moment_closure",
        ),
    )
    if higher_moment_closure != angular_closure_defaults[angular_name]:
        raise ValueError(
            "physics.angular_scattering model and higher_moment_closure mismatch"
        )
    moment_table: MomentTableConfig | None = None
    if angular_name == "moment_table":
        _reject_unknown_fields(
            moment_table_raw,
            {"path", "format", "provenance", "extrapolation"},
            "physics.angular_scattering.moment_table",
        )
        if "path" not in moment_table_raw:
            raise ValueError("physics.angular_scattering.moment_table.path is required")
        moment_table = MomentTableConfig(
            path=_as_path(moment_table_raw["path"], base)
            or Path(moment_table_raw["path"]),
            format=cast(
                Literal["normalized_legendre_moments"],
                _validate_literal(
                    str(
                        moment_table_raw.get(
                            "format", "normalized_legendre_moments"
                        )
                    ),
                    {"normalized_legendre_moments"},
                    "physics.angular_scattering.moment_table.format",
                ),
            ),
            provenance=cast(
                Literal["dcs_derived", "model_derived", "unknown"],
                _validate_literal(
                    str(moment_table_raw.get("provenance", "unknown")),
                    {"dcs_derived", "model_derived", "unknown"},
                    "physics.angular_scattering.moment_table.provenance",
                ),
            ),
            extrapolation=cast(
                Literal["error"],
                _validate_literal(
                    str(moment_table_raw.get("extrapolation", "error")),
                    {"error"},
                    "physics.angular_scattering.moment_table.extrapolation",
                ),
            ),
        )
    ee_enabled = _bool_field(ee_raw, "enabled", False, "physics.electron_electron.enabled")
    ee_model = str(ee_raw.get("model", "none")).lower()
    if ee_enabled and ee_model not in {"relaxation_postprocess", "fp_energy"}:
        raise ValueError(
            "physics.electron_electron.model must be relaxation_postprocess or fp_energy when enabled"
        )
    if not ee_enabled and ee_model != "none":
        raise ValueError(
            "physics.electron_electron.model must be none when electron_electron is disabled"
        )
    ee_strength_model = cast(
        Literal["simple_relaxation", "density_based"],
        _validate_literal(
            str(ee_raw.get("strength_model", "simple_relaxation")).lower(),
            {"simple_relaxation", "density_based"},
            "physics.electron_electron.strength_model",
        ),
    )
    ee_relaxation_fraction = float(ee_raw.get("relaxation_fraction", 0.05))
    if not (0.0 <= ee_relaxation_fraction <= 1.0):
        raise ValueError("physics.electron_electron.relaxation_fraction must be in [0, 1]")
    ee_fallback_temperature = float(ee_raw.get("fallback_temperature_eV", 2.0))
    if ee_fallback_temperature <= 0.0:
        raise ValueError(
            "physics.electron_electron.fallback_temperature_eV must be positive"
        )
    ionization_energy_sharing = cast(
        Literal["equal", "primary_secondary", "loss_only"],
        _validate_literal(
            str(ion_raw.get("energy_sharing", "equal")),
            {"equal", "primary_secondary", "loss_only"},
            "physics.ionization.energy_sharing",
        ),
    )
    secondary_electron_energy_eV = float(
        ion_raw.get("secondary_electron_energy_eV", 0.0)
    )
    if secondary_electron_energy_eV < 0.0:
        raise ValueError("physics.ionization.secondary_electron_energy_eV must be >= 0")
    if ionization_energy_sharing != "primary_secondary" and secondary_electron_energy_eV != 0.0:
        raise ValueError(
            "physics.ionization.secondary_electron_energy_eV is only used with primary_secondary"
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
    magnetic_B_T = _float_field(
        magnetic_raw, "B_T", 0.0, "physics.field.magnetic_field.B_T"
    )
    magnetic_angle = _float_field(
        magnetic_raw,
        "angle_EB_deg",
        0.0,
        "physics.field.magnetic_field.angle_EB_deg",
    )
    if not np_isfinite(magnetic_B_T) or magnetic_B_T < 0.0:
        raise ValueError("physics.field.magnetic_field.B_T must be finite and >= 0")
    if not np_isfinite(magnetic_angle):
        raise ValueError("physics.field.magnetic_field.angle_EB_deg must be finite")
    if magnetic_angle < 0.0 or magnetic_angle > 180.0:
        raise ValueError(
            "physics.field.magnetic_field.angle_EB_deg must be in [0, 180]"
        )
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
                B_T=magnetic_B_T,
                angle_EB_deg=magnetic_angle,
            ),
        ),
        angular_scattering=AngularScatteringConfig(
            model=cast(
                Literal["isotropic", "momentum_power", "maxent_p1", "moment_table"],
                angular_name,
            ),
            higher_moment_closure=higher_moment_closure,
            moment_table=moment_table,
        ),
        electron_electron=ElectronElectronConfig(
            enabled=ee_enabled,
            model=cast(Literal["none", "relaxation_postprocess", "fp_energy"], ee_model),
            strength_model=ee_strength_model,
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
            energy_sharing=ionization_energy_sharing,
            secondary_electron_energy_eV=secondary_electron_energy_eV,
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
        finite_k=_parse_finite_k(finite_k_raw),
    )


def _parse_solvers(raw: dict[str, Any], base: Path) -> SolversConfig:
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
    removed_mc = sorted(set(mc_raw) & REMOVED_MONTE_CARLO_SOLVER_FIELDS)
    if removed_mc:
        raise ValueError(
            f"{MIGRATION_ERROR}; solvers.monte_carlo is the internal product "
            f"backend only, remove fields {removed_mc}"
        )
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
        Literal["pn_closure_direct", "pn_dcs"],
        _validate_literal(
            str(mt_raw.get("method", "pn_closure_direct")),
            {"pn_closure_direct", "pn_dcs"},
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

    mc_population_model = cast(
        Literal["fixed_particle_single_daughter", "weighted_branching"],
        _validate_literal(
            str(mc_raw.get("population_model", "fixed_particle_single_daughter")),
            {"fixed_particle_single_daughter", "weighted_branching"},
            "solvers.monte_carlo.population_model",
        ),
    )
    mc_particles = (
        int(mc_raw["particles"]) if mc_raw.get("particles") is not None else None
    )
    mc_max_collisions = (
        int(mc_raw["max_collisions"])
        if mc_raw.get("max_collisions") is not None
        else None
    )
    mc_warmup_collisions = (
        int(mc_raw["warmup_collisions"])
        if mc_raw.get("warmup_collisions") is not None
        else None
    )
    if mc_particles is not None and mc_particles <= 0:
        raise ValueError("solvers.monte_carlo.particles must be positive")
    if mc_max_collisions is not None and mc_max_collisions <= 0:
        raise ValueError("solvers.monte_carlo.max_collisions must be positive")
    if mc_warmup_collisions is not None and mc_warmup_collisions < 0:
        raise ValueError("solvers.monte_carlo.warmup_collisions must be nonnegative")
    mc_seed = int(mc_raw["seed"]) if mc_raw.get("seed") is not None else None

    monte_carlo = MonteCarloProductConfig(
        population_model=mc_population_model,
        seed=mc_seed,
        particles=mc_particles,
        warmup_collisions=mc_warmup_collisions,
        max_collisions=mc_max_collisions,
    )
    return SolversConfig(two_term, multi_term, monte_carlo)


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
    if "allow_unsupported_fallback" in policy_raw:
        raise ValueError(
            f"{MIGRATION_ERROR}; remove feature_policy.allow_unsupported_fallback"
        )
    _reject_unknown_fields(
        policy_raw,
        {"unsupported", "degraded"},
        "feature_policy",
    )
    return FeaturePolicyConfig(
        unsupported=cast(
            UnsupportedPolicy,
            _validate_literal(
                str(policy_raw.get("unsupported", "fail")),
                {"fail", "skip_solver"},
                "feature_policy.unsupported",
            ),
        ),
        degraded=cast(
            DegradedPolicy,
            _validate_literal(
                str(policy_raw.get("degraded", "record")),
                {"fail", "record"},
                "feature_policy.degraded",
            ),
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


def _load_config_from_raw(raw: dict[str, Any], cfg_path: str | Path) -> SwarmConfig:
    cfg_path = Path(cfg_path).resolve()
    base = cfg_path.parent
    _validate_schema(raw)
    run = _parse_run(raw)
    conditions = _parse_conditions(raw)
    cross_sections = _parse_cross_sections(raw, base)
    physics = _parse_physics(raw, base)
    solvers = _parse_solvers(raw, base)
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
    )


def load_config(path: str | Path) -> SwarmConfig:
    """Load and validate a schema v2 product YAML configuration file."""

    cfg_path = Path(path).resolve()
    return _load_config_from_raw(_read_mapping(cfg_path), cfg_path)
