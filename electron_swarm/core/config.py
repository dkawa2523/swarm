"""Product YAML configuration model and validation.

Product mode uses schema v2 only. Legacy v1 runner-mode YAML is rejected at
load time with an explicit migration message.
"""

from __future__ import annotations

from dataclasses import dataclass, field as dc_field
from pathlib import Path
from typing import Any, Literal

from electron_swarm.core.legacy_schema import MIGRATION_ERROR
from electron_swarm.core.solver_registry import CANONICAL_SOLVER_IDS

SolverId = Literal["two_term", "multi_term", "monte_carlo"]
UnsupportedPolicy = Literal["fail", "skip_solver"]
DegradedPolicy = Literal["fail", "record"]
TwoTermBackend = Literal["native_sg"]


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
    backend: TwoTermBackend = "native_sg"
    nonconservative_model: Literal["growth", "ignore"] = "growth"
    min_momentum_cross_section_m2: float = 1.0e-24


@dataclass(slots=True)
class MultiTermProductConfig:
    method: Literal["pn_closure_direct", "pn_dcs"] = "pn_closure_direct"
    lmax: int = 4


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


def _load_config_from_raw(raw: dict[str, Any], cfg_path: str | Path) -> SwarmConfig:
    """Compatibility wrapper for development tooling that parses raw mappings."""

    from electron_swarm.core.config_parser import load_config_from_raw

    return load_config_from_raw(raw, cfg_path)


def load_config(path: str | Path) -> SwarmConfig:
    """Load and validate a schema v2 product YAML configuration file."""

    from electron_swarm.core.config_parser import read_mapping

    cfg_path = Path(path).resolve()
    return _load_config_from_raw(read_mapping(cfg_path), cfg_path)
