"""YAML configuration model and validation.

The schema separates solver-neutral physics/input from solver-specific controls.
Monte Carlo and Boltzmann two-term controls are intentionally namespaced so a
single YAML file can run either solver or both without ambiguous variables.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal, cast

import yaml

from electron_swarm.core.solver_registry import PRIMARY_SOLVERS, RUN_MODES

RunMode = Literal[
    "monte_carlo",
    "boltzmann_two_term",
    "multiterm_boltzmann",
    "both",
    "all",
]
BoltzmannBackend = Literal["auto", "native_bolsig", "internal", "bolos"]
PrimarySolver = Literal["monte_carlo", "boltzmann_two_term", "multiterm_boltzmann"]


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
    gas_mixture: list[GasComponent] = field(default_factory=list)


@dataclass(slots=True)
class CrossSectionFileConfig:
    path: Path
    species: str | None = None
    format: str | None = None


@dataclass(slots=True)
class CrossSectionsConfig:
    format: str = "csv"
    files: list[CrossSectionFileConfig] = field(default_factory=list)
    high_energy_extrapolation: Literal["zero", "hold", "error"] = "zero"


@dataclass(slots=True)
class EnergyGridConfig:
    min_eV: float = 1.0e-4
    max_eV: float = 100.0
    n: int = 600
    spacing: Literal["linear", "quadratic", "log"] = "quadratic"


@dataclass(slots=True)
class MultiTermEnergyGridConfig:
    min_eV: float = 1.0e-3
    max_eV: float = 200.0
    n: int = 220
    spacing: Literal["linear", "log", "log_linear"] = "log_linear"
    linear_until_eV: float = 2.0


@dataclass(slots=True)
class AdaptiveGridConfig:
    """Automatic upper-energy-domain control for Boltzmann solves."""

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
    backend: BoltzmannBackend = "native_bolsig"
    energy_grid: EnergyGridConfig = field(default_factory=EnergyGridConfig)
    adaptive_grid: AdaptiveGridConfig = field(default_factory=AdaptiveGridConfig)
    convergence: ConvergenceConfig = field(default_factory=ConvergenceConfig)
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
    lmax: int = 3
    method: Literal["hybrid", "operator", "moment_closure"] = "moment_closure"
    hydrodynamic: bool = False
    dense_threshold: int = 280
    eedf_shape: Literal["maxwellian", "druyvesteyn"] = "maxwellian"
    lmax_convergence_tolerance: float = 0.03
    field_coupling_scale: float = 1.0
    energy_grid: MultiTermEnergyGridConfig = field(
        default_factory=MultiTermEnergyGridConfig
    )


@dataclass(slots=True)
class MonteCarloAdapterConfig:
    """Controls for delegating to the existing particle Monte Carlo code."""

    enabled: bool = True
    command: str | None = None
    python_api: str | None = None  # module:function
    working_directory: Path | None = None
    timeout_s: float | None = None
    environment: dict[str, str] = field(default_factory=dict)
    output_summary_csv: Path | None = None
    output_eedf_csv: Path | None = None
    passthrough: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class RunConfig:
    mode: RunMode = "boltzmann_two_term"
    e_over_n_Td: list[float] = field(default_factory=lambda: [100.0])
    case_prefix: str = "case"


@dataclass(slots=True)
class CompatibilityOutputConfig:
    write_legacy_tables: bool = True
    primary_solver: PrimarySolver = "monte_carlo"


@dataclass(slots=True)
class OutputConfig:
    directory: Path = Path("outputs")
    base_name: str = "swarm"
    write_plots: bool = True
    write_eedf: bool = True
    write_rates: bool = True
    float_format: str = "%.10e"
    compatibility: CompatibilityOutputConfig = field(
        default_factory=CompatibilityOutputConfig
    )


@dataclass(slots=True)
class SwarmConfig:
    run: RunConfig = field(default_factory=RunConfig)
    conditions: ConditionsConfig = field(default_factory=ConditionsConfig)
    cross_sections: CrossSectionsConfig = field(default_factory=CrossSectionsConfig)
    monte_carlo: MonteCarloAdapterConfig = field(
        default_factory=MonteCarloAdapterConfig
    )
    boltzmann_two_term: BoltzmannTwoTermConfig = field(
        default_factory=BoltzmannTwoTermConfig
    )
    multiterm_boltzmann: MultiTermBoltzmannConfig = field(
        default_factory=MultiTermBoltzmannConfig
    )
    output: OutputConfig = field(default_factory=OutputConfig)
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


def load_config(path: str | Path) -> SwarmConfig:
    """Load and validate a YAML configuration file."""

    cfg_path = Path(path).resolve()
    base = cfg_path.parent
    raw = _read_mapping(cfg_path)

    run_raw = raw.get("run", {}) or {}
    run = RunConfig(
        mode=_validate_literal(
            str(run_raw.get("mode", "boltzmann_two_term")),
            set(RUN_MODES),
            "run.mode",
        ),
        e_over_n_Td=_list_float(
            run_raw.get("e_over_n_Td", run_raw.get("E_over_N_Td", [100.0]))
        ),
        case_prefix=str(run_raw.get("case_prefix", "case")),
    )
    if any(v <= 0 for v in run.e_over_n_Td):
        raise ValueError("All run.e_over_n_Td values must be positive")

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
    cross_sections = CrossSectionsConfig(
        format=xs_format,
        files=xs_files,
        high_energy_extrapolation=high_energy_extrapolation,
    )

    mc_raw = raw.get("monte_carlo", {}) or {}
    monte_carlo = MonteCarloAdapterConfig(
        enabled=bool(mc_raw.get("enabled", True)),
        command=mc_raw.get("command"),
        python_api=mc_raw.get("python_api"),
        working_directory=_as_path(mc_raw.get("working_directory"), base),
        timeout_s=(
            float(mc_raw["timeout_s"]) if mc_raw.get("timeout_s") is not None else None
        ),
        environment={
            str(k): str(v) for k, v in (mc_raw.get("environment", {}) or {}).items()
        },
        output_summary_csv=_as_path(mc_raw.get("output_summary_csv"), base),
        output_eedf_csv=_as_path(mc_raw.get("output_eedf_csv"), base),
        passthrough=mc_raw.get("passthrough", {}) or {},
    )

    b_raw = raw.get("boltzmann_two_term", {}) or {}
    g_raw = b_raw.get("energy_grid", {}) or {}
    a_raw = b_raw.get("adaptive_grid", {}) or {}
    c_raw = b_raw.get("convergence", {}) or {}
    boltzmann = BoltzmannTwoTermConfig(
        enabled=bool(b_raw.get("enabled", True)),
        backend=_validate_literal(
            str(b_raw.get("backend", "native_bolsig")),
            {"auto", "native_bolsig", "internal", "bolos"},
            "boltzmann_two_term.backend",
        ),
        energy_grid=EnergyGridConfig(
            min_eV=float(g_raw.get("min_eV", 1.0e-4)),
            max_eV=float(g_raw.get("max_eV", 100.0)),
            n=int(g_raw.get("n", 600)),
            spacing=_validate_literal(
                str(g_raw.get("spacing", "quadratic")),
                {"linear", "quadratic", "log"},
                "boltzmann_two_term.energy_grid.spacing",
            ),
        ),
        adaptive_grid=AdaptiveGridConfig(
            enabled=bool(a_raw.get("enabled", True)),
            max_cycles=int(a_raw.get("max_cycles", 4)),
            mean_energy_multiplier=float(a_raw.get("mean_energy_multiplier", 15.0)),
            min_max_eV=float(a_raw.get("min_max_eV", 20.0)),
            max_max_eV=float(a_raw.get("max_max_eV", 2000.0)),
            tail_probability=float(a_raw.get("tail_probability", 1.0e-8)),
            tail_cells_fraction=float(a_raw.get("tail_cells_fraction", 0.05)),
            edge_to_peak=float(a_raw.get("edge_to_peak", 1.0e-10)),
        ),
        convergence=ConvergenceConfig(
            max_iterations=int(c_raw.get("max_iterations", 120)),
            tolerance=float(c_raw.get("tolerance", 1.0e-8)),
            eigenvalue_tolerance=float(c_raw.get("eigenvalue_tolerance", 1.0e-8)),
            residual_tolerance=float(c_raw.get("residual_tolerance", 1.0e-7)),
            relaxation=float(c_raw.get("relaxation", 0.7)),
            clip_negative=bool(c_raw.get("clip_negative", True)),
        ),
        initial_electron_temperature_eV=float(
            b_raw.get("initial_electron_temperature_eV", 2.0)
        ),
        nonconservative_model=_validate_literal(
            str(b_raw.get("nonconservative_model", "growth")),
            {"growth", "ignore"},
            "boltzmann_two_term.nonconservative_model",
        ),
        ionization_energy_sharing=_validate_literal(
            str(b_raw.get("ionization_energy_sharing", "equal")),
            {"equal", "primary_secondary", "loss_only"},
            "boltzmann_two_term.ionization_energy_sharing",
        ),
        secondary_electron_energy_eV=float(
            b_raw.get("secondary_electron_energy_eV", 0.0)
        ),
        min_momentum_cross_section_m2=float(
            b_raw.get("min_momentum_cross_section_m2", 1.0e-24)
        ),
    )
    if boltzmann.energy_grid.n < 8:
        raise ValueError("boltzmann_two_term.energy_grid.n must be >= 8")
    if (
        boltzmann.energy_grid.min_eV < 0
        or boltzmann.energy_grid.max_eV <= boltzmann.energy_grid.min_eV
    ):
        raise ValueError("Invalid Boltzmann energy grid bounds")
    if not (0.0 < boltzmann.adaptive_grid.tail_cells_fraction < 1.0):
        raise ValueError(
            "boltzmann_two_term.adaptive_grid.tail_cells_fraction must be between 0 and 1"
        )
    if not (0.0 < boltzmann.convergence.relaxation <= 1.0):
        raise ValueError(
            "boltzmann_two_term.convergence.relaxation must be in (0, 1]"
        )

    mt_raw = raw.get("multiterm_boltzmann", {}) or {}
    mt_g_raw = mt_raw.get("energy_grid", {}) or {}
    multiterm = MultiTermBoltzmannConfig(
        enabled=bool(mt_raw.get("enabled", True)),
        lmax=int(mt_raw.get("lmax", 3)),
        method=_validate_literal(
            str(mt_raw.get("method", "moment_closure")),
            {"hybrid", "operator", "moment_closure"},
            "multiterm_boltzmann.method",
        ),
        hydrodynamic=bool(mt_raw.get("hydrodynamic", False)),
        dense_threshold=int(mt_raw.get("dense_threshold", 280)),
        eedf_shape=_validate_literal(
            str(mt_raw.get("eedf_shape", "maxwellian")),
            {"maxwellian", "druyvesteyn"},
            "multiterm_boltzmann.eedf_shape",
        ),
        lmax_convergence_tolerance=float(
            mt_raw.get("lmax_convergence_tolerance", 0.03)
        ),
        field_coupling_scale=float(mt_raw.get("field_coupling_scale", 1.0)),
        energy_grid=MultiTermEnergyGridConfig(
            min_eV=float(mt_g_raw.get("min_eV", 1.0e-3)),
            max_eV=float(mt_g_raw.get("max_eV", 200.0)),
            n=int(mt_g_raw.get("n", 220)),
            spacing=_validate_literal(
                str(mt_g_raw.get("spacing", "log_linear")),
                {"linear", "log", "log_linear"},
                "multiterm_boltzmann.energy_grid.spacing",
            ),
            linear_until_eV=float(mt_g_raw.get("linear_until_eV", 2.0)),
        ),
    )
    if multiterm.lmax < 1:
        raise ValueError("multiterm_boltzmann.lmax must be >= 1")
    if multiterm.energy_grid.n < 8:
        raise ValueError("multiterm_boltzmann.energy_grid.n must be >= 8")
    if (
        multiterm.energy_grid.min_eV < 0
        or multiterm.energy_grid.max_eV <= multiterm.energy_grid.min_eV
    ):
        raise ValueError("Invalid multiterm_boltzmann energy grid bounds")
    if multiterm.energy_grid.spacing in {"log", "log_linear"} and multiterm.energy_grid.min_eV <= 0:
        raise ValueError(
            "multiterm_boltzmann log grids require energy_grid.min_eV > 0"
        )

    out_raw = raw.get("output", {}) or {}
    compat_raw = out_raw.get("compatibility", {}) or {}
    output = OutputConfig(
        directory=_as_path(out_raw.get("directory", "outputs"), base)
        or (base / "outputs"),
        base_name=str(out_raw.get("base_name", "swarm")),
        write_plots=bool(out_raw.get("write_plots", True)),
        write_eedf=bool(out_raw.get("write_eedf", True)),
        write_rates=bool(out_raw.get("write_rates", True)),
        float_format=str(out_raw.get("float_format", "%.10e")),
        compatibility=CompatibilityOutputConfig(
            write_legacy_tables=bool(
                compat_raw.get("write_legacy_tables", True)
            ),
            primary_solver=_validate_literal(
                str(compat_raw.get("primary_solver", "monte_carlo")),
                set(PRIMARY_SOLVERS),
                "output.compatibility.primary_solver",
            ),
        ),
    )

    return SwarmConfig(
        run=run,
        conditions=conditions,
        cross_sections=cross_sections,
        monte_carlo=monte_carlo,
        boltzmann_two_term=boltzmann,
        multiterm_boltzmann=multiterm,
        output=output,
        source_path=cfg_path,
    )
