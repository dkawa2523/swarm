"""Validated configuration for external swarm workflow sweeps."""

from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
import math
from pathlib import Path
from typing import Any, Literal

import yaml

from electron_swarm import SwarmConfig, load_config

from ..selection import file_sha256, physical_context, read_selection
from ..quality.monte_carlo.policy import (
    MonteCarloConvergencePolicy,
    MonteCarloPolicyError,
    SamplingPlanEntry,
    parse_convergence_policy,
    validate_sampling_budget,
)
from ..quality.policy import QualityThresholds, parse_quality_thresholds


__all__ = [
    "DeterministicExecutionConfig",
    "MeanEnergySupportConfig",
    "MixtureSpec",
    "WorkflowConfig",
    "load_workflow",
]


@dataclass(frozen=True, slots=True)
class MixtureSpec:
    mixture_id: int
    fractions: dict[str, float]


@dataclass(frozen=True, slots=True)
class DeterministicExecutionConfig:
    workers: int = 1
    global_memory_budget_mb: int | None = None
    configured: bool = False


@dataclass(frozen=True, slots=True)
class MeanEnergySupportConfig:
    """Bounded deterministic continuation requested by a downstream model."""

    required_max_mean_energy_eV: float
    relative_guard: float = 0.1
    maximum_field_step_factor: float = 4.0
    maximum_steps: int = 2
    maximum_e_over_n_Td: float | None = None


@dataclass(frozen=True, slots=True)
class WorkflowConfig:
    path: Path
    base_config_path: Path
    database_path: Path
    e_over_n_Td: tuple[float, ...]
    mixtures: tuple[MixtureSpec, ...]
    mc_enabled: bool = False
    mc_e_over_n_Td: tuple[float, ...] | None = None
    mc_replicas: int = 1
    mc_base_seed: int | None = None
    mc_workers: int = 1
    mc_reuse_database: Path | None = None
    mc_sampling_plan: tuple[SamplingPlanEntry, ...] = ()
    mc_convergence: MonteCarloConvergencePolicy | None = None
    mc_previous_decision: Path | None = None
    quality: QualityThresholds = field(default_factory=QualityThresholds)
    deterministic_execution: DeterministicExecutionConfig = field(
        default_factory=DeterministicExecutionConfig
    )
    mean_energy_support: MeanEnergySupportConfig | None = None


def _read_mapping(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as fp:
        data = yaml.safe_load(fp) or {}
    if not isinstance(data, dict):
        raise ValueError("workflow YAML root must be a mapping")
    return data


def _resolve_path(value: object, base: Path, field_name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ValueError(f"{field_name} must be a path")
    path = Path(value)
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def _positive_float_sequence(value: object, field_name: str) -> tuple[float, ...]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"{field_name} must be a non-empty list")
    items = tuple(float(item) for item in value)
    if any(not math.isfinite(item) or item <= 0.0 for item in items):
        raise ValueError(f"{field_name} values must be finite and positive")
    return items


def _positive_int(value: object, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{field_name} must be a positive integer")
    return int(value)


def _nonnegative_int(value: object, field_name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{field_name} must be a nonnegative integer")
    return int(value)


def _parse_deterministic_execution(
    raw: object,
    *,
    deterministic_enabled: bool,
) -> DeterministicExecutionConfig:
    if raw is None:
        return DeterministicExecutionConfig()
    if not isinstance(raw, dict):
        raise ValueError("workflow execution must be a mapping")
    unknown = set(raw) - {"deterministic"}
    if unknown:
        raise ValueError(f"Unsupported workflow execution fields: {sorted(unknown)}")
    deterministic_raw = raw.get("deterministic", {}) or {}
    if not isinstance(deterministic_raw, dict):
        raise ValueError("workflow execution.deterministic must be a mapping")
    deterministic_unknown = set(deterministic_raw) - {
        "workers",
        "global_memory_budget_mb",
    }
    if deterministic_unknown:
        raise ValueError(
            "Unsupported workflow execution.deterministic fields: "
            f"{sorted(deterministic_unknown)}"
        )
    workers = _positive_int(
        deterministic_raw.get("workers", 1),
        "workflow execution.deterministic.workers",
    )
    raw_budget = deterministic_raw.get("global_memory_budget_mb")
    memory_budget = (
        _positive_int(
            raw_budget,
            "workflow execution.deterministic.global_memory_budget_mb",
        )
        if raw_budget is not None
        else None
    )
    if deterministic_enabled and workers > 1 and memory_budget is None:
        raise ValueError(
            "workflow execution.deterministic.global_memory_budget_mb is "
            "required when deterministic workers is greater than one"
        )
    return DeterministicExecutionConfig(
        workers=workers,
        global_memory_budget_mb=memory_budget,
        configured=bool(deterministic_raw),
    )


def _parse_mean_energy_support(
    raw: object,
    *,
    e_over_n_Td: tuple[float, ...],
    mixtures: tuple[MixtureSpec, ...],
    enabled_solvers: tuple[str, ...],
) -> MeanEnergySupportConfig | None:
    if raw is None:
        return None
    if not isinstance(raw, dict):
        raise ValueError("workflow mean_energy_support must be a mapping")
    allowed = {
        "required_max_mean_energy_eV",
        "relative_guard",
        "maximum_field_step_factor",
        "maximum_steps",
        "maximum_e_over_n_Td",
    }
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise ValueError(f"Unsupported workflow mean_energy_support fields: {unknown}")
    if len(mixtures) != 1 or len(enabled_solvers) != 1:
        raise ValueError(
            "workflow mean_energy_support requires one mixture and one enabled "
            "deterministic solver"
        )
    if enabled_solvers[0] == "monte_carlo":
        raise ValueError(
            "workflow mean_energy_support cannot extend monte_carlo; use its "
            "statistical campaign"
        )
    if any(right <= left for left, right in zip(e_over_n_Td, e_over_n_Td[1:])):
        raise ValueError(
            "workflow mean_energy_support requires strictly increasing E/N anchors"
        )

    def finite_number(name: str, default: object) -> float:
        value = raw.get(name, default)
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"workflow mean_energy_support.{name} must be finite")
        result = float(value)
        if not math.isfinite(result):
            raise ValueError(f"workflow mean_energy_support.{name} must be finite")
        return result

    required = finite_number("required_max_mean_energy_eV", None)
    guard = finite_number("relative_guard", 0.1)
    step_factor = finite_number("maximum_field_step_factor", 4.0)
    maximum_steps = raw.get("maximum_steps", 2)
    if required <= 0.0:
        raise ValueError(
            "workflow mean_energy_support.required_max_mean_energy_eV must be positive"
        )
    if guard < 0.0:
        raise ValueError(
            "workflow mean_energy_support.relative_guard must be nonnegative"
        )
    if step_factor <= 1.0:
        raise ValueError(
            "workflow mean_energy_support.maximum_field_step_factor must exceed one"
        )
    if (
        isinstance(maximum_steps, bool)
        or not isinstance(maximum_steps, int)
        or maximum_steps <= 0
    ):
        raise ValueError(
            "workflow mean_energy_support.maximum_steps must be a positive integer"
        )
    raw_field_limit = raw.get("maximum_e_over_n_Td")
    field_limit = (
        finite_number("maximum_e_over_n_Td", None)
        if raw_field_limit is not None
        else None
    )
    if field_limit is not None and field_limit <= e_over_n_Td[-1]:
        raise ValueError(
            "workflow mean_energy_support.maximum_e_over_n_Td must exceed the "
            "configured upper E/N anchor"
        )
    return MeanEnergySupportConfig(
        required_max_mean_energy_eV=required,
        relative_guard=guard,
        maximum_field_step_factor=step_factor,
        maximum_steps=maximum_steps,
        maximum_e_over_n_Td=field_limit,
    )


def _transport_correlation_lag(value: object, field_name: str) -> int:
    lag = _positive_int(value, field_name)
    if lag < 4 or lag & (lag - 1):
        raise ValueError(
            f"{field_name} must be a power of two greater than or equal to 4"
        )
    return lag


def _mc_transport_estimator(
    value: object,
    field_name: str,
) -> Literal["single_field", "paired_field_parity"]:
    if value not in {"single_field", "paired_field_parity"}:
        raise ValueError(f"{field_name} must be single_field or paired_field_parity")
    return value


def _parse_mc_sampling_plan(
    raw: object,
    *,
    anchors: tuple[float, ...],
    base_config: SwarmConfig,
    default_replicas: int,
) -> tuple[SamplingPlanEntry, ...]:
    """Resolve sparse workflow overrides into one canonical row per anchor."""

    if raw is None:
        rows: list[object] = []
    elif isinstance(raw, list):
        rows = raw
    else:
        raise ValueError("workflow mc.sampling_plan must be a list")

    allowed = {
        "e_over_n_Td",
        "particles",
        "warmup_collisions",
        "max_collisions",
        "tail_max_collisions",
        "replicas",
        "transport_correlation_lag_barriers",
        "transport_estimator",
    }
    required = {"e_over_n_Td"}
    anchor_set = set(anchors)
    overrides: dict[float, SamplingPlanEntry] = {}
    mc_defaults = base_config.solvers.monte_carlo
    default_particles = int(mc_defaults.particles or 256)
    default_warmup = int(mc_defaults.warmup_collisions or 0)
    default_production = int(mc_defaults.max_collisions or 100)
    default_tail = (
        int(mc_defaults.tail_max_collisions)
        if mc_defaults.tail_max_collisions is not None
        else None
    )
    default_lag = _transport_correlation_lag(
        getattr(mc_defaults, "transport_correlation_lag_barriers", 64),
        "solvers.monte_carlo.transport_correlation_lag_barriers",
    )
    default_transport_estimator = _mc_transport_estimator(
        getattr(mc_defaults, "transport_estimator", "single_field"),
        "solvers.monte_carlo.transport_estimator",
    )
    for index, item in enumerate(rows):
        field_name = f"workflow mc.sampling_plan[{index}]"
        if not isinstance(item, dict):
            raise ValueError(f"{field_name} must be a mapping")
        unknown = sorted(set(item) - allowed)
        if unknown:
            raise ValueError(f"{field_name} has unsupported fields: {unknown}")
        missing = sorted(required - set(item))
        if missing:
            raise ValueError(f"{field_name} is missing required fields: {missing}")

        e_over_n = float(item["e_over_n_Td"])
        if not math.isfinite(e_over_n) or e_over_n <= 0.0:
            raise ValueError(f"{field_name}.e_over_n_Td must be finite and positive")
        if e_over_n not in anchor_set:
            raise ValueError(
                f"{field_name}.e_over_n_Td={e_over_n:g} is not an MC anchor"
            )
        if e_over_n in overrides:
            raise ValueError(
                f"workflow mc.sampling_plan has duplicate E/N anchor {e_over_n:g} Td"
            )
        production_collisions = _positive_int(
            item.get("max_collisions", default_production),
            f"{field_name}.max_collisions",
        )
        raw_tail_collisions = item.get("tail_max_collisions", default_tail)
        tail_collisions = (
            0
            if raw_tail_collisions is None
            else _positive_int(
                raw_tail_collisions,
                f"{field_name}.tail_max_collisions",
            )
        )
        overrides[e_over_n] = SamplingPlanEntry(
            e_over_n_Td=e_over_n,
            particles=_positive_int(
                item.get("particles", default_particles),
                f"{field_name}.particles",
            ),
            warmup_collisions=_nonnegative_int(
                item.get("warmup_collisions", default_warmup),
                f"{field_name}.warmup_collisions",
            ),
            max_collisions=production_collisions,
            tail_max_collisions=tail_collisions,
            replicas=_positive_int(
                item.get("replicas", default_replicas),
                f"{field_name}.replicas",
            ),
            transport_correlation_lag_barriers=_transport_correlation_lag(
                item.get("transport_correlation_lag_barriers", default_lag),
                f"{field_name}.transport_correlation_lag_barriers",
            ),
            transport_estimator=_mc_transport_estimator(
                item.get("transport_estimator", default_transport_estimator),
                f"{field_name}.transport_estimator",
            ),
        )

    return tuple(
        overrides.get(
            e_over_n,
            SamplingPlanEntry(
                e_over_n_Td=e_over_n,
                particles=default_particles,
                warmup_collisions=default_warmup,
                max_collisions=default_production,
                tail_max_collisions=(default_tail if default_tail is not None else 0),
                replicas=default_replicas,
                transport_correlation_lag_barriers=default_lag,
                transport_estimator=default_transport_estimator,
            ),
        )
        for e_over_n in anchors
    )


def _normalize_mixture(
    raw: object,
    *,
    mixture_id: int,
    base_species: tuple[str, ...],
) -> MixtureSpec:
    if not isinstance(raw, dict) or not raw:
        raise ValueError("workflow mixtures entries must be non-empty mappings")
    unknown = sorted(set(str(key) for key in raw) - set(base_species))
    if unknown:
        raise ValueError(
            f"workflow mixture {mixture_id} has unknown species: {unknown}"
        )
    values: dict[str, float] = {}
    for species in base_species:
        fraction = float(raw.get(species, 0.0))
        if not math.isfinite(fraction) or fraction < 0.0:
            raise ValueError(
                f"workflow mixture {mixture_id} has invalid fraction for {species}"
            )
        values[species] = fraction
    total = sum(values.values())
    if total <= 0.0:
        raise ValueError(f"workflow mixture {mixture_id} fractions are all zero")
    return MixtureSpec(
        mixture_id=mixture_id,
        fractions={species: value / total for species, value in values.items()},
    )


def _load_workflow_raw(path: Path) -> tuple[dict[str, Any], Path]:
    workflow_path = path.resolve()
    return _read_mapping(workflow_path), workflow_path


def load_workflow(path: str | Path) -> WorkflowConfig:
    raw, workflow_path = _load_workflow_raw(Path(path))
    base_dir = workflow_path.parent
    allowed = {
        "base_config",
        "database",
        "e_over_n_Td",
        "mixtures",
        "mc",
        "quality",
        "execution",
        "mean_energy_support",
    }
    unknown = set(raw) - allowed
    if unknown:
        raise ValueError(f"Unsupported workflow fields: {sorted(unknown)}")

    base_config_path = _resolve_path(raw.get("base_config"), base_dir, "base_config")
    base_config = load_config(base_config_path)
    base_species = tuple(gas.species for gas in base_config.conditions.gas_mixture)
    if len(set(base_species)) != len(base_species):
        raise ValueError("base_config gas_mixture species must be unique")

    mixtures_raw = raw.get("mixtures")
    if not isinstance(mixtures_raw, list) or not mixtures_raw:
        raise ValueError("workflow mixtures must be a non-empty list")
    mixtures = tuple(
        _normalize_mixture(
            item,
            mixture_id=index,
            base_species=base_species,
        )
        for index, item in enumerate(mixtures_raw)
    )
    e_over_n_Td = _positive_float_sequence(raw.get("e_over_n_Td"), "e_over_n_Td")
    deterministic_enabled = any(
        item.enabled and item.id != "monte_carlo" for item in base_config.run.solvers
    )
    deterministic_execution = _parse_deterministic_execution(
        raw.get("execution"),
        deterministic_enabled=deterministic_enabled,
    )
    enabled_solvers = tuple(
        str(item.id) for item in base_config.run.solvers if item.enabled
    )
    mean_energy_support = _parse_mean_energy_support(
        raw.get("mean_energy_support"),
        e_over_n_Td=e_over_n_Td,
        mixtures=mixtures,
        enabled_solvers=enabled_solvers,
    )

    mc_raw = raw.get("mc", {}) or {}
    if not isinstance(mc_raw, dict):
        raise ValueError("workflow mc must be a mapping")
    mc_unknown = set(mc_raw) - {
        "replicas",
        "base_seed",
        "e_over_n_Td",
        "workers",
        "reuse_database",
        "sampling_plan",
        "convergence",
        "previous_decision",
        "previous_decision_sha256",
    }
    if mc_unknown:
        raise ValueError(f"Unsupported workflow mc fields: {sorted(mc_unknown)}")
    mc_replicas = _positive_int(mc_raw.get("replicas", 1), "workflow mc.replicas")
    mc_enabled = any(
        item.enabled and item.id == "monte_carlo" for item in base_config.run.solvers
    )
    raw_base_seed = mc_raw.get("base_seed")
    if mc_enabled and (
        "base_seed" not in mc_raw
        or isinstance(raw_base_seed, bool)
        or not isinstance(raw_base_seed, int)
    ):
        raise ValueError(
            "workflow mc.base_seed is required and must be an integer "
            "when monte_carlo is enabled"
        )
    if raw_base_seed is not None and (
        isinstance(raw_base_seed, bool) or not isinstance(raw_base_seed, int)
    ):
        raise ValueError("workflow mc.base_seed must be an integer")
    mc_base_seed = int(raw_base_seed) if raw_base_seed is not None else None
    mc_workers = _positive_int(mc_raw.get("workers", 1), "workflow mc.workers")
    mc_reuse_database = (
        _resolve_path(
            mc_raw["reuse_database"],
            base_dir,
            "workflow mc.reuse_database",
        )
        if mc_raw.get("reuse_database") is not None
        else None
    )
    mc_e_over_n_Td = (
        _positive_float_sequence(mc_raw["e_over_n_Td"], "mc.e_over_n_Td")
        if mc_raw.get("e_over_n_Td") is not None
        else None
    )
    mc_anchors = mc_e_over_n_Td or e_over_n_Td
    if len(set(mc_anchors)) != len(mc_anchors):
        raise ValueError("workflow MC E/N anchors must be unique")
    mc_sampling_plan = _parse_mc_sampling_plan(
        mc_raw.get("sampling_plan"),
        anchors=mc_anchors,
        base_config=base_config,
        default_replicas=mc_replicas,
    )
    convergence = (
        parse_convergence_policy(mc_raw["convergence"])
        if "convergence" in mc_raw
        else None
    )
    previous_decision = (
        _resolve_path(
            mc_raw["previous_decision"],
            base_dir,
            "mc.previous_decision",
        )
        if mc_raw.get("previous_decision") is not None
        else None
    )
    if previous_decision is not None and convergence is None:
        raise MonteCarloPolicyError("mc.previous_decision requires mc.convergence")
    if previous_decision is None and "previous_decision_sha256" in mc_raw:
        raise MonteCarloPolicyError(
            "mc.previous_decision_sha256 requires mc.previous_decision"
        )
    if previous_decision is not None and mc_raw.get(
        "previous_decision_sha256"
    ) != file_sha256(previous_decision):
        raise MonteCarloPolicyError(
            "previous MC decision is missing its hash or changed after workflow generation"
        )
    if convergence is not None:
        if not mc_enabled or len(mixtures) != 1:
            raise MonteCarloPolicyError(
                "a bounded MC campaign requires one mixture and an enabled "
                "monte_carlo solver"
            )
        validate_sampling_budget(
            mc_sampling_plan,
            convergence,
            previous_decision=previous_decision,
        )
        if previous_decision is not None:
            previous = read_selection(previous_decision)
            identity = previous.get("mc_run_identity", {})
            if (
                identity.get("base_config_sha256") != file_sha256(base_config_path)
                or identity.get("base_seed") != mc_base_seed
            ):
                raise MonteCarloPolicyError(
                    "MC config or base seed changed after the previous decision"
                )
            if previous.get("physical_context") != physical_context(base_config):
                raise MonteCarloPolicyError(
                    "MC physical context changed after the previous decision"
                )
            if identity.get("cross_sections_sha256") != _combined_hash(
                _cross_section_file_hashes(base_config)
            ):
                raise MonteCarloPolicyError(
                    "MC cross sections changed after the previous decision"
                )
            expected_species = [
                {
                    "species": gas.species,
                    "fraction": mixtures[0].fractions[gas.species],
                    "mass_amu": gas.mass_amu,
                }
                for gas in sorted(
                    base_config.conditions.gas_mixture,
                    key=lambda gas: gas.species,
                )
            ]
            if identity.get("mixture", {}).get("species") != expected_species:
                raise MonteCarloPolicyError(
                    "MC mixture changed after the previous decision"
                )

    return WorkflowConfig(
        path=workflow_path,
        base_config_path=base_config_path,
        database_path=_resolve_path(raw.get("database"), base_dir, "database"),
        e_over_n_Td=e_over_n_Td,
        mixtures=mixtures,
        mc_enabled=mc_enabled,
        mc_e_over_n_Td=mc_e_over_n_Td,
        mc_replicas=mc_replicas,
        mc_base_seed=mc_base_seed,
        mc_workers=mc_workers,
        mc_reuse_database=mc_reuse_database,
        mc_sampling_plan=mc_sampling_plan,
        mc_convergence=convergence,
        mc_previous_decision=previous_decision,
        quality=parse_quality_thresholds(raw.get("quality")),
        deterministic_execution=deterministic_execution,
        mean_energy_support=mean_energy_support,
    )


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as fp:
        for chunk in iter(lambda: fp.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _cross_section_file_hashes(config: SwarmConfig) -> dict[str, str]:
    return {
        str(file_cfg.path): _sha256_file(Path(file_cfg.path))
        for file_cfg in config.cross_sections.files
    }


def _combined_hash(values: dict[str, str]) -> str:
    digest = sha256()
    for key, value in sorted(values.items()):
        digest.update(key.encode("utf-8"))
        digest.update(b"\0")
        digest.update(value.encode("ascii"))
        digest.update(b"\0")
    return digest.hexdigest()
