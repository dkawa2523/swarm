"""External E/N and gas-mixture workflow sweeps."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any

import yaml

from electron_swarm import (
    GasComponent,
    RequestedSolverConfig,
    SolverId,
    SwarmConfig,
    load_config,
    run,
)

from .aggregate import (
    QualityThresholds,
    aggregate_workflow_results,
    parse_quality_thresholds,
    stable_mc_seed,
)
from .store import WorkflowStore


@dataclass(frozen=True, slots=True)
class MixtureSpec:
    mixture_id: int
    fractions: dict[str, float]


@dataclass(frozen=True, slots=True)
class WorkflowConfig:
    path: Path
    base_config_path: Path
    database_path: Path
    e_over_n_Td: tuple[float, ...]
    mixtures: tuple[MixtureSpec, ...]
    mc_replicas: int = 1
    mc_base_seed: int | None = None
    quality: QualityThresholds = field(default_factory=QualityThresholds)


@dataclass(frozen=True, slots=True)
class SweepSummary:
    database_path: Path
    mixtures: int
    e_over_n_values: int
    cases_written: int


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
        raise ValueError(f"workflow mixture {mixture_id} has unknown species: {unknown}")
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
    allowed = {"base_config", "database", "e_over_n_Td", "mixtures", "mc", "quality"}
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

    mc_raw = raw.get("mc", {}) or {}
    if not isinstance(mc_raw, dict):
        raise ValueError("workflow mc must be a mapping")
    mc_unknown = set(mc_raw) - {"replicas", "base_seed"}
    if mc_unknown:
        raise ValueError(f"Unsupported workflow mc fields: {sorted(mc_unknown)}")
    mc_replicas = int(mc_raw.get("replicas", 1))
    if mc_replicas <= 0:
        raise ValueError("workflow mc.replicas must be positive")
    mc_base_seed = (
        int(mc_raw["base_seed"]) if mc_raw.get("base_seed") is not None else None
    )

    return WorkflowConfig(
        path=workflow_path,
        base_config_path=base_config_path,
        database_path=_resolve_path(raw.get("database"), base_dir, "database"),
        e_over_n_Td=_positive_float_sequence(raw.get("e_over_n_Td"), "e_over_n_Td"),
        mixtures=mixtures,
        mc_replicas=mc_replicas,
        mc_base_seed=mc_base_seed,
        quality=parse_quality_thresholds(raw.get("quality")),
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


def _enabled_solver_ids(config: SwarmConfig) -> list[SolverId]:
    return [item.id for item in config.run.solvers if item.enabled]


def _case_prefix(base: SwarmConfig, mixture: MixtureSpec, replicate: int) -> str:
    return f"{base.run.case_prefix}_m{mixture.mixture_id:04d}_r{replicate:04d}"


def _clone_config(
    base: SwarmConfig,
    *,
    mixture: MixtureSpec,
    e_over_n_values: tuple[float, ...],
    solver_ids: list[SolverId],
    replicate: int,
    seed: int | None = None,
) -> SwarmConfig:
    config = deepcopy(base)
    by_species = {gas.species: gas for gas in base.conditions.gas_mixture}
    config.conditions.gas_mixture = [
        GasComponent(
            species=species,
            fraction=mixture.fractions[species],
            mass_amu=by_species[species].mass_amu,
        )
        for species in by_species
    ]
    config.run.e_over_n_Td = list(e_over_n_values)
    config.run.solvers = [
        RequestedSolverConfig(id=solver_id, enabled=True) for solver_id in solver_ids
    ]
    config.run.case_prefix = _case_prefix(base, mixture, replicate)
    if "monte_carlo" in solver_ids and seed is not None:
        config.solvers.monte_carlo.seed = seed
    return config


def _mixture_species_rows(
    base: SwarmConfig,
    mixture: MixtureSpec,
) -> list[tuple[str, float, float]]:
    mass_by_species = {gas.species: gas.mass_amu for gas in base.conditions.gas_mixture}
    return [
        (species, mixture.fractions[species], mass_by_species[species])
        for species in mass_by_species
    ]


def run_sweep(path: str | Path) -> SweepSummary:
    workflow = load_workflow(path)
    base_config = load_config(workflow.base_config_path)
    solver_ids = _enabled_solver_ids(base_config)
    if not solver_ids:
        raise ValueError("base_config run.solvers has no enabled solvers")

    xs_hashes = _cross_section_file_hashes(base_config)
    provenance = {
        "base_config_path": str(workflow.base_config_path),
        "base_config_sha256": _sha256_file(workflow.base_config_path),
        "cross_sections_sha256": _combined_hash(xs_hashes),
        "cross_section_files_json": json.dumps(xs_hashes, sort_keys=True),
    }

    cases_written = 0
    with WorkflowStore(workflow.database_path) as store:
        store.set_provenance(provenance)
        for mixture in workflow.mixtures:
            store.write_mixture(
                mixture.mixture_id,
                _mixture_species_rows(base_config, mixture),
            )

            grouped_solvers = [
                solver_id for solver_id in solver_ids if solver_id != "monte_carlo"
            ]
            if grouped_solvers:
                config = _clone_config(
                    base_config,
                    mixture=mixture,
                    e_over_n_values=workflow.e_over_n_Td,
                    solver_ids=grouped_solvers,
                    replicate=0,
                )
                for case in run(config, write=False).cases:
                    store.write_case(
                        mixture_id=mixture.mixture_id,
                        replicate=0,
                        case=case,
                    )
                    cases_written += 1

            if "monte_carlo" in solver_ids:
                for e_over_n in workflow.e_over_n_Td:
                    for replicate in range(workflow.mc_replicas):
                        config = _clone_config(
                            base_config,
                            mixture=mixture,
                            e_over_n_values=(e_over_n,),
                            solver_ids=["monte_carlo"],
                            replicate=replicate,
                            seed=stable_mc_seed(
                                base_seed=workflow.mc_base_seed,
                                mixture_id=mixture.mixture_id,
                                e_over_n_Td=e_over_n,
                                replicate=replicate,
                                solver="monte_carlo",
                            ),
                        )
                        for case in run(config, write=False).cases:
                            store.write_case(
                                mixture_id=mixture.mixture_id,
                                replicate=replicate,
                                case=case,
                            )
                            cases_written += 1
        aggregate_workflow_results(store.connection, workflow.quality)

    return SweepSummary(
        database_path=workflow.database_path,
        mixtures=len(workflow.mixtures),
        e_over_n_values=len(workflow.e_over_n_Td),
        cases_written=cases_written,
    )
