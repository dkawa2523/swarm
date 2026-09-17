"""Schema-v2 YAML loading and product configuration composition."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.config_sections.common import (
    as_mapping_section,
    reject_unknown_fields,
)
from electron_swarm.core.config_sections.conditions import (
    parse_conditions,
    parse_cross_sections,
)
from electron_swarm.core.config_sections.physics import parse_physics
from electron_swarm.core.config_sections.run_output import (
    parse_comparison,
    parse_feature_policy,
    parse_output,
    parse_run,
)
from electron_swarm.core.config_sections.solvers import parse_solvers
from electron_swarm.core.config_validation import validate_config
from electron_swarm.core.legacy_schema import reject_removed_schema

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
OUTPUT_FIELDS = {"directory", "base_name", "float_format"}


def read_mapping(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as fp:
        data = yaml.safe_load(fp)
    if data is None:
        data = {}
    if not isinstance(data, dict):
        raise ValueError("YAML root must be a mapping")
    return data


def _validate_document_shape(raw: dict[str, Any]) -> None:
    reject_removed_schema(raw)
    reject_unknown_fields(raw, TOP_LEVEL_FIELDS, "top-level")
    out_raw = as_mapping_section(raw, "output", "output")
    reject_unknown_fields(out_raw, OUTPUT_FIELDS, "output")


def load_config_from_raw(raw: dict[str, Any], cfg_path: str | Path) -> SwarmConfig:
    cfg_path = Path(cfg_path).resolve()
    base = cfg_path.parent
    _validate_document_shape(raw)
    config = SwarmConfig(
        schema_version=2,
        run=parse_run(raw),
        conditions=parse_conditions(raw),
        cross_sections=parse_cross_sections(raw, base),
        physics=parse_physics(raw, base),
        solvers=parse_solvers(raw),
        comparison=parse_comparison(raw),
        feature_policy=parse_feature_policy(raw),
        output=parse_output(raw, base),
        source_path=cfg_path,
    )
    validate_config(config)
    return config


def load_config(path: str | Path) -> SwarmConfig:
    """Load and validate a schema-v2 product configuration."""

    cfg_path = Path(path).resolve()
    return load_config_from_raw(read_mapping(cfg_path), cfg_path)
