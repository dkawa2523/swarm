"""One input document for calculation conditions and optional campaign controls."""

from copy import deepcopy
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.config_parser import load_config_from_raw
from electron_swarm.core.yaml_input import read_yaml_mapping


@dataclass(frozen=True, slots=True)
class CalculationDocument:
    path: Path
    template: SwarmConfig
    mixtures: tuple[dict[str, float], ...]
    campaign: dict[str, Any] | None

    def calculation_sha256(self) -> str:
        payload = asdict(self.template)
        # Publication and the document's location do not change a calculation.
        for key in ("source_path", "output", "comparison"):
            payload.pop(key)
        payload["mixtures"] = self.mixtures
        encoded = json.dumps(payload, default=str, sort_keys=True, allow_nan=False)
        return sha256(encoded.encode("utf-8")).hexdigest()


def read_calculation(path: str | Path) -> CalculationDocument:
    """Resolve aligned composition columns; each row is a complete gas mixture."""
    path = Path(path).resolve()
    raw = read_yaml_mapping(path)
    if "base_config" in raw:
        raise ValueError(
            "workflow migration required: inline the schema-v2 calculation, "
            "put execution controls under campaign, and remove base_config"
        )
    has_campaign = "campaign" in raw
    campaign = raw.pop("campaign", None)
    if has_campaign and not isinstance(campaign, dict):
        raise ValueError("campaign must be a mapping")
    conditions = raw.get("conditions")
    if not isinstance(conditions, dict):
        return CalculationDocument(path, load_config_from_raw(raw, path), (), campaign)
    gases = conditions.get("gas_mixture", [])
    if not isinstance(gases, list) or any(not isinstance(g, dict) for g in gases):
        # The product parser owns the ordinary shape/type error.
        return CalculationDocument(path, load_config_from_raw(raw, path), (), campaign)
    columns = [g.get("fraction", 1.0) for g in gases]
    lengths = {len(value) for value in columns if isinstance(value, list)}
    if lengths and (len(lengths) != 1 or 0 in lengths):
        raise ValueError("gas_mixture fraction columns must have the same nonzero length")
    if lengths and any(not isinstance(value, list) for value in columns):
        raise ValueError("use fraction arrays for every species in a mixture sweep")
    count = next(iter(lengths), 1)
    mixtures = []
    configurations = []
    for index in range(count):
        item = deepcopy(raw)
        for gas, column in zip(item.get("conditions", {}).get("gas_mixture", []), columns):
            gas["fraction"] = column[index] if isinstance(column, list) else column
        config = load_config_from_raw(item, path)
        configurations.append(config)
        mixtures.append({g.species: g.fraction for g in config.conditions.gas_mixture})
    if len({tuple(sorted(m.items())) for m in mixtures}) != len(mixtures):
        raise ValueError("gas_mixture fraction columns contain duplicate mixtures")
    return CalculationDocument(path, configurations[0], tuple(mixtures), campaign)


def load_case_template(path: str | Path) -> SwarmConfig:
    """Load the first resolved calculation used as the campaign's case template."""
    return read_calculation(path).template


def calculation_sha256(path: str | Path) -> str:
    return read_calculation(path).calculation_sha256()
