"""Small helpers for optional collision-extension YAML sections."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def raw_config_from_source(config: object) -> dict[str, Any]:
    path = getattr(config, "source_path", None)
    if path is None:
        return {}
    try:
        data = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    except OSError:
        return {}
    return data if isinstance(data, dict) else {}

