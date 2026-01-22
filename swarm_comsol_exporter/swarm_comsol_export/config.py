from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import yaml


def load_yaml(path: Union[str, Path]) -> Dict[str, Any]:
    path = Path(path)
    with path.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping/dict: {path}")
    return data


def deep_get(d: Dict[str, Any], keys: List[str], default: Any = None) -> Any:
    cur: Any = d
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


def deep_set(d: Dict[str, Any], keys: List[str], value: Any) -> None:
    cur: Any = d
    for k in keys[:-1]:
        if k not in cur or not isinstance(cur[k], dict):
            cur[k] = {}
        cur = cur[k]
    cur[keys[-1]] = value


@dataclass
class ComsolExportConfig:
    enabled: bool
    run_dir: Path
    output_dir: Optional[Path]
    input_cfg: Dict[str, Any]
    processing: Dict[str, Any]
    output: Dict[str, Any]

    @staticmethod
    def from_dict(cfg: Dict[str, Any], run_dir_override: Optional[Union[str, Path]] = None) -> "ComsolExportConfig":
        root = cfg.get("comsol_export", cfg)  # allow either top-level or nested
        enabled = bool(root.get("enabled", False))

        run_dir = Path(root.get("run_dir", "."))
        if run_dir_override is not None:
            run_dir = Path(run_dir_override)

        output_dir = root.get("output_dir", None)
        output_dir = None if output_dir in (None, "null") else Path(output_dir)

        input_cfg = dict(root.get("input", {}))
        processing = dict(root.get("processing", {}))
        output = dict(root.get("output", {}))

        return ComsolExportConfig(
            enabled=enabled,
            run_dir=run_dir,
            output_dir=output_dir,
            input_cfg=input_cfg,
            processing=processing,
            output=output,
        )
