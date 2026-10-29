"""Unambiguous YAML input shared by product and workflow entry points."""

from pathlib import Path
from typing import Any

import yaml


class _UniqueKeyLoader(yaml.SafeLoader):
    def construct_mapping(self, node, deep=False):
        result = {}
        for key_node, value_node in node.value:
            if key_node.tag == "tag:yaml.org,2002:merge":
                raise ValueError("YAML merge keys are ambiguous; write explicit fields")
            key = self.construct_object(key_node, deep=deep)
            if not isinstance(key, str):
                raise ValueError("YAML configuration mapping keys must be strings")
            if key in result:
                raise ValueError(
                    f"duplicate YAML key {key!r} at line {key_node.start_mark.line + 1}"
                )
            result[key] = self.construct_object(value_node, deep=deep)
        return result


def read_yaml_mapping(path: str | Path) -> dict[str, Any]:
    source = Path(path)
    try:
        data = yaml.load(source.read_text(encoding="utf-8"), Loader=_UniqueKeyLoader)
    except (yaml.YAMLError, ValueError) as exc:
        raise ValueError(f"invalid YAML input {source}: {exc}") from exc
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping: {source}")
    return data
