"""Load the repository GEC-CCP YAML mapping contract."""

from __future__ import annotations

from pathlib import Path

import yaml

from .config_values import (
    boolean_value,
    mapping_value,
    reject_unknown_keys,
    resolved_path,
)
from .contracts import GecCcpMapping, GecCcpWorkflowError, GecPathSpec
from .mapping_sections import (
    parse_bundle,
    parse_closure,
    parse_model,
    parse_reactions,
    parse_result,
    parse_run,
)


ROOT_KEYS = {
    "schema_version",
    "model",
    "bundle",
    "reactions",
    "closure",
    "initialization",
    "run",
    "results",
    "logs",
}


def load_gec_ccp_mapping(
    mapping_path: str | Path,
    *,
    bundle_path: str | Path | None = None,
) -> GecCcpMapping:
    path = Path(mapping_path).resolve()
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict):
        raise GecCcpWorkflowError("GEC CCP mapping root must be a mapping")
    reject_unknown_keys(raw, ROOT_KEYS, "mapping")
    if int(raw.get("schema_version", 0)) != 2:
        raise GecCcpWorkflowError(
            "GEC CCP mapping requires schema_version: 2; migrate run.powers_W "
            "to the single direct-solve setting run.power_W"
        )
    if "initialization" in raw:
        raise GecCcpWorkflowError(
            "initialization is obsolete; GEC CCP runs use native COMSOL "
            "Physics Initial Values without solution promotion or remapping"
        )

    root = path.parent
    run_raw = mapping_value(raw, "run")
    include_builtin_reference = boolean_value(
        run_raw, "include_builtin_reference", default=False
    )
    result = parse_result(mapping_value(raw, "results"), root)
    bundle = parse_bundle(
        mapping_value(raw, "bundle"),
        root,
        bundle_path=bundle_path,
        result_role=result.role,
    )
    closure = parse_closure(
        mapping_value(raw, "closure"),
        expected_source=bundle.expected_source,
        result_role=result.role,
    )
    logs_raw = mapping_value(raw, "logs")
    reject_unknown_keys(logs_raw, {"path"}, "logs")
    return GecCcpMapping(
        path=path,
        root=root,
        model=parse_model(
            mapping_value(raw, "model"),
            root,
            include_builtin_reference=include_builtin_reference,
        ),
        bundle=bundle,
        reactions=parse_reactions(raw.get("reactions")),
        closure=closure,
        run=parse_run(
            run_raw,
            root,
            include_builtin_reference=include_builtin_reference,
            source=bundle.expected_source,
            result_role=result.role,
            closure=closure,
        ),
        results=result,
        logs=GecPathSpec(
            resolved_path(root, logs_raw.get("path"), "logs.path")
        ),
    )


__all__ = ["load_gec_ccp_mapping"]
