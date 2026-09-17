"""Load the strict run mapping layered over the GEC-ICP model contract."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import yaml

from .mapping import load_gec_icp_mapping
from .run_contracts import (
    GecIcpBundleSpec,
    GecIcpClosureSpec,
    GecIcpConvergenceSpec,
    GecIcpFunctionEedfSpec,
    GecIcpRunMapping,
    GecIcpRunSpec,
    GecIcpWorkflowError,
)


_ROOT_FIELDS = {
    "schema_version",
    "model_mapping",
    "bundle",
    "closure",
    "run",
    "results",
    "logs",
}
_BUNDLE_FIELDS = {
    "path",
    "expected_source",
    "expected_field_type",
    "expected_transport_definition",
    "expected_mc_qualification_profile",
}
_CLOSURE_FIELDS = {
    "source_field",
    "electron_transport",
    "reaction_model",
    "elastic_energy_loss_model",
    "chemistry_owner",
    "rf_owner",
    "function_eedf",
}
_FUNCTION_EEDF_FIELDS = {
    "table",
    "function_tag",
    "interpolation",
    "extrapolation",
}
_RUN_FIELDS = {
    "output_mph",
    "power_W",
    "frequency_Hz",
    "gas_temperature_K",
    "pressure_Pa",
    "final_time_s",
    "output_points_per_decade",
    "clear_saved_solution",
    "convergence",
}
_CONVERGENCE_FIELDS = {
    "electron_inventory_relative_change",
    "ion_inventory_relative_change",
    "metastable_inventory_relative_change",
    "mean_energy_relative_change",
    "absorbed_power_relative_change",
    "coil_power_relative_error",
    "mean_energy_support_relative_guard",
}
_PATH_FIELDS = {"path"}
_SOURCES = {"two_term", "monte_carlo", "propagator", "composite"}
_FUNCTION_TABLE = "eedf_f0_comsol_2d.csv"
_FUNCTION_INTERPOLATION = "structured_spreadsheet_linear_projection"


def load_gec_icp_run_mapping(
    mapping_path: str | Path,
    *,
    bundle_path: str | Path | None = None,
) -> GecIcpRunMapping:
    """Load a source-specific ICP execution map without legacy aliases."""

    path = Path(mapping_path).resolve()
    raw = _read_yaml(path)
    _reject_unknown(raw, _ROOT_FIELDS, "mapping")
    if isinstance(raw.get("schema_version"), bool) or raw.get("schema_version") != 2:
        raise GecIcpWorkflowError("GEC ICP run mapping requires schema_version: 2")
    root = path.parent
    model_mapping_path = _resolved_path(
        _required_string(raw, "model_mapping", "mapping"), root
    )
    model_mapping = load_gec_icp_mapping(model_mapping_path, validate_files=True)

    bundle_raw = _required_mapping(raw, "bundle", "mapping")
    closure_raw = _required_mapping(raw, "closure", "mapping")
    run_raw = _required_mapping(raw, "run", "mapping")
    results_raw = _required_mapping(raw, "results", "mapping")
    logs_raw = _required_mapping(raw, "logs", "mapping")
    _reject_unknown(bundle_raw, _BUNDLE_FIELDS, "bundle")
    _reject_unknown(closure_raw, _CLOSURE_FIELDS, "closure")
    _reject_unknown(run_raw, _RUN_FIELDS, "run")
    _reject_unknown(results_raw, _PATH_FIELDS, "results")
    _reject_unknown(logs_raw, _PATH_FIELDS, "logs")

    source = _required_choice(bundle_raw, "expected_source", "bundle", _SOURCES)
    selected_bundle = (
        Path(bundle_path).resolve()
        if bundle_path is not None
        else _resolved_path(_required_string(bundle_raw, "path", "bundle"), root)
    )
    mc_profile = bundle_raw.get("expected_mc_qualification_profile")
    if source == "monte_carlo":
        if mc_profile != "function_eedf_restricted_lmea":
            raise GecIcpWorkflowError(
                "Monte Carlo ICP mapping requires the "
                "function_eedf_restricted_lmea qualification profile"
            )
    elif source == "composite":
        if mc_profile != "function_eedf_restricted_lmea":
            raise GecIcpWorkflowError(
                "composite ICP mapping requires the "
                "function_eedf_restricted_lmea qualification profile"
            )
    elif mc_profile is not None:
        raise GecIcpWorkflowError(
            "expected_mc_qualification_profile is only valid for "
            "monte_carlo or composite"
        )

    function_raw = _required_mapping(closure_raw, "function_eedf", "closure")
    _reject_unknown(function_raw, _FUNCTION_EEDF_FIELDS, "closure.function_eedf")
    table = _required_string(function_raw, "table", "closure.function_eedf")
    interpolation = _required_string(
        function_raw, "interpolation", "closure.function_eedf"
    )
    if (table, interpolation) != (_FUNCTION_TABLE, _FUNCTION_INTERPOLATION):
        raise GecIcpWorkflowError(
            f"Function-EEDF must use {_FUNCTION_TABLE} with {_FUNCTION_INTERPOLATION}"
        )

    convergence_raw = _required_mapping(run_raw, "convergence", "run")
    _reject_unknown(convergence_raw, _CONVERGENCE_FIELDS, "run.convergence")
    output_mph = _resolved_path(_required_string(run_raw, "output_mph", "run"), root)
    if output_mph.suffix.lower() != ".mph":
        raise GecIcpWorkflowError("run.output_mph must name an .mph file")
    if output_mph == model_mapping.model.input_mph:
        raise GecIcpWorkflowError("run.output_mph must not overwrite the input MPH")

    return GecIcpRunMapping(
        path=path,
        root=root,
        model_mapping_path=model_mapping_path,
        model_mapping=model_mapping,
        bundle=GecIcpBundleSpec(
            path=selected_bundle,
            expected_source=source,  # type: ignore[arg-type]
            expected_field_type=_required_exact(
                bundle_raw, "expected_field_type", "bundle", "dc"
            ),
            expected_transport_definition=_required_string(
                bundle_raw, "expected_transport_definition", "bundle"
            ),
            expected_mc_qualification_profile=(
                str(mc_profile) if mc_profile is not None else None
            ),
        ),
        closure=GecIcpClosureSpec(
            source_field=_required_exact(
                closure_raw, "source_field", "closure", "steady_dc"
            ),
            electron_transport=_required_exact(
                closure_raw,
                "electron_transport",
                "closure",
                "swarm_mobility_comsol_einstein",
            ),
            reaction_model=_required_exact(
                closure_raw, "reaction_model", "closure", "function_eedf"
            ),
            elastic_energy_loss_model=_required_exact(
                closure_raw,
                "elastic_energy_loss_model",
                "closure",
                "comsol_cross_section_integral",
            ),
            chemistry_owner=_required_exact(
                closure_raw,
                "chemistry_owner",
                "closure",
                "comsol_embedded_cross_sections",
            ),
            rf_owner=_required_exact(
                closure_raw,
                "rf_owner",
                "closure",
                "comsol_frequency_transient",
            ),
            function_eedf=GecIcpFunctionEedfSpec(
                table=table,
                function_tag=_required_string(
                    function_raw, "function_tag", "closure.function_eedf"
                ),
                interpolation=interpolation,  # type: ignore[arg-type]
                extrapolation=_required_exact(
                    function_raw,
                    "extrapolation",
                    "closure.function_eedf",
                    "constant",
                ),
            ),
        ),
        run=GecIcpRunSpec(
            output_mph=output_mph,
            power_W=_positive_number(run_raw, "power_W", "run"),
            frequency_Hz=_positive_number(run_raw, "frequency_Hz", "run"),
            gas_temperature_K=_positive_number(run_raw, "gas_temperature_K", "run"),
            pressure_Pa=_positive_number(run_raw, "pressure_Pa", "run"),
            final_time_s=_positive_number(run_raw, "final_time_s", "run"),
            output_points_per_decade=_bounded_integer(
                run_raw, "output_points_per_decade", "run", 1, 20
            ),
            clear_saved_solution=_required_true(run_raw, "clear_saved_solution", "run"),
            convergence=GecIcpConvergenceSpec(
                **{
                    name: _fraction(convergence_raw, name, "run.convergence")
                    for name in sorted(_CONVERGENCE_FIELDS)
                }
            ),
        ),
        output_directory=_resolved_path(
            _required_string(results_raw, "path", "results"), root
        ),
        log_path=_resolved_path(_required_string(logs_raw, "path", "logs"), root),
    )


def _read_yaml(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise GecIcpWorkflowError(f"GEC ICP run mapping does not exist: {path}")
    try:
        raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise GecIcpWorkflowError(f"cannot read GEC ICP run mapping: {path}") from exc
    if not isinstance(raw, dict):
        raise GecIcpWorkflowError("GEC ICP run mapping root must be a mapping")
    return raw


def _required_mapping(raw: dict[str, Any], name: str, owner: str) -> dict[str, Any]:
    value = raw.get(name)
    if not isinstance(value, dict):
        raise GecIcpWorkflowError(f"{owner}.{name} must be a mapping")
    return value


def _required_string(raw: dict[str, Any], name: str, owner: str) -> str:
    value = raw.get(name)
    if not isinstance(value, str) or not value.strip():
        raise GecIcpWorkflowError(f"{owner}.{name} must be a nonempty string")
    return value.strip()


def _required_choice(
    raw: dict[str, Any], name: str, owner: str, choices: set[str]
) -> str:
    value = _required_string(raw, name, owner)
    if value not in choices:
        raise GecIcpWorkflowError(
            f"{owner}.{name} must be one of {sorted(choices)}, got {value!r}"
        )
    return value


def _required_exact(raw: dict[str, Any], name: str, owner: str, expected: str) -> Any:
    value = raw.get(name)
    if value != expected:
        raise GecIcpWorkflowError(f"{owner}.{name} must be {expected!r}, got {value!r}")
    return value


def _positive_number(raw: dict[str, Any], name: str, owner: str) -> float:
    value = raw.get(name)
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(float(value))
        or value <= 0
    ):
        raise GecIcpWorkflowError(f"{owner}.{name} must be a positive number")
    return float(value)


def _fraction(raw: dict[str, Any], name: str, owner: str) -> float:
    value = _positive_number(raw, name, owner)
    if value >= 1.0:
        raise GecIcpWorkflowError(f"{owner}.{name} must be less than 1")
    return value


def _bounded_integer(
    raw: dict[str, Any], name: str, owner: str, minimum: int, maximum: int
) -> int:
    value = raw.get(name)
    if isinstance(value, bool) or not isinstance(value, int):
        raise GecIcpWorkflowError(f"{owner}.{name} must be an integer")
    if not minimum <= value <= maximum:
        raise GecIcpWorkflowError(
            f"{owner}.{name} must be between {minimum} and {maximum}"
        )
    return value


def _required_true(raw: dict[str, Any], name: str, owner: str) -> bool:
    value = raw.get(name)
    if value is not True:
        raise GecIcpWorkflowError(f"{owner}.{name} must be true")
    return True


def _resolved_path(value: str, root: Path) -> Path:
    path = Path(value).expanduser()
    return (path if path.is_absolute() else root / path).resolve()


def _reject_unknown(raw: dict[str, Any], allowed: set[str], owner: str) -> None:
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise GecIcpWorkflowError(
            f"unsupported {owner} field(s): " + ", ".join(unknown)
        )


__all__ = ["load_gec_icp_run_mapping"]
