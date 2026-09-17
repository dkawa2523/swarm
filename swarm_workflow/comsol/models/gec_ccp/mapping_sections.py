"""Parse and validate the independent sections of a GEC-CCP mapping."""

from __future__ import annotations

import math
from pathlib import Path
import re
from typing import Any

import electron_swarm.solvers.monte_carlo.evidence as _mc_evidence
from electron_swarm.solvers.two_term.transport import (
    TWO_TERM_TEMPORAL_GROWTH_TRANSPORT_SCHEMA_VERSION,
)

from swarm_workflow.tables.contracts import (
    MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA,
)

from .config_values import (
    boolean_value,
    expected_dc_field_type,
    reject_unknown_keys,
    required_text,
    resolved_path,
)
from .contracts import (
    FUNCTION_EEDF_PREINTEGRATED_INELASTIC,
    FUNCTION_EEDF_INTERPOLATION,
    FUNCTION_EEDF_TABLE,
    GEC_COMSOL_ELASTIC_ENERGY_LOSS,
    GEC_EXTERNAL_ELASTIC_ENERGY_LOSS,
    GEC_LOOKUP_INTERPOLATION,
    GEC_RESTRICTED_TRANSPORT_CLOSURE,
    GEC_RESULT_ROLES,
    GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE,
    GecBundleSpec,
    GecCcpWorkflowError,
    GecClosureSpec,
    GecFunctionEedfSpec,
    GecModelSpec,
    GecReactionSpec,
    GecResultSpec,
    GecRunSpec,
)


MODEL_KEYS = {
    "input_mph",
    "baseline_output_mph",
    "external_output_mph",
    "external_initialization_mph",
    "component",
    "physics",
    "plasma_feature",
    "time_periodic_study",
    "time_periodic_feature",
    "external_time_periodic_study",
    "external_time_periodic_feature",
    "external_solution",
    "conversion_study",
    "expected_original_eedf",
}
BUNDLE_KEYS = {
    "path",
    "expected_source",
    "expected_field_type",
    "expected_transport_definition",
    "expected_mc_transport_estimator_schema",
    "expected_mc_tail_estimator_schema",
    "expected_mc_qualification_profile",
    "expected_two_term_transport_kernel_schema",
    "require_solver_selection",
}
CLOSURE_KEYS = {
    "mode",
    "source_field",
    "electron_transport",
    "reaction_model",
    "external_rate_processes",
    "lookup_jacobian",
    "interpolation",
    "function_eedf",
    "zero_field_isotropization_Td",
    "thermal_diffusion_model",
    "gradient_response_policy",
    "elastic_energy_loss_model",
}
RUN_KEYS = {
    "power_parameter",
    "power_W",
    "powers_W",
    "energy_damping",
    "nonlinear_method",
    "nonlinear_globalization",
    "initial_values",
    "source_stabilization",
    "reaction_source_stabilization",
    "axis_dataset",
    "radial_dataset",
    "period_dataset",
    "external_period_dataset",
    "phase_dataset",
    "waveform_dataset",
    "baseline_waveform_dataset",
    "external_waveform_dataset",
    "include_builtin_reference",
    "support_policy",
    "low_energy_guard_bundle",
    "validation_mean_energy_floor_eV",
}


def parse_result(raw: dict[str, Any], root: Path) -> GecResultSpec:
    reject_unknown_keys(raw, {"output_directory", "role"}, "results")
    role = required_text(raw, "role")
    if role not in GEC_RESULT_ROLES:
        raise GecCcpWorkflowError(
            "results.role must be physical_target, diagnostic_control, "
            "ablation, or historical_validation"
        )
    return GecResultSpec(
        output_directory=resolved_path(
            root, raw.get("output_directory"), "results.output_directory"
        ),
        role=role,
    )


def parse_model(
    raw: dict[str, Any], root: Path, *, include_builtin_reference: bool
) -> GecModelSpec:
    reject_unknown_keys(raw, MODEL_KEYS, "model")
    baseline_raw = raw.get("baseline_output_mph")
    if include_builtin_reference:
        baseline_output = resolved_path(root, baseline_raw, "model.baseline_output_mph")
    elif baseline_raw is not None:
        raise GecCcpWorkflowError(
            "model.baseline_output_mph is only valid when "
            "run.include_builtin_reference is true"
        )
    else:
        baseline_output = None
    if "external_initialization_mph" in raw:
        raise GecCcpWorkflowError(
            "model.external_initialization_mph is obsolete; external input "
            "is applied to a copy of model.input_mph"
        )
    return GecModelSpec(
        input_mph=resolved_path(root, raw.get("input_mph"), "model.input_mph"),
        baseline_output_mph=baseline_output,
        output_mph=resolved_path(
            root, raw.get("external_output_mph"), "model.external_output_mph"
        ),
        component=required_text(raw, "component"),
        physics=required_text(raw, "physics"),
        plasma_feature=required_text(raw, "plasma_feature"),
        study=required_text(raw, "time_periodic_study"),
        time_periodic_feature=required_text(raw, "time_periodic_feature"),
        external_study=required_text(raw, "external_time_periodic_study"),
        external_time_periodic_feature=required_text(
            raw, "external_time_periodic_feature"
        ),
        external_solution=required_text(raw, "external_solution"),
        conversion_study=required_text(raw, "conversion_study"),
        expected_original_eedf=required_text(raw, "expected_original_eedf"),
    )


def parse_reactions(raw: object) -> tuple[GecReactionSpec, ...]:
    if not isinstance(raw, list) or not raw:
        raise GecCcpWorkflowError("reactions must be a non-empty list")
    reactions = tuple(
        GecReactionSpec(
            name=required_text(item, "name"),
            feature=required_text(item, "feature"),
            process_type=required_text(item, "process_type"),
        )
        for item in raw
        if isinstance(item, dict)
    )
    if len(reactions) != len(raw):
        raise GecCcpWorkflowError("every reactions entry must be a mapping")
    for index, item in enumerate(raw):
        reject_unknown_keys(
            item, {"name", "feature", "process_type"}, f"reactions[{index}]"
        )
    actual = tuple((reaction.feature, reaction.process_type) for reaction in reactions)
    if actual != (
        ("eir1", "elastic"),
        ("eir2", "excitation"),
        ("eir3", "ionization"),
    ):
        raise GecCcpWorkflowError(
            "the argon GEC-CCP model requires reactions in the fixed order "
            "eir1:elastic, eir2:excitation, eir3:ionization"
        )
    return reactions


def _optional_text(raw: dict[str, Any], name: str) -> str | None:
    value = raw.get(name)
    return None if value is None else str(value).strip()


def _validate_mc_bundle_contract(
    *,
    source: str,
    result_role: str,
    estimator_schema: str | None,
    tail_schema: str | None,
    qualification_profile: str | None,
) -> None:
    supported_estimators = {
        _mc_evidence.DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
        _mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
    }
    if source == "monte_carlo" and result_role == "physical_target":
        if estimator_schema not in supported_estimators:
            raise GecCcpWorkflowError(
                "a Monte Carlo physical target must declare bundle."
                "expected_mc_transport_estimator_schema as a supported "
                "direct MC estimator"
            )
        if qualification_profile is not None and qualification_profile not in {
            "full_transport",
            MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA,
        }:
            raise GecCcpWorkflowError(
                "bundle.expected_mc_qualification_profile is unsupported"
            )
        if (
            tail_schema is not None
            and tail_schema != _mc_evidence.WEIGHTED_MC_TAIL_ESTIMATOR_SCHEMA_VERSION
        ):
            raise GecCcpWorkflowError(
                "bundle.expected_mc_tail_estimator_schema is unsupported"
            )
    elif source != "monte_carlo" and estimator_schema is not None:
        raise GecCcpWorkflowError(
            "bundle.expected_mc_transport_estimator_schema is only valid "
            "for a Monte Carlo source"
        )
    if source != "monte_carlo" and qualification_profile is not None:
        raise GecCcpWorkflowError(
            "bundle.expected_mc_qualification_profile is only valid for a "
            "Monte Carlo source"
        )
    if source != "monte_carlo" and tail_schema is not None:
        raise GecCcpWorkflowError(
            "bundle.expected_mc_tail_estimator_schema is only valid for a "
            "Monte Carlo source"
        )


def parse_bundle(
    raw: dict[str, Any],
    root: Path,
    *,
    bundle_path: str | Path | None,
    result_role: str,
) -> GecBundleSpec:
    reject_unknown_keys(raw, BUNDLE_KEYS, "bundle")
    source = required_text(raw, "expected_source")
    if source not in {"two_term", "propagator", "monte_carlo"}:
        raise GecCcpWorkflowError(
            "bundle.expected_source must be two_term, propagator, or monte_carlo"
        )
    estimator_schema = _optional_text(raw, "expected_mc_transport_estimator_schema")
    tail_schema = _optional_text(raw, "expected_mc_tail_estimator_schema")
    qualification_profile = _optional_text(raw, "expected_mc_qualification_profile")
    two_term_schema = _optional_text(raw, "expected_two_term_transport_kernel_schema")
    _validate_mc_bundle_contract(
        source=source,
        result_role=result_role,
        estimator_schema=estimator_schema,
        tail_schema=tail_schema,
        qualification_profile=qualification_profile,
    )
    if two_term_schema is not None:
        if source != "two_term":
            raise GecCcpWorkflowError(
                "bundle.expected_two_term_transport_kernel_schema is only "
                "valid for a two_term source"
            )
        if two_term_schema != TWO_TERM_TEMPORAL_GROWTH_TRANSPORT_SCHEMA_VERSION:
            raise GecCcpWorkflowError(
                f"unsupported two-term transport-kernel schema: {two_term_schema}"
            )
    return GecBundleSpec(
        path=(
            Path(bundle_path).resolve()
            if bundle_path is not None
            else resolved_path(root, raw.get("path"), "bundle.path")
        ),
        expected_source=source,
        expected_field_type=expected_dc_field_type(raw),
        expected_transport_definition=required_text(
            raw, "expected_transport_definition"
        ),
        expected_mc_transport_estimator_schema=estimator_schema,
        expected_mc_tail_estimator_schema=tail_schema,
        expected_mc_qualification_profile=qualification_profile,
        expected_two_term_transport_kernel_schema=two_term_schema,
        require_solver_selection=boolean_value(
            raw, "require_solver_selection", default=False
        ),
    )


def _parse_transport(
    raw: dict[str, Any], expected_source: str
) -> tuple[str, str, str, float | None]:
    if "mode" in raw:
        raise GecCcpWorkflowError(
            "closure.mode is obsolete; migrate independently to "
            "closure.electron_transport and closure.reaction_model"
        )
    transport = required_text(raw, "electron_transport")
    if transport == "swarm_field_aligned_full":
        raise GecCcpWorkflowError(
            "closure.electron_transport=swarm_field_aligned_full is obsolete; "
            "migrate to comsol_specify_all_restricted and explicitly select "
            "thermal_diffusion_model"
        )
    if transport == "all_swarm_coefficients":
        raise GecCcpWorkflowError(
            "closure.electron_transport=all_swarm_coefficients is obsolete; "
            "use comsol_specify_all_restricted and explicitly declare its "
            "gradient-response assumptions"
        )
    if transport not in {
        "comsol",
        "swarm_mobility_einstein",
        GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE,
        GEC_RESTRICTED_TRANSPORT_CLOSURE,
    }:
        raise GecCcpWorkflowError(
            "closure.electron_transport must be comsol, "
            "swarm_mobility_einstein, swarm_hybrid_einstein_de, or "
            "comsol_specify_all_restricted"
        )
    thermal = required_text(raw, "thermal_diffusion_model")
    if thermal not in {
        "off_restricted_diagonal",
        "comsol_grad_diffusivity",
    }:
        raise GecCcpWorkflowError(
            "closure.thermal_diffusion_model must be "
            "off_restricted_diagonal or comsol_grad_diffusivity"
        )
    if thermal == "comsol_grad_diffusivity" and (
        transport != GEC_RESTRICTED_TRANSPORT_CLOSURE
    ):
        raise GecCcpWorkflowError(
            "closure.thermal_diffusion_model=comsol_grad_diffusivity is "
            "only valid with electron_transport=comsol_specify_all_restricted"
        )
    gradient = required_text(raw, "gradient_response_policy")
    if gradient not in {"standard_local_energy", "require_full"}:
        raise GecCcpWorkflowError(
            "closure.gradient_response_policy must be "
            "standard_local_energy or require_full"
        )
    if gradient == "require_full":
        raise GecCcpWorkflowError(
            "gradient_response_policy=require_full is unavailable: no "
            "current GEC external closure identifies the independent 2x2 "
            "density/energy-gradient response matrix"
        )
    zero_raw = raw.get("zero_field_isotropization_Td")
    if (
        transport
        in {
            GEC_RESTRICTED_TRANSPORT_CLOSURE,
            GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE,
        }
        and expected_source == "monte_carlo"
    ):
        try:
            zero_field = float(zero_raw)
        except (TypeError, ValueError) as exc:
            raise GecCcpWorkflowError(
                "closure.zero_field_isotropization_Td must be a positive "
                "number for Monte Carlo field-aligned energy transport"
            ) from exc
        if not math.isfinite(zero_field) or zero_field <= 0.0:
            raise GecCcpWorkflowError(
                "closure.zero_field_isotropization_Td must be positive and finite"
            )
    elif zero_raw is not None:
        raise GecCcpWorkflowError(
            "closure.zero_field_isotropization_Td is only valid for a "
            "monte_carlo restricted or hybrid field-aligned tensor"
        )
    else:
        zero_field = None
    return transport, thermal, gradient, zero_field


def _parse_external_rate_processes(
    raw: object, *, reaction_model: str, result_role: str
) -> tuple[str, ...]:
    if raw is None:
        return (
            ("elastic", "excitation", "ionization")
            if reaction_model == "external_rates"
            else ()
        )
    if reaction_model != "external_rates":
        raise GecCcpWorkflowError(
            "closure.external_rate_processes is only valid with "
            "reaction_model=external_rates"
        )
    if not isinstance(raw, list) or not all(
        isinstance(item, str) and item.strip() for item in raw
    ):
        raise GecCcpWorkflowError(
            "closure.external_rate_processes must be a non-empty list"
        )
    processes = tuple(item.strip() for item in raw)
    if (
        not processes
        or len(set(processes)) != len(processes)
        or not set(processes) <= {"elastic", "excitation", "ionization"}
    ):
        raise GecCcpWorkflowError(
            "closure.external_rate_processes must contain unique elastic, "
            "excitation, or ionization process ids"
        )
    if set(processes) != {"elastic", "excitation", "ionization"} and (
        result_role != "ablation"
    ):
        raise GecCcpWorkflowError(
            "a partial external-rate selection is restricted to results.role=ablation"
        )
    return processes


def _parse_elastic_energy_loss(
    raw: dict[str, Any],
    *,
    reaction_model: str,
    electron_transport: str,
    expected_source: str,
) -> str:
    model = str(
        raw.get("elastic_energy_loss_model", GEC_COMSOL_ELASTIC_ENERGY_LOSS)
    ).strip()
    if model not in {
        GEC_COMSOL_ELASTIC_ENERGY_LOSS,
        GEC_EXTERNAL_ELASTIC_ENERGY_LOSS,
    }:
        raise GecCcpWorkflowError(
            "closure.elastic_energy_loss_model must be comsol_mratio or "
            "external_solver_native"
        )
    if model != GEC_EXTERNAL_ELASTIC_ENERGY_LOSS:
        return model
    incompatibilities: list[str] = []
    if reaction_model not in {
        "external_rates",
        "function_eedf",
        FUNCTION_EEDF_PREINTEGRATED_INELASTIC,
    }:
        incompatibilities.append(
            "reaction_model must use the external Swarm EEDF or rates"
        )
    if electron_transport not in {
        "swarm_mobility_einstein",
        GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE,
        GEC_RESTRICTED_TRANSPORT_CLOSURE,
    }:
        incompatibilities.append(
            "electron_transport must consume external Swarm transport"
        )
    if expected_source not in {"two_term", "propagator", "monte_carlo"}:
        incompatibilities.append(
            "bundle source must be two_term, propagator, or monte_carlo"
        )
    if incompatibilities:
        raise GecCcpWorkflowError(
            "external_solver_native elastic energy loss is incompatible: "
            + "; ".join(incompatibilities)
        )
    return model


def _parse_function_eedf(
    raw: object, *, reaction_model: str, expected_source: str
) -> GecFunctionEedfSpec | None:
    if reaction_model not in {
        "function_eedf",
        FUNCTION_EEDF_PREINTEGRATED_INELASTIC,
    }:
        if raw is not None:
            raise GecCcpWorkflowError(
                "closure.function_eedf is only valid for "
                "closure.reaction_model=function_eedf or "
                "function_eedf_preintegrated_inelastic"
            )
        return None
    if not isinstance(raw, dict):
        raise GecCcpWorkflowError(
            "closure.function_eedf must be a mapping for a Function-EEDF reaction model"
        )
    reject_unknown_keys(
        raw,
        {"table", "function_tag", "interpolation", "extrapolation"},
        "closure.function_eedf",
    )
    table = required_text(raw, "table")
    table_path = Path(table)
    if table_path.is_absolute() or ".." in table_path.parts:
        raise GecCcpWorkflowError(
            "closure.function_eedf.table must be a relative path inside "
            "the COMSOL bundle"
        )
    if table != FUNCTION_EEDF_TABLE:
        raise GecCcpWorkflowError(
            "closure.function_eedf.table must select the canonical artifact "
            f"{FUNCTION_EEDF_TABLE}; the wide C1 table "
            "is offline evidence only"
        )
    function_tag = required_text(raw, "function_tag")
    if re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", function_tag) is None:
        raise GecCcpWorkflowError(
            "closure.function_eedf.function_tag must be a COMSOL-safe tag"
        )
    interpolation = str(raw.get("interpolation", "")).strip()
    if interpolation != FUNCTION_EEDF_INTERPOLATION:
        raise GecCcpWorkflowError(
            f"closure.function_eedf.interpolation must be {FUNCTION_EEDF_INTERPOLATION}"
        )
    extrapolation = str(raw.get("extrapolation", "constant")).strip()
    if extrapolation != "constant":
        raise GecCcpWorkflowError(
            "closure.function_eedf.extrapolation must be constant"
        )
    return GecFunctionEedfSpec(
        table=table,
        function_tag=function_tag,
        interpolation=FUNCTION_EEDF_INTERPOLATION,
        extrapolation="constant",
    )


def parse_closure(
    raw: dict[str, Any], *, expected_source: str, result_role: str
) -> GecClosureSpec:
    reject_unknown_keys(raw, CLOSURE_KEYS, "closure")
    transport, thermal, gradient, zero_field = _parse_transport(raw, expected_source)
    reaction_model = required_text(raw, "reaction_model")
    if reaction_model not in {
        "comsol_eedf",
        "external_rates",
        "function_eedf",
        FUNCTION_EEDF_PREINTEGRATED_INELASTIC,
    }:
        raise GecCcpWorkflowError(
            "closure.reaction_model must be comsol_eedf, external_rates, "
            "function_eedf, or function_eedf_preintegrated_inelastic"
        )
    processes = _parse_external_rate_processes(
        raw.get("external_rate_processes"),
        reaction_model=reaction_model,
        result_role=result_role,
    )
    elastic_model = _parse_elastic_energy_loss(
        raw,
        reaction_model=reaction_model,
        electron_transport=transport,
        expected_source=expected_source,
    )
    source_field = str(raw.get("source_field", "steady_dc")).strip()
    if source_field not in {"steady_dc", "rf_time_periodic"}:
        raise GecCcpWorkflowError(
            "closure.source_field must be steady_dc or rf_time_periodic"
        )
    if source_field != "steady_dc":
        raise GecCcpWorkflowError(
            "closure.source_field=rf_time_periodic is not supported by the "
            "single-valued mean-energy GEC closure; use steady_dc"
        )
    lookup_jacobian = required_text(raw, "lookup_jacobian")
    if lookup_jacobian != "exact":
        raise GecCcpWorkflowError(
            "closure.lookup_jacobian must be exact; lagged/nojac closure "
            "coupling is not a supported product mode"
        )
    interpolation = str(raw.get("interpolation", GEC_LOOKUP_INTERPOLATION)).strip()
    if interpolation in {"piecewise_cubic", "differentiable_log_cubic_spline"}:
        raise GecCcpWorkflowError(
            f"closure.interpolation={interpolation} is obsolete; use "
            f"{GEC_LOOKUP_INTERPOLATION}"
        )
    if interpolation != GEC_LOOKUP_INTERPOLATION:
        raise GecCcpWorkflowError(
            f"closure.interpolation must be {GEC_LOOKUP_INTERPOLATION}"
        )
    closure = GecClosureSpec(
        source_field=source_field,
        electron_transport=transport,
        zero_field_isotropization_Td=zero_field,
        thermal_diffusion_model=thermal,
        gradient_response_policy=gradient,
        reaction_model=reaction_model,
        external_rate_processes=processes,
        elastic_energy_loss_model=elastic_model,
        lookup_jacobian=lookup_jacobian,
        interpolation=interpolation,
        function_eedf=_parse_function_eedf(
            raw.get("function_eedf"),
            reaction_model=reaction_model,
            expected_source=expected_source,
        ),
    )
    if expected_source == "propagator" and (
        transport != "swarm_mobility_einstein"
        or reaction_model != "function_eedf"
        or elastic_model != GEC_EXTERNAL_ELASTIC_ENERGY_LOSS
    ):
        raise GecCcpWorkflowError(
            "propagator P1 is restricted to Function-EEDF, external native "
            "elastic energy loss, and swarm_mobility_einstein transport"
        )
    return closure


def _parse_positive_power(raw: dict[str, Any]) -> float:
    if "powers_W" in raw:
        raise GecCcpWorkflowError(
            "run.powers_W is obsolete; use one direct-solve value in run.power_W"
        )
    try:
        value = float(raw["power_W"])
    except (KeyError, TypeError, ValueError) as exc:
        raise GecCcpWorkflowError("run.power_W must be a positive number") from exc
    if not math.isfinite(value) or value <= 0.0:
        raise GecCcpWorkflowError("run.power_W must be a positive finite value")
    return value


def _parse_support_policy(
    raw: dict[str, Any],
    root: Path,
    *,
    source: str,
    result_role: str,
    closure: GecClosureSpec,
) -> tuple[str, Path | None]:
    policy = str(raw.get("support_policy", "strict")).strip()
    if policy not in {"strict", "qualified_guard_sensitivity"}:
        raise GecCcpWorkflowError(
            "run.support_policy must be strict or qualified_guard_sensitivity"
        )
    guard = (
        resolved_path(
            root,
            raw.get("low_energy_guard_bundle"),
            "run.low_energy_guard_bundle",
        )
        if raw.get("low_energy_guard_bundle") is not None
        else None
    )
    if (policy == "qualified_guard_sensitivity") != (guard is not None):
        raise GecCcpWorkflowError(
            "run.support_policy=qualified_guard_sensitivity requires exactly "
            "one run.low_energy_guard_bundle"
        )
    production_target = bool(
        source == "monte_carlo"
        and closure.reaction_model == "external_rates"
        and closure.electron_transport == "swarm_mobility_einstein"
        and closure.elastic_energy_loss_model == GEC_EXTERNAL_ELASTIC_ENERGY_LOSS
    )
    ablation_target = bool(
        source == "monte_carlo"
        and result_role == "ablation"
        and (
            closure.reaction_model == "external_rates"
            or closure.electron_transport == "swarm_mobility_einstein"
            or closure.elastic_energy_loss_model == GEC_EXTERNAL_ELASTIC_ENERGY_LOSS
        )
    )
    if policy == "qualified_guard_sensitivity" and not (
        production_target or ablation_target
    ):
        raise GecCcpWorkflowError(
            "qualified_guard_sensitivity is restricted to the Monte Carlo "
            "external-rate mobility-Einstein restricted-LMEA target and "
            "its explicit one-factor ablations"
        )
    return policy, guard


def _parse_validation_floor(
    raw: dict[str, Any], *, result_role: str, source: str
) -> float | None:
    value = raw.get("validation_mean_energy_floor_eV")
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise GecCcpWorkflowError(
            "run.validation_mean_energy_floor_eV must be a positive number"
        ) from exc
    if not math.isfinite(result) or result <= 0.0:
        raise GecCcpWorkflowError(
            "run.validation_mean_energy_floor_eV must be positive and finite"
        )
    if result_role != "historical_validation" or source != "monte_carlo":
        raise GecCcpWorkflowError(
            "run.validation_mean_energy_floor_eV is restricted to a "
            "Monte Carlo historical_validation run"
        )
    return result


def parse_run(
    raw: dict[str, Any],
    root: Path,
    *,
    include_builtin_reference: bool,
    source: str,
    result_role: str,
    closure: GecClosureSpec,
) -> GecRunSpec:
    reject_unknown_keys(raw, RUN_KEYS, "run")
    nonlinear = str(raw.get("nonlinear_globalization", "model_defined")).strip()
    if nonlinear not in {
        "model_defined",
        "automatic_newton_no_recovery",
        "double_dogleg",
    }:
        raise GecCcpWorkflowError(
            "run.nonlinear_globalization must be model_defined, "
            "automatic_newton_no_recovery, or double_dogleg"
        )
    policy, guard = _parse_support_policy(
        raw,
        root,
        source=source,
        result_role=result_role,
        closure=closure,
    )
    obsolete = sorted(
        field
        for field in ("energy_damping", "nonlinear_method", "initial_values")
        if field in raw
    )
    if obsolete:
        raise GecCcpWorkflowError(
            "obsolete GEC run settings: "
            + ", ".join(obsolete)
            + "; use native Physics Initial Values and the model-defined "
            "solver settings"
        )
    if "waveform_dataset" in raw:
        raise GecCcpWorkflowError(
            "run.waveform_dataset is ambiguous; use run.external_waveform_dataset"
        )
    baseline_waveform_raw = raw.get("baseline_waveform_dataset")
    baseline_period_raw = raw.get("period_dataset")
    if include_builtin_reference:
        baseline_period = required_text(raw, "period_dataset")
        baseline_waveform = required_text(raw, "baseline_waveform_dataset")
    elif baseline_waveform_raw is not None or baseline_period_raw is not None:
        raise GecCcpWorkflowError(
            "run.period_dataset and run.baseline_waveform_dataset are only "
            "valid when run.include_builtin_reference is true"
        )
    else:
        baseline_period = None
        baseline_waveform = None
    run = GecRunSpec(
        power_parameter=required_text(raw, "power_parameter"),
        power_W=_parse_positive_power(raw),
        include_builtin_reference=include_builtin_reference,
        nonlinear_globalization=nonlinear,
        source_stabilization=boolean_value(raw, "source_stabilization", default=False),
        reaction_source_stabilization=boolean_value(
            raw, "reaction_source_stabilization", default=False
        ),
        axis_dataset=required_text(raw, "axis_dataset"),
        radial_dataset=required_text(raw, "radial_dataset"),
        period_dataset=baseline_period,
        external_period_dataset=required_text(raw, "external_period_dataset"),
        phase_dataset=required_text(raw, "phase_dataset"),
        baseline_waveform_dataset=baseline_waveform,
        external_waveform_dataset=required_text(raw, "external_waveform_dataset"),
        support_policy=policy,
        low_energy_guard_bundle=guard,
        validation_mean_energy_floor_eV=_parse_validation_floor(
            raw, result_role=result_role, source=source
        ),
    )
    valid_external_waveforms = {"dset5", run.phase_dataset}
    if run.baseline_waveform_dataset is not None:
        valid_external_waveforms.add(run.baseline_waveform_dataset)
    if run.external_waveform_dataset not in valid_external_waveforms:
        raise GecCcpWorkflowError(
            "run.external_waveform_dataset must be the native periodic "
            "waveform dataset, dset5 for a cloned study, or the external "
            "phase_dataset for direct periodic-solution export"
        )
    return run
