"""Generate COMSOL Java sources for the canonical GEC-CCP workflow."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from swarm_workflow.comsol.java import java_path, java_string
import swarm_workflow.comsol.models.gec_ccp.closure_arguments as closure_arguments
import swarm_workflow.comsol.models.gec_ccp.java_apply_sections as java_apply_sections
import swarm_workflow.comsol.models.gec_ccp.java_export_sections as java_export_sections
import swarm_workflow.comsol.models.gec_ccp.java_support as java_support
from swarm_workflow.comsol.models.gec_ccp.closure import (
    _reaction_uses_preintegrated_rate,
    _uses_external_elastic_energy_loss,
    _uses_external_rates,
    _uses_transport,
)
from swarm_workflow.comsol.models.gec_ccp.contracts import (
    GecCcpMapping,
    GecCcpWorkflowError,
)
from swarm_workflow.comsol.models.gec_ccp.validation.joint_consistency import (
    _active_closure_mean_energy_support,
)


APPLY_CLASS = "SwarmGecCcpApply"
BASELINE_RUN_CLASS = "SwarmGecCcpBaselineRun"
EXTERNAL_RUN_CLASS = "SwarmGecCcpExternalRun"
BASELINE_EXPORT_CLASS = "SwarmGecCcpBaselineExport"
EXTERNAL_EXPORT_CLASS = "SwarmGecCcpExternalExport"
NATIVE_EEDF_AUDIT_CLASS = "SwarmGecCcpNativeEedfAudit"


def generate_apply_java(
    mapping: GecCcpMapping,
    *,
    input_mph: Path | None = None,
    preintegrated_rate_rows: list[dict[str, Any]] | None = None,
) -> str:
    context = java_apply_sections.build_apply_java_context(
        mapping,
        input_mph=input_mph,
        preintegrated_rate_rows=preintegrated_rate_rows,
    )
    lines = java_support.java_header(APPLY_CLASS, context.input_mph, "swarmGecApply")
    renderers = (
        java_apply_sections.render_stabilization,
        java_apply_sections.render_transport_mode,
        java_apply_sections.render_closure_arguments,
        java_apply_sections.render_transport_functions,
        java_apply_sections.render_elastic_interpolation,
        java_apply_sections.render_transport_bindings,
        java_apply_sections.render_reaction_bindings,
        java_apply_sections.render_external_elastic_owner,
        java_apply_sections.render_function_eedf,
        java_apply_sections.render_apply_footer,
    )
    for render in renderers:
        lines.extend(render(context))
    return "\n".join(lines)


def generate_run_java(
    mapping: GecCcpMapping,
    *,
    run_role: Literal["baseline", "external"],
    class_name: str,
    input_mph: Path,
    output_mph: Path,
    study: str | None = None,
    solution: str = "sol1",
    conversion_source_study: str | None = None,
    time_periodic_feature: str | None = None,
    convert_periodic_solution: bool = True,
) -> str:
    study = study or mapping.model.study
    conversion_source_study = conversion_source_study or study
    if time_periodic_feature is None:
        time_periodic_feature = (
            mapping.model.external_time_periodic_feature
            if run_role == "external"
            else mapping.model.time_periodic_feature
        )
    lines = java_support.java_header(class_name, input_mph, class_name)
    if run_role == "external":
        lines.extend(java_export_sections.render_elastic_domain_assertion(mapping))
    direct_solve_settings = [
        "    // Make both optional source stabilizations explicit. OFF/OFF",
        "    // is the schema default and was sufficient in the cold audit.",
        f"    model.component({java_string(mapping.model.component)})"
        f".physics({java_string(mapping.model.physics)})"
        '.prop("Stabilization").set("SourceStabilization", '
        + java_support.java_boolean(mapping.run.source_stabilization)
        + ");",
        f"    model.component({java_string(mapping.model.component)})"
        f".physics({java_string(mapping.model.physics)})"
        '.prop("Stabilization").set("ReactionSourceStabilization", '
        + java_support.java_boolean(mapping.run.reaction_source_stabilization)
        + ");",
    ]
    if mapping.run.nonlinear_globalization == "double_dogleg":
        nonlinear_solver = (
            f'model.sol({java_string(solution)}).feature("s1").feature("fc1")'
        )
        direct_solve_settings.extend(
            [
                "    // Apply the same equation-preserving trust-region",
                "    // globalization to baseline and external cold solves.",
                "    // No damping schedule, continuation, or source term is",
                "    // introduced; COMSOL retains its Double-Dogleg defaults.",
                f'    {nonlinear_solver}.set("dtech", "ddog");',
                f'    {nonlinear_solver}.set("resscale", "scalefieldwise");',
            ]
        )
    elif mapping.run.nonlinear_globalization == "automatic_newton_no_recovery":
        nonlinear_solver = (
            f'model.sol({java_string(solution)}).feature("s1").feature("fc1")'
        )
        direct_solve_settings.extend(
            [
                "    // Diagnostic, equation-preserving recovery ablation.",
                "    // Keep COMSOL's model-defined Newton method, damping,",
                "    // iteration limit, equations, and initial values.",
                f'    {nonlinear_solver}.set("useminsteprecovery", "off");',
            ]
        )
    uses_exact_periodic_locks = solution in {"sol1", mapping.model.external_solution}
    if uses_exact_periodic_locks:
        direct_solve_settings.extend(
            [
                "    // Do not change any further nonlinear solver settings.",
                "    // Keep COMSOL's native heavy-species equations, mass-",
                "    // fraction coordinate, boundary row, and initial values.",
                "    // The native cold-direct control converges without a",
                "    // heavy-species transformation or weak-form override.",
                "    // Every heavy species in this pure-Ar model has the",
                "    // same 0.04 kg/mol mass. Lock the exact algebraic",
                "    // identity to remove (1-exp(W))+exp(W) cancellation",
                "    // at nonphysical Newton trial points; no physical term",
                "    // or accepted-state value is changed.",
                f"    model.component({java_string(mapping.model.component)})"
                f".physics({java_string(mapping.model.physics)})"
                f".feature({java_string(mapping.model.plasma_feature)})"
                '.featureInfo("info").set("ptp.Mn", '
                'new String[]{"0.04[kg/mol]"});',
                "    // COMSOL's time-periodic unknowns are logarithmic:",
                "    // ne=exp(Ne_per), en=exp(En_per). Preserve the exact",
                "    // mean-energy identities without evaluating exp/exp,",
                "    // which otherwise becomes 0/0 at damped Newton trials.",
                f"    model.component({java_string(mapping.model.component)})"
                f".physics({java_string(mapping.model.physics)})"
                f".feature({java_string(mapping.model.plasma_feature)})"
                '.featureInfo("info").set("ptp.ebar", '
                'new String[]{"exp(En_per-Ne_per)*1[V]"});',
                f"    model.component({java_string(mapping.model.component)})"
                f".physics({java_string(mapping.model.physics)})"
                f".feature({java_string(mapping.model.plasma_feature)})"
                '.featureInfo("info").set("ptp.Te", '
                'new String[]{"2*ptp.ebar/3"});',
            ]
        )
    periodic = (
        f"model.study({java_string(study)})"
        f".feature({java_string(time_periodic_feature)})"
    )
    variables = f'model.sol({java_string(solution)}).feature("v1")'
    direct_solve_settings.extend(
        [
            "    // Cold start from the native physical initial state.",
            "    // No saved vector, manual field, or continuation is used.",
            f'    {periodic}.set("useinitsol", false);',
            f'    {periodic}.set("initmethod", "init");',
            f'    {periodic}.set("initstudy", "zero");',
            f'    {periodic}.set("initsol", "current");',
            f'    {variables}.set("initmethod", "init");',
            f'    {variables}.set("initsol", "zero");',
            f"    model.sol({java_string(solution)}).clearSolutionData();",
        ]
    )
    restore_native_equation_view: list[str] = []
    if uses_exact_periodic_locks:
        feature_info = (
            f"model.component({java_string(mapping.model.component)})"
            f".physics({java_string(mapping.model.physics)})"
            f".feature({java_string(mapping.model.plasma_feature)})"
            '.featureInfo("info")'
        )
        restore_native_equation_view = [
            "    // The exact algebraic locks are periodic-solver guards only.",
            "    // Restore COMSOL's native definitions before std2 and before",
            "    // persisting the model, so physical-time datasets use En/Ne.",
        ]
        for identifier in ("ptp.Mn", "ptp.ebar", "ptp.Te"):
            restore_native_equation_view.extend(
                [
                    f"    if ({feature_info}.removeLock("
                    f"{java_string(identifier)}) == null || "
                    f"{feature_info}.isLocked({java_string(identifier)})) {{",
                    "      throw new IllegalStateException("
                    + java_string(
                        "Failed to restore native Equation View expression: "
                        + identifier
                    )
                    + ");",
                    "    }",
                ]
            )
    restore_physical_closure_arguments = (
        _restore_physical_closure_argument_lines(mapping)
        if run_role == "external"
        else []
    )
    refresh_saved_solution_definitions = (
        [
            "    // Refresh stored solution definitions after switching the",
            "    // closure argument to ptp.ebar. This does not re-solve or",
            "    // provide an initial value; it makes existing periodic data",
            "    // evaluate the same constitutive functions after conversion.",
            "    for (String solutionTag : model.sol().tags()) {",
            "      model.sol(solutionTag).updateSolution();",
            "    }",
        ]
        if restore_physical_closure_arguments
        else []
    )
    lines.extend(
        [
            "    // Match the original model: one solve at its configured 1 W",
            "    // operating point, with no power or coefficient continuation.",
            *direct_solve_settings,
            "      model.param().set("
            + java_string(mapping.run.power_parameter)
            + f', "{mapping.run.power_W:.17g}[W]");',
            f"    model.study({java_string(study)}).run();",
            *restore_native_equation_view,
            *restore_physical_closure_arguments,
            *refresh_saved_solution_definitions,
            "    model.save(" + java_string(java_path(output_mph)) + ");",
        ]
    )
    if convert_periodic_solution:
        lines.extend(
            [
                "    // Convert the converged periodic solution to physical RF time.",
                f"    model.study({java_string(mapping.model.conversion_study)})"
                '.feature("tptd").set("notstudy", '
                f"{java_string(conversion_source_study)});",
                f"    model.study({java_string(mapping.model.conversion_study)})"
                '.feature("tptd").set("notstudystep", '
                f"{java_string(time_periodic_feature)});",
                f"    model.study({java_string(mapping.model.conversion_study)}).run();",
                "    model.save(" + java_string(java_path(output_mph)) + ");",
            ]
        )
    lines.extend(
        [
            "    ModelUtil.remove(" + java_string(class_name) + ");",
            "  }",
            "}",
            "",
        ]
    )
    return "\n".join(lines)


def generate_export_java(
    mapping: GecCcpMapping,
    *,
    class_name: str,
    input_mph: Path,
    output_dir: Path,
    period_dataset: str,
    waveform_dataset: str,
    include_external_closure_audit: bool = False,
) -> str:
    context = java_export_sections.build_export_java_context(
        mapping,
        output_directory=output_dir,
        period_dataset=period_dataset,
        waveform_dataset=waveform_dataset,
        include_external_closure_audit=include_external_closure_audit,
    )
    lines = java_support.java_header(class_name, input_mph, class_name)
    if include_external_closure_audit:
        lines.extend(java_export_sections.render_elastic_domain_assertion(mapping))
    lines.extend(java_export_sections.render_field_exports(context))
    lines.extend(java_export_sections.render_conservation_exports(context))
    lines.extend(java_export_sections.render_export_footer(class_name))
    return "\n".join(lines)


def _restore_physical_closure_argument_lines(
    mapping: GecCcpMapping,
) -> list[str]:
    """Make a converged closure evaluable on native physical-time datasets."""

    if not (
        _uses_transport(mapping.closure)
        or _uses_external_rates(mapping.closure)
        or _uses_external_elastic_energy_loss(mapping.closure)
    ):
        return []
    support = _active_closure_mean_energy_support(mapping)
    if not support.get("passed"):
        raise GecCcpWorkflowError(
            "cannot restore physical closure arguments without valid support"
        )
    component = java_string(mapping.model.component)
    target = f'model.component({component}).variable("swClosureVars")'
    base = "log(ptp.ebar/1[V])"
    expressions: dict[str, str] = {}
    if mapping.bundle.expected_source == "monte_carlo":
        transport_support = support.get("transport_intersection_mean_energy_eV")
        if _uses_transport(mapping.closure):
            if not isinstance(transport_support, list):
                raise GecCcpWorkflowError(
                    "Monte Carlo physical transport support is unavailable"
                )
            transport_minimum, transport_maximum = (
                closure_arguments.closure_argument_range(mapping, transport_support)
            )
            expressions["sw_logeps_transport"] = (
                closure_arguments.smooth_log_energy_argument(
                    base,
                    minimum_eV=transport_minimum,
                    maximum_eV=transport_maximum,
                )
            )
        for reaction in mapping.reactions:
            if not _reaction_uses_preintegrated_rate(mapping.closure, reaction):
                continue
            rate_support = support.get("rate_supports_mean_energy_eV", {}).get(
                reaction.process_type
            )
            if not isinstance(rate_support, list):
                raise GecCcpWorkflowError(
                    "Monte Carlo physical rate support is unavailable for "
                    f"{reaction.process_type}"
                )
            rate_minimum, rate_maximum = closure_arguments.closure_argument_range(
                mapping, rate_support
            )
            expressions[f"sw_logeps_rate_{reaction.process_type}"] = (
                closure_arguments.smooth_log_energy_argument(
                    base,
                    minimum_eV=rate_minimum,
                    maximum_eV=rate_maximum,
                )
            )
        if _uses_external_elastic_energy_loss(mapping.closure):
            elastic_support = support.get("elastic_energy_loss_support_mean_energy_eV")
            if not isinstance(elastic_support, list):
                raise GecCcpWorkflowError(
                    "Monte Carlo physical elastic-loss support is unavailable"
                )
            elastic_minimum, elastic_maximum = closure_arguments.closure_argument_range(
                mapping, elastic_support
            )
            expressions["sw_logeps_el"] = closure_arguments.smooth_log_energy_argument(
                base,
                minimum_eV=elastic_minimum,
                maximum_eV=elastic_maximum,
            )
    else:
        minimum_eV, maximum_eV = support["common_intersection_mean_energy_eV"]
        expressions["sw_logeps"] = closure_arguments.smooth_log_energy_argument(
            base,
            minimum_eV=float(minimum_eV),
            maximum_eV=float(maximum_eV),
        )
    return [
        "    // The nonlinear solve used En_per-Ne_per directly. After",
        "    // convergence, use COMSOL's native mean-energy variable so",
        "    // the same closure remains defined on std2 physical-time data.",
        *[
            f"    {target}.set({java_string(name)}, {java_string(expression)});"
            for name, expression in expressions.items()
        ],
    ]
